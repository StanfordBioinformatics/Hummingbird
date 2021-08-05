#!/usr/bin/env python3

import argparse
import json
import logging
import sys
from builtins import input
from pprint import pformat

from Hummingbird import validator
from Hummingbird.downsample import Downsample
from Hummingbird.hummingbird_utils import DOWNSAMPLE, PROFILING, PLATFORM, Predictor, bcolors, speedup_efficiency, \
    cost_efficiency
from Hummingbird.instance import Instance, GCPInstance
from Hummingbird.profiler import Profiler


def main():
    """The main pipeline."""
    parser = argparse.ArgumentParser()
    parser.add_argument('conf', help='Hummingbird JSON configuration')
    parser.add_argument('--fa_downsample', dest='downsample_tool',
                        choices=['seqtk', 'zless'],
                        default='seqtk',
                        help='the tool used for downsampling, default: %(default)s')
    parser.add_argument('-p', '--profiler', dest='profile_tool',
                        choices=['time', 'valgrind'],
                        default='time',
                        help='the tool used for profiling, default: %(default)s')
    args = parser.parse_args()

    with open(args.conf, 'r') as config_file:
        config = json.load(config_file)
        config[DOWNSAMPLE]['tool'] = args.downsample_tool

    logging.basicConfig(format='%(levelname)-8s :: %(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)

    try:
        validator.validate_config_file(config)
    except validator.ConfigurationError as e:
        logging.error(e)
        return

    logging.info('Preparing downsampling...')
    downsampler = Downsample(config)
    ds_dict = downsampler.subsample()
    logging.info('Downsampling done.')

    target = config[DOWNSAMPLE]['target']
    for i, workflow in enumerate(config[PROFILING]):
        # user_input = input('Do you want to continue? [y/N]: ')
        # if user_input != 'y' and user_input != 'Y':
        #     break
        if 'wdl_file' in workflow and 'backend_conf' in workflow:
            backend = 'cromwell'
        elif 'command' in workflow or 'script' in workflow:
            backend = 'bash'
        else:
            sys.exit('Invalid conf parameters.')

        logging.info('Preparing memory profiling...')
        logging.info(pformat(ds_dict))
        wf_conf = {PROFILING: workflow, DOWNSAMPLE: config[DOWNSAMPLE], PLATFORM: config[PLATFORM]}
        profiler = Profiler(backend, args.profile_tool, Profiler.mem_mode, wf_conf)
        profiling_dict = profiler.profile(ds_dict)
        logging.info('Memory profiling done.')
        logging.info(profiling_dict)

        all_valid = set()
        all_invalid = set()
        thread_list = workflow.get('thread', [2])
        predictor = Predictor(target, thread_list)
        for task in profiling_dict:
            print('==', task, '==')
            pf_dict = profiling_dict[task]
            predictions = predictor.extrapolate(pf_dict, task)
            for i, t in enumerate(thread_list):
                print('The memory usage for {:,} reads with {:>2} threads is predicted as {:,.0f} Kbytes.'.format(target, t, predictions[i]))
            reserved_mem = 2
            min_mem = [pred / 1000 / 1000 + reserved_mem for pred in predictions]
            valid, invalid = Instance.get_machine_types(wf_conf, min_mem)
            all_valid.update(valid)
            all_invalid.update(invalid)
            if valid:
                print('Sugguest to test on the following machine types:')
                for ins in valid:
                    print(bcolors.OKGREEN + ins.name + bcolors.ENDC, str(ins.mem) + 'GB')
                    ins.set_price()
            else:
                print('No pre-defined machine types found.')
            if invalid:
                print('Machine types might not have enouph memory and will be pruned:')
                for ins in invalid:
                    print(bcolors.FAIL + ins.name + bcolors.ENDC, str(ins.mem) + 'GB')
        all_valid.difference_update(all_invalid)
        cus = input('Do you want to include customized machine types? [y/N]: ')
        if cus.lower() == 'y' or cus.lower() == 'yes':
            cus_types = []
            while True:
                ins_name = input('instance types:')
                if ins_name:
                    cus_types.append(GCPInstance(name=ins_name))
                else:
                    break
            for ins in cus_types:
                ins.set_price()
            all_valid.update(cus_types)

        logging.info(all_valid)
        # if len(all_valid) <= 1:
        #     continue
        logging.info('Preparing runtime profiling...')
        profiler.mode = Profiler.time_mode

        multiplier = 1
        if wf_conf[DOWNSAMPLE].get('fullrun', False):
            ds_size = target
        else:
            ds_size = int(target * Downsample.runtime_frac * multiplier)
        while True:
            logging.info('Current downsample size: %s', str(ds_size))
            runtimes_dict = profiler.profile({ds_size: ds_dict[ds_size]}, all_valid)
            logging.info('Runtime profiling done.')
            logging.info(pformat(runtimes_dict))
            for task in runtimes_dict:
                print('==' + task + '==')
                runtimes = runtimes_dict[task][ds_size]
                zipped_runtimes = [(t, m) for t, m in zip(runtimes, all_valid) if t]
                runtimes, succeeded = [], []
                for t, m in zipped_runtimes:
                    runtimes.append(t)
                    succeeded.append(m)
                speedups = speedup_efficiency(zipped_runtimes)
                sorted_runtimes = sorted(zipped_runtimes, key=lambda x: x[0])
                prices = [ins.price for ins in succeeded]
                costs = [t * p for t, p in zip(runtimes, prices)]
                sorted_costs = sorted(zip(costs, succeeded))
                efficiencies = cost_efficiency(runtimes, costs)
                sorted_efficiencies = sorted(zip(efficiencies, succeeded), reverse=True)
                print(
                    'The fastest machine type: {}'.format(bcolors.OKGREEN + sorted_runtimes[0][1].name + bcolors.ENDC))
                print('The cheapest machine type: {}'.format(bcolors.OKGREEN + sorted_costs[0][1].name + bcolors.ENDC))
                print('The most cost-efficient machine type: {}'.format(
                    bcolors.OKGREEN + sorted_efficiencies[0][1].name + bcolors.ENDC))
            if ds_size >= int(target * 0.1) or input("Try larger downsample size? [Y/n]") == 'n':
                if profiler.output_dict is not None:  # CromwellProfiler returns None output_dict
                    ds_dict = profiler.output_dict  # update downsampled input as output from previous workflow
                break
            multiplier *= 10
            ds_size = int(target * Downsample.runtime_frac * multiplier)
            if ds_size not in ds_dict:
                logging.warning(
                    f"Unable to run profiling for downsample size of {ds_size} as it is not present in the downsample configuration.")
                break


if __name__ == "__main__":
    main()
