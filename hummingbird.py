from __future__ import division, print_function
from builtins import input
import argparse
import json
import logging
import math
import numpy as np
from sklearn import linear_model
from scipy.interpolate import Rbf, UnivariateSpline
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from downsample import Downsample
from profiling import Profiler
from instance import *
from hummingbird_utils import *

def extrapolate(method, known_data, target_value):
    """Use scipy interpolate module to extrapolate target value and plot."""
    x = known_data.keys()
    ys = np.array(known_data.values())
    result = []
    x_test = np.linspace(1000, target_value, num=100)
    plt.figure(figsize=(18, 10))
    for i in range(len(ys[0])):
        y = ys[:,i]
        if method == 'spline':
            f = UnivariateSpline(x, y, k=1, ext='extrapolate')
        elif method == 'rbf':
            f = Rbf(x, y)
        result.append(f(target_value))

        plt.subplot(len(ys[0]), 1, i + 1)
        plt.plot(x_test, f(x_test), 'b', x, y, 'ro')
    plt.show()
    return result

def regression(known_data, target_value, filename='plot/hummingbird.png'):
    """Use sklearn regression module to predict target value and plot."""
    X = np.array(known_data.keys())
    X = X.reshape((len(known_data), 1))
    ys = np.array(known_data.values())
    result = []
    x_test = np.linspace(1000, target_value, num=100)
    x_test = x_test.reshape(100, 1)
    plt.figure(figsize=(18, 10))
    for i in range(len(ys[0])):
        y = ys[:,i]
        regr = linear_model.LinearRegression()
        regr.fit(X, y)
        # Reshape data using array.reshape(-1, 1) if your data has a single feature
        target = np.array(target_value).reshape(-1, 1)
        res = max(np.max(y), np.asscalar(regr.predict(target)))
        result.append(res)

        plt.subplot(len(ys[0]), 1, i + 1)
        plt.plot(x_test, regr.predict(x_test), 'b', X, y, 'ro')
    plt.savefig(filename)
    return result

def main():
    """The main pipeline."""
    parser = argparse.ArgumentParser()
    parser.add_argument('conf', help='Hummingbird configuration')
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
        config['Downsample']['tool'] = args.downsample_tool

    logging.basicConfig(format='%(asctime)s: %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO)
    logging.info('Preparing downsampling...')
    downsampler = Downsample(config)
    ds_dict = downsampler.subsample()
    logging.info('Downsampling done.')
    logging.info(ds_dict)

    for i, workflow in enumerate(config['Profiling']):
        cont = input('Continue? (y/N): ')
        if cont != 'y':
            break
        if 'wdl_file' in workflow and 'backend_conf' in workflow:
            backend = 'cromwell'
        elif 'command' in workflow or 'script' in workflow:
            backend = 'bash'
        else:
            sys.exit('Invalid conf parameters.')

        logging.info('Preparing memory profiling...')
        wf_conf = {}
        wf_conf['Profiling'] = config['Profiling'][i]
        wf_conf['Downsample'] = config['Downsample']
        wf_conf['Platform'] = config['Platform']
        profiler = Profiler(backend, args.profile_tool, 'mem', wf_conf)
        profiling_dict = profiler.profile(ds_dict)
        logging.info('Memory profiling done.')
        logging.info(profiling_dict)

        all_valid = set()
        all_invalid = set()
        for task in profiling_dict:
            print('==', task, '==')
            pf_dict = profiling_dict[task]
            target = config['Downsample']['target']
            #print [np.asscalar(a) for a in extrapolate('spline', pf_dict, target)]
            predictions = regression(pf_dict, target, 'plot/' + task + '.png')
            thread_list = workflow.get('thread', [4])
            #thread_list = [t.strip() for t in thread_list]
            for i, t in enumerate(thread_list):
                print('The memory usage for {:,} reads with {:>2} threads is predicted as {:,.0f} Kbytes.'.format(target, t, predictions[i]))
                reserved_mem = 1
            min_mem = [pred/1000/1000 + reserved_mem for pred in predictions]
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
            print('Try customized machine type if necessary:')
            cus_types = []
            for i, t in enumerate(thread_list):
                print('min-core: {:2}\tmin-mem: {} GB'.format(t, min_mem[i]))
                cus_types.append(GCP_Instance('custom'+'-'+str(t), t, math.ceil(min_mem[i])))
            # cus = input('Do you want to include customized machine types? (Y/N): ')
            # if cus.lower() == 'y' or cus.lower() == 'yes':
            #     for ins in cus_types:
            #         price = input('Hourly price for {} core and {}GB memory:'.format(ins.cpu, ins.mem))
            #         ins.set_price(price)
                #valid += cus_types

        all_valid.difference_update(all_invalid)
        logging.info('Preparing runtime profiling...')
        profiler = Profiler(backend, args.profile_tool, 'time', wf_conf)
        # TODO: target unassinged if profiling_dict is empty
        if config['Downsample'].get('fullrun', False):
            ds_size = target
        else:
            ds_size = int(target * Downsample.runtime_size)
        runtimes_dict = profiler.profile({ds_size:ds_dict[ds_size]}, all_valid)
        logging.info('Runtime profiling done.')
        logging.info(runtimes_dict)
        for task in runtimes_dict:
            print('==' + task + '==')
            runtimes = runtimes_dict[task][ds_size]
            sorted_runtimes = sorted(zip(runtimes, all_valid))
            prices = [ins.price for ins in all_valid]
            costs = [t * p for t, p in zip(runtimes, prices)]
            #costs = [(t * ins.price, ins) for t, ins in sorted_runtimes]
            sorted_costs = sorted(zip(costs, all_valid))
            efficiencies = cost_efficiency(runtimes, costs)
            sorted_efficiencies = sorted(zip(efficiencies, all_valid), reverse=True)
            print('The fastest machine type: {}'.format(sorted_runtimes[0][1].name))
            print('The cheapest machine type: {}'.format(sorted_costs[0][1].name))
            print('The most cost-efficient machine type: {}'.format(sorted_efficiencies[0][1].name))

if __name__ == "__main__":
    main()
