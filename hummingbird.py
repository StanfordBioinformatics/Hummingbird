from __future__ import division, print_function
import argparse
import configparser
import logging
import numpy as np
from sklearn import linear_model
from scipy.interpolate import Rbf, UnivariateSpline
import matplotlib.pyplot as plt

from downsample import Downsample
from profiling import Profiler
from instance import Instance
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

def regression(known_data, target_value):
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
        result.append(np.asscalar(regr.predict(target_value)))

        plt.subplot(len(ys[0]), 1, i + 1)
        plt.plot(x_test, regr.predict(x_test), 'b', X, y, 'ro')
    #plt.show()
    return result

def main():
    """The main pipeline."""
    parser = argparse.ArgumentParser(description='Process command line input')
    parser.add_argument('-c', '--conf', dest='config_file',
                        default='User_Provided_Input.conf')
    parser.add_argument('-d', '--downsample', dest='downsample_tool',
                        default='seqtk')
    parser.add_argument('-p', '--profiler', dest='profile_tool',
                        default='time')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_file)

    logging.basicConfig(format='%(asctime)s: %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.DEBUG)
    logging.info('Preparing downsampling...')
    downsampler = Downsample(args.downsample_tool, config)
    ds_dict = downsampler.subsample()
    logging.info('Downsampling done.')
    logging.debug(ds_dict)

    logging.info('Preparing memory profiling...')
    profiler = Profiler(args.profile_tool, 'mem', config)
    pf_dict = profiler.profile(ds_dict)
    logging.info('Memory profiling done.')
    logging.debug(pf_dict)

    target = int(config['Downsample']['target'])
    #print [np.asscalar(a) for a in extrapolate('spline', pf_dict, target)]
    predictions = regression(pf_dict, target)
    thread_list = config.get('Profiling', 'thread', fallback="4").split(',')
    thread_list = [t.strip() for t in thread_list]
    for i, t in enumerate(thread_list):
        print('The memory usage for {:,} reads with {:>2} threads is predicted as {:,.0f} Bytes.'.format(target, t, predictions[i]))

    reserved_mem = 8
    min_mem = [pred/1000/1000/1000 + reserved_mem for pred in predictions]
    valid, invalid = Instance.get_machine_types(config, min_mem)
    if valid:
        print('Sugguest to test on the following machine types:')
        for ins in valid:
            print(bcolors.OKGREEN + ins.name + bcolors.ENDC, str(ins.mem) + 'GB')
    else:
        print('No pre-defined machine types found.')
    if invalid:
        print('Machine types might not have enouph memory and will be pruned:')
        for ins in invalid:
            print(bcolors.FAIL + ins.name + bcolors.ENDC, str(ins.mem) + 'GB')
    print('Try customized machine type if necessary:')
    for i, t in enumerate(thread_list):
        print('min-core: {:2}\tmin-mem: {} GB'.format(t, min_mem[i]))

    logging.info('Preparing runtime profiling...')
    profiler = Profiler(args.profile_tool, 'time', config)
    ds_size = int(target * 0.01)
    rt_dict = profiler.profile({ds_size:ds_dict[ds_size]}, valid)
    logging.info('Runtime profiling done.')
    logging.debug(rt_dict)

if __name__ == "__main__":
    main()
