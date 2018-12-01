import argparse
import configparser
import numpy as np
from scipy.interpolate import Rbf, UnivariateSpline
from sklearn import linear_model
import matplotlib.pyplot as plt
from downsample import Downsample
from profiling import Profiler

def extrapolate(method, known_data, target_value):
    """Use scipy interpolate module to extrapolate target value and plot."""
    x = known_data.keys()
    ys = np.array(known_data.values())
    result = []
    x_test = np.linspace(1000, 100000000, num=100)
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
    x_test = np.linspace(1000, 100000000, num=100)
    x_test = x_test.reshape(100, 1)
    plt.figure(figsize=(18, 10))
    for i in range(len(ys[0])):
        y = ys[:,i]
        regr = linear_model.LinearRegression()
        regr.fit(X, y)
        result.append(regr.predict(target_value))

        plt.subplot(len(ys[0]), 1, i + 1)
        plt.plot(x_test, regr.predict(x_test), 'b', X, y, 'ro')
    plt.show()
    return result

def main():
    """The main pipeline."""
    parser = argparse.ArgumentParser(description='Process command line input')
    parser.add_argument('-c', '--conf', dest='config_file',
                        default='User_Provided_Input.conf')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_file)

    downsampler = Downsample('zless', config)
    ds_dict = downsampler.subsample()

    profiler = Profiler('time', config)
    pf_dict = profiler.profile(ds_dict)

    print pf_dict
    print extrapolate('spline', pf_dict, 700000000)
    print regression(pf_dict, 700000000)

if __name__ == "__main__":
    main()
