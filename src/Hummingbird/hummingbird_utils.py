import os
import sys
import numpy as np
from scipy.interpolate import Rbf, UnivariateSpline
from sklearn import linear_model
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

FA_EXT = ['.fa', '.fasta']
FQ_EXT = ['.fq', '.fastq']
SAM_EXT = ['.sam']
BAM_EXT = ['.bam', 'ubam', 'cram']
ZIP_EXT = ['.gz']
FA = 'fa'
FQ = 'fq'
SAM = 'sam'
BAM = 'bam'

PLATFORM = 'Platform'
DOWNSAMPLE = 'Downsample'
PROFILING = 'Profiling'

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def humanize(num):
    """A utility function to help generate human readable number string"""
    try:
        num = int(num)
    except:
        sys.exit("Unalbe to humanize input value.")
    for unit in ['', 'K', 'M']:
        if num % 1000:
            return '%d%s' % (num, unit)
        else:
            num /= 1000
    return "%d%s" % (num, 'G')

def cost_efficiency(run_times, costs):
    """Compute the cost efficiency as the metrics of instance."""
    speedups = np.reciprocal(run_times)
    return speedups / costs

def speedup_efficiency(zipped_runtimes):
    zipped_runtimes.sort(key=lambda t:(t[1].cpu, t[0]))
    print(zipped_runtimes)
    base_time = zipped_runtimes[0][0]
    base_core = zipped_runtimes[0][1].cpu
    for t, m in zipped_runtimes:
        ideal = float(m.cpu) / base_core
        speedup = base_time / t
        efficiency = speedup / ideal
        print(speedup, ideal, efficiency)

def spline(method, known_data, target_value):
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


class Predictor(object):
    def __init__(self, target, threads):
        self.target = target
        self.threads = threads

    def extrapolate(self, known_data, task):
        x = np.array(list(known_data.keys()))
        ys = np.array(list(known_data.values())).transpose()
        x = x.reshape(-1, 1)
        x_log = np.log(x)
        predictions = []
        x_test = np.linspace(1000, self.target, num=100)
        x_test = x_test.reshape(-1, 1)
        x_test_log = np.log(x_test)
        # Reshape data using array.reshape(-1, 1) if your data has a single feature
        target = np.array(self.target).reshape(-1, 1)
        target_log = np.log(target)

        def compute_model(x, y, target):
            regr = linear_model.LinearRegression()
            regr.fit(x, y)
            r2 = regr.score(x, y)
            pred = max(np.max(y), np.asscalar(regr.predict(target)))
            return (pred, r2)

        for i, y in enumerate(ys):
            pred, r2 = compute_model(x, y, target)
            pred_log, r2_log = compute_model(x_log, y, target_log)
            print('prediction:', pred, pred_log)
            print('R2:', r2, r2_log)
            predictions.append(pred_log)

        return predictions
