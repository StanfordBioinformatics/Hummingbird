import sys
import numpy as np
from sklearn import linear_model
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

colors = {
    'BWA-MEM': '#3498DB',
    'SortSam': '#E74C3C',
    'MarkDuplicate': 'pink',
    'BaseRecalibrator': 'green',
    'ApplyBQSR': 'orange',
    'HaplotypeCaller': 'purple',
    'Mutect2': 'blue',
    'Trim_adapter': '#1f77b4',
    'Bowtie2': '#ff7f0e',
    'Filter': '#2ca02c',
    'Bam2ta': '#d62728',
    'Xcor': '#9467bd',
    'Macs2': '#8c564b',
    'Spr': '#e377c2',
    'Overlap': '#7f7f7f',
    'Idr': '#bcbd22'
}

def extrapolate(known_data, actual_data, target, task):
    x = np.array(list(known_data.keys()), dtype=np.float64)
    ys = np.array(list(known_data.values())).transpose()
    ys /= 1000000.0
    actual_data = [truth / 1000000.0 for truth in actual_data]
    x = x.reshape(-1, 1)
    x /= target
    x *= 100
    x = np.log10(x)
    x_test = np.linspace(0.1, 100, num=10)
    x_test = x_test.reshape(-1, 1)
    x_test = np.log10(x_test)
    for i, y in enumerate(ys):
        regr = linear_model.LinearRegression()
        regr.fit(x, y)
        # Reshape data using array.reshape(-1, 1) if your data has a single feature
        target = np.array(100).reshape(-1, 1)
        target = np.log10(target)
        pred = max(np.max(y), np.asscalar(regr.predict(target)))
        #plt.subplot(len(ys), 1, i + 1)
        if len(ys) == 3:
            if i == 0:
                linestyle = ':'
                marker = 'x'
                lable = task + ' 8 vCPUs'
            elif i == 1:
                linestyle = '--'
                marker = '^'
                lable = task + ' 16 vCPUs'
            else:
                linestyle = '-'
                marker = 'o'
                lable = task + ' 32 vCPUs'
        elif len(ys) == 2:
            if i == 0:
                linestyle = ':'
                marker = 'x'
                lable = task + ' 2 vCPUs'
            else:
                linestyle = '-'
                marker = 'o'
                lable = task + ' 4 vCPUs'
        elements = []
        regression_line, = plt.plot(x_test, regr.predict(x_test), color=colors[task], linestyle=linestyle)
        datapoint_plot, = plt.plot(x, y, marker, color=colors[task], markersize=10)
        legend_handles.append((regression_line, datapoint_plot))
        legend_lables.append(lable)
        #plt.plot(target, actual_data[i]], 'ro')
        yerr = [[0], [pred - actual_data[i]]]
        plt.errorbar(target, actual_data[i], fmt=marker, markersize=10, color=colors[task], yerr=yerr, capsize=4, ecolor='black')
        print(plt.axis())
        xmin, xmax, ymin, ymax = plt.axis()
        offset = (ymax-ymin)*0.02*i
        plt.annotate("%.2f%%" % (100*(pred-actual_data[i])/actual_data[i]), xy=(target+0.1, (actual_data[i]+pred)/2+offset))
    ax = plt.gca()
    ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.1f}'))
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams.update({'font.size': 22})
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

def plot_haplotypecaller():
    plt.figure(figsize=(16, 9))

    # data = {70764612: [10537584.0, 15552648.0, 25266520.0], 7076461: [9060508.0, 10148120.0, 11261328.0], 707646: [6508520.0, 7113960.0, 7186892.0]}
    # truth = [11936674.67, 17565605.33, 26138992]
    # extrapolate(data, truth, 707646124, 'BWA-MEM')

    # data = {308272: [3584592.0, 3792768.0, 4197804.0], 30827265: [3644068.0, 3913304.0, 4427708.0], 3082726: [3611956.0, 3846216.0, 4300156.0]}
    # truth = [3676768, 3971904, 4544272]
    # extrapolate(data, truth, 308272658, 'blue')

    # data = {70764612: [5966824.0, 6064188.0], 7076461: [6106232.0, 6132024.0], 707646: [3536596.0, 3404308.0]}
    # truth = [6198496, 6183356]
    # extrapolate(data, truth, 707646124, 'SortSam')
    #
    # data = {70764612: [8725256.0, 8789432.0], 7076461: [6575720.0, 6625776.0], 707646: [5087004.0, 5155976.0]}
    # truth = [8755872, 8810644]
    # extrapolate(data, truth, 707646124, 'MarkDuplicate')
    #
    # data = {70764612: [19626.666666666668, 19570.666666666668], 7076461: [17226.666666666664, 17242.666666666664], 707646: [15429.333333333332, 15513.333333333336]}
    # truth = [23624, 23776]
    # extrapolate(data, truth, 707646124, 'BaseRecalibrator')
    #
    # data = {70764612: [2371590.6666666665, 2686521.3333333335], 7076461: [2653056.0, 2700041.333333333], 707646: [2284984.0, 2364088.0]}
    # truth = [2358020, 2455488]
    # extrapolate(data, truth, 707646124, 'ApplyBQSR')

    data = {1417289: [12013737.333333332, 6333733.333333334, 6597452.0], 14172891: [12020142.666666666, 12626956.0, 12860710.666666668], 141728916: [12651493.333333332, 13344673.333333336, 13693002.666666668]}
    truth = [9262152, 13120061.3333333, 13665805.3333333]
    extrapolate(data, truth, 1417289160, 'HaplotypeCaller')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('GATK HaplotypeCaller Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('haplotypecaller.pdf')

def plot_mutect():
    plt.figure(figsize=(16, 9))

    data = {811944: [5063106.666666667, 6336018.666666666, 11067504.0], 8119448: [5400526.666666666, 7497428.0, 10915470.666666666], 81194: [5161272.0, 6537250.666666667, 11535181.333333334]}
    truth = [6555428, 11945592, 14540496]
    extrapolate(data, truth, 81194486, 'Mutect2')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('GATK Mutect2 Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('mutect2.pdf')

def plot_trim():
    plt.figure(figsize=(16, 9))

    data = {308272: [13356.0, 13304.0, 13340.0], 30827265: [13344.0, 13304.0, 13320.0], 3082726: [13284.0, 13344.0, 13292.0]}
    truth = [13320, 13316, 13340]
    extrapolate(data, truth, 308272658, 'Trim_adapter')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('trim.pdf')

def plot_filter():
    plt.figure(figsize=(16, 9))

    data = {308272: [2068312.0, 2169028.0, 2339160.0], 30827265: [3667928.0, 4155688.0, 3868056.0], 3082726: [2981972.0, 3029340.0, 3026776.0]}
    truth = [4375888,4421372,4886696]
    extrapolate(data, truth, 308272658, 'Filter')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('ATAC-seq Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('filter.pdf')

def plot_atac1():
    plt.figure(figsize=(16, 9))

    data = {308272: [196456.0, 198480.0, 207344.0], 30827265: [2128712.0, 2132856.0, 2139284.0], 3082726: [1115436.0, 1115656.0, 1118124.0]}
    truth = [2383024,2386120,2394704]
    extrapolate(data, truth, 308272658, 'Bam2ta')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('ATAC-seq Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('atac1.pdf')

def plot_atac2():
    plt.figure(figsize=(16, 9))

    data = {308272: [45408.0, 44964.0, 44956.0], 30827265: [3050292.0, 3050272.0, 3050292.0], 3082726: [314820.0, 314868.0, 314760.0]}
    truth = [4298268,4298364,4298256]
    extrapolate(data, truth, 308272658, 'Xcor')

    data = {308272: [141496.0, 141496.0], 30827265: [9636716.0, 12046012.0], 3082726: [1111720.0, 1389480.0]}
    truth = [11698440,11706940]
    extrapolate(data, truth, 308272658, 'Macs2')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('ATAC-seq Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('atac2.pdf')

def plot_spr():
    plt.figure(figsize=(16, 9))

    data = {308272: [16204.0, 16196.0], 30827265: [1419480.0, 1419468.0], 3082726: [147184.0, 147296.0]}
    truth = [11752640,11752540]
    extrapolate(data, truth, 308272658, 'Spr')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('ATAC-seq Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('spr.pdf')

def plot_atac3():
    plt.figure(figsize=(16, 9))

    data = {308272: [81656.0, 81668.0], 30827265: [182228.0, 182276.0], 3082726: [225696.0, 225780.0]}
    truth = [222320,222312]
    extrapolate(data, truth, 308272658, 'Overlap')

    data = {308272: [231160.0, 231004.0], 30827265: [671076.0, 671388.0], 3082726: [513376.0, 512628.0]}
    truth = [887312,887536]
    extrapolate(data, truth, 308272658, 'Idr')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('ATAC-seq Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('atac3.pdf')

def plot_bowtie():
    plt.figure(figsize=(16, 9))

    data = {308272: [3584592.0, 3792768.0, 4197804.0], 30827265: [3644068.0, 3913304.0, 4427708.0], 3082726: [3611956.0, 3846216.0, 4300156.0]}
    truth = [3676768,3971904,4544272]
    extrapolate(data, truth, 308272658, 'Bowtie2')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('Bowtie2 Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('bowtie.pdf')

def plot_bwa():
    plt.figure(figsize=(16, 9))

    data = {70764612: [10537584.0, 15552648.0, 25266520.0], 7076461: [9060508.0, 10148120.0, 11261328.0], 707646: [6508520.0, 7113960.0, 7186892.0]}
    truth = [11936674.67, 17565605.33, 26138992]
    extrapolate(data, truth, 707646124, 'BWA-MEM')

    plt.legend(legend_handles, legend_lables)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('BWA-MEM Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('bwa.pdf')

def plot_preprocess():
    plt.figure(figsize=(16, 9))

    data = {70764612: [5966824.0, 6064188.0], 7076461: [6106232.0, 6132024.0], 707646: [3536596.0, 3404308.0]}
    truth = [6198496, 6183356]
    extrapolate(data, truth, 707646124, 'SortSam')

    data = {70764612: [8725256.0, 8789432.0], 7076461: [6575720.0, 6625776.0], 707646: [5087004.0, 5155976.0]}
    truth = [8755872, 8810644]
    extrapolate(data, truth, 707646124, 'MarkDuplicate')

    data = {70764612: [19626.666666666668, 19570.666666666668], 7076461: [17226.666666666664, 17242.666666666664], 707646: [15429.333333333332, 15513.333333333336]}
    truth = [23624, 23776]
    extrapolate(data, truth, 707646124, 'BaseRecalibrator')

    data = {70764612: [2371590.6666666665, 2686521.3333333335], 7076461: [2653056.0, 2700041.333333333], 707646: [2284984.0, 2364088.0]}
    truth = [2358020, 2455488]
    extrapolate(data, truth, 707646124, 'ApplyBQSR')

    plt.legend(legend_handles, legend_lables, fontsize=15)
    plt.xticks([-1, 0, 1, 2], ('0.1%', '1%', '10%', '100%'))
    plt.title('GATK Data Pre-processing Memory Profiling and Prediction')
    plt.xlabel('Downsample Percentage')
    plt.ylabel('Memory Usage (GB)')
    plt.savefig('preprocess.pdf')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        exit('Please provide pipeline name')
    legend_handles = []
    legend_lables = []
    if sys.argv[1] == 'haplotypecaller':
        plot_haplotypecaller()
    elif sys.argv[1] == 'mutect2':
        plot_mutect()
    elif sys.argv[1] == 'atac':
        plot_trim()
        legend_handles = []
        legend_lables = []
        plot_filter()
        legend_handles = []
        legend_lables = []
        plot_atac1()
        legend_handles = []
        legend_lables = []
        plot_atac2()
        legend_handles = []
        legend_lables = []
        plot_spr()
        legend_handles = []
        legend_lables = []
        plot_atac3()
    elif sys.argv[1] == 'bwa':
        plot_bwa()
    elif sys.argv[1] == 'preprocess':
        plot_preprocess()
    elif sys.argv[1] == 'bowtie':
        plot_bowtie()
    else:
        exit('Invalid name')
