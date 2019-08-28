import numpy as np

FA_EXT = ['.fa', '.fasta']
FQ_EXT = ['.fq', '.fastq']
SAM_EXT = ['.sam']
BAM_EXT = ['.bam', 'ubam', 'cram']
ZIP_EXT = ['.gz']
FA = 'fa'
FQ = 'fq'
SAM = 'sam'
BAM = 'bam'

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
