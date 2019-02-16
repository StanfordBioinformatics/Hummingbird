FA_EXT = ['.fa', '.fasta', '.fq', '.fastq']
ZIP_EXT = ['.gz']

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
    if not isinstance(num, int) and num.isdigit():
        num = int(num)
    else:
        return num
    for unit in ['', 'K', 'M']:
        if num % 1000:
            return '%d%s' % (num, unit)
        else:
            num /= 1000
    return "%d%s" % (num, 'G')
