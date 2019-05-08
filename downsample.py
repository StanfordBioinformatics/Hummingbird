import os
import csv
import sys
from collections import defaultdict
from scheduler import Scheduler
from hummingbird_utils import *

# COMMAND = '''
# dsub \
# --provider google-v2 \
# --project gbsc-gcp-project-cba \
# --regions us-west1 \
# --image xingziye/seqtk:latest \
# --command 'seqtk sample ${INPUT_FILE} ${COUNT} > ${OUTPUT_FILE}' \
# --tasks dsub.tsv
# '''
# TEST_COMMAND = 'dsub --help'
MACHINE_TYPE = 'n1-standard-4'
DISK_SIZE = '500'
SIZE_RUNTIME = 0.01

class Downsample(object):
    """Generate subsamples for input data.

    Attributes:
        tool: A string describing downsample methods.
        conf: A dictionary of configrations.
    """
    def __init__(self, tool, conf):
        self.tool = tool
        self.conf = conf

    def subsample(self):
        """Generate subsamples for input data according to the input format."""
        input_dict = self.conf['Downsample']['input']
        fa_list = list()
        for path in input_dict.values():
            path = path.strip()
            base, extension = os.path.splitext(path)
            while extension in ZIP_EXT:
                base, extension = os.path.splitext(base)
            if extension.lower() in FA_EXT:
                fa_list.append(path)
            else:
                sys.exit("Unsupported input format.")
        return self.downsample_fa(fa_list)

    def downsample_fa(self, filenames):
        """Generate FASTA/FASTQ downsampled files.

        Each file specified in the filenames will be subsampled by the sizes
            specified in conf.

        Args:
            filenames: A list of paths which are the input files that
                need to be subsampled.

        Returns:
            A dict mapping the number of downsampled size to a list of downsampled
                filenames.
        """
        entry_counts = self.conf['Downsample'].get('size', [0.001, 0.01, 0.1])
        if any([count < 1 for count in entry_counts]):
            self.tool = 'seqtk'
        # 1% downsample for runtime profiling
        if self.tool == 'seqtk' and SIZE_RUNTIME not in entry_counts:
            entry_counts.append(SIZE_RUNTIME)
        elif self.tool == 'zless' and self.conf['Downsample']['target'] * SIZE_RUNTIME not in entry_counts:
            entry_counts.append(self.conf['Downsample']['target'] * SIZE_RUNTIME)
        bucket_path = 'gs://' + self.conf['Platform']['bucket']
        output_path = self.conf['Downsample']['output'].strip('/')
        log_path = bucket_path + '/' + self.conf['Downsample']['logging']
        downsampled = defaultdict(list)

        with open('downsample.tsv', 'w') as dsub_tsv:
            tsv_writer = csv.writer(dsub_tsv, delimiter='\t')
            tsv_writer.writerow(['--env COUNT', '--input INPUT_FILE', '--output OUTPUT_FILE'])

            for filename in filenames:
                base, extension = os.path.splitext(filename)
                while extension in ZIP_EXT:
                    base, extension = os.path.splitext(base)
                for entry_count in entry_counts:
                    if isinstance(entry_count, float):
                        count_int = int(self.conf['Downsample']['target'] * entry_count)
                    else:
                        count_int = entry_count
                    target_file = os.path.basename(base) + '_' + self.tool + '_' + humanize(count_int) + extension
                    target_path = '/'.join([bucket_path, output_path, target_file])
                    downsampled[count_int].append(target_path)
                    if self.tool == 'seqtk':
                        tsv_writer.writerow([entry_count, filename, target_path])
                    elif self.tool == 'zless':
                        tsv_writer.writerow([entry_count*4, filename, target_path])
        scheduler = Scheduler('dsub', self.conf)
        if self.tool == 'seqtk':
            scheduler.add_argument('--image', 'xingziye/seqtk:latest')
            scheduler.add_argument('--command', "'seqtk sample ${INPUT_FILE} ${COUNT} > ${OUTPUT_FILE}'")
        elif self.tool == 'zless':
            scheduler.add_argument('--image', 'xingziye/seqtk:latest')
            scheduler.add_argument('--command', "'zless ${INPUT_FILE} | head -n ${COUNT} > ${OUTPUT_FILE}'")
        scheduler.add_argument('--tasks', 'downsample.tsv')
        scheduler.add_argument('--logging', log_path)
        scheduler.add_argument('--machine-type', MACHINE_TYPE)
        scheduler.add_argument('--disk-size', DISK_SIZE)
        scheduler.add_argument('--wait')
        scheduler.add_argument('--skip')
        p = scheduler.run()
        p.wait()
        return downsampled
