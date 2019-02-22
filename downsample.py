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
        filenames = self.conf['Downsample']['input'].splitlines()
        file_dict = defaultdict(list)
        for filename in filenames:
            filename = filename.strip()
            base, extension = os.path.splitext(filename)
            while extension in ZIP_EXT:
                base, extension = os.path.splitext(base)
            if extension.lower() in FA_EXT:
                file_dict['fa'].append(filename)
        return self.downsample_fa(file_dict['fa'])

    def downsample_fa(self, filenames):
        """Generate FASTA/FASTQ downsampled files.

        Each file specified in the filenames will be subsampled by the sizes
            specified in conf.

        Args:
            filenames: A list of strings which are the filenames of input that
                need to be subsampled.

        Returns:
            A dict mapping a string of downsampled size to a list of downsampled
                filenames.
        """
        entry_counts = self.conf['Downsample']['size']
        if entry_counts == 'auto':
            entry_counts = ['0.001', '0.01', '0.1']
            self.tool = 'seqtk'
        else:
            entry_counts = entry_counts.split(',')
        entry_counts = [entry.strip() for entry in entry_counts]
        # 1% downsample for runtime profiling
        if '0.01' not in entry_counts:
            entry_counts.append('0.01')
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
                    if not entry_count.isdigit():
                        try:
                            target_size = int(self.conf['Downsample']['target'])
                            entry_count = int(float(entry_count) * target_size)
                            # numerical value exists according to equivalent fraction value
                            if str(entry_count) in entry_counts:
                                continue
                        except ValueError:
                            sys.exit("The downsample input size is invalid.")
                    if self.tool == 'seqtk':
                        target_file = os.path.basename(base) + '_seqtk_' + humanize(entry_count) + extension
                    elif self.tool == 'zless':
                        target_file = os.path.basename(base) + '_zless_' + humanize(entry_count) + extension
                    target_path = '/'.join([bucket_path, output_path, target_file])
                    downsampled[entry_count].append(target_path)
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
