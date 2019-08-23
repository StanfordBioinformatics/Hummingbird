import os
import csv
import sys
from collections import defaultdict
from scheduler import Scheduler
from hummingbird_utils import *

class Downsample(object):
    """Generate subsamples for input data.

    Attributes:
        conf: A dictionary of configrations.
        tool: A string describing downsample methods for input files.
        target: An integer indicating the size of whole input.
        sizes: An array of downsample size parameters for tool.
        counts: An array of actual downsampled sample size.
    """
    runtime_size = 0.01
    default_instance = 'n1-highmem-8'
    default_disksize = '500'

    def __init__(self, conf):
        self.conf = conf
        self.tool = conf['Downsample']['tool']
        self.target = conf['Downsample']['target']
        self.sizes = conf['Downsample'].get('size', [])
        fullrun = conf['Downsample'].get('fullrun', False)
        if fullrun:
            self.sizes += [self.target]
        if not self.sizes: # when neither size nor fullrun specified
            self.sizes = [0.001, 0.01, 0.1]
        if any([s < 1 for s in self.sizes]):
            self.tool = 'seqtk'
        if not fullrun: # TODO pickard
            if self.tool == 'seqtk' and Downsample.runtime_size not in self.sizes:
                self.sizes.append(Downsample.runtime_size)
            elif self.tool == 'zless' and int(self.target * Downsample.runtime_size) not in self.sizes:
                self.sizes.append(int(self.target * Downsample.runtime_size))
        self.counts = []
        for s in self.sizes:
            self.counts.append(int(self.target * s) if s < 1 else s)
        conf['Downsample']['size'] = self.sizes
        conf['Downsample']['count'] = self.counts

    def subsample(self):
        """Generate subsamples for input data according to the input format."""
        input_dict = self.conf['Downsample']['input']
        type_dict = defaultdict(dict)
        for key in input_dict:
            path = input_dict[key]

            base, extension = os.path.splitext(path)
            while extension in ZIP_EXT:
                base, extension = os.path.splitext(base)
            if extension.lower() in FA_EXT:
                type_dict[FA][key] = path
            elif extension.lower() in SAM_EXT:
                type_dict[SAM][key] = path
            else:
                sys.exit("Unsupported input format.")

        downsampled = {}
        for type in type_dict:
            downsampled.update(self.downsample_by_type(type_dict[type], type))
        return downsampled

    def downsample_by_type(self, filenames, type):
        """Execute downsampling job by types.

        Each file specified in the filenames will be subsampled by the sizes
            specified in conf.

        Args:
            filenames: A dict mapping key name to the path of the file.
            type: A string tells which method to use to downsample.

        Returns:
            A dict mapping the downsampling size to the input dict with the
                corresponding downsampled files' path.
        """
        # This condition avoids empty dsub tsv input.
        if len(self.counts) == 1 and self.counts[0] == self.conf['Downsample']['target']:
            return {self.conf['Downsample']['target']: filenames}
        if type == SAM:
            self.tool = 'samtools'
        bucket_dir = 'gs://' + self.conf['Platform']['bucket']
        output_path = self.conf['Downsample']['output'].strip('/')
        log_path = bucket_dir + '/' + self.conf['Downsample']['logging']
        downsampled = defaultdict(dict)

        with open('downsample.tsv', 'w') as dsub_tsv:
            tsv_writer = csv.writer(dsub_tsv, delimiter='\t')
            tsv_writer.writerow(['--env COUNT', '--input INPUT_FILE', '--output OUTPUT_FILE'])

            for key in filenames:
                filename = filenames[key]
                base, extension = os.path.splitext(filename)
                while extension in ZIP_EXT:
                    base, extension = os.path.splitext(base)
                for size, count_int in zip(self.sizes, self.counts):
                    if count_int == self.conf['Downsample']['target']:
                        downsampled[count_int][key] = filename
                        continue
                    target_file = os.path.basename(base) + '_' + self.tool + '_' + humanize(count_int) + extension
                    target_path = '/'.join([bucket_dir, output_path, target_file])
                    downsampled[count_int][key] = target_path
                    if self.tool in ['picard', 'samtools', 'seqtk']:
                        tsv_writer.writerow([size, filename, target_path])
                    elif self.tool == 'zless':
                        tsv_writer.writerow([size * 4, filename, target_path])
                    else:
                        sys.exit()

        scheduler = Scheduler('dsub', self.conf)
        if self.tool == 'picard':
            scheduler.add_argument('--image', 'xingziye/seqdownsample:latest')
            scheduler.add_argument('--command', "'java -jar /app/picard.jar DownsampleSam I=${INPUT_FILE} O=${OUTPUT_FILE} STRATEGY=Chained P=${COUNT} ACCURACY=0.0001'")
        elif self.tool == 'samtools':
            scheduler.add_argument('--image', 'xingziye/seqdownsample:latest')
            scheduler.add_argument('--command', "'samtools view -bs ${COUNT} ${INPUT_FILE} > ${OUTPUT_FILE}'")
        elif self.tool == 'seqtk':
            scheduler.add_argument('--image', 'xingziye/seqtk:latest')
            scheduler.add_argument('--command', "'seqtk sample ${INPUT_FILE} ${COUNT} > ${OUTPUT_FILE}'")
        elif self.tool == 'zless':
            scheduler.add_argument('--image', 'xingziye/seqtk:latest')
            scheduler.add_argument('--command', "'zless ${INPUT_FILE} | head -n ${COUNT} > ${OUTPUT_FILE}'")
        scheduler.add_argument('--tasks', 'downsample.tsv')
        scheduler.add_argument('--logging', log_path)
        scheduler.add_argument('--machine-type', Downsample.default_instance)
        scheduler.add_argument('--disk-size', Downsample.default_disksize)
        # scheduler.add_argument('--boot-disk-size', '50')
        scheduler.add_argument('--wait')
        scheduler.add_argument('--skip')
        p = scheduler.run()
        p.wait()
        return downsampled
