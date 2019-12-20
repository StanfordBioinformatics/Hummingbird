import os
import csv
import sys
import tempfile
from collections import defaultdict
from scheduler import Scheduler
from hummingbird_utils import *

class Downsample(object):
    """Generate subsamples for input data.

    Attributes:
        conf: A dictionary of configrations.
        target: An integer indicating the size of whole input.
        fractions: An array of downsample fraction parameters.
        counts: An array of actual downsampled sample size.
        fullrun: Wheather bypass downsampling and run the whole input.
        index: Wheather to index the inputs.
    """
    runtime_frac = 0.001
    default_instance = 'n1-highmem-8'
    default_disk_size = '500'
    default_sam_tool = 'samtools'
    default_fa_tool = 'seqtk'

    def __init__(self, conf):
        self.conf = conf
        self.target = conf[DOWNSAMPLE]['target']
        self.fractions = conf[DOWNSAMPLE].get('fractions', [0.001, 0.01, 0.1])
        self.fullrun = conf[DOWNSAMPLE].get('fullrun', False)
        self.index = conf[DOWNSAMPLE].get('index', False)

        if Downsample.runtime_frac not in self.fractions:
            self.fractions.append(Downsample.runtime_frac)
        if self.fullrun:
            self.fractions = [1]
        self.counts = [int(self.target * f) for f in self.fractions]
        conf[DOWNSAMPLE]['count'] = self.counts

    def subsample(self):
        """Generate subsamples for input data according to the input format."""
        input_dict = self.conf[DOWNSAMPLE]['input']
        type_dict = defaultdict(dict)
        for key in input_dict:
            path = input_dict[key]

            base, extension = os.path.splitext(path)
            while extension in ZIP_EXT:
                base, extension = os.path.splitext(base)
            ext = extension.lower()
            if ext in FA_EXT:
                type_dict[FA][key] = path
            elif ext in FQ_EXT:
                type_dict[FQ][key] = path
            elif ext in SAM_EXT:
                type_dict[SAM][key] = path
            elif ext in BAM_EXT:
                type_dict[BAM][key] = path
            else:
                sys.exit("Unsupported input format.")

        downsampled = {}
        for type in type_dict:
            downsampled.update(self.downsample_by_type(type_dict[type], type))
        return downsampled

    def downsample_by_type(self, key_file_dict, type):
        """Execute downsampling job by types.

        Each file specified in the key_file_dict will be subsampled by the sizes
            specified in conf.

        Args:
            key_file_dict: A dict mapping key name to the path of the file.
            type: A string tells which method to use to downsample.

        Returns:
            A dict mapping the downsampling size to the input dict with the
                corresponding downsampled files' path.
        """
        # This condition avoids empty dsub tsv input.
        if self.fullrun and not self.index:
            return {self.conf[DOWNSAMPLE]['target']: key_file_dict}

        if type == SAM or type == BAM:
            tool = Downsample.default_sam_tool
        else:
            tool = Downsample.default_fa_tool

        bucket_dir = 'gs://' + self.conf['Platform']['bucket']
        output_path = self.conf[DOWNSAMPLE]['output'].strip('/')
        log_path = bucket_dir + '/' + self.conf[DOWNSAMPLE]['logging']
        downsampled = defaultdict(dict)

        dsub_tsv = tempfile.NamedTemporaryFile()
        tsv_writer = csv.writer(dsub_tsv, delimiter='\t')
        row = ['--env COUNT', '--input INPUT_FILE']
        if not self.fullrun:
            row.append('--output OUTPUT_FILE')
        if self.index:
            row.append('--output OUTPUT_INDEX')
        tsv_writer.writerow(row)

        for key in key_file_dict:
            input_path = key_file_dict[key]
            base, extension = os.path.splitext(input_path)
            while extension in ZIP_EXT:
                base, extension = os.path.splitext(base)
            for frac, count_int in zip(self.fractions, self.counts):
                target_file = os.path.basename(base) + '_' + tool + '_' + humanize(count_int) + extension
                target_path = '/'.join([bucket_dir, output_path, target_file])
                if self.fullrun:
                    target_path = input_path
                downsampled[count_int][key] = target_path
                if tool in ['picard', 'samtools', 'seqtk']:
                    row = [frac, input_path]
                else: # zless
                    row = [frac * 4, input_path]
                if not self.fullrun:
                    row.append(target_path)
                if self.index:
                    if type == BAM:
                        ext = '.bai'
                    else: # fasta and fastq
                        ext = '.fai'
                    idx_key = key + '_IDX'
                    index_path = target_path + ext
                    downsampled[count_int][idx_key] = index_path
                    row.append(index_path)
                tsv_writer.writerow(row)

        scheduler = Scheduler('dsub', self.conf)
        if not self.fullrun:
            if tool == 'picard':
                scheduler.add_argument('--image', 'xingziye/seqdownsample:latest')
                command = "'java -jar /app/picard.jar DownsampleSam I=${INPUT_FILE} O=${OUTPUT_FILE} STRATEGY=Chained P=${COUNT} ACCURACY=0.0001'"
            elif tool == 'samtools':
                scheduler.add_argument('--image', 'xingziye/seqdownsample:latest')
                command = "'samtools view -bs ${COUNT} ${INPUT_FILE} > ${OUTPUT_FILE}'"
            elif tool == 'seqtk':
                scheduler.add_argument('--image', 'xingziye/seqtk:latest')
                command = "'seqtk sample ${INPUT_FILE} ${COUNT} > ${OUTPUT_FILE}'"
            elif tool == 'zless':
                scheduler.add_argument('--image', 'xingziye/seqtk:latest')
                command = "'zless ${INPUT_FILE} | head -n ${COUNT} > ${OUTPUT_FILE}'"
            if self.index:
                command = command.rstrip("'")
                if type == FA:
                    command += "; samtools faidx ${OUTPUT_FILE} -o ${OUTPUT_INDEX}'"
                elif type == FQ:
                    command += "; samtools fqidx ${OUTPUT_FILE} -o ${OUTPUT_INDEX}'"
                else:
                    command += "; samtools index ${OUTPUT_FILE} ${OUTPUT_INDEX}'"
        else:
            scheduler.add_argument('--image', 'xingziye/seqdownsample:latest')
            if type == FA:
                command = "'samtools faidx ${INPUT_FILE} -o ${OUTPUT_INDEX}'"
            elif type == FQ:
                command = "'samtools fqidx ${INPUT_FILE} -o ${OUTPUT_INDEX}'"
            else:
                command = "'samtools index ${INPUT_FILE} ${OUTPUT_INDEX}'"
        scheduler.add_argument('--command', command)
        scheduler.add_argument('--tasks', dsub_tsv.name)
        scheduler.add_argument('--logging', log_path)
        scheduler.add_argument('--machine-type', Downsample.default_instance)
        scheduler.add_argument('--disk-size', Downsample.default_disk_size)
        # scheduler.add_argument('--boot-disk-size', '50')
        scheduler.add_argument('--wait')
        if not self.conf[DOWNSAMPLE].get('force', False):
            scheduler.add_argument('--skip')
        dsub_tsv.seek(0)
        try:
            p = scheduler.run()
            p.wait()
        finally:
            dsub_tsv.close()
        return downsampled
