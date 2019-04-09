import os
import csv
from collections import defaultdict
from google.cloud import storage
from instance import *
from scheduler import Scheduler
from hummingbird_utils import *

# COMMAND = '''
# dsub \
# --provider google-v2 \
# --project gbsc-gcp-project-cba \
# --regions us-west1 \
# --image vandhanak/broadcloudgatk36gatk37-cgs:v1 \
# --command 'seqtk sample ${INPUT_FILE} ${COUNT} > ${OUTPUT_FILE}' \
# --tasks dsub.tsv
# '''
MACHINE_TYPE_PREFIX = 'n1-highmem-'

class Profiler(object):
    """Profile and collect the results for subsampled input data.

    Attributes:
        tool: A string describing downsample methods.
        conf: A dictionary of configrations.
        client: Cloud storage client instance.
    """
    def __init__(self, tool, mode, conf):
        self.tool = tool
        self.mode = mode
        self.conf = conf
        self.client = storage.Client(project=conf['Platform']['project'])

    def profile(self, input_dict, machines=None):
        """Perform profiling on given input and collect result from cloud
            storage bucket.

        Args:
            input_dict: A dict mapping a string of downsampled size to a list of
                downsampled filenames.

        Returns:
            A dict mapping a string of downsampled size to a list of
                memory/time usage sizes in the same order of threads specified in
                config file.
        """
        if self.tool == 'time':
            if self.mode == 'mem':
                result_dict = self.gnutime_profile('"%M"', input_dict, machines)
            elif self.mode == 'time':
                result_dict = self.gnutime_profile('"%e"', input_dict, machines)

        bucket = self.client.get_bucket(self.conf['Platform']['bucket'])
        profiling_dict = defaultdict(list)
        for entry_count in result_dict:
            files_list = result_dict[entry_count]
            for result_file in files_list:
                blob = bucket.blob(result_file)
                res_set = float(blob.download_as_string())
                profiling_dict[entry_count].append(res_set)
        return profiling_dict

    def gnutime_profile(self, format, input_dict, machines=None):
        """Perform profiling on given input.

        Args:
            input_dict: A dict mapping a string of downsampled size to a list of
                downsampled filenames.

        Returns:
            A dict mapping a string of downsampled size to a list of
                paths of the results of memory/time usage profiling.
        """
        if not machines:
            machines = list()
            thread_list = self.conf.get('Profiling', 'thread', fallback="4").split(',')
            for thread in thread_list:
                thread = thread.strip()
                machines.append(GCP_Instance(MACHINE_TYPE_PREFIX + thread, thread, None))
        dir_input = self.conf.get('Profiling', 'reference', fallback=None)
        extra_input = self.conf.get('Profiling', 'extra_input', fallback="").splitlines()
        extra_output = self.conf.get('Profiling', 'output', fallback="").splitlines()
        output_base = 'gs://' + self.conf['Platform']['bucket'] + '/'
        log_path = output_base + self.conf['Downsample']['logging']
        script_name = 'profiling.sh'
        result_dict = defaultdict(list)

        # Prepare profiling script
        if self.conf.get('Profiling', 'script', fallback=None):
            if self.conf['Profiling']['script'] == script_name:
                _, ext = os.path.splitext(script_name)
                script_name = os.path.basename(script_name) + '_' + ext
            with open(self.conf['Profiling']['script'], 'r') as ori, open(script_name, 'w') as copy:
                copy.write('apt-get -qq update\n')
                copy.write('apt-get -qq install time\n')
                for line in ori:
                    copy.write('/usr/bin/time -f {} -o ${{RESULT_FILE}} -a '.format(format) + line + '\n')
        else:
            with open(script_name, 'w') as script:
                script.write('apt-get -qq update\n')
                script.write('apt-get -qq install time\n')
                for line in self.conf['Profiling']['command'].splitlines():
                    script.write('/usr/bin/time -f {} -o ${{RESULT_FILE}} -a '.format(format) + line + '\n')

        # Prepare dsub tsv file and run
        # dsub only one type of machine at a time
        # each machine type will have an individual tsv
        procs = []
        for machine in machines:
            #thread = thread.strip()
            scheduler = Scheduler('dsub', self.conf)
            scheduler.add_argument('--script', script_name)
            tsv_filename = 'profiling_' + machine.name + '.tsv'
            with open(tsv_filename, 'w') as dsub_tsv:
                tsv_writer = csv.writer(dsub_tsv, delimiter='\t')
                headline = ['--env THREAD', '--output RESULT_FILE']
                for i in range(1, len(input_dict.values()[0]) + 1):
                    headline.append('--input INPUT_FILE' + str(i))
                if dir_input:
                    headline.append('--input-recursive REF_FILE')
                # Add header for extra inputs other than downsampled input
                for i, extra in enumerate(extra_input, start=1):
                    if extra:
                        headline.append('--input EXTRA_FILE' + str(i))
                for i, output in enumerate(extra_output, start=1):
                    if output:
                        headline.append('--output OUTPUT_FILE' + str(i))
                tsv_writer.writerow(headline)

                for entry_count in input_dict:
                    result_path = self.conf['Profiling']['result'] + '/' + self.mode + '/' + machine.name + '_' + humanize(entry_count) + '.out'
                    result_dict[entry_count].append(result_path)
                    result_addr = output_base + result_path
                    row = [machine.get_core(), result_addr] + input_dict[entry_count]
                    if dir_input:
                        row.append(dir_input)
                    for extra in extra_input:
                        extra = extra.strip()
                        if extra:
                            row.append(extra)
                    for output in extra_output:
                        output = output.strip()
                        if output:
                            output = output_base + output
                            full_ext = ''
                            while True:
                                output, ext = os.path.splitext(output)
                                if ext:
                                    full_ext = ext + full_ext
                                else:
                                    break
                            output = output + '_' + machine.get_core() + "_" + humanize(entry_count) + full_ext
                            row.append(output)
                    tsv_writer.writerow(row)

            scheduler.add_argument('--image', self.conf['Profiling']['image'])
            scheduler.add_argument('--tasks', tsv_filename)
            scheduler.add_argument('--logging', log_path)
            scheduler.add_argument('--machine-type', machine.name)
            scheduler.add_argument('--wait')
            scheduler.add_argument('--skip')
            proc = scheduler.run()
            procs.append(proc)

        for proc in procs:
            proc.wait()
        return result_dict

# COMMAND = '''
# dsub \
# --provider google-v2 \
# --project gbsc-gcp-project-cba \
# --regions us-west1 \
# --image quay.io/encode-dcc/atac-seq-pipeline:v1.1.3 \
# --script cromwell.sh

# '''
