import os
import csv
from collections import defaultdict
from google.cloud import storage
from scheduler import Scheduler, humanize

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
    def __init__(self, tool, conf):
        self.tool = tool
        self.conf = conf
        self.client = storage.Client(project=conf['Platform']['project'])

    def profile(self, input_dict):
        """Perform profiling on given input and collect result from cloud
            storage bucket.

        Args:
            input_dict: A dict mapping a string of downsampled size to a list of
                downsampled filenames.

        Returns:
            A dict mapping a string of downsampled size to a list of
                memory usage sizes in the same order of threads specified in
                config file.
        """
        if self.tool == 'time':
            output_dict = self.time_profile(input_dict)

        bucket = self.client.get_bucket(self.conf['Program']['bucket'])
        profiling_dict = defaultdict(list)
        for entry_count in output_dict:
            files_list = output_dict[entry_count]
            entry_count_int = int(entry_count)
            for output_file in files_list:
                blob = bucket.blob(output_file)
                max_res_set = int(blob.download_as_string())
                profiling_dict[entry_count_int].append(max_res_set)
        return profiling_dict

    def time_profile(self, input_dict):
        """Perform profiling on given input.

        Args:
            input_dict: A dict mapping a string of downsampled size to a list of
                downsampled filenames.

        Returns:
            A dict mapping a string of downsampled size to a list of
                paths of the results of memory usage profiling.
        """
        thread_list = self.conf['Program']['thread'].split(',')
        output_base = 'gs://' + self.conf['Program']['bucket'] + '/'
        script_name = 'profiling.sh'
        output_dict = defaultdict(list)

        if self.conf.get('Program', 'script', fallback=None):
            if self.conf['Program']['script'] == script_name:
                _, ext = os.path.splitext(script_name)
                script_name = os.path.basename(script_name) + '_' + ext
            with open(self.conf['Program']['script'], 'r') as ori, open(script_name, 'w') as copy:
                copy.write('apt-get -qq update\n')
                copy.write('apt-get -qq install time\n')
                for line in ori:
                    copy.write('/usr/bin/time -f "%M" -o ${OUTPUT_FILE} -a ' + line + '\n')
        else:
            with open(script_name, 'w') as script:
                script.write('apt-get -qq update\n')
                script.write('apt-get -qq install time\n')
                for line in self.conf['Program']['command'].splitlines():
                    script.write('/usr/bin/time -f "%M" -o ${OUTPUT_FILE} -a ' + line + '\n')

        procs = []
        for thread in thread_list:
            thread = thread.strip()
            scheduler = Scheduler('dsub', self.conf)
            scheduler.add_argument('--script', script_name)
            tsv_filename = 'profiling_' + thread + '.tsv'
            with open(tsv_filename, 'w') as dsub_tsv:
                tsv_writer = csv.writer(dsub_tsv, delimiter='\t')
                headline = ['--env THREAD', '--output OUTPUT_FILE', '--input-recursive REF_FILE']
                for i in range(1, len(input_dict.values()[0]) + 1):
                    headline.append('--input INPUT_FILE' + str(i))
                tsv_writer.writerow(headline)

                ref_path = self.conf['Program']['reference']
                for entry_count in input_dict:
                    output_path = self.conf['Program']['output'] + '_' + thread + '_' + humanize(entry_count)
                    output_dict[entry_count].append(output_path)
                    output_addr = output_base + output_path
                    row = [thread, output_addr, ref_path] + input_dict[entry_count]
                    tsv_writer.writerow(row)

            scheduler.add_argument('--image', self.conf['Program']['image'])
            scheduler.add_argument('--tasks', tsv_filename)
            scheduler.add_argument('--logging', self.conf['Program']['logging'])
            scheduler.add_argument('--machine-type', MACHINE_TYPE_PREFIX + thread)
            scheduler.add_argument('--wait')
            scheduler.add_argument('--skip')
            proc = scheduler.run_after()
            procs.append(proc)

        for proc in procs:
            proc.wait()
        return output_dict
