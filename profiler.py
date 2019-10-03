import os
import csv
import time
import tempfile
import shutil
from collections import defaultdict
from google.cloud import storage
from instance import *
from scheduler import Scheduler
from hummingbird_utils import *
try:
    from urllib import unquote # python2
except ImportError:
    from urllib.parse import unquote # python3

PROFILING = 'Profiling'
MACHINE_TYPE_PREFIX = 'n1-highmem-'
DEFAULT_THREAD = 8

class Profiler(object):
    """Profile and collect the results for subsampled input data.

    Attributes:
        tool: A string describing downsample methods.
        conf: A dictionary of configrations.
        mode: Option for memory of runtime profiling.
        client: Cloud storage client instance.
    """
    default_boot_disk_size = '50'
    default_disk_size = '1000'
    mem_mode = 'mem'
    time_mode = 'time'

    def __init__(self, backend, tool, mode, conf):
        if backend == 'cromwell':
            self.profiler = CromwellProfiler(mode, conf)
        else:
            self.profiler = BashProfiler(mode, conf)
        self.tool = tool
        self.mode = mode
        self.conf = conf
        self.client = storage.Client(project=conf['Platform']['project'])

    def profile(self, input_dict, machines=None):
        """Perform profiling on given input and collect result from cloud
            storage bucket.

        Args:
            input_dict: A dict mapping the number of downsampled size to a list of
                downsampled filenames.

        Returns:
            A dict mapping the number of downsampled size to a list of
                memory/time usage sizes in the same order of threads specified in
                config file.
        """
        if machines is None: # do nothing with empty list
            machines = list()
            thread_list = self.conf[PROFILING].get('thread', [DEFAULT_THREAD])
            for thread in thread_list:
                #thread = thread.strip()
                machines.append(GCP_Instance(MACHINE_TYPE_PREFIX + str(thread), thread, None))
        if self.tool == 'time':
            result_dict = self.profiler.profile(input_dict, machines)

        bucket = self.client.get_bucket(self.conf['Platform']['bucket'])
        profiling_dict = dict()
        empty_list = []
        for entry_count in result_dict:
            dir_list = result_dict[entry_count]
            for i, dir_prefix in enumerate(dir_list):
                blobs = list(bucket.list_blobs(prefix=dir_prefix))
                if not blobs:
                    empty_list.append((entry_count, i))
                for blob in blobs:
                    res_set = float(blob.download_as_string())
                    basename = os.path.basename(unquote(blob.path))
                    taskname, _ = os.path.splitext(basename)
                    if taskname not in profiling_dict:
                        profiling_dict[taskname] = defaultdict(list)
                    profiling_dict[taskname][entry_count].append(res_set)
        for taskname in profiling_dict:
            for entry_count, i in empty_list:
                profiling_dict[taskname][entry_count].insert(i, None)
        return profiling_dict

class BaseProfiler(object):
    def __init__(self, mode, conf):
        self.mode = mode
        self.conf = conf

class CromwellProfiler(BaseProfiler):
    def profile(self, input_dict, machines):
        result_dict = defaultdict(list)
        dsub_script = tempfile.NamedTemporaryFile()
        CMD = '''echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
apt-get -qq update
apt-get -qq install time wget
wget -nv -c https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
echo '{"final_workflow_outputs_dir": "/mnt/data/final_outputs"}' > options.json
'''
        if 'imports' in self.conf[PROFILING]:
            CMD += '''unzip ${IMPORT_ZIP}'''
        CMD += '''
java -Dconfig.file=${BACKEND_CONF} -jar cromwell-34.jar run ${WDL_FILE} --inputs ${INPUT_JSON} -m meta.json -o options.json
mkdir /mnt/data/call_logs
for path in $(grep -Po '(?<="callRoot": ").*(?=")' meta.json)
do
call=$(grep -Po '(?<=\/call-)[^\/]*' <<< $path)
grep -Po '%s' $path/execution/stderr.background > /mnt/data/call_logs/$call.txt
done
cp -r -T /mnt/data/call_logs ${RESULT_DIR}
find /mnt/data/final_outputs -type f -execdir cp {} ${OUTPUT_DIR} \;
'''
        if self.mode == Profiler.mem_mode:
            CMD = CMD % '^(\d*\.)?\d+'
        else:
            CMD = CMD % '(\d*\.)?\d+$'
        dsub_script.write(CMD)
        bucket_base = 'gs://' + self.conf['Platform']['bucket'] + '/'
        backend_conf = bucket_base + self.conf[PROFILING]['backend_conf']
        wdl_file = bucket_base + self.conf[PROFILING]['wdl_file']
        log_path = bucket_base + self.conf[PROFILING]['logging']
        procs = []
        temp = tempfile.mkdtemp()
        for machine in machines:
            scheduler = Scheduler('dsub', self.conf)
            scheduler.add_argument('--script', dsub_script.name)
            tsv_filename = os.path.join(temp, 'cromwell_profiling_' + machine.name + '.tsv')
            thread_list = self.conf[PROFILING].get('thread', [DEFAULT_THREAD])
            json_inputs = self.conf[PROFILING]['json_input'][thread_list.index(machine.cpu)]
            with open(tsv_filename, 'w') as dsub_tsv:
                tsv_writer = csv.writer(dsub_tsv, delimiter='\t')
                headline = ['--input BACKEND_CONF',
                            '--input WDL_FILE',
                            '--input INPUT_JSON',
                            '--output-recursive RESULT_DIR',
                            '--output-recursive OUTPUT_DIR']
                if 'imports' in self.conf[PROFILING]:
                    headline += ['--input IMPORT_ZIP']
                    import_zip = bucket_base + self.conf[PROFILING]['imports']
                else:
                    import_zip = None
                tsv_writer.writerow(headline)
                for entry_count, json_input in zip(self.conf['Downsample']['count'], json_inputs):
                    # Enforce the entry_count has the same order as json inputs
                    if entry_count not in input_dict:
                        continue
                    json_input = bucket_base + json_input
                    result_prefix = self.conf[PROFILING]['result'] + '/' + humanize(str(entry_count)) + '/' + machine.name + '/' + self.mode + '/'
                    result_dict[entry_count].append(result_prefix)
                    result_dir = bucket_base + result_prefix
                    output_dir = bucket_base + self.conf[PROFILING]['output'] + '/' + humanize(str(entry_count)) + '/' + machine.name + '/'
                    row = [backend_conf, wdl_file, json_input, result_dir, output_dir]
                    if import_zip:
                        row += [import_zip]
                    tsv_writer.writerow(row)
            scheduler.add_argument('--image', self.conf[PROFILING]['image'])
            scheduler.add_argument('--tasks', tsv_filename)
            scheduler.add_argument('--logging', log_path)
            scheduler.add_argument('--machine-type', machine.name)
            scheduler.add_argument('--boot-disk-size', Profiler.default_boot_disk_size) # google providers put the Docker container's /tmp directory on the boot disk, 10 GB by defualt
            scheduler.add_argument('--disk-size', Profiler.default_disk_size)
            scheduler.add_argument('--wait')
            if not self.conf[PROFILING].get('force', False):
                scheduler.add_argument('--skip')
            dsub_script.seek(0)
            proc = scheduler.run()
            procs.append(proc)
            time.sleep(1)

        try:
            for proc in procs:
                proc.wait()
        finally:
            dsub_script.close()
            shutil.rmtree(temp)
        return result_dict

class BashProfiler(BaseProfiler):
    def profile(self, input_dict, machines):
        """Perform profiling on given input.

        Args:
            input_dict: A dict mapping the number of downsampled size to a list of
                downsampled filenames.

        Returns:
            A dict mapping the number of downsampled size to a list of
                paths of the results of memory/time usage profiling.
        """
        url_base = 'gs://' + self.conf['Platform']['bucket'] + '/'
        log_path = url_base + self.conf[PROFILING]['logging']
        dsub_script = tempfile.NamedTemporaryFile()
        result_dict = defaultdict(list)

        # Prepare profiling script
        dsub_script.write('apt-get -qq update\n')
        dsub_script.write('apt-get -qq install time\n')
        if self.conf[PROFILING].get('script'):
            with open(self.conf[PROFILING]['script'], 'r') as ori:
                for line in ori:
                    dsub_script.write('/usr/bin/time -a -f "%e %M" -o result.txt ' + line + '\n')
        else:
            for line in self.conf[PROFILING]['command'].splitlines():
                dsub_script.write('/usr/bin/time -a -f "%e %M" -o result.txt ' + line + '\n')
        dsub_script.write("awk '{print $1}' result.txt > ${RUNTIME_RESULT}\n")
        dsub_script.write("awk '{print $2}' result.txt > ${RSS_RESULT}\n")

        # Prepare dsub tsv file and run
        # dsub only one type of machine at a time
        # each machine type will have an individual tsv
        procs = []
        temp = tempfile.mkdtemp()
        def add_headline(flag):
            if flag in self.conf[PROFILING]:
                for key in self.conf[PROFILING][flag]:
                    headline.append('--' + flag + ' ' + key)
        for machine in machines:
            scheduler = Scheduler('dsub', self.conf)
            scheduler.add_argument('--script', dsub_script.name)
            tsv_filename = os.path.join(temp, 'profiling_' + machine.name + '.tsv')
            with open(tsv_filename, 'w') as dsub_tsv:
                tsv_writer = csv.writer(dsub_tsv, delimiter='\t')
                headline = ['--env THREAD', '--output RUNTIME_RESULT', '--output RSS_RESULT']
                any_count, any_input_dict = input_dict.popitem() # Pick an arbitrary element to get keys including index and put item back
                input_dict[any_count] = any_input_dict
                for key in any_input_dict.keys():
                    headline.append('--input ' + key)
                add_headline('input') # additional inputs for profiling stage
                add_headline('input-recursive')
                add_headline('output')
                tsv_writer.writerow(headline)

                for entry_count in input_dict:
                    mem_path = self.conf[PROFILING]['result'] + '/' + humanize(str(entry_count)) + '/' + machine.name + '/' + Profiler.mem_mode + '/'
                    time_path = self.conf[PROFILING]['result'] + '/' + humanize(str(entry_count)) + '/' + machine.name + '/' + Profiler.time_mode + '/'
                    result_path = mem_path if self.mode == Profiler.mem_mode else time_path
                    result_dict[entry_count].append(result_path)
                    mem_addr = url_base + mem_path + 'script.txt'
                    time_addr = url_base + time_path + 'script.txt'
                    row = [str(machine.get_core()), time_addr, mem_addr] + input_dict[entry_count].values()
                    if 'input' in self.conf[PROFILING]:
                        for path in self.conf[PROFILING]['input'].values():
                            if path.startswith("gs://"):
                                row.append(path)
                            else:
                                row.append(url_base + path)
                    if 'input-recursive' in self.conf[PROFILING]:
                        for path in self.conf[PROFILING]['input-recursive'].values():
                            if path.startswith("gs://"):
                                row.append(path)
                            else:
                                row.append(url_base + path)
                    if 'output' in self.conf[PROFILING]:
                        for path in self.conf[PROFILING]['output'].values():
                            extension = ""
                            basename, ext = os.path.splitext(path)
                            while ext:
                                extension = ext + extension
                                basename, ext = os.path.splitext(basename)
                            path = basename + '_' + str(machine.get_core()) + '_' + humanize(str(entry_count)) + extension
                            row.append(url_base + path)
                    tsv_writer.writerow(row)

            scheduler.add_argument('--image', self.conf[PROFILING]['image'])
            scheduler.add_argument('--tasks', tsv_filename)
            scheduler.add_argument('--logging', log_path)
            scheduler.add_argument('--machine-type', machine.name)
            scheduler.add_argument('--boot-disk-size', Profiler.default_boot_disk_size) # 10 GB by defualt
            scheduler.add_argument('--disk-size', Profiler.default_disk_size)
            scheduler.add_argument('--wait')
            if not self.conf[PROFILING].get('force', False):
                scheduler.add_argument('--skip')
            dsub_script.seek(0)
            proc = scheduler.run()
            procs.append(proc)
            time.sleep(1)

        try:
            for proc in procs:
                proc.wait()
        finally:
            dsub_script.close()
            shutil.rmtree(temp)
        return result_dict
