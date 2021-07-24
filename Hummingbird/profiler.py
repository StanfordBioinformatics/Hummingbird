# !/usr/bin/env python3

import copy
import csv
import json
import logging
import os
import shutil
import sys
import tempfile
import time
from collections import defaultdict
from multiprocessing import Pool
from string import Template

from retry import retry

from .hummingbird_utils import PLATFORM, PROFILING, humanize, DOWNSAMPLE, ZIP_EXT
from .instance import GCPInstance, AWSInstance, AzureInstance
from .scheduler import AWSBatchScheduler, AzureBatchScheduler, Scheduler

try:
    from urllib import unquote  # python2
except ImportError:
    from urllib.parse import unquote  # python3

DEFAULT_VCPU = 2  # profiler's thread == vCPU in the current implementation


class Profiler(object):
    """Profile and collect the results for subsampled input data.

    Attributes:
        tool: A string describing downsample methods.
        conf: A dictionary of configrations.
        mode: Option for memory of runtime profiling.
        client: Cloud storage client instance.
    """
    default_boot_disk_size = '50'
    default_disk_size = '500'
    mem_mode = 'mem'
    time_mode = 'time'

    def __init__(self, backend, tool, mode, conf):
        self.backend = backend
        self.tool = tool
        self.mode = mode
        self.conf = conf
        self.service = conf[PLATFORM]['service']
        if self.service == 'aws':
            import boto3
            self.client = boto3.client('s3', region_name=conf[PLATFORM]['regions'])
        elif self.service == 'gcp':
            from google.cloud import storage
            self.client = storage.Client(project=conf[PLATFORM]['project'])
        elif self.service == 'azure':
            import logging
            logging.getLogger('azure').setLevel(logging.WARNING)
            from azure.storage.blob import BlobServiceClient
            self.client = BlobServiceClient.from_connection_string(conf[PLATFORM]['storage_connection_string'])
            self.container_client = self.client.get_container_client(container=conf[PLATFORM]['storage_container'])
        self.output_dict = None

    def _set_profiler(self, conf):
        if self.backend == 'cromwell':
            self.profiler = CromwellProfiler(self.mode, conf)
        else:
            self.profiler = BashProfiler(self.mode, conf)

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
        if self.service == 'aws':
            profiling_dict = dict()
            profiling_dict['script'] = defaultdict(list)
            result_dict = self.aws_batch_profile(input_dict, machines)
            bucket_name = self.conf[PLATFORM]['bucket']
            tries = self.conf[PROFILING].get('tries', 1)
            for entry_count in result_dict:
                dir_list = result_dict[entry_count]
                for dir_prefix in dir_list:
                    total = 0.0
                    for t in range(tries):
                        obj = self.client.get_object(Bucket=bucket_name, Key=dir_prefix + 'try' + str(t) + '.txt')
                        value = json.loads(obj['Body'].read())
                        print(t, value)
                        total += value
                    profiling_dict['script'][entry_count].append(total / tries)
            return profiling_dict
        elif self.service == 'azure':
            profiling_dict = dict()
            profiling_dict['script'] = defaultdict(list)
            result_dict, jobs_info = self.az_batch_profile(input_dict, machines)
            container_name = self.conf[PLATFORM]['storage_container']
            if jobs_info:
                for entry_count in result_dict:
                    dir_list = result_dict[entry_count]
                    for dir_prefix in dir_list:
                        matching_jobs = [job_info for job_info in jobs_info if job_info['result_path'] == dir_prefix]
                        total = 0.0
                        for job_info in matching_jobs:
                            path = os.path.join(dir_prefix, 'try' + job_info['task_id'] + '.txt')
                            blob_client = self.client.get_blob_client(container=container_name, blob=path)

                            from azure.core.exceptions import ResourceNotFoundError
                            try:
                                value = self.get_azure_blob_value(blob_client)
                                value = json.loads(value)
                                total += value
                            except ResourceNotFoundError:
                                pass

                        profiling_dict['script'][entry_count].append(total / len(matching_jobs))
            return profiling_dict

        if machines is None:  # do nothing with empty list
            machines = list()
            thread_list = self.conf[PROFILING].get('thread', [DEFAULT_VCPU])
            for thread in thread_list:
                machines.append(GCPInstance('n1-highmem-' + str(thread)))
        tries = self.conf[PROFILING].get('tries', 1)
        if tries > 1:
            # Apply parallel processing
            pool = Pool(tries)
            results = []
            for i in range(tries):
                new_conf = copy.deepcopy(self.conf)
                new_conf[PROFILING]['result'] += '/try' + str(i + 1)
                self._set_profiler(new_conf)
                results.append(pool.apply_async(self.profiler, (input_dict, machines)))
            try:
                result_dicts = [r.get() for r in results]  # Collecting results
            except:
                pool.terminate()
                pool.join()
        else:
            self._set_profiler(self.conf)
            result_dicts = [self.profiler.profile(input_dict, machines)]
        if self.mode == Profiler.mem_mode:
            self.output_dict = self.profiler.output_dict

        bucket = self.client.get_bucket(self.conf[PLATFORM]['bucket'])
        profiling_dict = dict()
        for result_dict in result_dicts:
            print(result_dict)
            for entry_count in result_dict:
                dir_list = result_dict[entry_count]
                for i, dir_prefix in enumerate(dir_list):
                    blobs = list(bucket.list_blobs(prefix=dir_prefix))
                    print(blobs)
                    for blob in blobs:  # blobs are named as task.txt
                        res_set = float(blob.download_as_string())
                        basename = os.path.basename(unquote(blob.path))
                        taskname, _ = os.path.splitext(basename)
                        if taskname not in profiling_dict:
                            profiling_dict[taskname] = defaultdict(lambda: [0] * len(dir_list))
                        profiling_dict[taskname][entry_count][i] += res_set / tries
                        print(taskname, entry_count, res_set)
        return profiling_dict

    @retry(tries=5, delay=1, max_delay=5, logger=None)
    def get_azure_blob_value(self, blob_client):
        return blob_client.download_blob().readall()

    def aws_batch_profile(self, input_dict, machines):
        if machines is None:  # do nothing with empty list
            machines = list()
            thread_list = self.conf[PROFILING].get('thread', [DEFAULT_VCPU])
            for thread in thread_list:
                if thread not in AWSInstance.thread_suffix:
                    logging.warning("Skipping thread '%d' for AWS as it is not supported", thread)
                    continue
                machines.append(AWSInstance('r5' + AWSInstance.thread_suffix[thread]))
        tries = self.conf[PROFILING].get('tries', 1)

        url_base = 's3://' + self.conf[PLATFORM]['bucket'] + '/'
        result_dict = defaultdict(list)
        output_dict = defaultdict(dict)
        prev_mem = 0  # used to get output from the largest instances
        jobs = []

        for machine in machines:
            for entry_count in input_dict:
                job_script = tempfile.NamedTemporaryFile(mode='w')
                job_script.write('apt-get -qq update && apt-get -qq install time\n')
                local_name_dict = {'THREAD': machine.get_core()}
                for input_key, input_s3_path in input_dict[entry_count].items():
                    local_name = os.path.basename(input_s3_path)
                    job_script.write(' '.join(['aws', 's3', 'cp', input_s3_path, local_name]) + '\n')
                    local_name_dict[input_key] = local_name
                if 'input' in self.conf[PROFILING]:
                    for key, path in self.conf[PROFILING]['input'].items():
                        local_name_dict[key] = path
                        job_script.write(' '.join(['aws', 's3', 'cp', url_base + path, path]) + '\n')
                if 'input-recursive' in self.conf[PROFILING]:
                    for key, path in self.conf[PROFILING]['input-recursive'].items():
                        local_name_dict[key] = path
                        job_script.write(' '.join(['aws', 's3', 'cp', url_base + path, path, '--recursive']) + '\n')
                if 'output' in self.conf[PROFILING]:
                    for key in self.conf[PROFILING]['output'].keys():
                        path = self.conf[PROFILING]['output'][key]
                        local_name = os.path.basename(path)
                        local_name_dict[key] = local_name
                if 'output-recursive' in self.conf[PROFILING]:
                    for key in self.conf[PROFILING]['output-recursive'].keys():
                        path = self.conf[PROFILING]['output-recursive'][key]
                        local_name = os.path.basename(path)
                        local_name_dict[key] = local_name
                if self.conf[PROFILING].get('script'):
                    with open(self.conf[PROFILING]['script'], 'r') as ori:
                        for line in ori:
                            sub_line = Template(line).safe_substitute(local_name_dict)
                            job_script.write('/usr/bin/time -a -f "%e %M" -o result.txt ' + sub_line + '\n')
                else:
                    for line in self.conf[PROFILING]['command'].splitlines():
                        sub_line = Template(line).safe_substitute(local_name_dict)

                        job_script.write('/usr/bin/time -a -f "%e %M" -o result.txt ' + sub_line + '\n')
                job_script.write("awk '{print $1}' result.txt > time_result.txt\n")
                job_script.write("awk '{print $2}' result.txt > mem_result.txt\n")

                mem_path = self.conf[PROFILING]['result'] + '/' + humanize(
                    str(entry_count)) + '/' + machine.name + '/' + Profiler.mem_mode + '/'
                time_path = self.conf[PROFILING]['result'] + '/' + humanize(
                    str(entry_count)) + '/' + machine.name + '/' + Profiler.time_mode + '/'
                result_path = mem_path if self.mode == Profiler.mem_mode else time_path
                result_dict[entry_count].append(result_path)
                mem_addr = url_base + mem_path + 'try$AWS_BATCH_JOB_ARRAY_INDEX.txt'
                time_addr = url_base + time_path + 'try$AWS_BATCH_JOB_ARRAY_INDEX.txt'
                job_script.write(' '.join(['aws', 's3', 'cp', 'time_result.txt', time_addr]) + '\n')
                job_script.write(' '.join(['aws', 's3', 'cp', 'mem_result.txt', mem_addr]) + '\n')

                if 'output' in self.conf[PROFILING]:
                    for key in self.conf[PROFILING]['output'].keys():
                        path = self.conf[PROFILING]['output'][key]
                        basename, ext = os.path.splitext(path)
                        path = basename + '_' + str(machine.get_core()) + '_' + humanize(str(entry_count)) + ext
                        job_script.write(' '.join(['aws', 's3', 'cp', local_name_dict[key], url_base + path]) + '\n')
                        if machine.mem > prev_mem:
                            output_dict[entry_count][key] = url_base + path
                if 'output-recursive' in self.conf[PROFILING]:
                    for key in self.conf[PROFILING]['output-recursive'].keys():
                        path = self.conf[PROFILING]['output-recursive'][key]
                        path += '_' + str(machine.get_core()) + '_' + humanize(str(entry_count))
                        job_script.write(
                            ' '.join(['aws', 's3', 'cp', local_name_dict[key], url_base + path, '--recursive']) + '\n')
                        if machine.mem > prev_mem:
                            output_dict[entry_count][key] = url_base + path

                response = self.client.list_objects_v2(
                    Bucket=self.conf[PLATFORM]['bucket'],
                    Prefix=result_path + 'try')
                if len(response.get('Contents', [])) >= tries and not self.conf[PROFILING].get('force', False):
                    job_script.close()
                    logging.info('Result files exist. Job skipped.')
                else:
                    job_script.seek(0)
                    disk_size = self.conf[PROFILING].get('disk', Profiler.default_disk_size)
                    image = self.conf[PROFILING].get('image')
                    scheduler = AWSBatchScheduler(self.conf, machine, disk_size, job_script.name, image=image)
                    job = scheduler.submit_job(tries)
                    jobs.append(job)
                    time.sleep(1)
            prev_mem = machine.mem
        if jobs:
            AWSBatchScheduler(self.conf, None, None, None).wait_jobs(jobs)
        self.output_dict = output_dict
        return result_dict

    def az_batch_profile(self, input_dict, machines):
        if machines is None:  # do nothing with empty list
            machines = list()
            thread_list = self.conf[PROFILING].get('thread', [DEFAULT_VCPU])
            for thread in thread_list:
                machines.append(AzureInstance(self.conf, name=AzureInstance.machine_thread_mapping[int(thread)][-1]))
        tries = self.conf[PROFILING].get('tries', 1)

        result_dict = defaultdict(list)
        output_dict = defaultdict(dict)
        prev_mem = 0  # used to get output from the largest instances
        jobs = []
        metadata = 'pool=$AZ_BATCH_POOL_ID job=$AZ_BATCH_JOB_ID task=$AZ_BATCH_TASK_ID'

        storage_account = self.conf[PLATFORM]['storage_account']
        storage_container = self.conf[PLATFORM]['storage_container']
        for machine in machines:
            for entry_count in input_dict:
                job_script = tempfile.NamedTemporaryFile(mode='w')
                job_script.write('apt-get -qq update && apt-get -qqy install time\n')
                local_name_dict = {'THREAD': machine.get_core()}
                for input_key, key_space in input_dict[entry_count].items():
                    local_name = os.path.basename(key_space)
                    job_script.write(' '.join(['az', 'storage', 'blob', 'download',
                                               '--container-name', storage_container, '--name', key_space,
                                               '--file', local_name, '--account-name', storage_account]) + '\n')
                    local_name_dict[input_key] = local_name
                if 'input' in self.conf[PROFILING]:
                    for key, path in self.conf[PROFILING]['input'].items():
                        local_name_dict[key] = path
                        job_script.write(' '.join(['az', 'storage', 'blob', 'download',
                                                   '--container-name', storage_container, '--name', path, '--file',
                                                   path,
                                                   '--account-name', storage_account]) + '\n')
                if 'input-recursive' in self.conf[PROFILING]:
                    for key, path in self.conf[PROFILING]['input-recursive'].items():
                        local_name_dict[key] = path
                        job_script.write(' '.join(['mkdir', '-p', path]) + '\n')
                        job_script.write(' '.join(['az', 'storage', 'blob', 'download-batch',
                                                   '--source', storage_container, '--destination', '.',
                                                   '--pattern', f'{path}/*',
                                                   '--account-name', storage_account]) + '\n')
                if 'output' in self.conf[PROFILING]:
                    for key in self.conf[PROFILING]['output'].keys():
                        path = self.conf[PROFILING]['output'][key]
                        local_name = os.path.basename(path)
                        local_name_dict[key] = local_name
                if 'output-recursive' in self.conf[PROFILING]:
                    for key in self.conf[PROFILING]['output-recursive'].keys():
                        path = self.conf[PROFILING]['output-recursive'][key]
                        local_name = os.path.basename(path)
                        local_name_dict[key] = local_name
                if self.conf[PROFILING].get('script'):
                    with open(self.conf[PROFILING]['script'], 'r') as ori:
                        for line in ori:
                            sub_line = Template(line).safe_substitute(local_name_dict)
                            job_script.write('/usr/bin/time -a -f "%e %M" -o result.txt ' + sub_line + '\n')
                else:
                    for line in self.conf[PROFILING]['command'].splitlines():
                        sub_line = Template(line).safe_substitute(local_name_dict)
                        job_script.write('/usr/bin/time -a -f "%e %M" -o result.txt ' + sub_line + '\n')
                job_script.write("awk '{print $1}' result.txt > time_result.txt\n")
                job_script.write("awk '{print $2}' result.txt > mem_result.txt\n")

                mem_path = self.conf[PROFILING]['result'] + '/' + humanize(
                    str(entry_count)) + '/' + machine.name + '/' + Profiler.mem_mode + '/'
                time_path = self.conf[PROFILING]['result'] + '/' + humanize(
                    str(entry_count)) + '/' + machine.name + '/' + Profiler.time_mode + '/'
                result_path = mem_path if self.mode == Profiler.mem_mode else time_path
                result_dict[entry_count].append(result_path)
                mem_addr = mem_path + 'try$AZ_BATCH_TASK_ID.txt'
                time_addr = time_path + 'try$AZ_BATCH_TASK_ID.txt'
                job_script.write(' '.join([
                    'az', 'storage', 'blob', 'upload', '--container-name', storage_container, '--name', time_addr,
                    '--file', 'time_result.txt', '--account-name', storage_account, '--metadata', metadata]) + '\n')
                job_script.write(' '.join([
                    'az', 'storage', 'blob', 'upload', '--container-name', storage_container, '--name', mem_addr,
                    '--file', 'mem_result.txt', '--account-name', storage_account, '--metadata', metadata]) + '\n')

                if 'output' in self.conf[PROFILING]:
                    for key in self.conf[PROFILING]['output'].keys():
                        path = self.conf[PROFILING]['output'][key]
                        basename, ext = os.path.splitext(path)
                        path = basename + '_' + str(machine.get_core()) + '_' + humanize(str(entry_count)) + ext
                        job_script.write(' '.join([
                            'az', 'storage', 'blob', 'upload', '--container-name', storage_container,
                            '--file', local_name_dict[key], '--name', path, '--account-name', storage_account,
                            '--metadata', metadata]) + '\n')
                        if machine.mem > prev_mem:
                            output_dict[entry_count][key] = path
                if 'output-recursive' in self.conf[PROFILING]:
                    for key in self.conf[PROFILING]['output-recursive'].keys():
                        path = self.conf[PROFILING]['output-recursive'][key]
                        path += '_' + str(machine.get_core()) + '_' + humanize(str(entry_count))
                        job_script.write(' '.join([
                            'az', 'storage', 'blob', 'upload-batch', '--destination', storage_container,
                            '--source', local_name_dict[key], '--destination-path', path, '--account-name',
                            storage_account,
                            '--metadata', metadata]) + '\n')
                        if machine.mem > prev_mem:
                            output_dict[entry_count][key] = path

                blobs = list(self.container_client.list_blobs(name_starts_with=result_path + 'try', include='metadata'))
                if blobs and len(blobs) >= tries and not self.conf[PROFILING].get('force', False):
                    job_script.close()
                    for blob in blobs:
                        jobs.append({'pool_id': blob.metadata['pool'], 'job_id': blob.metadata['job'],
                                     'task_id': blob.metadata['task'], 'result_path': result_path})
                    print("Result files exist. Job skipped.")
                else:
                    job_script.seek(0)
                    disk_size = self.conf[PROFILING].get('disk', Profiler.default_disk_size)
                    image = self.conf[PROFILING].get('image')
                    scheduler = AzureBatchScheduler(self.conf, machine, disk_size, job_script.name, image=image)
                    job_info = scheduler.submit_job(tries=tries)
                    job_info['result_path'] = result_path
                    jobs.append(job_info)
                    time.sleep(1)
            prev_mem = machine.mem
        if jobs:
            scheduler = AzureBatchScheduler(self.conf, None, None, None)
            scheduler.wait_for_tasks_to_complete([info['job_id'] for info in jobs])

        self.output_dict = output_dict
        return result_dict, jobs


class BaseProfiler(object):
    def __init__(self, mode, conf):
        self.mode = mode
        self.conf = conf
        self.output_dict = None


class CromwellProfiler(BaseProfiler):
    def profile(self, input_dict, machines):
        result_dict = defaultdict(list)
        dsub_script = tempfile.NamedTemporaryFile(mode='w')
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
mkdir /mnt/data/mem_logs
mkdir /mnt/data/time_logs
for path in $(grep -Po '(?<="callRoot": ").*(?=")' meta.json)
do
call=$(grep -Po '(?<=\/call-)[^\/]*' <<< $path)
grep -Po '^(\d*\.)?\d+' $path/execution/stderr.background > /mnt/data/mem_logs/$call.txt
grep -Po '(\d*\.)?\d+$' $path/execution/stderr.background > /mnt/data/time_logs/$call.txt
done
cp -r -T /mnt/data/mem_logs ${MEM_RESULT}
cp -r -T /mnt/data/time_logs ${TIME_RESULT}
find /mnt/data/final_outputs -type f -execdir cp {} ${OUTPUT_DIR} \;
'''
        # if self.mode == Profiler.mem_mode:
        #     CMD = CMD % '^(\d*\.)?\d+'
        # else:
        #     CMD = CMD % '(\d*\.)?\d+$'
        dsub_script.write(CMD)
        bucket_base = 'gs://' + self.conf[PLATFORM]['bucket'] + '/'
        backend_conf = bucket_base + self.conf[PROFILING]['backend_conf']
        wdl_file = bucket_base + self.conf[PROFILING]['wdl_file']
        log_path = bucket_base + self.conf[PROFILING]['logging']
        procs = []
        temp = tempfile.mkdtemp()
        for machine in machines:
            scheduler = Scheduler('dsub', self.conf)
            scheduler.add_argument('--script', dsub_script.name)
            tsv_filename = os.path.join(temp, 'cromwell_profiling_' + machine.name + '.tsv')
            thread_list = self.conf[PROFILING].get('thread', [DEFAULT_VCPU])
            json_inputs = self.conf[PROFILING]['json_input'][thread_list.index(machine.cpu)]
            with open(tsv_filename, 'w') as dsub_tsv:
                tsv_writer = csv.writer(dsub_tsv, delimiter='\t')
                headline = ['--input BACKEND_CONF',
                            '--input WDL_FILE',
                            '--input INPUT_JSON',
                            '--output-recursive MEM_RESULT',
                            '--output-recursive TIME_RESULT',
                            '--output-recursive OUTPUT_DIR']
                if 'imports' in self.conf[PROFILING]:
                    headline += ['--input IMPORT_ZIP']
                    import_zip = bucket_base + self.conf[PROFILING]['imports']
                else:
                    import_zip = None
                tsv_writer.writerow(headline)
                for entry_count, json_input in zip(self.conf[DOWNSAMPLE]['count'], json_inputs):
                    # Enforce the entry_count has the same order as json inputs
                    # if entry_count not in input_dict:
                    #     continue
                    json_input = bucket_base + json_input
                    mem_path = self.conf[PROFILING]['result'] + '/' + humanize(
                        str(entry_count)) + '/' + machine.name + '/' + Profiler.mem_mode + '/'
                    time_path = self.conf[PROFILING]['result'] + '/' + humanize(
                        str(entry_count)) + '/' + machine.name + '/' + Profiler.time_mode + '/'
                    result_path = mem_path if self.mode == Profiler.mem_mode else time_path
                    result_dict[entry_count].append(result_path)
                    mem_result = bucket_base + mem_path
                    time_result = bucket_base + time_path
                    output_dir = bucket_base + self.conf[PROFILING]['output'] + '/' + humanize(
                        str(entry_count)) + '/' + machine.name + '/'
                    row = [backend_conf, wdl_file, json_input, mem_result, time_result, output_dir]
                    if import_zip:
                        row += [import_zip]
                    tsv_writer.writerow(row)
            scheduler.add_argument('--image', self.conf[PROFILING]['image'])
            scheduler.add_argument('--tasks', tsv_filename)
            scheduler.add_argument('--logging', log_path)
            scheduler.add_argument('--machine-type', machine.name)
            scheduler.add_argument('--boot-disk-size',
                                   Profiler.default_boot_disk_size)  # google providers put the Docker container's /tmp directory on the boot disk, 10 GB by defualt
            scheduler.add_argument('--disk-size', self.conf[PROFILING].get('disk', Profiler.default_disk_size))
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
    # A hook to allow pickling object for multiprocessing
    def __call__(self, input_dict, machines):
        return self.profile(input_dict, machines)

    def profile(self, input_dict, machines):
        """Perform profiling on given input.

        Args:
            input_dict: A dict mapping the number of downsampled size to a list of
                downsampled filenames.

        Returns:
            A dict mapping the number of downsampled size to a list of
                paths of the results of memory/time usage profiling.
        """
        url_base = 'gs://' + self.conf[PLATFORM]['bucket'] + '/'
        log_path = url_base + self.conf[PROFILING]['logging']
        dsub_script = tempfile.NamedTemporaryFile(mode='w')  # 'w' mode for python3 csv.writer
        result_dict = defaultdict(list)
        output_dict = defaultdict(dict)
        prev_mem = 0  # used to get output from the largest instances

        # Prepare profiling script
        # dsub_script.write('echo "deb [check-valid-until=no] http://cdn-fastly.deb.debian.org/debian jessie main" > /etc/apt/sources.list.d/jessie.list\n')
        # dsub_script.write('echo "deb [check-valid-until=no] http://archive.debian.org/debian jessie-backports main" > /etc/apt/sources.list.d/jessie-backports.list\n')
        # dsub_script.write('sed -i "/deb http:\/\/deb.debian.org\/debian jessie-updates main/d" /etc/apt/sources.list\n')
        # dsub_script.write('apt-get -o Acquire::Check-Valid-Until=false update\n')
        # dsub_script.write('apt-get -qq install time\n')
        dsub_script.write('apt-get -qq update && apt-get -qq install time\n')
        dsub_script.write('df --output=used /mnt/data | tail -1\n')
        if self.conf[PROFILING].get('script'):
            with open(self.conf[PROFILING]['script'], 'r') as ori:
                for line in ori:
                    dsub_script.write('/usr/bin/time -a -f "%e %M" -o result.txt ' + line + '\n')
        else:
            for line in self.conf[PROFILING]['command'].splitlines():
                dsub_script.write('/usr/bin/time -a -f "%e %M" -o result.txt ' + line + '\n')
        dsub_script.write('df --output=used /mnt/data | tail -1\n')
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
                any_count, any_input_dict = input_dict.popitem()  # Pick an arbitrary element to get keys including index and put item back
                input_dict[any_count] = any_input_dict
                for key in any_input_dict.keys():
                    headline.append('--input ' + key)
                add_headline('input')  # additional inputs for profiling stage
                add_headline('input-recursive')
                add_headline('output')
                add_headline('output-recursive')
                tsv_writer.writerow(headline)

                for entry_count in input_dict:
                    mem_path = self.conf[PROFILING]['result'] + '/' + humanize(
                        str(entry_count)) + '/' + machine.name + '/' + Profiler.mem_mode + '/'
                    time_path = self.conf[PROFILING]['result'] + '/' + humanize(
                        str(entry_count)) + '/' + machine.name + '/' + Profiler.time_mode + '/'
                    result_path = mem_path if self.mode == Profiler.mem_mode else time_path
                    result_dict[entry_count].append(result_path)
                    mem_addr = url_base + mem_path + 'script.txt'
                    time_addr = url_base + time_path + 'script.txt'
                    row = [str(machine.get_core()), time_addr, mem_addr] + list(input_dict[entry_count].values())
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
                        for key in self.conf[PROFILING]['output'].keys():
                            path = self.conf[PROFILING]['output'][key]
                            basename, ext = os.path.splitext(path)
                            while ext in ZIP_EXT:
                                basename, ext = os.path.splitext(basename)
                            path = basename + '_' + str(machine.get_core()) + '_' + humanize(str(entry_count)) + ext
                            row.append(url_base + path)
                            if machine.mem > prev_mem:
                                output_dict[entry_count][key] = url_base + path
                    if 'output-recursive' in self.conf[PROFILING]:
                        for key in self.conf[PROFILING]['output-recursive'].keys():
                            path = self.conf[PROFILING]['output-recursive'][key]
                            path += '_' + str(machine.get_core()) + '_' + humanize(str(entry_count))
                            row.append(url_base + path)
                            if machine.mem > prev_mem:
                                output_dict[entry_count][key] = url_base + path
                    tsv_writer.writerow(row)
                prev_mem = machine.mem

            scheduler.add_argument('--image', self.conf[PROFILING]['image'])
            scheduler.add_argument('--tasks', tsv_filename)
            scheduler.add_argument('--logging', log_path)
            scheduler.add_argument('--machine-type', machine.name)
            scheduler.add_argument('--boot-disk-size', Profiler.default_boot_disk_size)  # 10 GB by defualt
            scheduler.add_argument('--disk-size', self.conf[PROFILING].get('disk', Profiler.default_disk_size))
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
        except:
            sys.exit("Exit. Unfinished dsub jobs will still be runing.")
        finally:
            dsub_script.close()
            shutil.rmtree(temp)
        self.output_dict = output_dict
        return result_dict
