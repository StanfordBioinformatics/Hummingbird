import sys
from datetime import datetime, timedelta
from typing import List

from retry import retry

try:
    from email import encoders
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMEText import MIMEText
except ImportError:
    from email import encoders
    from email.mime.multipart import MIMEMultipart
    from email.mime.text import MIMEText
from string import Template
import boto3
import json
import subprocess
import time
import os
import logging

from Hummingbird.hummingbird_utils import *
from Hummingbird.instance import *

class Scheduler(object):
    """Dsub scheduler construction and execution."""
    def __init__(self, tool, conf):
        self.tool = tool
        if self.tool == 'dsub':
            self.cmd = 'dsub'
            self.add_argument('--provider', 'google-v2')
            self.add_argument('--project', conf['Platform']['project'])
            if 'zones' in conf['Platform']:
                self.add_argument('--zones', conf['Platform']['zones'])
            else:
                self.add_argument('--regions', conf['Platform']['regions'])

    def add_argument(self, argname, value=None):
        """Add one argument to cmd string."""
        if value:
            self.cmd += ' ' + argname + ' ' + value
        else:
            self.cmd += ' ' + argname

    def run(self):
        """Run cmd as subprocess and return the Popen object."""
        logging.debug(self.cmd)
        return subprocess.Popen(self.cmd, shell=True)


class BaseBatchSchduler(object):
    job_def_name = 'hummingbird-job'
    compute_env_prefix = 'hummingbird-env-'
    startup_script_template = '''#cloud-boothook
#!/bin/bash
cloud-init-per once docker_options echo 'OPTIONS="$${OPTIONS} --storage-opt dm.basesize=$basesize"' >> /etc/sysconfig/docker
'''


class AWSBatchScheduler(BaseBatchSchduler):
    def __init__(self, conf, machine, disk_size, script):
        self.conf = conf
        self.machine = machine
        self.disk_size = disk_size
        self.script = script
        self.batch_client = boto3.client('batch')
        super(AWSBatchScheduler, self).__init__()

    def update_laungch_template(self):
        basesize = str(self.disk_size) + 'G'
        startup_script = Template(self.startup_script_template).substitute(basesize=basesize)
        payload = MIMEText(startup_script, 'cloud-boothook')
        mime = MIMEMultipart()
        mime.attach(payload)
        user_data = MIMEText(mime.as_string(), 'plain')
        encoders.encode_base64(user_data)
        user_data_base64 = user_data.get_payload()

        with open('AWS/launch-template-data.json') as f:
            data = json.load(f)
            data['LaunchTemplateData']['BlockDeviceMappings'][0]['Ebs']['VolumeSize'] = self.disk_size
            data['LaunchTemplateData']['UserData'] = user_data_base64

        #subprocess.call(['aws', 'ec2', 'delete-launch-template', '--launch-template-name', 'hummingbird_disk_launch_template'])
        subprocess.call(['aws', 'ec2', '--region', 'us-west-2', 'create-launch-template-version', '--cli-input-json', json.dumps(data)])

    def create_compute_environment(self):
        env_name = self.compute_env_prefix + self.machine.name.replace('.', '_') + '-' + str(self.disk_size)
        output = subprocess.check_output(['aws', 'batch', 'describe-compute-environments', '--compute-environments', env_name])
        desc_json = json.loads(output)
        if desc_json['computeEnvironments']: # Create if not exist
            return env_name

        with open('AWS/compute_environment.json') as f:
            data = json.load(f)
            data['computeEnvironmentName'] = env_name
            data['computeResources']['instanceTypes'].append(self.machine.name)
        subprocess.call(['aws', 'batch', 'create-compute-environment', '--cli-input-json', json.dumps(data)])

        while True:
            time.sleep(1)
            output = subprocess.check_output(['aws', 'batch', 'describe-compute-environments', '--compute-environments', env_name])
            desc_json = json.loads(output)
            if desc_json['computeEnvironments']:
                break
        return env_name

    def update_job_queue(self, env_name):
        job_queue_name = env_name + '-queue'
        output = subprocess.check_output(['aws', 'batch', 'describe-job-queues', '--job-queues', job_queue_name])
        desc_json = json.loads(output)
        env = {"order": 1, "computeEnvironment": env_name}
        data = {"computeEnvironmentOrder": [env]}
        if desc_json['jobQueues']: # Create if not exist
            data["jobQueue"] = job_queue_name
            subprocess.call(['aws', 'batch', 'update-job-queue', '--cli-input-json', json.dumps(data)])
        else:
            data["jobQueueName"] = job_queue_name
            data['state'] = 'ENABLED'
            data['priority'] = 100
            subprocess.call(['aws', 'batch', 'create-job-queue', '--cli-input-json', json.dumps(data)])
        time.sleep(1)
        return job_queue_name

    def reg_job_def(self):
        # subprocess.call(['aws', 'batch', 'deregister-job-definition', '--job-definition', 'hummingbird-job:1'])
        with open('AWS/job-definition.json') as f:
            data = json.load(f)
            data['containerProperties']['vcpus'] = self.machine.cpu
            data['containerProperties']['memory'] = int(self.machine.mem) * 1024
        subprocess.call(['aws', 'batch', 'register-job-definition', '--cli-input-json', json.dumps(data)])

    def submit_job(self, tries=1):
        #self.update_laungch_template()
        job_queue_name = self.update_job_queue(self.create_compute_environment())

        jobname = os.path.basename(self.script)
        s3_path = 's3://' + self.conf[PLATFORM]['bucket'] + '/script/' + jobname + '.sh'
        subprocess.call(['aws', 's3', 'cp', self.script, s3_path])
        data = dict()
        data['vcpus'] = self.machine.cpu
        data['memory'] = int(self.machine.mem * 1024 * 0.9)
        data['command'] = [jobname + '.sh']
        data['environment'] = [{"name": "BATCH_FILE_TYPE", "value": "script"},
                               {"name": "BATCH_FILE_S3_URL", "value": s3_path}]
        arguments = ['aws', 'batch', 'submit-job', '--job-name', jobname,
                     '--job-queue', job_queue_name,
                     '--job-definition', self.job_def_name,
                     '--container-overrides', json.dumps(data)]
        if tries > 1:
            array_properties = {"size": tries}
            arguments += ['--array-properties', json.dumps(array_properties)]
        output = subprocess.check_output(arguments)
        desc_json = json.loads(output)
        return desc_json['jobId']

    @staticmethod
    def wait_jobs(jobs_list):
        client = boto3.client('batch')
        while True:
            finished = True
            response = client.describe_jobs(jobs=jobs_list)
            for job in response['jobs']:
                if job['status'] != 'SUCCEEDED' and job['status'] != 'FAILED':
                    finished = False
                    break
            if finished:
                break
            time.sleep(60)


class AzureBatchScheduler(BaseBatchSchduler):
    def __init__(self, conf, machine, disk_size, script):
        self.conf = conf
        self.machine = machine
        self.disk_size = disk_size
        self.script = script
        self.script_target_name = os.path.basename(self.script) + '.sh' if script else None
        self.task_definition = self._get_task_definition()
        self.image = self.task_definition['image']
        self.batch_client = self._get_azure_batch_client(conf)
        self.container_client = self._get_azure_container_client(conf)
        super(AzureBatchScheduler, self).__init__()

    @staticmethod
    def _get_azure_batch_client(conf):
        from azure.batch import batch_auth, BatchServiceClient
        creds = batch_auth.SharedKeyCredentials(conf[PLATFORM]['batch_account'], conf[PLATFORM]['batch_key'])
        batch_url = f"https://{conf[PLATFORM]['batch_account']}.{conf[PLATFORM]['location']}.batch.azure.com"
        return BatchServiceClient(creds, batch_url=batch_url)

    @staticmethod
    def _get_azure_container_client(conf):
        from azure.storage.blob import BlobServiceClient
        client = BlobServiceClient.from_connection_string(conf[PLATFORM]['storage_connection_string'])
        container_client = client.get_container_client(container=conf[PLATFORM]['storage_container'])
        return container_client

    @staticmethod
    def _get_task_definition():
        path = os.path.join(os.path.dirname(__file__), 'Azure/task.json')
        with open(path, 'r') as task:
            return json.load(task)[0]

    @retry(tries=3, delay=1)
    def create_pool(self):
        from azure.batch import models as batchmodels

        pool_id = self.compute_env_prefix + self.machine.name + '-' + str(self.disk_size)

        pool = self.get_pool(pool_id)
        if pool is not None:
            return pool_id

        sku_to_use, image_ref_to_use = self.select_latest_verified_vm_image_with_node_agent_sku()

        container_configuration = batchmodels.ContainerConfiguration(container_image_names=[self.image])

        config = batchmodels.VirtualMachineConfiguration(
            image_reference=image_ref_to_use,
            node_agent_sku_id=sku_to_use,
            data_disks=[batchmodels.DataDisk(disk_size_gb=self.disk_size, lun=1)],
            container_configuration=container_configuration,
        )

        pool = batchmodels.PoolAddParameter(
            id=pool_id,
            display_name=pool_id,
            virtual_machine_configuration=config,
            vm_size=self.machine.name,
        )

        if self.conf[PLATFORM].get('low_priority', False):
            pool.target_low_priority_nodes = 1
        else:
            pool.target_dedicated_nodes = 1

        self.batch_client.pool.add(pool)

        while self.get_pool(pool_id) is None:
            time.sleep(1)

        return pool_id

    @retry(tries=3, delay=1)
    def get_pool(self, name: str):
        from azure.batch.models import BatchErrorException
        try:
            pool = self.batch_client.pool.get(name)
            if pool and getattr(pool, 'id') == name:
                return pool
        except BatchErrorException:
            pool = None

        return pool

    @retry(tries=3, delay=1)
    def select_latest_verified_vm_image_with_node_agent_sku(
            self, publisher='microsoft-azure-batch', offer='ubuntu-server-container', sku_starts_with='16-04'):
        """Select the latest verified image that Azure Batch supports given
        a publisher, offer and sku (starts with filter).
        :param batch_client: The batch client to use.
        :type batch_client: `batchserviceclient.BatchServiceClient`
        :param str publisher: vm image publisher
        :param str offer: vm image offer
        :param str sku_starts_with: vm sku starts with filter
        :rtype: tuple
        :return: (node agent sku id to use, vm image ref to use)
        """
        # get verified vm image list and node agent sku ids from service
        from azure.batch import models as batchmodels

        options = batchmodels.AccountListSupportedImagesOptions(filter="verificationType eq 'verified'")
        images = self.batch_client.account.list_supported_images(account_list_supported_images_options=options)

        # pick the latest supported sku
        skus_to_use = [
            (image.node_agent_sku_id, image.image_reference) for image in images
            if image.image_reference.publisher.lower() == publisher.lower() and
               image.image_reference.offer.lower() == offer.lower() and
               image.image_reference.sku.startswith(sku_starts_with)
        ]

        # pick first
        agent_sku_id, image_ref_to_use = skus_to_use[0]
        return agent_sku_id, image_ref_to_use

    @retry(tries=3, delay=1)
    def create_job(self, pool_id: str):
        from azure.batch import models as batchmodels

        job_queue_name = pool_id + '-queue'
        job = batchmodels.JobAddParameter(
            id=job_queue_name,
            display_name=job_queue_name,
            pool_info=batchmodels.PoolInformation(pool_id=pool_id)
        )

        try:
            self.batch_client.job.add(job)
        except batchmodels.BatchErrorException as err:
            if err.error.code != "JobExists":
                raise
            else:
                logging.info("Job {!r} already exists".format(job_queue_name))

        return job

    @retry(tries=3, delay=1)
    def add_task(self, job_id: str, default_max_tries=None):
        """
        Adds a task for each input file in the collection to the specified job.
        :param str job_id: The ID of the job to which to add the tasks.
         created for each input file.
        :param int default_max_tries: Fallback max tries.
        :output task: Azure Batch task
        """
        from azure.batch import models as batchmodels

        if 'id' in self.task_definition:
            task_id = self.task_definition.get('id')
        else:
            task_id = os.path.basename(self.script)
        display_name = self.task_definition.get('displayName', task_id)

        logging.info('Adding {} tasks to job [{}]...'.format(task_id, job_id))

        container_settings = batchmodels.TaskContainerSettings(
            image_name=self.image,
            container_run_options='--rm'
        )

        platform = self.conf[PLATFORM]
        environment_settings = [
            batchmodels.EnvironmentSetting(name='AZURE_SUBSCRIPTION_ID', value=platform['subscription']),
            batchmodels.EnvironmentSetting(name='AZURE_STORAGE_ACCOUNT', value=platform['storage_account']),
            batchmodels.EnvironmentSetting(name='AZURE_STORAGE_CONTAINER', value=platform['storage_container']),
            batchmodels.EnvironmentSetting(name='AZURE_STORAGE_CONNECTION_STRING', value=platform['storage_connection_string']),
            batchmodels.EnvironmentSetting(name='BLOB_NAME', value=self.script_target_name),
        ]

        if 'environmentSettings' in self.task_definition and self.task_definition['environmentSettings'] is not None:
            environment_settings.extend([
                batchmodels.EnvironmentSetting(**setting) for setting in self.task_definition['environmentSettings']
            ])

        constraints = None
        if 'constraints' in self.task_definition and self.task_definition['constraints']:
            constraints = batchmodels.TaskConstraints(
                max_wall_clock_time=self.task_definition['constraints'].get('maxWallClockTime', "P1D"),
                max_task_retry_count=self.task_definition['constraints'].get('maxTaskRetryCount', default_max_tries),
                retention_time=self.task_definition['constraints'].get('retentionTime', "P1D"),
            ),

        user_identity = batchmodels.UserIdentity(
            auto_user=batchmodels.AutoUserSpecification(
                scope=batchmodels.AutoUserScope.pool,
                elevation_level=batchmodels.ElevationLevel.admin
            )
        )

        task = batchmodels.TaskAddParameter(
            id=task_id,
            display_name=display_name,
            command_line=self.task_definition['commandLine'],
            constraints=constraints[0],
            container_settings=container_settings,
            environment_settings=environment_settings,
            user_identity=user_identity,
        )

        for validation in task.validate():
            logging.info(validation)

        self.batch_client.task.add(job_id=job_id, task=task)

        return task

    @retry(tries=10, delay=1, backoff=2, max_delay=10)
    def wait_for_tasks_to_complete(self, job_ids: List[str], timeout=timedelta(hours=24)):
        """
        Returns when all tasks in the specified job reach the Completed state.
        :param str job_ids: The id of the jobs whose tasks should be monitored.
        :param timedelta timeout: The duration to wait for task completion. If all
        tasks in the specified job do not reach Completed state within this time
        period, an exception will be raised.
        """
        from azure.batch import models as batchmodels

        timeout_expiration = datetime.now() + timeout

        print("Monitoring all tasks for 'Completed' state, timeout in {}...".format(timeout), end='')

        while datetime.now() < timeout_expiration:
            completed_jobs = 0
            for job_id in job_ids:
                print('.', end='')
                sys.stdout.flush()
                tasks = self.batch_client.task.list(job_id)

                incomplete_tasks = [task for task in tasks if task.state != batchmodels.TaskState.completed]
                if not incomplete_tasks:
                    completed_jobs += 1

            if len(job_ids) == completed_jobs:
                print()
                return True
            else:
                time.sleep(5)

        print()
        raise RuntimeError("ERROR: Tasks did not reach 'Completed' state within "
                           "timeout period of " + str(timeout))

    @retry(tries=3, delay=1)
    def upload_script(self):
        if not self.script:
            return

        with open(self.script, 'rb') as data:
            self.container_client.upload_blob(
                name=self.script_target_name,
                data=data,
            )

    def submit_job(self, tries=1):
        pool_id = self.create_pool()
        job = self.create_job(pool_id)
        self.upload_script()
        task = self.add_task(job.id, default_max_tries=tries)
        return {
            'pool_id': pool_id,
            'task_id': task.id,
            'job_id': job.id,
        }
