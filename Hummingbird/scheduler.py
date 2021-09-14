#!/usr/bin/env python3

import json
import logging
import os
import subprocess
import sys
import time
from datetime import datetime, timedelta
from typing import List

from retry import retry

from .errors import SchedulerException
from .hummingbird_utils import PLATFORM


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


class AWSBatchScheduler(BaseBatchSchduler):
    def __init__(self, conf, machine, disk_size, script, **kwargs):
        self.conf = conf
        self.machine = machine
        self.disk_size = disk_size
        self.script = script
        self.image = kwargs.get('image')
        self.region = conf[PLATFORM]['regions']
        import boto3
        self.batch_client = boto3.client('batch', region_name=self.region)
        self.ec2_client = boto3.client('ec2', region_name=self.region)
        self.s3_bucket = boto3.resource('s3').Bucket(self.conf[PLATFORM]['bucket'])
        self.cf_client = boto3.client('cloudformation', region_name=self.region)
        self.cf_stack_name = conf[PLATFORM]['cloudformation_stack_name']
        super(AWSBatchScheduler, self).__init__()

    def create_or_update_launch_template(self):
        with open('AWS/launch-template-data.json') as f:
            data = json.load(f)
            data['LaunchTemplateData']['BlockDeviceMappings'][-1]['Ebs']['VolumeSize'] = int(self.disk_size)

        from botocore.exceptions import ClientError
        try:
            response = self.ec2_client.describe_launch_templates(LaunchTemplateNames=[data['LaunchTemplateName']])
        except ClientError:
            response = {}

        if not response.get('LaunchTemplates'):
            logging.info('Creating launch template %s as it does not exist', data['LaunchTemplateName'])
            self.ec2_client.create_launch_template(**data)
        else:
            logging.info('Creating a new version for launch template %s', data['LaunchTemplateName'])
            self.ec2_client.create_launch_template_version(**data)

    def create_or_update_compute_environment(self, cf_output):
        with open('AWS/compute_environment.json') as f:
            data = json.load(f)

            compute_env_name = self.cf_stack_name + '-' + self.machine.name.replace('.', '_') + '-' + str(self.disk_size)
            desc_json = self.batch_client.describe_compute_environments(computeEnvironments=[compute_env_name])
            if desc_json['computeEnvironments']:
                logging.info('Skipping creation of AWS Batch Compute environment %s as it already exists', compute_env_name)
                return compute_env_name

            compute_resources = data['computeResources']
            data['computeEnvironmentName'] = compute_env_name
            compute_resources['instanceTypes'].append(self.machine.name)
            if 'EC2KeyPair' in cf_output and cf_output['EC2KeyPair']:
                compute_resources['ec2KeyPair'] = cf_output

            data['serviceRole'] = cf_output['BatchServiceRoleARN']
            compute_resources['subnets'] = [cf_output['PrivateSubnet1'], cf_output['PrivateSubnet2']]
            compute_resources['securityGroupIds'] = [cf_output['BatchEC2SecurityGroup']]
            compute_resources['instanceRole'] = cf_output['ECSInstanceProfileRoleARN']

        data['tags'] = {'Name': compute_env_name}
        logging.info('Attempting to create AWS Batch Compute environment: %s', compute_env_name)
        self.batch_client.create_compute_environment(**data)

        import botocore.waiter
        try:
            logging.info('Waiting for AWS Batch Compute environment %s to provision...', compute_env_name)
            waiter = self.get_compute_environment_waiter(compute_env_name)
            waiter.wait(computeEnvironments=[compute_env_name])
        except botocore.waiter.WaiterError as e:
            msg = f"There was an error with the AWS Batch Compute Environment: {compute_env_name}"
            logging.exception(msg)
            raise SchedulerException(msg)

        logging.info('Successfully created AWS Batch Compute environment: %s', compute_env_name)
        return compute_env_name

    def get_compute_environment_waiter(self, waiter_id):
        from botocore.waiter import WaiterModel
        model = WaiterModel({
            'version': 2,
            'waiters': {
                waiter_id: {
                    'delay': 1,
                    'operation': 'DescribeComputeEnvironments',
                    'maxAttempts': 20,
                    'acceptors': [
                        {
                            'expected': 'VALID',
                            'matcher': 'pathAll',
                            'state': 'success',
                            'argument': 'computeEnvironments[].status'
                        },
                        {
                            'expected': 'INVALID',
                            'matcher': 'pathAny',
                            'state': 'failure',
                            'argument': 'computeEnvironments[].status'
                        }
                    ]
                }
            }
        })
        from botocore import waiter
        return waiter.create_waiter_with_client(waiter_id, model, self.batch_client)

    def create_or_update_job_queue(self, env_name):
        job_queue_name = env_name + '-queue'
        desc_json = self.batch_client.describe_job_queues(jobQueues=[job_queue_name])
        env = {"order": 1, "computeEnvironment": env_name}
        data = {'computeEnvironmentOrder': [env]}
        if desc_json['jobQueues']:  # Create if not exist
            data["jobQueue"] = job_queue_name
            logging.info('Attempting to update AWS Batch Job Queue: %s', job_queue_name)
            self.batch_client.update_job_queue(**data)
        else:
            data['jobQueueName'] = job_queue_name
            data['state'] = 'ENABLED'
            data['priority'] = 100
            data['tags'] = {'Name': job_queue_name, 'ComputeEnvironment': env_name}
            logging.info('Attempting to create AWS Batch Job Queue: %s', job_queue_name)
            self.batch_client.create_job_queue(**data)

        from botocore.waiter import WaiterError
        try:
            logging.info('Ensuring AWS Batch Job Queue %s is valid...', job_queue_name)
            job_queue_waiter = self.get_compute_job_queue_waiter(job_queue_name)
            job_queue_waiter.wait(jobQueues=[job_queue_name])
            logging.info('AWS Batch Job Queue %s is valid', job_queue_name)
        except WaiterError as e:
            msg = f"There was an error with the AWS Batch Job Queue: {job_queue_name}"
            logging.exception(msg)
            raise SchedulerException(msg)

        return job_queue_name

    def register_job_definition(self, cf_output, compute_env_name, job_queue_name):
        with open('AWS/job-definition.json') as f:
            data = json.load(f)
            data['containerProperties']['vcpus'] = self.machine.cpu
            data['containerProperties']['memory'] = int(self.machine.mem) * 1024
            data['containerProperties']['jobRoleArn'] = cf_output['ECSTaskExecutionRoleARN']
            if self.image:
                data['containerProperties']['image'] = self.image
        job_definition_name = data.get('jobDefinitionName', self.job_def_name)
        data.setdefault('tags', {})
        data['tags'].update({'Name': job_definition_name, 'ComputeEnvironment': compute_env_name, 'JobQueue': job_queue_name})
        self.batch_client.register_job_definition(**data)
        logging.info('Successfully registered AWS Batch Job Definition: %s', job_definition_name)
        return job_definition_name

    def get_cf_stack_output(self):
        logging.info('Attempting to query Cloudformation Stack: %s', self.cf_stack_name)
        response = self.cf_client.describe_stacks(StackName=self.cf_stack_name)
        stacks = response['Stacks']
        if not stacks or 'Outputs' not in stacks[0] or not stacks[0]['Outputs']:
            msg = f"Unable to query Cloudformation Stack {self.cf_stack_name}"
            logging.exception(msg)
            raise SchedulerException(msg)

        cf_output = {}
        for key in ['PrivateSubnet1', 'PrivateSubnet2', 'BatchEC2SecurityGroup', 'ECSInstanceProfileRoleARN', 'ECSTaskExecutionRoleARN', 'BatchServiceRoleARN']:
            for kv in stacks[0]['Outputs']:
                if kv['OutputKey'] == key:
                    cf_output[key] = kv['OutputValue']

            if key not in cf_output:
                msg = f"Cloudformation stack {self.cf_stack_name} is missing required output: {key}"
                logging.exception(msg)
                raise SchedulerException(msg)

        logging.info('Successfully queried Cloudformation Stack: %s', self.cf_stack_name)
        return cf_output

    def submit_job(self, tries=1):
        cf_output = self.get_cf_stack_output()
        self.create_or_update_launch_template()
        compute_env_name = self.create_or_update_compute_environment(cf_output)
        job_queue_name = self.create_or_update_job_queue(compute_env_name)
        job_definition_name = self.register_job_definition(cf_output, compute_env_name, job_queue_name)

        jobname = os.path.basename(self.script)
        s3_path = 'script/' + jobname + '.sh'
        self.s3_bucket.upload_file(self.script, s3_path)
        data = dict()
        data['vcpus'] = self.machine.cpu
        data['memory'] = int(self.machine.mem * 1024 * 0.9)
        data['command'] = [jobname + '.sh']
        data['environment'] = [
            {"name": "BATCH_FILE_TYPE", "value": "script"},
            {"name": "BATCH_FILE_S3_URL", "value": "s3://{}/{}".format(self.conf[PLATFORM]['bucket'], s3_path)}
        ]
        arguments = {
            'jobName': jobname,
            'jobQueue': job_queue_name,
            'jobDefinition': job_definition_name,
            'containerOverrides': data,
            'tags': {
                'Name': jobname,
                'ComputeEnvironment': compute_env_name,
                'JobQueue': job_queue_name,
                'JobDefinitionName': job_definition_name
            },
            'propagateTags': True
        }
        if tries > 1:
            arguments['arrayProperties'] = {'size': tries}

        desc_json = self.batch_client.submit_job(**arguments)
        job_id = desc_json['jobId']

        logging.info('You can observe the job status via AWS Console: '
                     'https://console.aws.amazon.com/batch/home?region=%s#jobs/%s/%s',
                     self.region, 'array-job' if tries > 1 else 'detail', job_id)
        return job_id

    def wait_jobs(self, jobs_list):
        from botocore.waiter import WaiterError
        waiter_id = '_'.join(jobs_list)
        logging.info('Waiting for AWS Batch Jobs %s to finish...', jobs_list)
        try:
            job_waiter = self.get_compute_job_waiter(waiter_id)
            job_waiter.wait(jobs=jobs_list)
        except WaiterError as e:
            msg = f"There was an error with AWS Batch Jobs {jobs_list}"
            logging.exception(msg)
            raise SchedulerException(msg)

        logging.info('AWS Batch Jobs %s have completed', jobs_list)

    def get_compute_job_waiter(self, waiter_id):
        from botocore.waiter import WaiterModel, create_waiter_with_client
        model = WaiterModel({
            'version': 2,
            'waiters': {
                waiter_id: {
                    'delay': 60,
                    'operation': 'DescribeJobs',
                    'maxAttempts': 24 * 60 * 2,  # timeout of 2 days
                    'acceptors': [
                        {
                            'expected': 'SUCCEEDED',
                            'matcher': 'pathAll',
                            'state': 'success',
                            'argument': 'jobs[].status'
                        },
                        {
                            'expected': 'FAILED',
                            'matcher': 'pathAny',
                            'state': 'failure',
                            'argument': 'jobs[].status'
                        }
                    ]
                }
            }
        })
        return create_waiter_with_client(waiter_id, model, self.batch_client)

    def get_compute_job_queue_waiter(self, waiter_id):
        from botocore.waiter import WaiterModel
        model = WaiterModel({
            'version': 2,
            'waiters': {
                waiter_id: {
                    'delay': 10,
                    'operation': 'DescribeJobQueues',
                    'maxAttempts': 20,
                    'acceptors': [
                        {
                            'expected': 'VALID',
                            'matcher': 'pathAll',
                            'state': 'success',
                            'argument': 'jobQueues[].status'
                        },
                        {
                            'expected': 'INVALID',
                            'matcher': 'pathAny',
                            'state': 'failure',
                            'argument': 'jobQueues[].status'
                        }
                    ]
                }
            }
        })
        from botocore import waiter
        return waiter.create_waiter_with_client(waiter_id, model, self.batch_client)


class AzureBatchScheduler(BaseBatchSchduler):
    def __init__(self, conf, machine, disk_size, script, **kwargs):
        self.conf = conf
        self.machine = machine
        self.disk_size = disk_size
        self.script = script
        self.script_target_name = os.path.basename(self.script) + '.sh' if script else None
        self.task_definition = self._get_task_definition()
        self.image = kwargs.get('image', self.task_definition['image'])
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
        skus_to_use = []
        for image in images:
            if image.image_reference.publisher.lower() == publisher.lower() \
                    and image.image_reference.offer.lower() == offer.lower() \
                    and image.image_reference.sku.startswith(sku_starts_with):
                skus_to_use.append((image.node_agent_sku_id, image.image_reference))

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
                raise SchedulerException(f"Unable to create job {job_queue_name}")
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
            batchmodels.EnvironmentSetting(name='AZURE_STORAGE_CONNECTION_STRING',
                                           value=platform['storage_connection_string']),
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
        raise SchedulerException("ERROR: Tasks did not reach 'Completed' state within timeout period of " + str(timeout))

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
