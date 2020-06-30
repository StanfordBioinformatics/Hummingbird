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

from hummingbird_utils import *
from instance import *

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

class BatchScheduler(object):
    job_def_name = 'hummingbird-job'
    compute_env_prefix = 'hummingbird-env-'
    startup_script_template = '''#cloud-boothook
#!/bin/bash
cloud-init-per once docker_options echo 'OPTIONS="$${OPTIONS} --storage-opt dm.basesize=$basesize"' >> /etc/sysconfig/docker
'''

    def __init__(self, conf, machine, disk_size, script):
        self.conf = conf
        self.machine = machine
        self.disk_size = disk_size
        self.script = script
        self.batch_client = boto3.client('batch')

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
        with open('AWS/job-defination.json') as f:
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

if __name__ == "__main__":
    ins = AWS_Instance('r4.xlarge')
    scheduler = BatchScheduler(ins, 75, "s3://")
    scheduler.update_job_queue(scheduler.create_compute_environment())
    scheduler.reg_job_def()
