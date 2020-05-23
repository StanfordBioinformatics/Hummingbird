from email import encoders
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from string import Template
import json
import subprocess
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
    job_queue_name = 'hummingbird-job-queue'
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
        env_name = self.compute_env_prefix + self.machine.name.replace('.', '_')
        output = subprocess.check_output(['aws', 'batch', 'describe-compute-environments', '--compute-environments', env_name])
        desc_json = json.loads(output)
        if desc_json['computeEnvironments']: # Create if not exist
            return env_name

        with open('AWS/compute_environment.json') as f:
            data = json.load(f)
            data['computeEnvironmentName'] = env_name
            data['computeResources']['instanceTypes'].append(self.machine.name)
        subprocess.call(['aws', 'batch', 'create-compute-environment', '--cli-input-json', json.dumps(data)])
        return env_name

    def update_job_queue(self, env_name):
        output = subprocess.check_output(['aws', 'batch', 'describe-job-queues', '--job-queues', self.job_queue_name])
        desc_json = json.loads(output)
        env = {"order": 1, "computeEnvironment": env_name}
        data = {"jobQueueName": self.job_queue_name, "computeEnvironmentOrder": [env]}
        if desc_json['jobQueues']: # Create if not exist
            subprocess.call(['aws', 'batch', 'update-job-queue', '--cli-input-json', json.dumps(data)])
        else:
            data['state'] = 'ENABLED'
            data['priority'] = 100
            subprocess.call(['aws', 'batch', 'create-job-queue', '--cli-input-json', json.dumps(data)])

    def reg_job_def(self):
        # subprocess.call(['aws', 'batch', 'deregister-job-definition', '--job-definition', 'hummingbird-job:1'])
        with open('AWS/job-defination.json') as f:
            data = json.load(f)
            data['containerProperties']['vcpus'] = self.machine.cpu
            data['containerProperties']['memory'] = int(self.machine.mem) * 1024
        subprocess.call(['aws', 'batch', 'register-job-definition', '--cli-input-json', json.dumps(data)])

    def submit_job(self):
        #self.update_laungch_template()
        self.update_job_queue(self.create_compute_environment())

        jobname = os.path.basename(self.script)
        s3_path = 's3://' + self.conf[PLATFORM]['bucket'] + '/script/' + jobname + '.sh'
        subprocess.call(['aws', 's3', 'cp', self.script, s3_path])
        data = dict()
        data['vcpus'] = self.machine.cpu
        data['memory'] = int(self.machine.mem) * 1024 - 1024
        data['command'] = [jobname + '.sh']
        data['environment'] = [{"name": "BATCH_FILE_TYPE", "value": "script"},
                               {"name": "BATCH_FILE_S3_URL", "value": s3_path}]
        subprocess.call(['aws', 'batch', 'submit-job', '--job-name', jobname,
                         '--job-queue', self.job_queue_name,
                         '--job-definition', self.job_def_name,
                         '--container-overrides', json.dumps(data)])

if __name__ == "__main__":
    ins = AWS_Instance('r4.xlarge')
    scheduler = BatchScheduler(ins, 75, "s3://")
    scheduler.update_job_queue(scheduler.create_compute_environment())
    scheduler.reg_job_def()
