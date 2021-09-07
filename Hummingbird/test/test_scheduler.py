import unittest

from botocore.exceptions import WaiterError, ClientError
from mock import patch, MagicMock, mock_open

from Hummingbird.errors import SchedulerException
from Hummingbird.hummingbird_utils import PLATFORM
from Hummingbird.scheduler import AWSBatchScheduler


class TestAWSScheduler(unittest.TestCase):
    conf = {PLATFORM: {'regions': 'us-west-2', 'bucket': 'local-bucket', 'cloudformation_stack_name': 'test'}}
    jobs = ['some-job-id']
    cf_stack_output = [
        {'OutputKey': 'PrivateSubnet1', 'OutputValue': 'subnet1'},
        {'OutputKey': 'PrivateSubnet2', 'OutputValue': 'subnet2'},
        {'OutputKey': 'BatchEC2SecurityGroup', 'OutputValue': 'sg-test'},
        {'OutputKey': 'ECSInstanceProfileRoleARN', 'OutputValue': 'ecsInstanceRole'},
        {'OutputKey': 'ECSTaskExecutionRoleARN', 'OutputValue': 'taskExecutionRole'},
        {'OutputKey': 'BatchServiceRoleARN', 'OutputValue': 'awsBatchServiceRole'}
    ]
    launch_template = """
{
  "LaunchTemplateName": "hummingbird_launch_template",
  "LaunchTemplateData": {
    "EbsOptimized": true,
    "BlockDeviceMappings": [
      {
        "Ebs": {
          "DeleteOnTermination": true,
          "VolumeType": "gp3",
          "VolumeSize": 100,
          "Encrypted": true
        },
        "DeviceName": "/dev/xvda"
      }
    ]
  }
}
    """

    def setUp(self):
        self.instance = AWSBatchScheduler(self.conf, None, 100, None)

    def test_instance_fields(self):
        instance = AWSBatchScheduler(self.conf, None, None, None)
        self.assertIsNotNone(instance.batch_client, 'batch_client field was not initialized')
        self.assertIsNotNone(instance.ec2_client, 'ec2_client field was not initialized')
        self.assertIsNotNone(instance.s3_bucket, 's3_bucket field was not initialized')

    @patch('botocore.waiter.create_waiter_with_client')
    def test_wait_jobs(self, create_waiter_with_client_mock):
        self.instance.wait_jobs(self.jobs)

        create_waiter_with_client_mock.return_value.wait.assert_called_once_with(jobs=self.jobs)

    @patch('logging.exception')
    @patch('botocore.waiter.create_waiter_with_client')
    def test_wait_jobs(self, create_waiter_with_client_mock, exception_mock):
        create_waiter_with_client_mock.return_value.wait.side_effect = WaiterError('', '', '')

        self.assertRaises(SchedulerException, self.instance.wait_jobs, self.jobs)
        exception_mock.assert_called_once()

    def test_get_compute_environment_waiter(self):
        waiter_id = 'some-waiter-id'

        compute_env_waiter = self.instance.get_compute_environment_waiter(waiter_id)

        self.assertEqual(waiter_id, compute_env_waiter.name)
        self.assertEqual(20, compute_env_waiter.config.max_attempts)
        self.assertEqual(1, compute_env_waiter.config.delay)

    def test_get_compute_job_queue_waiter(self):
        waiter_id = 'some-waiter-id'

        compute_env_waiter = self.instance.get_compute_job_queue_waiter(waiter_id)

        self.assertEqual(waiter_id, compute_env_waiter.name)
        self.assertEqual(20, compute_env_waiter.config.max_attempts)
        self.assertEqual(10, compute_env_waiter.config.delay)

    def test_get_compute_job_waiter(self):
        waiter_id = 'some-waiter-id'

        compute_env_waiter = self.instance.get_compute_job_waiter(waiter_id)

        self.assertEqual(waiter_id, compute_env_waiter.name)
        self.assertEqual(24 * 60 * 2, compute_env_waiter.config.max_attempts)
        self.assertEqual(60, compute_env_waiter.config.delay)

    @patch('boto3.client', return_value=MagicMock())
    @patch('builtins.open', new_callable=mock_open, read_data=launch_template)
    def test_create_or_update_launch_template_create(self, _, client_mock):
        self.instance.ec2_client = client_mock
        client_mock.describe_launch_templates.side_effect = ClientError({}, 'DescribeLaunchTemplate')

        self.instance.create_or_update_launch_template()

        client_mock.create_launch_template.assert_called_once()

    @patch('boto3.client', return_value=MagicMock())
    @patch('builtins.open', new_callable=mock_open, read_data=launch_template)
    def test_create_or_update_launch_template_create_version(self, _, client_mock):
        self.instance.ec2_client = client_mock
        client_mock.describe_launch_templates.return_value = {'LaunchTemplates': [self.launch_template]}

        self.instance.create_or_update_launch_template()

        client_mock.create_launch_template_version.assert_called_once()

    @patch('boto3.client', return_value=MagicMock())
    def test_get_cf_stack_output(self, client_mock):
        self.instance.cf_client = client_mock
        client_mock.describe_stacks.return_value = {'Stacks': [{'StackName': 'test', 'Outputs': self.cf_stack_output}]}

        self.instance.get_cf_stack_output()

        client_mock.describe_stacks.assert_called_once_with(StackName='test')

    @patch('boto3.client', return_value=MagicMock())
    @patch('logging.exception')
    def test_get_cf_stack_output_missing_key(self, _, client_mock):
        self.instance.cf_client = client_mock

        for kv in self.cf_stack_output:
            output = [item for item in self.cf_stack_output if item != kv]
            client_mock.describe_stacks.return_value = {'Stacks': [{'StackName': 'test', 'Outputs': output}]}

            self.assertRaises(SchedulerException, self.instance.get_cf_stack_output)
