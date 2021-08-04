import unittest

from botocore.exceptions import WaiterError
from mock import patch

from Hummingbird.hummingbird_utils import PLATFORM
from Hummingbird.scheduler import AWSBatchScheduler


class TestAWSScheduler(unittest.TestCase):
    conf = {PLATFORM: {'regions': 'us-west-2', 'bucket': 'local-bucket'}}
    jobs = ['some-job-id']
    instance = AWSBatchScheduler(conf, None, None, None)

    def test_instance_fields(self):
        self.assertIsNotNone(self.instance.batch_client, 'batch_client field was not initialized')
        self.assertIsNotNone(self.instance.ec2_client, 'ec2_client field was not initialized')
        self.assertIsNotNone(self.instance.s3_bucket, 's3_bucket field was not initialized')

    @patch('botocore.waiter.create_waiter_with_client')
    def test_wait_jobs(self, create_waiter_with_client_mock):
        self.instance.wait_jobs(self.jobs)

        create_waiter_with_client_mock.return_value.wait.assert_called_once_with(jobs=self.jobs)

    @patch('logging.error')
    @patch('botocore.waiter.create_waiter_with_client')
    def test_wait_jobs(self, create_waiter_with_client_mock, _):
        create_waiter_with_client_mock.return_value.wait.side_effect = WaiterError('', '', '')

        self.assertRaises(WaiterError, self.instance.wait_jobs, self.jobs)

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
