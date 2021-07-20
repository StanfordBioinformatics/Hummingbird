from Hummingbird.instance import AWSInstance
import unittest


class TestAWSInstance(unittest.TestCase):
    expected_supported_threads = {1, 2, 4, 8, 16, 32, 64}

    def test_valid_threads_are_supported(self):
        is_supported = AWSInstance.check_threads_supported(self.expected_supported_threads)

        self.assertTrue(is_supported, "Should support threads 1-64")

    def test_invalid_threads_are_not_supported(self):
        not_supported = AWSInstance.check_threads_supported([2048, 4096])

        self.assertFalse(not_supported, "Should not support invalid threads")

    def test_get_supported_threads(self):
        supported_threads = AWSInstance.get_supported_threads()

        self.assertTrue(self.expected_supported_threads.issubset(supported_threads), "Should supported threads 1-64")

    def test_instance_families(self):
        families = AWSInstance.get_instance_families()

        self.assertEqual({'m5', 'r5', 'c5', 'c5a', 'c5d', 'c5ad', 'i3'}, families)

    def test_is_instance_valid(self):
        self.assertTrue(AWSInstance.is_instance_valid('c5.xlarge'))

    def test_is_instance_not_valid(self):
        self.assertFalse(AWSInstance.is_instance_valid('c5.999xlarge'))
