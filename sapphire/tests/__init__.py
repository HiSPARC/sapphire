from unittest import defaultTestLoader, TestSuite, TextTestRunner
import os


def run_tests():
    """Collect and run all SAPPHiRE tests

    :return: `unittest.TextTestResult` object.

    """
    test_path = os.path.dirname(__file__)
    package_tests = defaultTestLoader.discover(start_dir=test_path)
    test_suite = TestSuite(tests=package_tests)
    test_result = TextTestRunner().run(test_suite)
    return test_result
