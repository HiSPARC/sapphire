"""Perform code tests.

This package contains all tests which verify the proper working of all
SAPPHiRE code. These tests can also be used to verify if SAPPHiRE was
installed correctly. Simply call the :func:`run_tests` function.

"""
from unittest import defaultTestLoader, TestSuite, TextTestRunner
import os


def run_tests():
    """Collect and run all SAPPHiRE tests

    :return: `unittest.TextTestResult` object containing the test results.

    """
    test_path = os.path.dirname(__file__)
    package_tests = defaultTestLoader.discover(start_dir=test_path)
    test_suite = TestSuite(tests=package_tests)
    test_result = TextTestRunner().run(test_suite)
    return test_result
