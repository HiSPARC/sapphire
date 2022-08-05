import argparse

from .esd_load_data import create_and_store_test_data


def main():
    descr = "Download and update test data for usage in tests."
    parser = argparse.ArgumentParser(description=descr)
    parser.parse_args()
    create_and_store_test_data()
