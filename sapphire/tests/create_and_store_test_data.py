import argparse

from .esd_load_data import create_and_store_test_data


def main():
    description = 'Download and update test data for usage in tests.'
    parser = argparse.ArgumentParser(description=description)
    parser.parse_args()
    create_and_store_test_data()
