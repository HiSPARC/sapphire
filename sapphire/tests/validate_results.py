"""Compare HDF5 output from a test to expected result"""

import tables

from numpy import all, array
from numpy.testing import assert_array_almost_equal


def validate_results(test, expected_path, actual_path):
    """Validate results by comparing in and output HDF5 files

    :param test: instance of the TestCase.
    :param expected_path: path to the reference data.
    :param actual_path: path to the output from the test.

    """
    with tables.open_file(expected_path, 'r') as expected_file, \
            tables.open_file(actual_path, 'r') as actual_file:
        for expected_node in expected_file.walk_nodes('/', 'Leaf'):
            try:
                actual_node = actual_file.get_node(expected_node._v_pathname)
            except tables.NoSuchNodeError:
                test.fail(f"Node '{expected_node._v_pathname}' does not exist in datafile")

            if type(expected_node) is tables.table.Table:
                validate_tables(test, expected_node, actual_node)
            elif type(expected_node) is tables.vlarray.VLArray:
                validate_vlarrays(test, expected_node, actual_node)
            elif type(expected_node) is tables.array.Array:
                validate_arrays(test, expected_node, actual_node)
            else:
                raise NotImplementedError
            validate_attributes(test, expected_node, actual_node)
        validate_attributes(test, expected_file.root, actual_file.root)


def validate_results_node(test, expected_path, actual_path, expected_node,
                          actual_node):
    """Validate results by comparing two specific nodes

    :param test: instance of the TestCase.
    :param expected_path: path to the reference data.
    :param actual_path: path to the output from the test.
    :param expected_node: path to the reference node.
    :param actual_node: path to the output node from the test.

    """
    with tables.open_file(expected_path, 'r') as expected_file, \
            tables.open_file(actual_path, 'r') as actual_file:
        expected = expected_file.get_node(expected_node)
        try:
            actual = actual_file.get_node(actual_node)
        except tables.NoSuchNodeError:
            test.fail(f"Node '{actual_node}' does not exist in datafile")

        if type(expected) is tables.table.Table:
            validate_tables(test, expected, actual)
        elif type(expected) is tables.vlarray.VLArray:
            validate_vlarrays(test, expected, actual)
        elif type(expected) is tables.array.Array:
            validate_arrays(test, expected, actual)
        else:
            raise NotImplementedError


def validate_tables(test, expected_node, actual_node):
    """Verify that two Tables are identical"""

    test.assertEqual(expected_node.nrows, actual_node.nrows,
                     f"Tables '{expected_node._v_pathname}' do not have the same length.")
    for colname in expected_node.colnames:
        test.assertIn(colname, actual_node.colnames,
                      f"Tables '{expected_node._v_pathname}' do not have the same columns.")
        expected_col = expected_node.col(colname)
        actual_col = actual_node.col(colname)
        assert_array_almost_equal(
            expected_col,
            actual_col,
            err_msg=f"Tables '{expected_node._v_pathname}' column '{colname}' do not match."
        )


def validate_vlarrays(test, expected_node, actual_node):
    """Verify that two VLArrays are identical"""

    test.assertEqual(expected_node.shape, actual_node.shape,
                     f"VLArrays '{expected_node._v_pathname}' do not have the same shape.")
    for expected_array, actual_array in zip(expected_node, actual_node):
        test.assertTrue(all(expected_array == actual_array),
                        f"VLArrays '{expected_node._v_pathname}' do not match.")


def validate_arrays(test, expected_node, actual_node):
    """Verify that two Arrays are identical"""

    test.assertEqual(expected_node.shape, actual_node.shape,
                     f"Arrays '{expected_node._v_pathname}' do not have the same shape.")
    test.assertTrue(all(array(expected_node.read()) == array(actual_node.read())),
                    f"Arrays '{expected_node._v_pathname}' do not match.")


def validate_attributes(test, expected_node, actual_node):
    """Verify that two nodes have the same user attributes"""

    test.assertEqual(expected_node._v_attrs._v_attrnamesuser,
                     actual_node._v_attrs._v_attrnamesuser)
