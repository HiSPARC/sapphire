"""Compare HDF5 output from a test to expected result"""

from numpy import array, all
from numpy.testing import assert_array_equal
import tables


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
                test.fail("Node '%s' does not exist in datafile" %
                          expected_node._v_pathname)
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
            test.fail("Node '%s' does not exist in datafile" % actual_node)
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
                     "Tables '%s' do not have the same length." %
                     expected_node._v_pathname)
    for colname in expected_node.colnames:
        test.assertIn(colname, actual_node.colnames,
                      "Tables '%s' do not have the same columns." %
                      expected_node._v_pathname)
        expected_col = expected_node.col(colname)
        actual_col = actual_node.col(colname)
        assert_array_equal(expected_col, actual_col,
                           "Tables '%s' column '%s' do not match." %
                           (expected_node._v_pathname, colname))


def validate_vlarrays(test, expected_node, actual_node):
    """Verify that two VLArrays are identical"""

    test.assertEqual(expected_node.shape, actual_node.shape,
                     "VLArrays '%s' do not have the same shape." %
                     expected_node._v_pathname)
    for i, j in zip(expected_node, actual_node):
        test.assertTrue(all(i == j),
                        "VLArrays '%s' do not match." %
                        expected_node._v_pathname)


def validate_arrays(test, expected_node, actual_node):
    """Verify that two Arrays are identical"""

    test.assertEqual(expected_node.shape, actual_node.shape,
                     "Arrays '%s' do not have the same shape." %
                     expected_node._v_pathname)
    test.assertTrue(all(array(expected_node.read()) ==
                        array(actual_node.read())),
                    "Arrays '%s' do not match." %
                    expected_node._v_pathname)


def validate_attributes(test, expected_node, actual_node):
    """Verify that two nodes have the same user attributes"""

    test.assertEqual(expected_node._v_attrs._v_attrnamesuser,
                     actual_node._v_attrs._v_attrnamesuser)
