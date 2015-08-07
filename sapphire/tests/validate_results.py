"""Compare HDF5 output from a test to expected result"""

from numpy import array
import tables


def validate_results(test, expected_path, actual_path):
    """Validate simulation results

    :param test: instance of the TestCase.
    :param expected_path: path to the reference data.
    :param actual_path: path to the output from the test.

    """
    with tables.open_file(expected_path) as expected_file, \
            tables.open_file(actual_path) as actual_file:
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


def validate_tables(test, expected_node, actual_node):
    """Verify that two Tables are identical"""

    if not expected_node.nrows == actual_node.nrows:
        test.fail("Tables '%s' do not have the same length." %
                  expected_node._v_pathname)
    for colname in expected_node.colnames:
        if colname not in actual_node.colnames:
            test.fail("Tables '%s' do not have the same columns." %
                      expected_node._v_pathname)
        expected_col = expected_node.col(colname)
        actual_col = actual_node.col(colname)
        test.assertTrue((expected_col == actual_col).all())


def validate_vlarrays(test, expected_node, actual_node):
    """Verify that two VLArrays are identical"""

    if expected_node.shape == actual_node.shape:
        test.assertTrue((array(expected_node.read()) ==
                         array(actual_node.read())).all())
    else:
        test.fail("VLArrays '%s' do not have the same length." %
                  expected_node._v_pathname)


def validate_arrays(test, expected_node, actual_node):
    """Verify that two Arrays are identical"""

    if expected_node.shape == actual_node.shape:
        test.assertTrue((array(expected_node.read()) ==
                         array(actual_node.read())).all())
    else:
        test.fail("Arrays '%s' do not have the same length." %
                  expected_node._v_pathname)
