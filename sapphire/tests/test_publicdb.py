import unittest
from datetime import datetime
import logging
import shutil
import os

from mock import sentinel, Mock, patch
import tables

from sapphire import publicdb
from sapphire.tests.validate_results import validate_results

from esd_load_data import create_tempfile_path


self_path = os.path.dirname(__file__)
test_data_src_path = os.path.join(self_path, 'test_data/publicdb_src.h5')
test_data_path = os.path.join(self_path, 'test_data/publicdb.h5')


class DownloadDataTest(unittest.TestCase):

    def setUp(self):
        logging.disable(logging.CRITICAL)

    def tearDown(self):
        logging.disable(logging.NOTSET)

    @patch.object(publicdb, '_store_data')
    @patch.object(publicdb.urllib, 'urlretrieve')
    @patch.object(publicdb.xmlrpclib, 'ServerProxy')
    def test_download_data(self, mock_server, mock_retrieve, mock_store):
        start = datetime(2010, 1, 1, 11)
        end = datetime(2010, 1, 1, 13)
        file = Mock()
        mock_get_data_url = mock_server.return_value.hisparc.get_data_url
        mock_get_data_url.return_value = sentinel.url
        mock_retrieve.return_value = (sentinel.tmpdata, sentinel.headers)
        publicdb.download_data(file, sentinel.group, sentinel.station_id,
                               start, end, get_blobs=sentinel.blobs)
        mock_get_data_url.assert_called_once_with(sentinel.station_id, start,
                                                  sentinel.blobs)

        mock_get_data_url.side_effect = Exception("No data")
        publicdb.download_data(file, sentinel.group, sentinel.station_id,
                               start, end, get_blobs=sentinel.blobs)

        mock_get_data_url.side_effect = Exception("Unknown error")
        self.assertRaises(Exception, publicdb.download_data, file,
                          sentinel.group, sentinel.station_id, start,
                          end, get_blobs=sentinel.blobs)

    def test__store_data(self):
        # store data removes the source data when completed, so use a temp
        tmp_src_path = create_tempfile_path()
        shutil.copy(test_data_src_path, tmp_src_path)

        output_path = create_tempfile_path()
        start = datetime(2016, 4, 21)
        end = datetime(2016, 4, 21, 0, 1)
        filters = tables.Filters(complevel=1)
        with tables.open_file(output_path, 'w', filters=filters) as datafile:
            publicdb._store_data(datafile, '/s501', tmp_src_path, start, end)
        validate_results(self, test_data_path, output_path)
        os.remove(output_path)

    def test__store_data_no_end(self):
        # store data removes the source data when completed, so use a temp
        tmp_src_path = create_tempfile_path()
        shutil.copy(test_data_src_path, tmp_src_path)

        output_path = create_tempfile_path()
        start = datetime(2016, 4, 21)
        filters = tables.Filters(complevel=1)
        with tables.open_file(output_path, 'w', filters=filters) as datafile:
            publicdb._store_data(datafile, '/station_501', tmp_src_path, start,
                                 None)
        validate_results(self, test_data_src_path, output_path)
        os.remove(output_path)

    def test_datetimerange(self):
        combinations = [
            (datetime(2010, 1, 1, 11),
             datetime(2010, 1, 1, 13),
             [(datetime(2010, 1, 1, 11), datetime(2010, 1, 1, 13))]),
            (datetime(2010, 1, 1, 11),
             datetime(2010, 1, 2),
             [(datetime(2010, 1, 1, 11), None)]),
            (datetime(2010, 1, 1, 11),
             datetime(2010, 1, 2, 13),
             [(datetime(2010, 1, 1, 11), None),
              (datetime(2010, 1, 2), datetime(2010, 1, 2, 13))]),
            (datetime(2010, 1, 1, 11),
             datetime(2010, 1, 5, 13),
             [(datetime(2010, 1, 1, 11), None),
              (datetime(2010, 1, 2), None),
              (datetime(2010, 1, 3), None),
              (datetime(2010, 1, 4), None),
              (datetime(2010, 1, 5), datetime(2010, 1, 5, 13))])]
        for start, stop, result in combinations:
            self.assertEqual(list(publicdb.datetimerange(start, stop)), result)
            self.assertRaises(Exception, next,
                              publicdb.datetimerange(stop, start))

    def test__get_or_create_group(self):
        file = Mock()
        file.get_node.return_value = sentinel.file_group
        group = publicdb._get_or_create_group(file, sentinel.group)
        self.assertEqual(group, sentinel.file_group)

        file = Mock()
        file.get_node.side_effect = tables.NoSuchNodeError('no such node!')
        in_group = '/hisparc/station_501'
        out_group = publicdb._get_or_create_group(file, in_group)
        file.create_group.assert_called_once_with('/hisparc', 'station_501',
                                                  'Data group',
                                                  createparents=True)
        self.assertEqual(file.create_group.return_value, out_group)

    def test__get_or_create_node(self):
        file = Mock()
        src_node = Mock()
        file.get_node.return_value = sentinel.node

        node = publicdb._get_or_create_node(file, sentinel.group, src_node)
        file.get_node.assert_called_once_with(sentinel.group, src_node.name)
        self.assertEqual(node, sentinel.node)

        file.get_node.side_effect = tables.NoSuchNodeError('no such node!')
        # Raise exception because type of Mock src_node is not Table or VLArray
        self.assertRaises(Exception, publicdb._get_or_create_node, file,
                          sentinel.group, src_node)

        src_node = Mock(spec=tables.Table)
        src_node.description = sentinel.description
        node = publicdb._get_or_create_node(file, sentinel.group, src_node)
        file.create_table.assert_called_once_with(
            sentinel.group, src_node.name, src_node.description,
            src_node.title)

        src_node = Mock(spec=tables.VLArray)
        src_node.atom = sentinel.atom
        node = publicdb._get_or_create_node(file, sentinel.group, src_node)
        file.create_vlarray.assert_called_once_with(
            sentinel.group, src_node.name, src_node.atom, src_node.title)


if __name__ == '__main__':
    unittest.main()
