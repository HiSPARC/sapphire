import unittest

from mock import patch, sentinel, MagicMock
import tables

from sapphire import esd


class ESDTest(unittest.TestCase):

    @patch.object(esd.os.path, 'split')
    def test_create_table(self, mock_split):
        description = {'event_id': tables.UInt32Col(pos=0),
                       'timestamp': tables.Time32Col(pos=1),
                       'nanoseconds': tables.UInt32Col(pos=2),
                       'ext_timestamp': tables.UInt64Col(pos=3),
                       'pulseheights': tables.Int16Col(pos=4, shape=4),
                       'integrals': tables.Int32Col(pos=5, shape=4),
                       'n1': tables.Float32Col(pos=6),
                       'n2': tables.Float32Col(pos=7),
                       'n3': tables.Float32Col(pos=8),
                       'n4': tables.Float32Col(pos=9),
                       't1': tables.Float32Col(pos=10),
                       't2': tables.Float32Col(pos=11),
                       't3': tables.Float32Col(pos=12),
                       't4': tables.Float32Col(pos=13),
                       't_trigger': tables.Float32Col(pos=14)}
        mock_split.return_value = (sentinel.head, sentinel.tail)
        file = MagicMock()
        result = esd.create_table(file, sentinel.group)
        file.create_group.assert_called_once_with(sentinel.head, sentinel.tail)
        file.create_table.assert_called_once_with(sentinel.group, 'events',
                                                  description)
        self.assertEqual(result, file.create_table.return_value)

    def test_read_line_and_store_event(self):
        comment_input = ('#',)
        input = ('2014-01-27', '00:00:01', '1390780801', '637714403',
                '2', '139', '397', '2', '0', '1085', '5212', '0', '0.0',
                '0.3085', '1.5414', '0.0', '-999', '15.0', '12.5', '-999',
                '17.5')
        output =  [[0, 1390780801, 637714403, 1390780801637714403,
                   [2, 139, 397, 2], [0, 1085, 5212, 0],
                   0.0, 0.3085, 1.5414, 0.0, -999.0, 15.0, 12.5, -999.0, 17.5]]
        table = []
        timestamp = esd.read_line_and_store_event(comment_input, table)
        self.assertEqual(timestamp, 0.)
        self.assertEqual(len(table), 0)
        timestamp = esd.read_line_and_store_event(input, table)
        self.assertEqual(timestamp, int(input[2]))
        self.assertEqual(len(table), 1)
        self.assertEqual(table[0], output)


if __name__ == '__main__':
    unittest.main()
