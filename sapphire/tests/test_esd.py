import unittest

from sapphire import esd


class ESDTest(unittest.TestCase):

    def test_read_line_and_store_event(self):
        input = ('2014-01-27', '00:00:01', '1390780801', '637714403',
                '2', '139', '397', '2', '0', '1085', '5212', '0', '0.0',
                '0.3085', '1.5414', '0.0', '-999', '15.0', '12.5', '-999')
        output =  [[0, 1390780801, 637714403, 1390780801637714403,
                   [2, 139, 397, 2], [0, 1085, 5212, 0],
                   0.0, 0.3085, 1.5414, 0.0, -999.0, 15.0, 12.5, -999.0]]
        table = []
        timestamp = esd.read_line_and_store_event(input, table)
        self.assertEqual(timestamp, int(input[2]))
        self.assertEqual(len(table), 1)
        self.assertEqual(table[0], output)
