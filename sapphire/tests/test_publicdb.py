import unittest
from datetime import datetime

from sapphire import publicdb


class DownloadDataTest(unittest.TestCase):

    def test_datetimerange(self):
        combinations = [(datetime(2010, 1, 1, 11),
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
                          (datetime(2010, 1, 5), datetime(2010, 1, 5, 13))]),
                        ]
        for start, stop, result in combinations:
            self.assertEqual(list(publicdb.datetimerange(start, stop)), result)


if __name__ == '__main__':
    unittest.main()
