import unittest

from mock import patch, sentinel, MagicMock

from sapphire.analysis import process_traces


class TraceObservablesTests(unittest.TestCase):

    def setUp(self):
        self.traces = [[200] * 400 + [500] + [400] * 20 + [200] * 600]
        self.to = process_traces.TraceObservables(self.traces)

    def test_baselines(self):
        self.assertEqual(self.to.baselines, [200])

    def test_std_dev(self):
        self.assertEqual(self.to.std_dev, [0])

    def test_pulseheights(self):
        self.assertEqual(self.to.pulseheights, [300])

    def test_integrals(self):
        self.assertEqual(self.to.integrals, [300 + 200 * 20])


class MeanFilterTests(unittest.TestCase):

    def setUp(self):
        self.trace = [[200] * 400 + [500] + [400] * 20 + [200] * 600]
        self.traces = self.trace * 2
        self.mf = process_traces.MeanFilter()

    def test_init(self):
        self.mf = process_traces.MeanFilter(use_threshold=True,
                                            threshold=sentinel.threshold)
        self.assertEqual(self.mf.threshold, sentinel.threshold)
        self.assertEqual(self.mf.filter, self.mf.mean_filter_with_threshold)

        self.mf = process_traces.MeanFilter(use_threshold=False)
        self.assertRaises(AttributeError, lambda: self.mf.threshold)
        self.assertEqual(self.mf.filter, self.mf.mean_filter_without_threshold)

    @patch.object(process_traces.MeanFilter, 'filter_trace')
    def test_filter_traces(self, mock_filter_trace):
        mock_filter_trace.return_value = sentinel.filtered_trace
        self.assertEqual(self.mf.filter_traces(self.traces),
                         [sentinel.filtered_trace, sentinel.filtered_trace])

    def test_filter_trace(self):
        pass

    def test_mean_filter_with_threshold(self):
        pass

    def test_mean_filter_without_threshold(self):
        pass

if __name__ == '__main__':
    unittest.main()
