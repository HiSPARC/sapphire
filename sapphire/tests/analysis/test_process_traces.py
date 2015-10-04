import unittest
from itertools import cycle

from numpy import array
from mock import patch, sentinel, MagicMock

from sapphire.analysis import process_traces


class TraceObservablesTests(unittest.TestCase):

    def setUp(self):
        trace = ([200] * 400 + [500] + [510] + [400] * 10 + [200] * 600 +
                 [400] * 10 + [200])
        self.traces = array([trace, [0] * len(trace)]).T
        self.to = process_traces.TraceObservables(self.traces)

    def test_baselines(self):
        self.assertEqual(self.to.baselines, [200, 0, -1, -1])

    def test_std_dev(self):
        self.assertEqual(self.to.std_dev, [0, 0, -1, -1])

    def test_pulseheights(self):
        self.assertEqual(self.to.pulseheights, [310, 0, -1, -1])

    def test_integrals(self):
        self.assertEqual(self.to.integrals, [300 + 310 + 200 * 20, 0, -1, -1])

    def test_n_peaks(self):
        self.assertEqual(self.to.n_peaks, [2, 0, -1, -1])


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
        mock_filter = MagicMock()
        self.mf.filter = mock_filter
        mock_filter.side_effect = cycle([[sentinel.filtered_even] * 2,
                                         [sentinel.filtered_odd] * 2,
                                         [sentinel.filtered_recombined]])
        trace_segment = [sentinel.trace_even, sentinel.trace_odd]

        filtered_trace = self.mf.filter_trace(trace_segment * 4)

        self.assertEqual(filtered_trace, [sentinel.filtered_recombined])
        mock_filter.assert_any_call([sentinel.trace_even] * 4)
        mock_filter.assert_any_call([sentinel.trace_odd] * 4)
        mock_filter.assert_called_with([sentinel.filtered_even, sentinel.filtered_odd] * 2)

    def test_mean_filter_with_threshold(self):
        # Small deviations in first few elements
        # (199 + 201 + 199 + 201) / 4. = 200
        # |199 - 200| < threshold and |201 - 200| < threshold
        raw_trace = [199, 201, 199, 201]
        exp_trace = [200, 200, 200, 200]
        # mean/trace  m    m    m    m
        filtered_trace = self.mf.mean_filter_with_threshold(raw_trace)
        self.assertEqual(filtered_trace, exp_trace)

        # Large deviations in first few elements
        # (199 + 211 + 189 + 201) / 4. = 200
        # |211 - 200| > threshold or |189 - 200| > threshold
        raw_trace = [199, 211, 189, 201]
        exp_trace = [199, 211, 189, 201]
        # mean/trace  t    t    t    t
        filtered_trace = self.mf.mean_filter_with_threshold(raw_trace)
        self.assertEqual(filtered_trace, exp_trace)

        # Large jump in later elements
        # |236 - 201| > 2 * threshold
        raw_trace = [199, 201, 199, 201, 236]
        exp_trace = [200, 200, 200, 200, 236]
        # mean/trace  m    m    m    m    t
        filtered_trace = self.mf.mean_filter_with_threshold(raw_trace)
        self.assertEqual(filtered_trace, exp_trace)

        # Values on same side of local mean:
        # (201 + 199 + 201 + 202) / 4. = 200.75
        # 202 > 200.75 and 201 > 200.75
        raw_trace = [199, 201, 199, 201, 202]
        exp_trace = [200, 200, 200, 200, 202]
        # mean/trace  m    m    m    m    t
        filtered_trace = self.mf.mean_filter_with_threshold(raw_trace)
        self.assertEqual(filtered_trace, exp_trace)

        # Value far from local mean
        # 216 - (201 + 199 + 201 + 216) / 4. > threshold
        raw_trace = [199, 201, 199, 201, 216]
        exp_trace = [200, 200, 200, 200, 216]
        # mean/trace  m    m    m    m    t
        filtered_trace = self.mf.mean_filter_with_threshold(raw_trace)
        self.assertEqual(filtered_trace, exp_trace)

        # Value close to local mean:
        # 205 - (201 + 199 + 201 + 205) / 4. < threshold
        raw_trace = [199, 201, 199, 201, 205]
        exp_trace = [200, 200, 200, 200, 202]
        # mean/trace  m    m    m    m    m
        filtered_trace = self.mf.mean_filter_with_threshold(raw_trace)
        self.assertEqual(filtered_trace, exp_trace)

    def test_mean_filter_without_threshold(self):
        raw_trace = [199, 201, 199, 201, 216, 220, 219, 205, 200, 201]
        exp_trace = [200, 200, 200, 200, 204, 220, 219, 215, 200, 201]
        # mean/trace  m    m    m    m    m    t    t    m    t    t
        filtered_trace = self.mf.mean_filter_without_threshold(raw_trace)
        self.assertEqual(filtered_trace, exp_trace)


if __name__ == '__main__':
    unittest.main()
