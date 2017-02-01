import unittest
from itertools import cycle

from numpy import array
from mock import patch, sentinel, MagicMock

from sapphire.analysis import process_traces


class TraceObservablesTests(unittest.TestCase):

    def setUp(self):
        trace = ([200] * 400 + [500] + [510] + [400] * 10 + [200] * 600 +
                 [400] * 10 + [200])
        trace2 = ([203, 199] * 200 + [500] + [510] + [398, 402] * 5 +
                  [203, 199] * 300 + [400] * 10 + [200])
        self.traces = array([trace, trace2]).T
        self.to = process_traces.TraceObservables(self.traces)

    def test_baselines(self):
        self.assertEqual(self.to.baselines, [200, 201, -1, -1])

    def test_std_dev(self):
        self.assertEqual(self.to.std_dev, [0, 2000, -1, -1])

    def test_pulseheights(self):
        self.assertEqual(self.to.pulseheights, [310, 309, -1, -1])

    def test_integrals(self):
        self.assertEqual(self.to.integrals, [300 + 310 + 200 * 20,
                                             299 + 309 + 199 * 20, -1, -1])

    def test_n_peaks(self):
        self.assertEqual(self.to.n_peaks, [2, 2, -1, -1])


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

        # Check proper rounding of mean values ending in .5 to nearest even
        # (200 + 201 + 200 + 201) / 4. = 200.5 -> 200
        # (201 + 200 + 201 + 204) / 4. = 201.5 -> 202
        raw_trace = [200, 201, 200, 201, 204]
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

        # Check proper rounding of mean values ending in .5 to nearest even
        # (200 + 201 + 200 + 201) / 4. = 200.5 -> 200
        # (201 + 200 + 201 + 204) / 4. = 201.5 -> 202
        raw_trace = [200, 201, 200, 201, 204]
        exp_trace = [200, 200, 200, 200, 202]
        # mean/trace  m    m    m    m
        filtered_trace = self.mf.mean_filter_without_threshold(raw_trace)
        self.assertEqual(filtered_trace, exp_trace)


class DataReductionTests(unittest.TestCase):

    def setUp(self):
        self.dr = process_traces.DataReduction()

    def test_init(self):
        dr = process_traces.DataReduction(sentinel.threshold, sentinel.padding)
        self.assertEqual(dr.threshold, sentinel.threshold)
        self.assertEqual(dr.padding, sentinel.padding)

    def test_reduce_traces(self):
        pre = 400
        post = 300
        baseline = 200
        trace = ([baseline] * pre + [baseline + 50] + [baseline + 60] * 4 +
                 [baseline] * 600 + [baseline + 90] * 5 + [baseline] * post)
        traces = array([trace, [baseline] * len(trace)]).T
        reduced_traces = self.dr.reduce_traces(traces, [baseline] * 2)
        r_traces, left = self.dr.reduce_traces(traces, [baseline] * 2, True)
        r_traces_no_baseline = self.dr.reduce_traces(traces)
        self.assertTrue((reduced_traces == r_traces).all())
        self.assertTrue((reduced_traces == r_traces_no_baseline).all())
        self.assertEqual(len(reduced_traces),
                         len(trace) - pre - post + self.dr.padding * 2)
        self.assertEqual(left, pre - self.dr.padding)

        pre = 10
        post = 10
        baseline = 200
        trace = ([baseline] * pre + [baseline + 50] + [baseline + 60] * 4 +
                 [baseline] * 600 + [baseline + 90] * 5 + [baseline] * post)
        traces = array([trace, [baseline] * len(trace)]).T
        reduced_traces = self.dr.reduce_traces(traces, [baseline] * 2)
        r_traces, left = self.dr.reduce_traces(traces, [baseline] * 2, True)
        self.assertTrue((reduced_traces == traces).all())
        self.assertTrue((reduced_traces == r_traces).all())
        self.assertEqual(len(reduced_traces), len(trace))
        self.assertEqual(left, 0)

    def test_determine_cuts(self):
        pre = 400
        post = 300
        baseline = 200
        trace = ([baseline] * pre + [baseline + 50] + [baseline + 60] * 4 +
                 [baseline] * 600 + [baseline + 90] * 5 + [baseline] * post)
        traces = array([trace, [baseline] * len(trace)]).T
        left, right = self.dr.determine_cuts(traces, [baseline] * 2)
        self.assertEqual(left, pre)
        self.assertEqual(right, len(trace) - post)

        # No signal, return entire trace
        length = 400
        baseline = 200
        trace = [baseline] * length
        traces = array([trace, trace]).T
        left, right = self.dr.determine_cuts(traces, [baseline] * 2)
        self.assertEqual(left, 0)
        self.assertEqual(right, length)

    def test_add_padding(self):
        combinations = (((0, 20), (0, 46)),  # left at limit
                        ((4, 20), (0, 46)),  # left close to limit
                        ((50, 2400), (24, 2426)),  # left far from limit
                        ((50, 2400, 2400), (24, 2400)),  # right at limit
                        ((50, 2400, 2410), (24, 2410)),  # right close to limit
                        ((0, 200, 2400), (0, 226)),)  # right far from limit
        for input, expected in combinations:
            self.assertEqual(self.dr.add_padding(*input), expected)


if __name__ == '__main__':
    unittest.main()
