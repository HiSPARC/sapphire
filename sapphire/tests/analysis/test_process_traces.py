import unittest

from mock import patch, sentinel

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
