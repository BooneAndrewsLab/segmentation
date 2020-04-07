import unittest

import segmentation


class SegmentationTestCase(unittest.TestCase):
    def test_imported(self):
        _ = segmentation.mixture_model  # just faking it for now
        self.assertTrue(True)
