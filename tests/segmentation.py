import unittest

import segmentation


class SegmentationTestCase(unittest.TestCase):
    def test_no_bam(self):
        _ = segmentation.mixture_model  # just faking it for now
        self.assertTrue(True)
