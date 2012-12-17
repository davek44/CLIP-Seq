#!/usr/bin/env python
from optparse import OptionParser
import unittest

################################################################################
# test_clip_peaks.py
#
#
################################################################################


################################################################################
# convolute_lambda
################################################################################
class TestConvoluteLambda(unittest.TestCase):
    def setUp(self):
        self.fpkm_exon = 2
        self.fpkm_pre = 1
        self.window_size = 4

    def test1:
        gene = 'EEEEEEIIIEEIIIIIIEEE'
        #       ....
        junctions = [7,10,12,18]
        exon_counts = [4,4,4,3,2,1,1,2,2,2,1,0,0,0,1,2,3]

        true_lambda = [float(ec*self.fpkm_exon + (self.window_size-ec)*self.fpkm_pre)/self.window_size for ec in exon_counts]

        code_lambda = [0]*len(true_lambda)
        ji = 0
        for window_start in range(len(gene)-self.window_size+1):
            # update junctions index
            if ji < len(junctions) and junctions[ji] <= window_start:
                ji += 1

            # compute lambda
            code_lambda[window_start] = clip_peaks.convolute_lambda(window_start, window_start+self.window_size-1, self.fpkm_exon, self.fpkm_pre, junctions, 0)

        for i in range(len(true_lambda)):
            self.assertEqual(true_lambda[i], code_lambda[i])

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
