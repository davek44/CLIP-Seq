#!/usr/bin/env python
from optparse import OptionParser
import pdb, unittest
import clip_peaks

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
        self.fpkm_exon = 1
        self.fpkm_pre = 1
        self.window_size = 4

    def compute_lambdas(self, gene_len, junctions):
        code_lambda = []
        ji = 0
        for window_start in range(1,gene_len-self.window_size+2):
            # update junctions index
            if ji < len(junctions) and junctions[ji] <= window_start:
                ji += 1

            # compute lambda
            code_lambda.append(clip_peaks.convolute_lambda(window_start, window_start+self.window_size-1, self.fpkm_exon, self.fpkm_pre, junctions, ji))
            if False:
                pdb.runcall(clip_peaks.convolute_lambda, window_start, window_start+self.window_size-1, self.fpkm_exon, self.fpkm_pre, junctions, ji)

        return code_lambda

    def test0(self):
        # no junctions
        gene = 'EEEEEEEE'
        #           ....
        junctions = []
        exon_counts = [4]*5

        true_lambda = [float(ec*(self.fpkm_exon+self.fpkm_pre) + (self.window_size-ec)*self.fpkm_pre)/self.window_size/1000.0/1000000.0 for ec in exon_counts]
        code_lambda = self.compute_lambdas(len(gene),junctions)

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertEqual(true_lambda[i], code_lambda[i])

    def test1(self):
        # bunch of junctions
        gene = 'EEEEEEIIIEEIIIIIIEEE'
        #               ....
        junctions = [7,10,12,18]
        exon_counts = [4,4,4,3,2,1,1,2,2,2,1,0,0,0,1,2,3]

        true_lambda = [float(ec*(self.fpkm_exon+self.fpkm_pre) + (self.window_size-ec)*self.fpkm_pre)/self.window_size/1000.0/1000000.0 for ec in exon_counts]
        code_lambda = self.compute_lambdas(len(gene),junctions)

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertEqual(true_lambda[i], code_lambda[i])

    def test2(self):
        # bunch of junctions
        gene = 'EEEIIIIIIEEIIIEEEEEE'
        #                       ....
        junctions = [4,10,12,15]
        exon_counts = [3,2,1,0,0,0,1,2,2,2,1,1,2,3,4,4,4]

        true_lambda = [float(ec*(self.fpkm_exon+self.fpkm_pre) + (self.window_size-ec)*self.fpkm_pre)/self.window_size/1000.0/1000000.0 for ec in exon_counts]
        code_lambda = self.compute_lambdas(len(gene),junctions)

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertEqual(true_lambda[i], code_lambda[i])


################################################################################
# windows2peaks
################################################################################
class TestWindows2Peaks(unittest.TestCase):
    def setUp(self):
        self.txome_size = 8
        self.window_size = 4
        self.read_midpoints = [3,3,4,5,6,6]

        self.tx = clip_peaks.Gene('chr1','+',{})
        self.tx.add_exon(1,10)
        self.tx.fpkm_exon = 1
        self.tx.fpkm_pre = 0
        self.poisson_lambda = self.tx.fpkm_exon / 1000.0 / 1000000.0

    def test1(self):        
        window_counts = [3,4,6,4,2]
        window_stats = [(wc,.001) for wc in window_counts]

        true_peaks = [(3,6,6,clip_peaks.scan_stat_approx3(6, self.window_size, self.txome_size, self.poisson_lambda))]
        code_peaks = clip_peaks.windows2peaks(self.read_midpoints, [], window_stats, self.window_size, .05, self.tx, self.txome_size)

        self.assertEqual(len(true_peaks),len(code_peaks))
        for i in range(len(true_peaks)):
            self.assertEqual(true_peaks[i],code_peaks[i])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    unittest.main()
