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

        true_lambda = [float(ec*(self.fpkm_exon+self.fpkm_pre) + (self.window_size-ec)*self.fpkm_pre)/self.window_size for ec in exon_counts]
        code_lambda = self.compute_lambdas(len(gene),junctions)

        for i in range(len(true_lambda)):
            self.assertEqual(true_lambda[i], code_lambda[i])


    def test1(self):
        # bunch of junctions
        gene = 'EEEEEEIIIEEIIIIIIEEE'
        #               ....
        junctions = [7,10,12,18]
        exon_counts = [4,4,4,3,2,1,1,2,2,2,1,0,0,0,1,2,3]

        true_lambda = [float(ec*(self.fpkm_exon+self.fpkm_pre) + (self.window_size-ec)*self.fpkm_pre)/self.window_size for ec in exon_counts]
        code_lambda = self.compute_lambdas(len(gene),junctions)

        for i in range(len(true_lambda)):
            self.assertEqual(true_lambda[i], code_lambda[i])

    def test2(self):
        # bunch of junctions
        gene = 'EEEIIIIIIEEIIIEEEEEE'
        #                       ....
        junctions = [4,10,12,15]
        exon_counts = [3,2,1,0,0,0,1,2,2,2,1,1,2,3,4,4,4]

        true_lambda = [float(ec*(self.fpkm_exon+self.fpkm_pre) + (self.window_size-ec)*self.fpkm_pre)/self.window_size for ec in exon_counts]
        code_lambda = self.compute_lambdas(len(gene),junctions)

        for i in range(len(true_lambda)):
            self.assertEqual(true_lambda[i], code_lambda[i])

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    unittest.main()
