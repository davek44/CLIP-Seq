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
        self.isoform1 = clip_peaks.Gene('chr1', '+', {})
        self.isoform1.fpkm = 1

        self.pre_isoform1 = clip_peaks.Gene('chr1', '+', {})
        self.pre_isoform1.fpkm = 1
        self.pre_isoform1.junctions = []

        self.gene_transcripts = {'isoform1':self.isoform1, 'pre_isoform1':self.pre_isoform1}

        self.window_size = 4
        self.total_reads = 100

    def compute_lambdas(self, gene_len):
        code_lambda = []
        ji = {}
        for tid in self.gene_transcripts:
            ji[tid] = 0
        for window_start in range(1,gene_len-self.window_size+2):
            # update junctions index
            for tid in self.gene_transcripts:
                tjunctions = self.gene_transcripts[tid].junctions
                if ji[tid] < len(tjunctions) and tjunctions[ji[tid]] <= window_start:
                    ji[tid] += 1

            # compute lambda
            code_lambda.append(clip_peaks.convolute_lambda(window_start, window_start+self.window_size-1, self.gene_transcripts, ji, self.total_reads))
            if False:
                pdb.runcall(code_lambda.append(clip_peaks.convolute_lambda(window_start, window_start+self.window_size-1, self.gene_transcripts, ji, self.total_reads)))

        return code_lambda

    def test0(self):
        # no junctions
        isoform1 = 'EEEEEEEE'
        #           ....
        isoform1_exon_counts = [4]*5

        self.isoform1.junctions = []

        # iniitialize lambdas
        true_lambda = [0]*len(isoform1_exon_counts)

        # add isoform1 to lambda
        for i in range(len(isoform1_exon_counts)):
            ec = isoform1_exon_counts[i]
            true_lambda[i] += float(ec)/self.window_size*self.isoform1.fpkm

        # add pre_isoform1 to lambda
        for i in range(len(isoform1_exon_counts)):
            true_lambda[i] += self.pre_isoform1.fpkm

        # normalize to lambda
        true_lambda = [fpkm/1000.0*self.total_reads/1000000.0 for fpkm in true_lambda]

        code_lambda = self.compute_lambdas(len(isoform1))

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertEqual(true_lambda[i], code_lambda[i])

    def test1(self):
        # bunch of junctions
        isoform1 = 'EEEEEEIIIEEIIIIIIEEE'
        #               ....
        isoform1_exon_counts = [4,4,4,3,2,1,1,2,2,2,1,0,0,0,1,2,3]

        self.isoform1.junctions = [7,10,12,18]

        # iniitialize lambdas
        true_lambda = [0]*len(isoform1_exon_counts)

        # add isoform1 to lambda
        for i in range(len(isoform1_exon_counts)):
            ec = isoform1_exon_counts[i]
            true_lambda[i] += float(ec)/self.window_size*self.isoform1.fpkm

        # add pre_isoform1 to lambda
        for i in range(len(isoform1_exon_counts)):
            true_lambda[i] += self.pre_isoform1.fpkm

        # normalize to lambda
        true_lambda = [fpkm/1000.0*self.total_reads/1000000.0 for fpkm in true_lambda]

        code_lambda = self.compute_lambdas(len(isoform1))

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertTrue(abs(true_lambda[i] - code_lambda[i]) < 1e-9)

    def test2(self):
        # bunch of junctions
        isoform1 = 'EEEIIIIIIEEIIIEEEEEE'
        #                       ....        
        isoform1_exon_counts = [3,2,1,0,0,0,1,2,2,2,1,1,2,3,4,4,4]

        self.isoform1.junctions = [4,10,12,15]

        # iniitialize lambdas
        true_lambda = [0]*len(isoform1_exon_counts)

        # add isoform1 to lambda
        for i in range(len(isoform1_exon_counts)):
            ec = isoform1_exon_counts[i]
            true_lambda[i] += float(ec)/self.window_size*self.isoform1.fpkm

        # add pre_isoform1 to lambda
        for i in range(len(isoform1_exon_counts)):
            true_lambda[i] += self.pre_isoform1.fpkm

        # normalize to lambda
        true_lambda = [fpkm/1000.0*self.total_reads/1000000.0 for fpkm in true_lambda]

        code_lambda = self.compute_lambdas(len(isoform1))

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertTrue(abs(true_lambda[i] - code_lambda[i]) < 1e-9)


################################################################################
# windows2peaks
################################################################################
class TestWindows2Peaks(unittest.TestCase):
    def setUp(self):
        self.txome_size = 8
        self.total_reads = 100
        self.window_size = 4
        self.read_midpoints = [3,3,4,5,6,6]

        self.tx = clip_peaks.Gene('chr1','+',{})
        self.tx.add_exon(1,10)
        self.tx.fpkm = 1
        self.gene_transcripts = {'tx':self.tx}

        self.poisson_lambda = self.tx.fpkm / 1000.0 * self.total_reads / 1000000.0

    def test1(self):        
        window_counts = [3,4,6,4,2]
        window_stats = [(wc,.001) for wc in window_counts]

        true_peaks = [(3,6,6,clip_peaks.scan_stat_approx3(6, self.window_size, self.txome_size, self.poisson_lambda))]
        #code_peaks = clip_peaks.windows2peaks(self.read_midpoints, [], window_stats, self.window_size, .05, self.tx, self.txome_size)
        code_peaks = clip_peaks.windows2peaks(self.read_midpoints, self.gene_transcripts, 1, window_stats, self.window_size, .05, self.total_reads, self.txome_size)

        self.assertEqual(len(true_peaks),len(code_peaks))
        for i in range(len(true_peaks)):
            self.assertEqual(true_peaks[i],code_peaks[i])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    unittest.main()
