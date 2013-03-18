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


    ############################################################
    # compute_code_lambdas
    ############################################################
    def compute_code_lambdas(self):
        # initialize
        code_lambda = []
        ji = {}
        for tid in self.gene_transcripts:
            ji[tid] = 0
            
        # for each window
        gene_len = len(self.isoform1.labels)
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


    ############################################################
    # compute_true_lambdas
    #
    # Compute lambdas rigorously.
    ############################################################
    def compute_true_lambdas(self):
        isoforms = self.gene_transcripts.values()

        # initialize
        true_lambda = [0.0]*(len(isoforms[0].labels)-self.window_size+1)

        # consider each isoform
        for isoform in isoforms:
            exon_counts = self.count_exons(isoform.labels)
            for i in range(len(exon_counts)):
                true_lambda[i] += float(exon_counts[i])/self.window_size*isoform.fpkm

        # normalize to lambda
        true_lambda = [fpkm/1000.0*self.total_reads/1000000.0 for fpkm in true_lambda]

        return true_lambda


    ############################################################
    # count_exons
    # 
    # Count the number of exon bp in each window
    # across the transcript.
    ############################################################
    def count_exons(self, transcript):
        exon_counts = [0]*(len(transcript)-self.window_size+1)
        for i in range(len(exon_counts)):
            for j in range(i,i+self.window_size):
                exon_counts[i] += int(transcript[j] == 'E')
        return exon_counts


    ############################################################
    # find_junctions
    # 
    # Find the 1-based indexes of the exon-intron junctions
    # in the given transcript.
    ############################################################
    def find_junctions(self, transcript):
        junctions = [1]
        for i in range(1,len(transcript)):
            if transcript[i] != transcript[i-1]:
                junctions.append(i+1)
        junctions.append(len(transcript)+1)
        return junctions


    ############################################################
    def test0(self):
        # no junctions
        self.isoform1.labels = 'EEEEEEEE'
        self.isoform1.junctions = self.find_junctions(self.isoform1.labels)
        self.pre_isoform1.labels = 'E'*len(self.isoform1.labels)
        self.pre_isoform1.junctions = self.find_junctions(self.pre_isoform1.labels)

        true_lambda = self.compute_true_lambdas()
        code_lambda = self.compute_code_lambdas()

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertEqual(true_lambda[i], code_lambda[i])

    ############################################################
    def test1(self):
        # bunch of junctions
        self.isoform1.labels = 'EEEEEEIIIEEIIIIIIEEE'
        self.isoform1.junctions = self.find_junctions(self.isoform1.labels)
        self.pre_isoform1.labels = 'E'*len(self.isoform1.labels)
        self.pre_isoform1.junctions = self.find_junctions(self.pre_isoform1.labels)

        true_lambda = self.compute_true_lambdas()
        code_lambda = self.compute_code_lambdas()

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertTrue(abs(true_lambda[i] - code_lambda[i]) < 1e-9)

    ############################################################
    def test2(self):
        # bunch of junctions
        self.isoform1.labels = 'EEEIIIIIIEEIIIEEEEEE'
        self.isoform1.junctions = self.find_junctions(self.isoform1.labels)
        self.pre_isoform1.labels = 'E'*len(self.isoform1.labels)
        self.pre_isoform1.junctions = self.find_junctions(self.pre_isoform1.labels)

        true_lambda = self.compute_true_lambdas()
        code_lambda = self.compute_code_lambdas()

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertTrue(abs(true_lambda[i] - code_lambda[i]) < 1e-9)

    ############################################################
    def test3(self):
        self.isoform1.labels = 'EEEIIIIIIEEIIIEEEEEE'
        self.isoform1.junctions = self.find_junctions(self.isoform1.labels)

        self.isoform2 = clip_peaks.Gene('chr1', '+', {})
        self.isoform2.fpkm = 3
        self.isoform2.labels = 'EEEIIIIIEEIIIIEEEEEE'
        self.isoform2.junctions = self.find_junctions(self.isoform2.labels)

        self.pre_isoform1.labels = 'E'*len(self.isoform1.labels)
        self.pre_isoform1.junctions = self.find_junctions(self.pre_isoform1.labels)

        self.gene_transcripts = {'isoform1':self.isoform1, 'isoform2':self.isoform2, 'pre_isoform1':self.pre_isoform1}
        
        true_lambda = self.compute_true_lambdas()
        code_lambda = self.compute_code_lambdas()

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertTrue(abs(true_lambda[i] - code_lambda[i]) < 1e-9)

    ############################################################
    def test4(self):
        # multiple TSSs
        self.isoform1.labels = 'EEEEEEEEEE'
        self.isoform1.junctions = self.find_junctions(self.isoform1.labels)

        self.isoform2 = clip_peaks.Gene('chr1', '+', {})
        self.isoform2.fpkm = 3
        self.isoform2.labels = 'IIEEEEEEEE'
        self.isoform2.junctions = [3,11]

        self.pre_isoform1.labels = 'E'*len(self.isoform1.labels)
        self.pre_isoform1.junctions = self.find_junctions(self.pre_isoform1.labels)

        self.gene_transcripts = {'isoform1':self.isoform1, 'isoform2':self.isoform2, 'pre_isoform1':self.pre_isoform1}

        true_lambda = self.compute_true_lambdas()
        code_lambda = self.compute_code_lambdas()

        self.assertEqual(len(true_lambda),len(code_lambda))
        for i in range(len(true_lambda)):
            self.assertTrue(abs(true_lambda[i] - code_lambda[i]) < 1e-9)

    ############################################################
    def test5(self):
        # multiple endpoints
        self.isoform1.labels = 'EEEEEEEEEE'
        self.isoform1.junctions = self.find_junctions(self.isoform1.labels)

        self.isoform2 = clip_peaks.Gene('chr1', '+', {})
        self.isoform2.fpkm = 3
        self.isoform2.labels = 'EEEEEEEEII'
        self.isoform2.junctions = [1,9]

        self.pre_isoform1.labels = 'E'*len(self.isoform1.labels)
        self.pre_isoform1.junctions = self.find_junctions(self.pre_isoform1.labels)

        self.gene_transcripts = {'isoform1':self.isoform1, 'isoform2':self.isoform2, 'pre_isoform1':self.pre_isoform1}

        true_lambda = self.compute_true_lambdas()
        code_lambda = self.compute_code_lambdas()

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
        read_midpoints = [3,3,4,5,6,6]
        self.read_pos_weights = [(read_midpoints[i],1) for i in range(len(read_midpoints))]

        self.tx = clip_peaks.Gene('chr1','+',{})
        self.tx.add_exon(1,10)
        self.tx.fpkm = 1
        self.gene_transcripts = {'tx':self.tx}

        self.poisson_lambda = self.tx.fpkm / 1000.0 * self.total_reads / 1000000.0

    def test1(self):        
        window_counts = [3,4,6,4,2]
        window_stats = [(wc,.001) for wc in window_counts]

        true_peaks = [(3,6,6,clip_peaks.scan_stat_approx3(6, self.window_size, self.txome_size, self.poisson_lambda))]
        code_peaks = clip_peaks.windows2peaks(self.read_pos_weights, self.gene_transcripts, 1, window_stats, self.window_size, .05, self.total_reads, self.txome_size)

        self.assertEqual(len(true_peaks),len(code_peaks))
        for i in range(len(true_peaks)):
            self.assertEqual(true_peaks[i],code_peaks[i])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    unittest.main()
    #pdb.runcall(unittest.main())
