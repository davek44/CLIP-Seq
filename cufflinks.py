#!/usr/bin/env python
from optparse import OptionParser
from numpy import array, empty
from scipy.stats import norm
import math, sys
import stats

import rpy2
from rpy2.robjects.numpy2ri import numpy2ri
import rpy2.robjects as ro

ro.conversion.py2ri = numpy2ri

################################################################################
# cufflinks.py
#
# Code to support analysis of cufflinks output.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()


################################################################################
# fpkm_tracking
################################################################################
class fpkm_tracking:
    ############################################################################
    # Constructor
    #
    # Load the expression matrix 
    ############################################################################
    def __init__(self, fpkm_file='/Users/dk/research/common/data/lncrna/genes.fpkm_tracking', rand_zero=False):
        # obtain basic information
        fpkm_in = open(fpkm_file)
        headers = fpkm_in.readline().split()
        self.genes = []
        line = fpkm_in.readline()
        while line:
            a = line.split()
            self.genes.append(a[0])
            line = fpkm_in.readline()
        fpkm_in.close()

        self.gene_map = dict([(self.genes[i],i) for i in range(len(self.genes))])
        self.experiments = [h[:-5] for h in headers if h[-5:] == '_FPKM']

        # obtain expression
        self.expr = empty([len(self.genes), len(self.experiments)])
        g = 0

        fpkm_in = open(fpkm_file)
        line = fpkm_in.readline()
        line = fpkm_in.readline()
        while line:
            a = line.split('\t')
            a[-1] = a[-1].rstrip()
            e = 0
            for i in range(len(headers)):
                if headers[i][-5:] == '_FPKM':
                    self.expr[g,e] = float(a[i])
                    if rand_zero and self.expr[g,e] == 0:
                        self.expr[g,e] = abs(norm.rvs(scale=1e-6))
                    e += 1
            g += 1
            line = fpkm_in.readline()
        fpkm_in.close()

        print >> sys.stderr, 'Loaded expression of %d genes in %d experiments' % (g,e)


    ############################################################################
    # gene_entropy
    #
    # Return the entropy of the expression vector for the given gene.
    # Note:
    #  -Log creates negative values so it's not a distribution for sure.
    ############################################################################
    def gene_entropy(self, gene, log=False):
        gene_i = self.name_or_index(gene)

        gexpr = self.expr[gene_i,:]
        if log:
            gexpr = [math.log(e+1) for e in gexpr]
        gexpr = stats.normalize(gexpr)

        return stats.entropy(gexpr)


    ############################################################################
    # gene_expr
    #
    # Return an expression vector for the given gene.
    ############################################################################
    def gene_expr(self, gene):
        gene_i = self.name_or_index(gene)
        if gene_i:
            return self.expr[gene_i,:]
        else:
            return [0]*len(self.experiments)


    ############################################################################
    # gene_expr_exp
    #
    # Return FPKM for a given gene in a given experiment.
    ############################################################################
    def gene_expr_exp(self, gene, exp):
        gene_i = self.name_or_index(gene)
        if gene_i == None:
            print >> sys.stderr, '%s expression not found' % gene
            return 0
        else:
            for exp_i in range(len(self.experiments)):
                if self.experiments[exp_i] == exp:
                    return self.expr[gene_i,exp_i]
            print >> sys.stderr, '%s expression not found' % exp
            return 0


    ############################################################################
    # gene_expr_print
    #
    # Print expression data for the given gene.
    ############################################################################
    def gene_expr_print(self, gene):
        gene_i = self.name_or_index(gene)
            
        for j in range(len(self.experiments)):
            print '%-15s %8.3f' % (self.experiments[j], self.expr[gene_i,j])


    ############################################################################
    # gene_specificity
    #
    # Return tissue specificity for the given gene.
    ############################################################################
    def gene_specificity(self, gene, log=True):
        gene_i = self.name_or_index(gene)

        if gene_i == None:
            spec = 0
        else:
            gexpr = self.expr[gene_i,:]
            if log:
                gexpr = [math.log(e+1) for e in gexpr]

            if sum(gexpr) == 0:
                spec = 0
            else:
                gexpr = stats.normalize(gexpr)

                min_jsd = 1.0
                for j in range(len(self.experiments)):
                    q_j = [0]*len(self.experiments)
                    q_j[j] = 1.0

                    min_jsd = min(min_jsd, math.sqrt(stats.jsd(gexpr, q_j)))

                spec = 1.0 - min_jsd

        return spec


    ############################################################################
    # genes_jsd
    #
    # Jensen-Shannon divergence between two genes
    ############################################################################
    def genes_jsd(self, gene1, gene2, log=True):
        gene1_i = self.name_or_index(gene1)
        gene2_i = self.name_or_index(gene2)

        gexpr1 = self.expr[gene1_i,:]
        if log:
            gexpr1 = [math.log(e+1) for e in gexpr1]
        gexpr1 = stats.normalize(gexpr1)

        gexpr2 = self.expr[gene2_i,:]
        if log:
            gexpr2 = [math.log(e+1) for e in gexpr2]
        gexpr2 = stats.normalize(gexpr2)

        return stats.jsd(gexpr1, gexpr2)
            

    ############################################################################
    # name_or_index
    #
    # Given a name or index, return an index
    ############################################################################
    def name_or_index(self, gene):
        if type(gene) == str:
            if gene in self.gene_map:
                return self.gene_map[gene]
            else:
                print >> sys.stderr, 'Missing gene - %s' % gene
                return None
        elif type(gene) == int:
            return gene
        else:
            print >> sys.stderr, 'Bad gene input'
            return None

    
    ############################################################################
    # spearman
    #
    # Compute Spearman correlations for either all pairs of genes or one given
    # gene to all others.
    ############################################################################
    def spearman(self, gene=None):
        if type(gene) == str:
            cors = ro.r.cor(self.expr[self.gene_map[gene],:], self.expr.transpose(), method='spearman')
            return array(cors)[0]

        elif type(gene) == int:
            cors = ro.r.cor(self.expr[gene,:], self.expr.transpose(), method='spearman')
            return array(cors)[0]

        else:
            cors = ro.r.cor(self.expr.transpose(), method='spearman')
            return array(cors)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
