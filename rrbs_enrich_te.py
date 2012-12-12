#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import binom
import gzip, os, pdb, random, subprocess, sys
import pysam
import gff, stats

################################################################################
# rrbs_enrich_te.py
#
# Compute the density of CpG's and their methylation% across TEs, or
# alternatively across TEs in lincRNAs.
################################################################################

linc_gtf = '%s/research/common/data/lncrna/linc_catalog_trans.gtf' % os.environ['HOME']

cov_t = 1

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <rrbs files>'
    parser = OptionParser(usage)
    parser.add_option('-l', dest='lincrna', default=False, action='store_true', help='Only use lincRNA promoter-associated TEs [Default: %default]')
    parser.add_option('-r', dest='repeats_gff', default='%s/research/common/data/genomes/hg19/annotation/repeatmasker/hg19.fa.out.tp.gff' % os.environ['HOME'])
    (options,args) = parser.parse_args()

    if len(args) < 1:
        parser.error('Must provide a gff file for the feature of interest.')
    else:
        rrbs_files = args

    # count genomic bp
    genome_bp = count_hg19()

    # methylation genome-wide
    meth_pcts = [int(line.split()[10]) for rrbs_file in rrbs_files for line in open(rrbs_file) if not line.startswith('track') and int(line.split()[9]) > cov_t]

    if options.lincrna:
        tmp_gff = 'te_linc_%d.gff' % random.randint(0,1000000)
        p = subprocess.Popen('intersectBed -wa -u -a %s -b %s > %s' % (options.repeats_gff, linc_gtf, tmp_gff), shell=True)
        os.waitpid(p.pid,0)
        options.repeats_gff = tmp_gff

    # hash counted repeat genomic bp
    te_lengths = measure_te(options.repeats_gff)

    te_pct = {}
    for rrbs_file in rrbs_files:
        # intersect
        proc = subprocess.Popen('intersectBed -wo -a %s -b %s' % (rrbs_file,options.repeats_gff), shell=True, stdout=subprocess.PIPE)

        # hash stats by TE family
        line = proc.stdout.readline()
        while line:
            a = line.split('\t')

            if int(a[9]) > cov_t:
                te_kv = gff.gtf_kv(a[19])

                # save %
                te_pct.setdefault((te_kv['repeat'],te_kv['family']), []).append(int(a[10]))
                te_pct.setdefault(('*',te_kv['family']), []).append(int(a[10]))
                te_pct.setdefault(('*','*'), []).append(int(a[10]))

            line = proc.stdout.readline()
        proc.communicate()

    if options.lincrna:
        os.remove(tmp_gff)

    meth_mean = stats.mean(meth_pcts)
    print >> sys.stderr, 'Genome-wide methylation  %d  %.3f' % (len(meth_pcts),meth_mean/100)

    # compute stats, print table
    for te in te_pct:
        # binomial test
        # p = TE sum length / genome_bp
        # N = # genome CpGs
        # x = # TE CpGs

        p = float(te_lengths[te]) / genome_bp
        N = len(meth_pcts)
        x = len(te_pct[te])

        cpg_fold = x / float(p*N)

        if cpg_fold > 1:
            cpg_p = binom.sf(x-1, N, p)
        else:
            cpg_p = binom.cdf(x, N, p)

        # i've got samples from a distribution
        # i've got the distribution
        # are the samples from that distribution greater or lesser
        pct_fold = stats.mean(te_pct[te]) / stats.mean(meth_pcts)
        #pct_z, pct_p = stats.mannwhitneyu(te_pct[te], meth_pcts)

        #cols = (te[0], te[1], te_lengths[te], float(x)/te_lengths[te], cpg_fold, cpg_p, stats.mean(te_pct[te]), pct_fold, pct_z, pct_p)
        cols = (te[0], te[1], te_lengths[te], float(x)/te_lengths[te], cpg_fold, cpg_p, stats.mean(te_pct[te]), pct_fold)

        print '%-18s %-18s %10d %7.3f %7.3f %10.2e %7.3f %7.3f' % cols


################################################################################
# count_hg19
#
# Count the number of bp in hg19 where TEs could be.
################################################################################
def count_hg19():
    chrom_sizes_file = '%s/research/common/data/genomes/hg19/assembly/human.hg19.genome' % os.environ['HOME']
    gap_bed_file = '%s/research/common/data/genomes/hg19/assembly/hg19_gaps.bed' % os.environ['HOME']
    valid_chrs = ['chr%d' % c for c in range(1,23)] + ['chrX','chrY']

    genome_bp = 0
    for line in open(chrom_sizes_file):        
        a = line.split()
        if len(a) > 0 and a[0] in valid_chrs:
            genome_bp += int(a[1])

    for line in open(gap_bed_file):
        a = line.split()
        if a[0] in valid_chrs:
            genome_bp -= int(a[2])-int(a[1])

    return genome_bp


################################################################################
# measure_te
#
# Hash the number of bp covered by various repeats in the RepeatMasker gff file
# and the lincRNA gtf file.
################################################################################
def measure_te(rm_file):
    repeat_bp = {}
    for line in open(rm_file):
        a = line.split('\t')

        kv = gff.gtf_kv(a[8])
        rep = kv['repeat']
        family = kv['family']

        length = int(a[4]) - int(a[3]) + 1

        repeat_bp[(rep,family)] = repeat_bp.get((rep,family),0) + length
        repeat_bp[('*',family)] = repeat_bp.get(('*',family),0) + length
        repeat_bp[('*','*')] = repeat_bp.get(('*','*'),0) + length

    return repeat_bp


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
