#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import poisson
import copy, math, os, subprocess, sys, tempfile
import pysam
import gff

################################################################################
# clip_peaks.py
#
# Call peaks in CLIP-Seq data.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <clip_bam> <ref_gtf>'
    parser = OptionParser(usage)

    parser.add_option('-c', dest='control_bam', help='Control BAM file')

    parser.add_option('-w', dest='window_size', type='int', default=50, help='Window size for scan statistic [Default: %default]')
    parser.add_option('-p', dest='p_val', type='float', default=.01, help='P-value required of window scan statistic tests [Default: %default]')

    parser.add_option('--cuff_done', dest='cuff_done', action='store_true', default=False, help='A cufflinks run to estimate the model parameters is already done [Default: %default]')
    parser.add_option('-t', dest='threads', type='int', default=2, help='Number of threads to use [Default: %default]')

    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error(usage)
    else:
        clip_bam = args[0]
        ref_gtf = args[1]


    ############################################
    # parameterize
    ############################################
    if options.control_bam:
        # make a new gtf w/ unspliced RNAs
        pre_ref_gtf = prerna_gtf(ref_gtf)

        # run Cufflinks on new gtf file and control BAM
        if not options.cuff_done:
            subprocess.call('cufflinks -p %d -G %s %s' % (options.threads, pre_ref_gtf, options.control_bam), shell=True)

        # choose most expressed isoform and save spliced/unspliced fpkms
        isoform_fpkms = control_max_fpkm(pre_ref_gtf)

        # filter gtf for clip analysis
        clip_ref_gtf = isoform_gtf_filter(ref_gtf, isoform_fpkms)
    
    else:
        # make a new gtf file of only loci-spanning RNAs
        span_ref_gtf = span_gtf(ref_gtf)

        # run Cufflinks on new gtf file and CLIP BAM
        if not options.cuff_done:
            subprocess.call('cufflinks -p %d -G %s %s' % (options.threads, span_ref_gtf, clip_bam), shell=True)

        # not sure how to represent parameterization here...
        isoform_fpkms = {}

        # filter gtf for clip analysis
        clip_ref_gtf = span_ref_gtf


    ############################################
    # process genes
    ############################################
    # read genes
    transcripts = read_genes(clip_ref_gtf, key_id='transcript_id')

    # open clip-seq bam
    clip_in = pysam.Samfile(clip_bam, 'rb')

    # for each span
    for transcript_id in transcripts:
        tx = transcripts[transcript_id]

        # map reads to midpoints
        read_midpoints = map_midpoints(clip_in, tx.chrom, tx.exons[0].start, tx.exons[-1].end, tx.strand)

        # find splice junctions
        junctions = get_splice_junctions(tx)

        # count reads and compute p-values in windows
        window_stats = count_windows(tx.exons[0].start, tx.exons[-1].end, options.window_size, read_midpoints)

        # post-process windows to peaks
        allowed_sig_gap = 1
        peaks = windows2peaks(window_stats, options.p_val, allowed_sig_gap, tx.exons[0].start)

        # output peaks
        # ...

    clip_in.close()
        

################################################################################
# cigar_midpoint
# 
# Returned the aligned read's midpoint, considering the insertions and
# deletions in its CIGAR string (which includes splicing).
################################################################################
def cigar_midpoint(aligned_read):
    read_half = aligned_read.qlen / 2.0
    read_walked = 0
    genome_pos = aligned_read.pos

    for (operation,length) in aligned_read.cigar:
        # match
        if operation in [0,7,8]:
            if read_walked + length >= read_half:
                midpoint = genome_pos + (read_half - read_walked)
                break
            else:
                genome_pos += length
                read_walked += length

        # insertion
        elif operation in [1,3]:
            genome_pos += length

        # deletion
        elif operation == 2:
            read_walked += length

        else:
            print >> sys.stderr, 'Unknown CIGAR operation - %d, %s' % (operation, aligned_read.qname)

    return midpoint


################################################################################
# control_max_fpkm
#
# Choose the most expressed isoform for gene and return a hash mapping
# transcript_id's to tuples of the isoform's FPKM and its corresponding
# unspliced version's FPKM.
################################################################################
def control_max_fpkm(ref_gtf):
    # collect transcript spans to map spliced isoforms to their pre-RNA
    transcript_span = {}
    transcript_length = {}
    span_unspliced = {}
    transcripts = read_genes(ref_gtf, key_id='transcript_id')
    for tid in transcripts:
        tx = transcripts[tid]

        span_start = tx.exons[0].start
        span_end = tx.exons[-1].end

        transcript_span[tid] = (tx.chrom, span_start, span_end, tx.strand)
        transcript_length[tid] = 0
        for exon in tx.exons:
            transcript_length[tid] += exon.end-exon.start+1

        if len(tx.exons) == 1:
            span_unspliced[transcript_span[tid]] = tid

    # read FPKMS
    transcript_fpkm = {}
    fpkm_tracking_in = open('isoforms.fpkm_tracking')
    line = fpkm_tracking_in.readline()
    for line in fpkm_tracking_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        transcript_id = a[0]
        fpkm = float(a[9])

        transcript_fpkm[transcript_id] = fpkm
    fpkm_tracking_in.close()

    # choose isoform by FPKM and length
    isoform_fpkms = {}
    g2t = gff.g2t(ref_gtf)
    for gene_id in g2t:
        # focus on spliced transcripts
        spliced_transcripts = [tid for tid in g2t[gene_id] if not tid.startswith('UNSPLICED')]

        # verify abundance estimates
        for tid in spliced_transcripts:
            if not tid in transcript_fpkm:
                # this can happen if two transcripts are the same in the GTF file.
                # cufflinks will ignore one of them.
                print >> sys.stderr, 'WARNING: Missing FPKM for spliced transcript %s' % tid
                transcript_fpkm[tid] = 0

        # choose isoform
        max_tid = spliced_transcripts[0]
        for tid in spliced_transcripts[1:]:
            if transcript_fpkm[tid] > transcript_fpkm[max_tid]:
                max_tid = tid
            elif transcript_fpkm[tid] == transcript_fpkm[max_tid]:
                if transcript_length[tid] > transcript_length[max_tid]:
                    max_tid = tid

        # find unspliced FPKM
        unspliced_tid = span_unspliced[transcript_span[max_tid]]
        if unspliced_tid not in transcript_fpkm:
            # this can happen if two transcripts are the same except for differing strands
            # cufflinks will ignore one of them.
            print >> sys.stderr, 'WARNING: Missing FPKM for unspliced transcript %s' % unspliced_tid
            unspliced_fpkm = transcript_fpkm[max_tid] / 2.0
        else:
            unspliced_fpkm = transcript_fpkm[unspliced_tid]
        
        isoform_fpkms[max_tid] = (transcript_fpkm[max_tid], unspliced_fpkm)

    return isoform_fpkms


################################################################################
# count_windows
#
# Count the number of reads and compute the scan statistic p-value in each
# window through the gene.
################################################################################
def count_windows(gene_start, gene_end, window_size, read_midpoints):
    # set lambda using whole region
    poisson_lambda = float(len(read_midpoints)) / (gene_end - gene_start)

    midpoints_window_start = 0 # index of the first read_midpoint that fit in the window (except I'm allowing 0)
    midpoints_window_end = 0 # index of the first read_midpoint past the window

    window_stats = []

    for window_start in range(gene_start, gene_end-window_size+1):
        window_end = window_start + window_size

        # update midpoints start
        while midpoints_window_start < len(read_midpoints) and read_midpoints[midpoints_window_start] < window_start:
            midpoints_window_start += 1
        if midpoints_window_start >= len(read_midpoints):
            break

        # update midpoints end
        while midpoints_window_end < len(read_midpoints) and read_midpoints[midpoints_window_end] < window_end:
            midpoints_window_end += 1

        # count reads
        window_count = midpoints_window_end - midpoints_window_start

        # compute p-value
        if window_count > 2:
            p_val = scan_stat_approx3(window_count, window_size, gene_end-gene_start, poisson_lambda)
            window_stats.append((window_count,p_val))
        else:
            window_stats.append((window_count,1))

    return window_stats


################################################################################
# get_gene_regions
#
# Return a hash of gene_id's mapping to lists consisting of (chromosome, start,
# end, strand). Coordinates are GTF format.
################################################################################
def get_gene_regions(ref_gtf):
    gene_regions = {}

    transcripts = read_genes(ref_gtf, key_id='transcript_id')

    for tid in transcripts:
        tx = transcripts[tid]
        gid = tx.kv['gene_id']

        if not gid in gene_regions:
            gene_regions[gid] = [tx.chrom, tx.exons[0].start, tx.exons[-1].end, tx.strand]
        else:
            gene_regions[gid][1] = min(gene_regions[gid][1], tx.exons[0].start)
            gene_regions[gid][2] = max(gene_regions[gid][2], tx.exons[-1].end)

    return gene_regions


################################################################################
# get_splice_junctions
#
# Return a list of indexes mapping the splice junctions of the given
# transcript Gene object.
################################################################################
def get_splice_junctions(tx):
    junctions = []
    if len(tx.exons) > 1:
        junctions.append(tx.exons[0].end)
        for i in range(1,len(tx.exons)-1):
            junctions.append(tx.exons[i].start)
            junctions.append(tx.exons[i].end)
        junctions.append(tx.exons[-1].start)
    return junctions


################################################################################
# isoform_gtf_filter
#
# Filter a gtf file for only the isoforms in the given set.
################################################################################
def isoform_gtf_filter(ref_gtf, isoform_set):
    iso_ref_fd, iso_ref_gtf = tempfile.mkstemp()
    iso_ref_open = os.fdopen(iso_ref_fd, 'w')

    for line in open(ref_gtf):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        if kv['transcript_id'] in isoform_set:
            print >> iso_ref_open, line,

    iso_ref_open.close()

    return iso_ref_gtf


################################################################################
# map_midpoints
#
# Map reads to their alignment midpoints
################################################################################
def map_midpoints(clip_in, chrom, gene_start, gene_end, gene_strand):
    read_midpoints = []

    if chrom in clip_in.references:
        # for each read in span
        for aligned_read in clip_in.fetch(chrom, gene_start, gene_end):
            ar_strand = '+'
            if aligned_read.is_reverse:
                ar_strand = '-'

            # check strand and quality
            if gene_strand == ar_strand and aligned_read.mapq > 0:

                # map read to midpoint
                read_midpoints.append(cigar_midpoint(aligned_read))

        # in case of differing read alignment lengths
        read_midpoints.sort()

    return read_midpoints


################################################################################
# prerna_gtf
#
# Add unspliced preRNAs to the gtf file, focus on exons, and remove
# redundancies.
################################################################################
def prerna_gtf(ref_gtf):
    unspliced_index = 0
    unspliced_hash = set()

    transcripts = read_genes(ref_gtf, key_id='transcript_id')

    pre_ref_fd, pre_ref_gtf = tempfile.mkstemp()
    pre_ref_open = os.fdopen(pre_ref_fd, 'w')

    # add unspliced single exon transcripts to hash
    for tid in transcripts:
        tx = transcripts[tid]
        if len(tx.exons) == 1:
            tx_key = (tx.chrom, tx.exons[0].start, tx.exons[0].end, tx.strand)
            unspliced_hash.add(tx_key)
        
    # process transcripts
    for tid in transcripts:
        tx = transcripts[tid]
        pre_start = tx.exons[0].start
        pre_end = tx.exons[-1].end
        pre_key = (tx.chrom, pre_start, pre_end, tx.strand)

        for i in range(len(tx.exons)):
            cols = (tx.chrom, 'clip_peaks', 'exon', str(tx.exons[i].start), str(tx.exons[i].end), '.', tx.strand, '.', gff.kv_gtf(tx.kv))
            print >> pre_ref_open, '\t'.join(cols)

        if not pre_key in unspliced_hash:
            unspliced_hash.add(pre_key)
            pre_kv = copy.copy(tx.kv)
            pre_kv['transcript_id'] = 'UNSPLICED%d' % unspliced_index
            unspliced_index += 1
            cols = (tx.chrom, 'clip_peaks', 'exon', str(pre_start), str(pre_end), '.', tx.strand, '.', gff.kv_gtf(pre_kv))
            print >> pre_ref_open, '\t'.join(cols)

    pre_ref_open.close()

    return pre_ref_gtf


################################################################################
# read_genes
#
# Parse a gtf file and return a set of Gene objects in a hash keyed by the
# id given.
################################################################################
def read_genes(gtf_file, key_id='transcript_id', sort=True):
    genes = {}
    for line in open(gtf_file):
        a = line.split('\t')

        kv = gff.gtf_kv(a[8])
        if not kv[key_id] in genes:
            genes[kv[key_id]] = Gene(a[0], a[6], kv)

        if a[2] == 'exon':
            genes[kv[key_id]].add_exon(int(a[3]), int(a[4]))

    return genes


################################################################################
# scan_stat_approx3
#
# Approximation 3.3 to the unconditional, poisson distributed scan statistic
# defined on p.28 of Glaz, Naus, Wallenstein book.
################################################################################
def scan_stat_approx3(k, w, T, lambd):
    L = float(T)/w
    psi = float(lambd)*w
    sigma = (k-1.0)*(L-1.0)*poisson.pmf(k, psi)
    p_val = 1.0 - math.exp(-sigma)
    return p_val


################################################################################
# span_gtf
#
# Add unspliced preRNAs to the gtf file, focus on exons, and remove
# redundancies.
################################################################################
def span_gtf(ref_gtf):
    # obtain gene regions
    gene_regions = get_gene_regions(ref_gtf)

    # print
    span_ref_fd, span_ref_gtf = tempfile.mkstemp()
    span_ref_open = os.fdopen(span_ref_fd, 'w')

    for gid in gene_regions:
        g = gene_regions[gid]
        cols = [g[0], 'clip_peaks', 'exon', str(g[1]), str(g[2]), '.', g[3], '.', kv_gtf({'gene_id':gid, 'transcript_id':gid})]
        print >> span_ref_open, '\t'.join(cols)

    span_ref_open.close()

    return span_ref_gtf


################################################################################
# windows2peaks
#
# Convert window counts and p-values to peak calls.
################################################################################
def windows2peaks(window_stats, sig_p, allowed_sig_gap, gene_start):
    peaks = []
    window_peak_start = None
    insig_gap = 0

    for i in range(len(window_stats)):
        c, p = window_stats[i]

        if p < sig_p:
            if window_peak_start == None:
                window_peak_start = i
            insig_gap = 0
        elif window_peak_start != None:
            insig_gap += 1
            if insig_gap > allowed_sig_gap:
                # save peak
                peaks.append((gene_start+window_peak_start, gene_start+i-insig_gap))

                # reset
                window_peak_start = None
                insig_gap = 0
            else:
                # let it ride
                pass

    if window_peak_start != None:
        peaks.append((gene_start+window_peak_start, gene_start+len(window_stats)-1-insig_gap))

    return peaks


################################################################################
# Gene class
################################################################################
class Gene:
    def __init__(self, chrom, strand, kv):
        self.chrom = chrom
        self.strand = strand
        self.kv = kv
        self.exons = []

    def add_exon(self, start, end):
        self.exons.append(Exon(start,end))
        if len(self.exons) > 1 and self.exons[-2].end > start:
            self.exons.sort()

    def __str__(self):
        return '%s %s %s %s' % (self.chrom, self.strand, kv_gtf(self.kv), ','.join([ex.__str__() for ex in self.exons]))

################################################################################
# Exon class
################################################################################
class Exon:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __cmp__(self, x):
        if self.start < x.start:
            return -1
        elif self.start > x.start:
            return 1
        else:
            return 0

    def __str__(self):
        return 'exon(%d-%d)' % (self.start,self.end)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
