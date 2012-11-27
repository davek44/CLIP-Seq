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
    usage = 'usage: %prog [options] <ref_gtf> <clip_bam>'
    parser = OptionParser(usage)

    parser.add_option('-c', dest='control_bam', help='Control BAM file')

    parser.add_option('-w', dest='window', type='int', default=50, help='Window size for scan statistic [Default: %default]')
    parser.add_option('-p', dest='p_val', type='float', default=.01, help='P-value required of window scan statistic tests [Default: %default]')

    parser.add_option('-t', dest='threads', type='int', default=2, help='Number of threads to use [Default: %default]')

    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error(usage)
    else:
        ref_gtf = args[0]
        clip_bam = args[1]


    ############################################
    # parameterize
    ############################################
    if options.control_bam:
        # make a new gtf w/ unspliced RNAs
        pre_ref_gtf = prerna_gtf(ref_gtf)

        # run Cufflinks on new gtf file and control BAM
        subprocess.call('cufflinks -p %d -G %s %s' % (options.threads, pre_ref_gtf, options.control_bam), shell=True)

        # choose most expressed isoform and save spliced/unspliced fpkms
        isoform_fpkms = control_max_fpkm(ref_gtf)

        # filter gtf for clip analysis
        clip_ref_gtf = isoform_gtf_filter(ref_gtf, isoform_fpkms)
    
    else:
        # make a new gtf file of only loci-spanning RNAs
        span_ref_gtf = span_gtf(ref_gtf)

        # run Cufflinks on new gtf file and CLIP BAM
        subprocess.call('cufflinks -p %d -G %s %s' % (options.threads, span_ref_gtf, clip_bam), shell=True)

        # not sure how to represent parameterization here...
        isoform_fpkms = {}

        # filter gtf for clip analysis
        clip_ref_gtf = span_ref_gtf


    ############################################
    # process genes
    ############################################
    # get gene regions
    gene_regions = get_gene_regions(clip_ref_gtf)

    # open clip-seq bam
    clip_in = pysam.Samfile(clip_bam, 'rb')

    # for each span
    for gene_id in gene_regions:
        (chrom, gene_start, gene_end, gene_strand) = gene_regions[gene_id]

        ############################################
        # map reads to midpoints
        ############################################
        read_midpoints = []

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

        ############################################
        # process windows
        ############################################
        poisson_lambda = float(len(read_midpoints)) / (gene_end - gene_start)
        midpoints_window_start = 0 # index of the first read_midpoint that fit in the window (except I'm allowing 0)
        midpoints_window_end = 0 # index of the first read_midpoint past the window

        for window_start in range(gene_start, gene_end-options.window_size+1):
            window_end = window_start + options.window_size

            # update midpoints start
            while midpoints_window_start < len(read_midpoints) and read_midpoints[midpoints_window_start] < window_start:
                midpoints_window_start += 1
            if midpoints_window_start >= len(read_midpoints):
                break

            # update midpoints end
            while midpoints_window_end < len(read_midpoints) and read_midpoints[midpoints_window_end] < window_end:
                midpoints_window_end += 1

            # compute statistic
            window_count = midpoints_window_end - midpoints_window_start
            if window_count > 2:
                p_val = scan_stat_approx3(window_count, options.window_size, gene_end-gene_start, poisson_lambda)
                
                if p_val < options.p_val:
                    #cols = (chrom, window_start, window_end, 'peak', -math.log(p_val))
                    cols = (chrom, window_start, window_end, window_count, p_val)
                    # (score is supposed to be 1-1000)
                    print '\t'.join([str(c) for c in cols])

        ############################################
        # post-process windows to peaks
        ############################################
        
        # ...

        ############################################
        # output peaks
        ############################################

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
    span_unsplied = {}
    transcripts = gff.read_genes(ref_gtf, key_id='transcript_id')
    for tid in transcripts:
        tx = transcripts[tid]

        span_start = tx.exons[0].start
        span_end = tx.exons[-1].end

        transcript_span[tid] = (tx.chrom, span_start, span_end, tx.strand)
        transcript_length[tid] = 0
        for exon in tx.exons:
            transcript_length[tid] += exon.end-exon.start+1

        if tid.startswith('UNSPLICED'):
            span_unspliced[transcript_span[tid]] = tid

    # read FPKMS
    transcript_fpkm = {}
    for line in open('isoforms.fpkm_tracking'):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        transcript_id = a[0]
        fpkm = float(a[9])

        transcript_fpkm[transcript_id] = fpkm

    # choose isoform by FPKM and length
    isoform_fpkms = {}
    g2t = gff.g2t(ref_gtf)
    for gene_id in g2t:
        # focus on spliced transcripts
        spliced_transcripts = [tid for tid in g2t[gene_id] if not tid.startswith('UNSPLICED')]

        # choose isoform
        max_tid = spliced_transcripts[0]
        for tid in spliced_transcripts[1:]:
            if transcript_fpkm[tid] > transcript_fpkm[max_tid]:
                max_tid = tid
            elif transcript_fpkm[tid] == transcript_fpkm[max_tid]:
                if transcript_length[tid] > transcript_length[max_tid]:
                    max_tid = tid

        # find unspliced FPKM
        unspliced_tid = span_unspliced[transcript_span[tid]]
        isoform_fpkms[max_tid] = (transcript_fpkm[max_tid], transcript_fpkm[unspliced_tid])

    return isoform_fpkms


################################################################################
# get_gene_regions
#
# Return a hash of gene_id's mapping to lists consisting of (chromosome, start,
# end, strand). Coordinates are BED format.
################################################################################
def get_gene_regions(ref_gtf):
    gene_regions = {}

    transcripts = gff.read_genes(ref_gtf, key_id='transcript_id')

    for tid in transcripts:
        tx = transcripts[tid]
        gid = tx.kv['gene_id']

        if not gid in gene_regions:
            gene_regions[gid] = [tx.chrom, tx.exons[0].start-1, tx.exons[-1].end, tx.strand]
        else:
            gene_regions[gid][1] = min(gene_regions[gid][1], tx.exons[0].start-1)
            gene_regions[gid][2] = max(gene_regions[gid][2], tx.exons[-1].end)

    return gene_regions


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
# prerna_gtf
#
# Add unspliced preRNAs to the gtf file, focus on exons, and remove
# redundancies.
################################################################################
def prerna_gtf(ref_gtf):
    unspliced_index = 0
    unspliced_hash = set()

    transcripts = gff.read_genes(ref_gtf, key_id='transcript_id')

    pre_ref_fd, pre_ref_gtf = tempfile.mkstemp()
    pre_ref_open = os.fdopen(pre_ref_fd, 'w')

    for tid in transcripts:
        tx = transcripts[tid]
        pre_start = tx.exons[0].start
        pre_end = tx.exons[-1].end
        pre_key = (tx.chrom, pre_start, pre_end, tx.strand)

        for i in range(len(tx.exons)):
            cols = (tx.chrom, 'clip_peaks', 'exon', str(tx.exons[i].start), str(tx.exons[i].end), '.', tx.strand, '.', gff.kv_gtf(tx.kv))
            unspliced_hash.add((tx.chrom, tx.exons[i].start, tx.exons[i].end, tx.strand))
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
# __main__
################################################################################
if __name__ == '__main__':
    main()
