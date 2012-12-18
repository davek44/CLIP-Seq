#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import poisson
from bisect import bisect_left, bisect_right
import copy, math, os, subprocess, sys, tempfile
import pysam
import gff

################################################################################
# clip_peaks.py
#
# Call peaks in CLIP-Seq data.
#
# Notes on conventions:
# 1. All indexes are GFF-based. I.e. the first bp in a sequence is 1 and the
#    last is len(sequence). For annotations, the start marks the first bp of
#    the annotation and the end marks the last. The length is end-start+1.
#    
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

        # store transcripts
        transcripts = read_genes(ref_gtf, key_id='transcript_id')

        # set exon and preRNA FPKMs and filter for most expressed isoform
        set_fpkm_control(transcripts, pre_ref_gtf)
    
    else:
        # make a new gtf file of only loci-spanning RNAs
        span_ref_gtf = span_gtf(ref_gtf)

        # run Cufflinks on new gtf file and CLIP BAM
        if not options.cuff_done:
            subprocess.call('cufflinks -p %d -G %s %s' % (options.threads, span_ref_gtf, clip_bam), shell=True)

        # store span transcripts
        transcripts = read_genes(span_ref_gtf, key_id='transcript_id')

        # set "exon" FPKMs
        set_fpkm_span(transcripts)

    # compute # of tests we will perform
    txome_size = transcriptome_size(transcripts, options.window_size)


    ############################################
    # process genes
    ############################################
    # open clip-seq bam
    clip_in = pysam.Samfile(clip_bam, 'rb')

    # for each span
    for tid in transcripts:
        tx = transcripts[tid]

        # map reads to midpoints
        read_midpoints = map_midpoints(clip_in, tx.chrom, tx.exons[0].start, tx.exons[-1].end, tx.strand)

        # find splice junctions
        junctions = map_splice_junctions(tx)

        # count reads and compute p-values in windows
        window_stats = count_windows(clip_in, options.window_size, tx, read_midpoints, junctions, txome_size)

        # post-process windows to peaks
        peaks = windows2peaks(read_midpoints, junctions, window_stats, options.p_val, tx.exons[0].start, txome_size)

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
# convolute_lambda
#
# Determine the convoluted poisson lambda for the given window using the
# transcript's FPKM estimates.
#
# Recall that junctions contains the 1st bp of the next exon/intron.
################################################################################
def convolute_lambda(window_start, window_end, fpkm_exon, fpkm_pre, junctions, ji):
    # after junctions
    if ji >= len(junctions):
        fpkm_conv = fpkm_exon+fpkm_pre

    # next junction out of window
    elif window_end < junctions[ji]:
        if ji % 2 == 0:
            fpkm_conv = fpkm_exon+fpkm_pre
        else:
            fpkm_conv = fpkm_pre

    # junctions
    else:
        # window start to first junction
        if ji % 2 == 0: # exon
            fpkm_conv = (junctions[ji]-window_start)*(fpkm_exon+fpkm_pre)        
        else: # intron
            fpkm_conv = (junctions[ji]-window_start)*fpkm_pre

        # advance
        ji += 1

        # between junctions
        while ji < len(junctions) and junctions[ji] <= window_end:
            if ji % 2 == 0: # exon
                fpkm_conv += (junctions[ji]-junctions[ji-1])*(fpkm_exon+fpkm_pre)
            else: # intron
                fpkm_conv += (junctions[ji]-junctions[ji-1])*fpkm_pre

            ji += 1

        # back up
        ji -= 1

        # last junction to window end
        if ji % 2 == 0: # intron
            fpkm_conv += (window_end-junctions[ji]+1)*fpkm_pre
        else: # exon
            fpkm_conv += (window_end-junctions[ji]+1)*(fpkm_exon+fpkm_pre)

        # normalize
        fpkm_conv /= float(window_end-window_start+1)

    return fpkm_conv


################################################################################
# count_windows
#
# Count the number of reads and compute the scan statistic p-value in each
# window through the gene.
################################################################################
def count_windows(clip_in, window_size, tx, read_midpoints, junctions, txome_size):
    gene_start = tx.exons[0].start
    gene_end = tx.exons[-1].end

    # set lambda using whole region (some day, compare this to the cufflinks estimate)
    # poisson_lambda = float(len(read_midpoints)) / (gene_end - gene_start)

    midpoints_window_start = 0 # index of the first read_midpoint that fit in the window (except I'm allowing 0)
    midpoints_window_end = 0 # index of the first read_midpoint past the window

    junctions_i = 0 # index of the first junction ahead of the window start

    # to avoid redundant computation
    precomputed_pvals = {}

    window_stats = []

    for window_start in range(gene_start, gene_end-window_size+1):
        window_end = window_start + window_size - 1

        # update midpoints start
        while midpoints_window_start < len(read_midpoints) and read_midpoints[midpoints_window_start] < window_start:
            midpoints_window_start += 1
        if midpoints_window_start >= len(read_midpoints):
            break

        # update midpoints end
        while midpoints_window_end < len(read_midpoints) and read_midpoints[midpoints_window_end] <= window_end:
            midpoints_window_end += 1

        # count reads
        window_count = midpoints_window_end - midpoints_window_start

        # update junctions index (<= comparison because junctions holds the 1st bp of next exon/intron)
        while junctions_i < len(junctions) and junctions[junctions_i] <= window_start:
            junctions_i += 1

        # set lambda
        window_lambda = convolute_lambda(window_start, window_end, tx.fpkm_exon, tx.fpkm_pre, junctions, junctions_i)

        # compute p-value
        if window_count > 2:
            if (window_count,window_lambda) in precomputed_pvals:
                p_val = precomputed_pvals[(window_count,window_lambda)]
            else:
                p_val = scan_stat_approx3(window_count, window_size, txome_size, window_lambda)
                precomputed_pvals[(window_count,window_lambda)] = p_val
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
# map_splice_junctions
#
# Return a list of indexes mapping the splice junctions of the given
# transcript Gene object.
#
# For each junction, save the first bp of the next exon/intron.
################################################################################
def map_splice_junctions(tx):
    junctions = []
    if len(tx.exons) > 1:
        junctions.append(tx.exons[0].end+1)
        for i in range(1,len(tx.exons)-1):
            junctions.append(tx.exons[i].start)
            junctions.append(tx.exons[i].end+1)
        junctions.append(tx.exons[-1].start)
    return junctions


################################################################################
# map_midpoints
#
# Map reads to their alignment midpoints, filtering for strand and quality.
# Return a sorted list of midpoints.
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
# merge_windows
#
# Merge adjacent significant windows and save index tuples.
################################################################################
def merge_windows(window_stats, sig_p, gene_start, allowed_sig_gap = 1):
    merged_windows = []
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
                # save window
                merged_windows.append((gene_start+window_peak_start, gene_start+i-insig_gap))

                # reset
                window_peak_start = None
                insig_gap = 0
            else:
                # let it ride
                pass

    if window_peak_start != None:
        merged_windows.append((gene_start+window_peak_start, gene_start+len(window_stats)-1-insig_gap))

    return merged_windows


################################################################################
# peak_stats
#
# Compute a new p-value for the final peak.
################################################################################
def peak_stats(windows_counts, junctions, txome_size):
    peaks = []
    for wstart, wend, wcount in windows:
        junctions_i = bisect_left(junctions, wstart)
        peak_lambda = convolute_lambda(wstart, wend, fpkm_exon, fpkm_pre, junctions, junctions_i)
        p_val = scan_stat_approx3(wcount, wend-wstart+1, txome_size, peak_lambda)
        peaks.append((wstart,wend,wcount,p_val))
    return peaks


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
# set_fpkm_control
#
# For each gene:
# 1. Choose the most expressed isoform.
# 2. Set it's exonic and preRNA FPKMs.
# 3. Filter the unchosen isoform out of 'ref_transcripts'
################################################################################
def set_fpkm_control(ref_transcripts, add_gtf):
    # collect transcript spans to map spliced isoforms to their pre-RNA
    transcript_span = {}
    span_unspliced = {}
    add_transcripts = read_genes(add_gtf, key_id='transcript_id')
    for tid in add_transcripts:
        tx = add_transcripts[tid]

        span_start = tx.exons[0].start
        span_end = tx.exons[-1].end

        transcript_span[tid] = (tx.chrom, span_start, span_end, tx.strand)

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

    # choose isoform by FPKM
    g2t = gff.g2t(add_gtf)
    for gene_id in g2t:
        # focus on processed transcripts
        processed_transcripts = [tid for tid in g2t[gene_id] if not tid.startswith('UNSPLICED')]

        # verify abundance estimates
        for tid in processed_transcripts:
            if not tid in transcript_fpkm:
                # this can happen if two transcripts are the same in the GTF file.
                # cufflinks will ignore one of them.
                print >> sys.stderr, 'WARNING: Missing FPKM for spliced transcript %s' % tid
                transcript_fpkm[tid] = 0

        # choose isoform
        max_tid = processed_transcripts[0]
        for tid in processed_transcripts[1:]:
            if transcript_fpkm[tid] > transcript_fpkm[max_tid]:
                max_tid = tid

        # set exonic transcript FPKM
        ref_transcripts[max_tid].fpkm_exon = transcript_fpkm[max_tid]

        # set preRNA FPKM
        if len(transcripts[max_tid].exons) == 1:
            pass # leave unset
        else:
            # find unspliced
            unspliced_tid = span_unspliced[transcript_span[max_tid]]
            if unspliced_tid not in transcript_fpkm:
                # this can happen if two transcripts are the same except for differing strands
                # cufflinks will ignore one of them.
                print >> sys.stderr, 'WARNING: Missing FPKM for unspliced transcript %s' % unspliced_tid
                ref_transcripts[max_tid].fpkm_pre = transcript_fpkm[max_tid] / 2.0
            else:
                ref_transcripts[max_tid].fpkm_pre = transcript_fpkm[unspliced_tid]

    # remove unset transcripts
    for tid in ref_transcripts.keys():
        if ref_transcripts[tid].fpkm_exon == None:
            del ref_transcripts[tid]


################################################################################
# set_fpkm_span
#
# Set the "exonic" FPKMs for each gene span.
################################################################################
def set_fpkm_control(ref_transcripts):
    # read FPKMS
    fpkm_tracking_in = open('isoforms.fpkm_tracking')
    line = fpkm_tracking_in.readline()
    for line in fpkm_tracking_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        transcript_id = a[0]
        fpkm = float(a[9])

        ref_transcripts[transcript_id].fpkm_exon = fpkm
    fpkm_tracking_in.close()

    # remove unset transcripts
    for tid in ref_transcripts.keys():
        if ref_transcripts[tid].fpkm_exon == None:
            # this can happen if two transcripts are the same except for differing strands
            # cufflinks will ignore one of them.
            print >> sys.stderr, 'WARNING: Missing FPKM for gene span %s' % tid

            del ref_transcripts[tid]


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
# transcriptome_size
#
# Compute the number of window tests we will perform by considering the size of
# the transcriptome with window_size's subtracted.
################################################################################
def transcriptome_size(transcripts, window_size):
    txome_size = 0
    for tid in transcripts:
        tx = transcripts[tid]
        txome_size += tx.exons[-1].end - tx.exons[0].start - window_size + 1
    return txome_size


################################################################################
# trim_windows_count
#
# Trim each window to be tight around read midpoints.
################################################################################
def trim_windows_count(windows, read_midpoints):
    trimmed_windows = []
    for wstart, wend in windows:
        trim_start_i = read_midpoints[bisect_left(read_midpoints, wstart)]
        trim_end_i = read_midpoints[bisect_right(read_midpoints, wend)]
        read_count = trim_end_i - trim_start_i + 1
        trimmed_windows.append((read_midpoints[trim_start_i], read_midpoints[trim_end_i], read_count))
    return trimmed_windows


################################################################################
# windows2peaks
#
# Convert window counts and p-values to peak calls.
################################################################################
def windows2peaks(read_midpoints, junctions, window_stats, sig_p, gene_start, txome_size):
    merged_windows = merge_windows(window_stats, sig_p, gene_start)
    trimmed_windows_counts = trim_windows_count(merged_windows, read_midpoints, gene_start)
    peaks = peak_stats(trimmed_windows_counts, junctions, txome_size)
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

        self.fpkm_exon = None
        self.fpkm_pre = None

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
