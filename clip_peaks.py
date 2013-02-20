#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import poisson
from bisect import bisect_left, bisect_right
import copy, math, os, pdb, subprocess, sys, tempfile
import pybedtools, pysam
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

    #parser.add_option('--min_control_fpkm', dest='min_control_fpkm', type='float', default=0.25, help='Minimum FPKM to allow the transcripts covering a window to take [Default: %default]')

    parser.add_option('-o', dest='out_dir', default='peaks', help='Output directory [Default: %default]')

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

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    ############################################
    # parameterize
    ############################################
    if options.control_bam:

        # make a new gtf w/ unspliced RNAs
        pre_ref_gtf = prerna_gtf(ref_gtf, options.out_dir)

        # store transcripts
        transcripts = read_genes(pre_ref_gtf, key_id='transcript_id')
        g2t = gff.g2t(pre_ref_gtf)

        # set junctions
        set_transcript_junctions(transcripts)

        # run Cufflinks on new gtf file and control BAM
        if not options.cuff_done:
            subprocess.call('cufflinks -o %s -p %d -G %s %s' % (options.out_dir, options.threads, pre_ref_gtf, options.control_bam), shell=True)

        # set FPKMs
        set_transcript_fpkms(transcripts, options.out_dir)

        # set exon and preRNA FPKMs and filter for most expressed isoform
        #set_fpkm_control(transcripts, pre_ref_gtf)
    
    else:
        # make a new gtf file of only loci-spanning RNAs
        span_ref_gtf = span_gtf(ref_gtf, options.out_dir)

        # run Cufflinks on new gtf file and CLIP BAM
        if not options.cuff_done:
            subprocess.call('cufflinks -o %s -p %d -G %s %s' % (options.out_dir, options.threads, span_ref_gtf, clip_bam), shell=True)

        # store span transcripts
        transcripts = read_genes(span_ref_gtf, key_id='transcript_id')

        # set "exon" FPKMs
        set_fpkm_span(transcripts)

    # count transcriptome CLIP reads
    total_reads = int(subprocess.check_output('intersectBed -bed -u -s -abam %s -b %s/transcripts.gtf | cut -f4 | sort -u | wc -l' % (clip_bam, options.out_dir), shell=True))

    # compute # of tests we will perform
    txome_size = transcriptome_size(transcripts, options.window_size)


    ############################################
    # process genes
    ############################################
    # open clip-seq bam
    clip_in = pysam.Samfile(clip_bam, 'rb')
    
    # open peak output gff
    peaks_out = open('%s/peaks.gff' % options.out_dir, 'w')
    peak_id = 1

    # for each gene
    for gene_id in g2t:

        # make a more focused transcript hash for this gene
        gene_transcripts = {}
        for tid in g2t[gene_id]:
            gene_transcripts[tid] = transcripts[tid]

        # obtain basic gene attributes
        (gchrom, gstrand, gstart, gend) = gene_attrs(gene_transcripts)

        # map reads to midpoints
        read_midpoints = map_midpoints(clip_in, gchrom, gstart, gend, gstrand)

        # find splice junctions
        #junctions = map_splice_junctions(tx)

        # count reads and compute p-values in windows
        window_stats = count_windows(clip_in, options.window_size, read_midpoints, gene_transcripts, gstart, gend, total_reads, txome_size)

        # post-process windows to peaks
        peaks = windows2peaks(read_midpoints, gene_transcripts, gene_start, window_stats, options.window_size, options.p_val, total_reads, txome_size)

        # output peaks
        for pstart, pend, pcount, ppval in peaks:
            cols = [gchrom, 'clip_peaks', 'peak', str(pstart), str(pend), '.', gstrand, '.', 'id "PEAK%d"; gene_id "%s"; count "%d"; p "%.2e"' % (peak_id,gene_id,pcount,ppval)]
            print >> peaks_out, '\t'.join(cols)
            peak_id += 1

    clip_in.close()
    peaks_out.close()


################################################################################
# cigar_midpoint
# 
# Input
#  aligned_read: pysam AlignedRead object.
#
# Output
#  midpoint:     Midpoint of the alignment, considering the insertions and
#                 deletions in its CIGAR string (which includes splicing).
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
#
# Input
#  window_start:     Window start coordinate.
#  window_end:       Window end coordinate.
#  gene_transcripts: Hash mapping transcript_id to isoform Gene objects,
#                     containing only keys for a specific gene.
#  junctions_i:      Hash mapping transcript_id to the index of the first
#                     junction ahead of the window start in that transcript's
#                     junctions array.
#  total_reads:      Total number of reads aligned to the transcriptome.
# 
# Output
#  lambda:           Poisson lambda for this window.
################################################################################
def convolute_lambda(window_start, window_end, gene_transcripts, junctions_i, total_reads):
#def convolute_lambda(window_start, window_end, fpkm_exon, fpkm_pre, total_reads, junctions, ji):

    # initialize FPKM
    fpkm_conv = 0

    # for each isoform
    for tid in gene_transcripts:
        # shortcuts
        tx = gene_transcripts[tid]
        ji = junctions_i[tid]

        # determine transcript coefficient
        tcoef = 0.0

        # after junctions
        if ji >= len(tx.junctions):
            tcoef = 1.0

        # next junction out of window
        elif window_end < tx.junctions[ji]:
            if ji % 2 == 0: # in an exon
                tcoef = 1.0

        # junctions
        else:
            # window start to first junction
            if ji % 2 == 0: # exon
                tcoef = tx.junctions[ji] - window_start

            # advance
            ji += 1

            # between junctions
            while ji < len(tx.junctions) and tx.junctions[ji] <= window_end:
                if ji % 2 == 0: # exon
                    tcoef += tx.junctions[ji] - tx.junctions[ji-1]
                ji += 1

            # back up
            ji -= 1

            # last junction to window end
            if ji % 2 == 1: # exon
                tcoef += window_end - tx.junctions[ji] + 1

            # normalize
            tcoef /= float(window_end-window_start+1)

        # add to fpkm
        fpkm_conv += tcoef * tx.fpkm

    # bump to min fpkm
    fpkm_conv = max(fpkm_conv, 0.1)

    # convert from fpkm to lambda
    return fpkm_conv / 1000.0*(total_reads/1000000.0)


################################################################################
# count_windows
#
# Count the number of reads and compute the scan statistic p-value in each
# window through the gene.
#
# Input
#  clip_in:          Open pysam BAM file for clip-seq alignments.
#  window_size:      Scan statistic window size.
#  read_midpoints:   Sorted list of read alignment midpoints.
#  gene_transcripts: Hash mapping transcript_id to isoform Gene objects,
#                     containing only keys for a specific gene.
#  gene_start:       Start of the gene's span.
#  gene_end:         End of the gene's span.
#  total_reads:      Total number of reads aligned to the transcriptome.
#  txome_size:       Total number of bp in the transcriptome.
#
# Output
#  window_stats:     List of tuples (alignment count, p value) for all windows.
################################################################################
def count_windows(clip_in, window_size, read_midpoints, gene_transcripts, gene_start, gene_end, total_reads, txome_size):
    # set lambda using whole region (some day, compare this to the cufflinks estimate)
    # poisson_lambda = float(len(read_midpoints)) / (gene_end - gene_start)

    midpoints_window_start = 0 # index of the first read_midpoint that fit in the window (except I'm allowing 0)
    midpoints_window_end = 0 # index of the first read_midpoint past the window

    # initialize index of the first junction ahead of the window start for each transcript
    junctions_i = {}
    for tid in gene_transcripts:
        junctions_i[tid] = 0

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

        # update junctions indexes (<= comparison because junctions holds the 1st bp of next exon/intron)
        for tid in gene_transcripts:
            tjunctions = gene_transcripts[tid].junctions
            while junctions_i[tid] < len(tjunctions) and tjunctions[junctions_i[tid]] <= window_start:
                junctions_i[tid] += 1

        # set lambda
        window_lambda = convolute_lambda(window_start, window_end, gene_transcripts, junctions_i, total_reads)

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

        # for debugging
        # print tx.exons[0].start+len(window_stats)-1, window_count, window_stats[-1][1], window_lambda

    return window_stats


################################################################################
# gene_attrs
#
# Input
#  gene_transcripts: Hash mapping transcript_id to isoform Gene objects,
#                      containing only keys for a specific gene.
#
# Output
#  gene_chrom:       Self explanatory...
#  gene_strand:
#  gene_start:
#  gene_end:
################################################################################
def gene_attrs(gene_transcripts):
    gene_start = None
    gene_end = None
    for tx in gene_transcripts.values():
        gene_chrom = tx.chrom
        gene_strand = tx.strand
        if gene_start == None:
            gene_start = tx.exons[0].start
        else:
            gene_start = min(gene_start, tx.exons[0].start)
        if gene_end == None:
            gene_end = tx.exons[-1].end
        else:
            gene_end = max(gene_end, tx.exons[-1].end)

    return gene_chrom, gene_strand, gene_start, gene_end


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
# map_midpoints
#
# Map reads from a gene of interest to their alignment midpoints, filtering for
# strand and quality.
#
# Input
#  clip_in:        Open pysam BAM file for clip-seq alignments.
#  gene_chrom:     Chromosome of the gene of interest.
#  gene_start:     Start coordinate of the gene of interest.
#  gene_end:       End coordinate of the gene of interest.
#  gene_strand:    Strand of the gene of interest.
#
# Output
#  read_midpoints: Sorted list of read alignment midpoints.
################################################################################
def map_midpoints(clip_in, gene_chrom, gene_start, gene_end, gene_strand):
    read_midpoints = []

    if gene_chrom in clip_in.references:
        # for each read in span
        for aligned_read in clip_in.fetch(gene_chrom, gene_start, gene_end):
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
# merge_peaks_count
#
# Merge the given trimmed windows with counts using bedtools and then recount
# the reads in the new intervals.
#
# Input
#  trimmed_windows: List of (start,end) tuples for significant windows, trimmed
#                    to be tight around read midpoints.
#  read_midpoints:  Sorted list of read alignment midpoints.
#
# Output
#  peaks:           List of (start,end,count) tuples for significant trimmed
#                    and merged windows
################################################################################
def merge_peaks_count(trimmed_windows, read_midpoints):
    # create BED intervals
    bed_intervals = []
    for wstart, wend in trimmed_windows:
        bed_a = ['chrFAKE', str(wstart-1), str(wend)]
        bed_intervals.append(pybedtools.create_interval_from_list(bed_a))
    bedtool = pybedtools.BedTool(bed_intervals)

    # merge BED intervals
    bedtool_merge = bedtool.merge(stream=True)

    # recount peaks
    peaks = []
    for bed_interval in bedtool_merge.features():
        pstart = bed_interval.start+1
        pend = bed_interval.end

        reads_start_i = bisect_left(read_midpoints, pstart)
        reads_end_i = bisect_right(read_midpoints, pend)
        read_count = reads_end_i - reads_start_i

        peaks.append((pstart, pend, read_count))

    return peaks


################################################################################
# merge_windows
#
# Merge adjacent significant windows and save index tuples.
#
# Input
#  window_stats:    Counts and p-values for each window.
#  window_size:     Scan statistic window size.
#  sig_p:           P-value at which to call window counts significant.
#  gene_start:      Start of the gene's span.
#  allowed_sig_gap: Max gap size between significant windows to perform a merge.
#
# Output
#  merged_windows:  List of (start,end) tuples for merged significant windows.
################################################################################
def merge_windows(window_stats, window_size, sig_p, gene_start, allowed_sig_gap = 1):
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
                merged_windows.append((gene_start+window_peak_start, gene_start+i-insig_gap+window_size-1))

                # reset
                window_peak_start = None
                insig_gap = 0
            else:
                # let it ride
                pass

    if window_peak_start != None:
        merged_windows.append((gene_start+window_peak_start, gene_start+len(window_stats)-1-insig_gap+window_size-1))

    return merged_windows


################################################################################
# peak_stats
#
# Compute a new p-value for the final peak.
#
# Input
#  windows_counts:   List of (start,end,count) tuples for significant trimmed
#                     and merged windows.
#  gene_transcripts: Hash mapping transcript_id to isoform Gene objects,
#                     containing only keys for a specific gene.
#  total_reads:      Total number of reads aligned to the transcriptome.
#  txome_size:       Total number of bp in the transcriptome.
#
# Output
#  peaks:            List of (start,end,count,p-val) tuples for peaks.
################################################################################
def peak_stats(windows_counts, gene_transcripts, total_reads, txome_size):
    peaks = []
    for wstart, wend, wcount in windows_counts:
        # set index of the first junction ahead of the window start for each transcript
        junctions_i = {}
        for tid in gene_transcripts:
            junctions_i[tid] = bisect_left(gene_transcripts[tid].junctions, wstart)

        peak_lambda = convolute_lambda(wstart, wend, gene_transcripts, junctions_i, total_reads)
        p_val = scan_stat_approx3(wcount, wend-wstart+1, txome_size, peak_lambda)
        peaks.append((wstart,wend,wcount,p_val))
    return peaks


################################################################################
# prerna_gtf
#
# Add unspliced preRNAs to the gtf file, focus on exons, and remove
# redundancies.
#
# Input
#  ref_gtf:     Reference GTF file to expand with unspliced preRNAs.
#  out_dir:     Directory in which to output the expanded GTF file.
#
# Output
#  prerna.gtf:  Expanded GTF file.
#  pre_ref_gtf: The filename of the expanded GTF file.
################################################################################
def prerna_gtf(ref_gtf, out_dir):
    unspliced_index = 0
    unspliced_hash = set()

    transcripts = read_genes(ref_gtf, key_id='transcript_id')

    pre_ref_gtf = '%s/prerna.gtf' % out_dir
    pre_ref_open = open(pre_ref_gtf, 'w')

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
#
# Input
#  gtf_file: GTF file from which to read genes.
#  key_id:   GTF key with which to hash the genes (transcript_id or gene_id)
#
# Output
#  genes:    Hash mapping key_id to Gene objects.
################################################################################
def read_genes(gtf_file, key_id='transcript_id'):
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
#
# k is the # of reads
# w is the window size
# T is the transcriptome size
# lambd is the reads/nt
#
# It seems to breakdown on the very low end, so I'm adding a check that
# that k is greater than psi.
################################################################################
def scan_stat_approx3(k, w, T, lambd):
    L = float(T)/w
    psi = float(lambd)*w
    if k < psi:
        p_val = 1.0
    else:
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
        if len(ref_transcripts[max_tid].exons) == 1:
            # irrelevant
            ref_transcripts[max_tid].fpkm_pre = 0
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
def set_fpkm_span(ref_transcripts):
    # read FPKMS
    fpkm_tracking_in = open('isoforms.fpkm_tracking')
    line = fpkm_tracking_in.readline()
    for line in fpkm_tracking_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        transcript_id = a[0]
        fpkm = float(a[9])

        ref_transcripts[transcript_id].fpkm_exon = fpkm
        ref_transcripts[transcript_id].fpkm_pre = 0
    fpkm_tracking_in.close()

    # remove unset transcripts
    for tid in ref_transcripts.keys():
        if ref_transcripts[tid].fpkm_exon == None:
            # this can happen if two transcripts are the same except for differing strands
            # cufflinks will ignore one of them.
            print >> sys.stderr, 'WARNING: Missing FPKM for gene span %s' % tid

            del ref_transcripts[tid]


################################################################################
# set_transcript_fpkms
#
# Input
#  transcripts: Hash mapping transcript_id to isoform Gene objects.
#  out_dir:     Directory in which to output the gene span GTF file.
# 
# Output
#  transcripts: Same hash, FPKM attribute set.
################################################################################
def set_transcript_fpkms(transcripts, out_dir):
    fpkm_in = open('%s/isoforms.fpkm_tracking' % out_dir)
    line = fpkm_in.readline()
    for line in fpkm_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        tid = a[0]
        fpkm = float(a[9])

        transcripts[tid].fpkm = fpkm


################################################################################
# set_transcript_junctions
#
# Input
#  transcripts: Hash mapping transcript_id to isoform Gene objects.
# 
# Output
#  transcripts:  Same hash, junctions attribute set.
#  tx.junctions: Sorted list of junction coordinates. For each junction, save
#                 the first bp of the next exon/intron.
################################################################################
def set_transcript_junctions(transcripts):
    for tid in transcripts:
        tx = transcripts[tid]
        if len(tx.exons) > 1:
            tx.junctions.append(tx.exons[0].end+1)
            for i in range(1,len(tx.exons)-1):
                tx.junctions.append(tx.exons[i].start)
                tx.junctions.append(tx.exons[i].end+1)
            tx.junctions.append(tx.exons[-1].start)


################################################################################
# span_gtf
#
# Input
#  ref_gtf:      Reference GTF file to convert to gene spans.
#  out_dir:      Directory in which to output the gene span GTF file.
#
# Output
#  span.gtf:     Gene span GTF file.
#  span_ref_gtf: The filename of the gene span GTF file.
################################################################################
def span_gtf(ref_gtf, out_dir):
    # obtain gene regions
    gene_regions = get_gene_regions(ref_gtf)

    # print
    span_ref_gtf = '%s/span.gtf' % out_dir
    span_ref_open = open(span_ref_gtf, 'w')

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
# trim_windows
#
# Input
#  windows:        List of (start,end) tuples for merged significant windows.
#  read_midpoints: Sorted list of read alignment midpoints.
#
# Output
#  trimmed_windows: List of (start,end) tuples for significant windows, trimmed
#                    to be tight around read midpoints.
################################################################################
def trim_windows(windows, read_midpoints):
    trimmed_windows = []
    for wstart, wend in windows:
        trim_start_i = bisect_left(read_midpoints, wstart)
        trim_end_i = bisect_right(read_midpoints, wend)
        trimmed_windows.append((int(read_midpoints[trim_start_i]), int(read_midpoints[trim_end_i-1]+0.5)))
    return trimmed_windows


################################################################################
# windows2peaks
#
# Convert window counts and p-values to peak calls.
#
# Input
#  read_midpoints:   Sorted list of read alignment midpoints.
#  gene_transcripts: Hash mapping transcript_id to isoform Gene objects,
#                     containing only keys for a specific gene.
#  gene_start:       Start of the gene's span.
#  window_stats:     Counts and p-values for each window.
#  window_size:      Scan statistic window size.
#  sig_p:            P-value at which to call window counts significant.
#  total_reads:      Total number of reads aligned to the transcriptome.
#  txome_size:       Total number of bp in the transcriptome.
#
# Output
#  peaks:            List of (start,end,count,p-val) tuples for peaks.
################################################################################
def windows2peaks(read_midpoints, gene_transcripts, gene_start, window_stats, window_size, sig_p, total_reads, txome_size):
    merged_windows = merge_windows(window_stats, window_size, sig_p, gene_start)
    trimmed_windows = trim_windows(merged_windows, read_midpoints)
    statless_peaks = merge_peaks_count(trimmed_windows, read_midpoints)
    peaks = peak_stats(statless_peaks, gene_transcripts, total_reads, txome_size)
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

        self.fpkm = None
        self.junctions = []

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
    #pdb.runcall(main)
