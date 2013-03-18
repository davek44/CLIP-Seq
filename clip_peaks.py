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
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <clip_bam> <ref_gtf>'
    parser = OptionParser(usage)

    # IO options
    parser.add_option('-c', dest='control_bam', help='Control BAM file')
    parser.add_option('-o', dest='out_dir', default='peaks', help='Output directory [Default: %default]')

    # peak calling options
    parser.add_option('-w', dest='window_size', type='int', default=50, help='Window size for scan statistic [Default: %default]')
    parser.add_option('-p', dest='p_val', type='float', default=.01, help='P-value required of window scan statistic tests [Default: %default]')

    # cufflinks options
    parser.add_option('--cuff_done', dest='cuff_done', action='store_true', default=False, help='A cufflinks run to estimate the model parameters is already done [Default: %default]')
    parser.add_option('-t', dest='threads', type='int', default=2, help='Number of threads to use [Default: %default]')

    # debug options
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='Verbose output [Default: %default]')
    parser.add_option('-g', '--gene', dest='gene_only', help='Call peaks on the specified gene only')
    parser.add_option('--print_windows', dest='print_windows', default=False, action='store_true', help='Print statistics for all windows [Default: %default]')

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
    if options.verbose:
        print >> sys.stderr, 'Estimating gene abundances...'

    if options.control_bam:
        # make a new gtf w/ unspliced RNAs
        update_ref_gtf = prerna_gtf(ref_gtf, options.out_dir)

        # run Cufflinks on new gtf file and control BAM
        if not options.cuff_done:
            subprocess.call('cufflinks -o %s -p %d -G %s %s' % (options.out_dir, options.threads, update_ref_gtf, options.control_bam), shell=True)
    
    else:
        # make a new gtf file of only loci-spanning RNAs
        update_ref_gtf = span_gtf(ref_gtf, options.out_dir)

        # run Cufflinks on new gtf file and CLIP BAM
        if not options.cuff_done:
            subprocess.call('cufflinks -o %s -p %d -G %s %s' % (options.out_dir, options.threads, update_ref_gtf, clip_bam), shell=True)

    # store transcripts
    transcripts = read_genes(update_ref_gtf, key_id='transcript_id')
    g2t = gff.g2t(update_ref_gtf)

    # set junctions
    set_transcript_junctions(transcripts)

    # set "exon" FPKMs
    set_transcript_fpkms(transcripts, options.out_dir, options.verbose)

    if options.verbose:
        print >> sys.stderr, 'Computing global statistics...'

    # count transcriptome CLIP reads (overestimates small RNA single ended reads by counting antisense)
    subprocess.call('intersectBed -abam %s -b %s/transcripts.gtf > %s/transcripts.bam' % (clip_bam, options.out_dir, options.out_dir), shell=True)
    total_reads = count_reads('%s/transcripts.bam' % options.out_dir)

    # compute # of tests we will perform
    txome_size = transcriptome_size(transcripts, options.window_size)


    ############################################
    # process genes
    ############################################
    # TODO: Can I convert to using transcripts.bam here? Does it affect performance given an indexing?
    # index
    subprocess.call('samtools index %s' % clip_bam, shell=True)

    # open clip-seq bam
    clip_in = pysam.Samfile(clip_bam, 'rb')
    
    # open peak output gff
    peaks_out = open('%s/peaks.gff' % options.out_dir, 'w')
    peak_id = 1

    # open window output
    windows_out = None
    if options.print_windows:
        windows_out = open('%s/window_stats.txt' % options.out_dir, 'w')

    # for each gene
    if options.gene_only:
        gene_ids = [options.gene_only]
    else:
        gene_ids = g2t.keys()

    for gene_id in gene_ids:
        if options.verbose:
            print >> sys.stderr, 'Processing %s...' % gene_id

        # make a more focused transcript hash for this gene
        gene_transcripts = {}
        for tid in g2t[gene_id]:
            gene_transcripts[tid] = transcripts[tid]

        # obtain basic gene attributes
        (gchrom, gstrand, gstart, gend) = gene_attrs(gene_transcripts)

        if options.verbose:
            print >> sys.stderr, '\tFetching alignments...'

        # choose a single event position and weight the reads
        read_pos_weights = position_reads(clip_in, gchrom, gstart, gend, gstrand)

        # find splice junctions
        #junctions = map_splice_junctions(tx)

        if options.verbose:
            print >> sys.stderr, '\tCounting and computing in windows...'

        # count reads and compute p-values in windows
        window_stats = count_windows(clip_in, options.window_size, read_pos_weights, gene_transcripts, gstart, gend, total_reads, txome_size, windows_out)

        if options.verbose:
            print >> sys.stderr, '\tRefining peaks...'

        # post-process windows to peaks
        peaks = windows2peaks(read_pos_weights, gene_transcripts, gstart, window_stats, options.window_size, options.p_val, total_reads, txome_size)        

        # output peaks
        for pstart, pend, pcount, ppval in peaks:
            if ppval > 0:
                peak_score = int(2000/math.pi*math.atan(-math.log(ppval,1000)))
            else:
                peak_score = 1000
            cols = [gchrom, 'clip_peaks', 'peak', str(pstart), str(pend), str(peak_score), gstrand, '.', 'id "PEAK%d"; gene_id "%s"; count "%.1f"; p "%.2e"' % (peak_id,gene_id,pcount,ppval)]
            print >> peaks_out, '\t'.join(cols)
            peak_id += 1

    clip_in.close()
    peaks_out.close()


################################################################################
# cigar_endpoint
# 
# Input
#  aligned_read: pysam AlignedRead object.
#
# Output
#  genome_pos:   Endpoint of the alignment, considering the insertions and
#                 deletions in its CIGAR string (which includes splicing).
################################################################################
def cigar_endpoint(aligned_read):
    genome_pos = aligned_read.pos

    for (operation,length) in aligned_read.cigar:
        # match
        if operation in [0,7,8]:
            genome_pos += length

        # insertion
        elif operation in [1,3]:
            genome_pos += length

        # deletion
        elif operation == 2:
            pass

        else:
            print >> sys.stderr, 'Unknown CIGAR operation - %d, %s' % (operation, aligned_read.qname)

    return genome_pos


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
# count_reads
#
# Input
#  bam_file:    BAM file with NH tags.
#
# Output
#  total_reads: Total number of reads aligned with positive mapping quality,
#                counting based on the NH tags and counting paired end reads
#                as half.
################################################################################
def count_reads(bam_file):
    total_reads = 0.0
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        if aligned_read.mapq > 0:
            if aligned_read.is_paired:
                total_reads += 0.5/aligned_read.opt('NH')
            else:
                total_reads += 1.0/aligned_read.opt('NH')
    return total_reads


################################################################################
# count_windows
#
# Count the number of reads and compute the scan statistic p-value in each
# window through the gene.
#
# Input
#  clip_in:          Open pysam BAM file for clip-seq alignments.
#  window_size:      Scan statistic window size.
#  read_pos_weights: Sorted list of read alignment positions and weights.
#  gene_transcripts: Hash mapping transcript_id to isoform Gene objects,
#                     containing only keys for a specific gene.
#  gene_start:       Start of the gene's span.
#  gene_end:         End of the gene's span.
#  total_reads:      Total number of reads aligned to the transcriptome.
#  txome_size:       Total number of bp in the transcriptome.
#  windows_out:      Open file if we should print window stats, or None.
#
# Output
#  window_stats:     List of tuples (alignment count, p value) for all windows.
################################################################################
def count_windows(clip_in, window_size, read_pos_weights, gene_transcripts, gene_start, gene_end, total_reads, txome_size, windows_out):
    # set lambda using whole region (some day, compare this to the cufflinks estimate)
    # poisson_lambda = float(len(read_pos_weights)) / (gene_end - gene_start)

    # get gene info
    tid0 = gene_transcripts.keys()[0]
    chrom = gene_transcripts[tid0].chrom
    gene_id = gene_transcripts[tid0].kv['gene_id']

    reads_window_start = 0 # index of the first read position that fits in the window (except I'm allowing 0)
    reads_window_end = 0 # index of the first read position past the window

    # initialize index of the first junction ahead of the window start for each transcript
    junctions_i = {}
    for tid in gene_transcripts:
        junctions_i[tid] = 0

    # to avoid redundant computation
    precomputed_pvals = {}

    window_stats = []

    for window_start in range(gene_start, gene_end-window_size+1):
        window_end = window_start + window_size - 1

        # update read_window_start
        while reads_window_start < len(read_pos_weights) and read_pos_weights[reads_window_start][0] < window_start:
            reads_window_start += 1
        if reads_window_start >= len(read_pos_weights):
            break

        # update read_window_end
        while reads_window_end < len(read_pos_weights) and read_pos_weights[reads_window_end][0] <= window_end:
            reads_window_end += 1

        # count reads
        #window_count = reads_window_end - reads_window_start
        window_count_float = sum([read_pos_weights[i][1] for i in range(reads_window_start,reads_window_end)])

        # round count
        window_count = int(window_count_float + 0.5)

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
        if windows_out:
            cols = (chrom, window_start, gene_id, window_count, window_stats[-1][1], window_lambda)
            print >> windows_out, '%-5s %9d %18s %5d %8.1e %8.2e' % cols

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
# Return a 
#
# Input:
#  transcripts:  Hash mapping transcript_id keys to Gene class instances.
#
# Output:
#  gene_regions: Hash mapping gene_id keys to lists consisting of (chromosome,
#                start, end, strand) tuples with coordinates in GTF format.
################################################################################
def get_gene_regions(transcripts):
    gene_regions = {}

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
# position_reads
#
# Map read alignments for a gene of interest to a single genomic position,
# filtering for strand and quality, and assign weights.
#
# Input
#  clip_in:          Open pysam BAM file for clip-seq alignments.
#  gene_chrom:       Chromosome of the gene of interest.
#  gene_start:       Start coordinate of the gene of interest.
#  gene_end:         End coordinate of the gene of interest.
#  gene_strand:      Strand of the gene of interest.
#
# Output
#  read_pos_weights: Sorted list of read alignment positions and weights.
################################################################################
def position_reads(clip_in, gene_chrom, gene_start, gene_end, gene_strand):
    read_pos_weights = []

    if gene_chrom in clip_in.references:
        # for each read in span
        for aligned_read in clip_in.fetch(gene_chrom, gene_start, gene_end):
            # assign strand
            try:
                ar_strand = aligned_read.opt('XS')
            except:
                ar_strand = '+'
                if aligned_read.is_reverse:
                    ar_strand = '-'

            # check strand and quality
            if gene_strand == ar_strand and aligned_read.mapq > 0:
                if aligned_read.is_paired:
                    # map read to endpoint (closer to fragment center)
                    if aligned_read.is_reverse:
                        read_pos_weights.append((aligned_read.pos, 0.5/aligned_read.opt('NH')))
                    else:
                        read_pos_weights.append((cigar_endpoint(aligned_read), 0.5/aligned_read.opt('NH')))
                else:
                    # map read to midpoint
                    read_pos_weights.append((cigar_midpoint(aligned_read), 1.0/aligned_read.opt('NH')))

        # in case of differing read alignment lengths
        read_pos_weights.sort()

    return read_pos_weights


################################################################################
# merge_peaks_count
#
# Merge the given trimmed windows with counts using bedtools and then recount
# the reads in the new intervals.
#
# Input
#  trimmed_windows:  List of (start,end) tuples for significant windows, trimmed
#                     to be tight around read midpoints.
#  read_pos_weights: Sorted list of read alignment positions and weights.
#
# Output
#  peaks:            List of (start,end,count) tuples for significant trimmed
#                     and merged windows
################################################################################
def merge_peaks_count(trimmed_windows, read_pos_weights):
    # create BED intervals
    bed_intervals = []
    for wstart, wend in trimmed_windows:
        bed_a = ['chrFAKE', str(wstart-1), str(wend)]
        bed_intervals.append(pybedtools.create_interval_from_list(bed_a))
    bedtool = pybedtools.BedTool(bed_intervals)

    # merge BED intervals
    bedtool_merge = bedtool.merge(stream=True)

    # recount peaks
    read_positions = [pos for (pos,w) in read_pos_weights]
    peaks = []
    for bed_interval in bedtool_merge.features():
        pstart = bed_interval.start+1
        pend = bed_interval.end

        reads_start_i = bisect_left(read_positions, pstart)
        reads_end_i = bisect_right(read_positions, pend)
        # TODO: Count using the weights
        #read_count = reads_end_i - reads_start_i
        read_count = sum([read_pos_weights[i][1] for i in range(reads_start_i,reads_end_i)])

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

        wcount_round = int(wcount + 0.5)
        p_val = scan_stat_approx3(wcount_round, wend-wstart+1, txome_size, peak_lambda)
        peaks.append((wstart,wend,wcount_round,p_val))
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
# set_transcript_fpkms
#
# Input
#  transcripts: Hash mapping transcript_id to isoform Gene objects.
#  out_dir:     Directory in which to output the gene span GTF file.
# 
# Output
#  transcripts: Same hash, FPKM attribute set.
################################################################################
def set_transcript_fpkms(transcripts, out_dir, verbose):
    # read from isoforms.fpkm_tracking
    fpkm_in = open('%s/isoforms.fpkm_tracking' % out_dir)
    line = fpkm_in.readline()
    for line in fpkm_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        tid = a[0]
        fpkm = float(a[9])

        transcripts[tid].fpkm = fpkm
    fpkm_in.close()

    # fill in those that go missing
    missing_fpkms = 0
    for tid in transcripts:
        if transcripts[tid].fpkm == None:
            if verbose:
                print >> sys.stderr, 'WARNING: Missing FPKM for %s' % tid
            missing_fpkms += 1
            transcripts[tid].fpkm = 0
    print >> sys.stderr, 'WARNING: %d genes missing FPKM, assigned 0.' % missing_fpkms


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
    transcripts = read_genes(ref_gtf, key_id='transcript_id')
    gene_regions = get_gene_regions(transcripts)

    # print
    span_ref_gtf = '%s/span.gtf' % out_dir
    span_ref_open = open(span_ref_gtf, 'w')

    for gid in gene_regions:
        g = gene_regions[gid]
        cols = [g[0], 'clip_peaks', 'exon', str(g[1]), str(g[2]), '.', g[3], '.', gff.kv_gtf({'gene_id':gid, 'transcript_id':gid})]
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
    gene_regions = get_gene_regions(transcripts)

    txome_size = 0
    for gene_id in gene_regions:
        (gchrom,gstart,gend,gstrand) = gene_regions[gene_id]
        txome_size += gend - gstart - window_size + 1

    return txome_size


################################################################################
# trim_windows
#
# Input
#  windows:          List of (start,end) tuples for merged significant windows.
#  read_pos_weights: Sorted list of read alignment positions and weights.
#
# Output
#  trimmed_windows:  List of (start,end) tuples for significant windows, trimmed
#                     to be tight around read midpoints.
################################################################################
def trim_windows(windows, read_pos_weights):
    read_positions = [pos for (pos,w) in read_pos_weights]
    trimmed_windows = []
    for wstart, wend in windows:
        trim_start_i = bisect_left(read_positions, wstart)
        trim_end_i = bisect_right(read_positions, wend)
        # TODO: When are these ever not integers?
        trimmed_windows.append((int(read_positions[trim_start_i]), int(read_positions[trim_end_i-1]+0.5)))
    return trimmed_windows


################################################################################
# windows2peaks
#
# Convert window counts and p-values to peak calls.
#
# Input
#  read_pos_weights: Sorted list of read alignment positions and weights.
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
def windows2peaks(read_pos_weights, gene_transcripts, gene_start, window_stats, window_size, sig_p, total_reads, txome_size):
    merged_windows = merge_windows(window_stats, window_size, sig_p, gene_start)
    trimmed_windows = trim_windows(merged_windows, read_pos_weights)
    statless_peaks = merge_peaks_count(trimmed_windows, read_pos_weights)
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
