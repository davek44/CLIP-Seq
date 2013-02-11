#!/usr/bin/env python
from optparse import OptionParser
import copy, gzip, os, subprocess
import pysam

################################################################################
# tophat_growing.py
#
# Align reads with tophat where we initially start with a 5' seed and grow
# outward from there only if the read multimaps.
#
# Assuming all reads are the same length.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bowtie index> <reads1[,reads2,...]>'
    parser = OptionParser(usage)

    # shrinking options
    parser.add_option('-s', '--initial_seed', dest='initial_seed', type='int', default=18, help='Seed length to initially consider for aligning the reads [Default: %default]')

    # tophat options
    parser.add_option('-p', '--num_threads', dest='num_threads', type='int', default=2, help='# of TopHat threads to launch [Default: %default]')
    parser.add_option('-G','--GTF', dest='gtf_file', help='Reference GTF file')
    parser.add_option('--transcriptome-index', dest='tx_index', default='txome', help='Transcriptome bowtie2 index [Default: %default]')
    #parser.add_option('--no-novel-juncs', dest='no_novel_juncs', type='bool', action='store_true', help='Do not search for novel splice junctions [Default: %default]')
    (options,args) = parser.parse_args()

    # parse required input
    if len(args) < 2:
        parser.error(usage)
    else:
        bowtie_index = args[0]
        fastq_files = args[1].split(',')

    # find read length
    full_read_length = fastq_read_length(fastq_files[0])
    print >> sys.stderr, 'Read length: %d' % full_read_length

    # initialize multimap set for first iteration
    multimap_set = None

    for read_len in range(options.initial_seed, full_read_length+1):
        # make a new fastq of only multimappers
        make_iter_fastq(fastq_files, multimap_set, read_len)

        # align
        subprocess.call('tophat -o thout%d -p %d -G %s --no-novel-juncs --transcriptome-index=%s %s iter.fq' % (read_len, options.num_threads, options.gtf_file, options.tx_index, bowtie_index), shell=True)

        # parse BAM to split unique and store aligned and multimapped
        aligned_set, new_multimap_set = parse_iter_bam(read_len)

        # split lost multimappers from previous iteration
        if read_len > options.initial_seed:
            split_lost_multi(read_len-1, aligned_set)

        # update multimap set
        multimap_set = new_multimap_set

        # for debug purposes for now
        #os.rename('iter.fq', 'thout%d/iter.fq' % read_len)

    # save remaining multimappers
    split_lost_multi(full_read_length, set())

    # clean up
    os.remove('iter.fq')

    # combine all alignments
    bam_files = []
    for read_len in range(options.initial_seed, full_read_length+1):
        bam_files.append('thout%d/unique.bam' % read_len)
        bam_files.append('thout%d/lost_multi.bam' % read_len)
    subprocess.call('samtools merge all.bam %s' % ' '.join(bam_files), shell=True)



################################################################################
# fastq_read_length
#
# Input
#  fastq_file: fastq file name
#
# Output
#  read length
################################################################################
def fastq_read_length(fastq_file):
    fastq_in = open(fastq_file)    
    header = fastq_in.readline()
    seq = fastq_in.readline()
    return len(seq)

    
################################################################################
# make_iter_fastq
#
# Input
#  fastq_files: List of fastq file names
#  reads_set:   Set containing read headers for reads we want (or None for all)
#  read_len:    Length to trim the reads to
#
# Output
#  iter.fq:     New fastq file containing the trimmed reads we want
################################################################################
def make_iter_fastq(fastq_files, reads_set, read_len):
    out_fq = open('iter.fq', 'w')

    for fq_file in fastq_files:
        if fq_file[-2:] == 'gz':
            fq_open = gzip.open(fq_file)
        else:
            fq_open = open(fq_file)

        header = fq_open.readline()
        while header:
            seq = fq_open.readline()
            mid = fq_open.readline()
            qual = fq_open.readline()

            if reads_set == None or header[1:].split()[0] in reads_set:
                print >> out_fq, header.rstrip()
                print >> out_fq, seq[:read_len].rstrip()
                print >> out_fq, mid.rstrip()
                print >> out_fq, qual[:read_len].rstrip()

            header = fq_open.readline()
        fq_open.close()

    out_fq.close()


################################################################################
# parse_iter_bam
#
# Input
#  read_len: Trimmed read length used to find filenames
#
# Output
#  unique.bam:   BAM file filtered for uniquely mapping reads
#  aligned_set:  Set of aligned read headers
#  multimap_set: Set of multimapping read headers
################################################################################
def parse_iter_bam(read_len):
    # original bam for header
    original_bam = pysam.Samfile('thout%d/accepted_hits.bam' % read_len, 'rb')

    # initialize uniquely mapped read BAM file
    unique_bam = pysam.Samfile('thout%d/unique.bam' % read_len, 'wb', template=original_bam)

    # initialize alignment sets
    aligned_set = set()
    multimap_set = set()

    for aligned_read in original_bam:
        # save aligned read header
        aligned_set.add(aligned_read.qname)

        if aligned_read.opt('NH') == 1:
            # unique
            unique_bam.write(aligned_read)

        else:
            # multimap
            multimap_set.add(aligned_read.qname)

    unique_bam.close()

    return aligned_set, multimap_set


################################################################################
# split_lost_multi
#
# Input
#  read_len:       Trimmed read length used to find filenames
#  aligned_set:    Set of aligned read headers to detect lost multimappers
#
# Output
#  lost_multi.bam: BAM file filtered for multimapping reads lost in the next iter
################################################################################
def split_lost_multi(read_len, aligned_set):
    # open original bam
    original_bam = pysam.Samfile('thout%d/accepted_hits.bam' % read_len, 'rb')

    # initialize lost multi mapped read BAM file
    lost_multi_bam = pysam.Samfile('thout%d/lost_multi.bam' % read_len, 'wb', template=original_bam)

    for aligned_read in original_bam:
        if aligned_read.opt('NH') > 1 and aligned_read.qname not in aligned_set:
            lost_multi_bam.write(aligned_read)

    lost_multi_bam.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
