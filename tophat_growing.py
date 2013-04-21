#!/usr/bin/env python
from optparse import OptionParser
import copy, gzip, os, shutil, subprocess, sys
from bitarray import bitarray
import pybloomfilter, pysam

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
    usage = 'usage: %prog [options] <bowtie index> <reads1> [<reads2> ...]'
    parser = OptionParser(usage)

    # shrinking options
    parser.add_option('-s', '--initial_seed', dest='initial_seed', type='int', default=18, help='Seed length to initially consider for aligning the reads [Default: %default]')

    # tophat options
    parser.add_option('-p', '--num_threads', dest='num_threads', type='int', default=2, help='# of TopHat threads to launch [Default: %default]')
    parser.add_option('-G','--GTF', dest='gtf_file', help='Reference GTF file')
    parser.add_option('--transcriptome-index', dest='tx_index', default='txome', help='Transcriptome bowtie2 index [Default: %default]')

    # output options
    parser.add_option('-o', dest='output_dir', default'.', help='Output directory [Default %default]')
    parser.add_option('--tmp', dest='keep_tmp', default=False, action='store_true', help='Keep temporary output files [Default: %default]')

    (options,args) = parser.parse_args()    

    ############################################
    # setup
    ############################################
    # parse required input
    if len(args) < 2:
        parser.error('Missing required input.')
    else:
        bowtie_index = os.path.abspath(args[0])
        fastq_files = [os.path.abspath(fq) for fq in args[1:]]

    # convert paths to absolute
    options.gtf_file = os.path.abspath(options.gtf_file)
    options.tx_index = os.path.abspath(options.tx_index)

    # change to output directory
    os.chdir(options.output_dir)

    # find read length
    full_read_length = fastq_read_length(fastq_files[0])

    # fastq bit array
    read_finalized = bitarray()

    # clean directories
    if os.path.isdir('tmp_sort'):
        shutil.rmtree('tmp_sort')
    for read_len in range(full_read_length+1):
        if os.path.isdir('thout%d' % read_len):
            shutil.rmtree('thout%d' % read_len)

    # make a tmp dir for sorting
    os.mkdir('tmp_sort')

    ############################################
    # initial iteration
    ############################################
    read_len = options.initial_seed

    # trim reads
    initial_fastq(fastq_files, read_len, read_finalized)

    # align fastq
    subprocess.call('tophat -o thout%d -p %d -G %s -M --no-novel-juncs --transcriptome-index=%s %s iter.fq' % (read_len, options.num_threads, options.gtf_file, options.tx_index, bowtie_index), shell=True)
    if options.keep_tmp:
        os.rename('iter.fq','thout%d/iter.fq' % read_len)

    ############################################
    # onward iterations
    ############################################
    read_len += 1
    while read_len <= full_read_length:
        # construct bloom filter for last iteration multimaps
        multimap_bf = construct_bloomfilter(read_len-1)

        # if no multimappers, quit
        if multimap_bf == None:
            break

        # split last iteration unique and multimappers
        split_iter_bam_bf(read_len-1, multimap_bf)

        # update fastq for multimappers and grow
        update_fastq(fastq_files, read_len, read_finalized, multimap_bf)

        # align iteration fastq
        subprocess.call('tophat -o thout%d -p %d -G %s -M --no-novel-juncs --transcriptome-index=%s %s iter.fq' % (read_len, options.num_threads, options.gtf_file, options.tx_index, bowtie_index), shell=True)
        if options.keep_tmp:
            os.rename('iter.fq','thout%d/iter.fq' % read_len)

        # split lost multimappers from previous iteration
        split_lost_multi(read_len-1)

        # grow read
        read_len += 1

    ############################################
    # conclusion
    ############################################
    # combine all alignments
    bam_files = ['thout%d/accepted_hits.bam' % (read_len-1)]
    for rl in range(options.initial_seed, read_len-1):
        bam_files.append('thout%d/unique.bam' % rl)
        bam_files.append('thout%d/lost_multi.bam' % rl)
    subprocess.call('samtools merge all.bam %s' % ' '.join(bam_files), shell=True)

    # clean up
    os.remove('multimap.txt')
    os.remove('multimap.bf')
    os.rmdir('tmp_sort')
    if not options.keep_tmp:
        os.remove('iter.fq')
        for rl in range(options.initial_seed, read_len):
            shutil.rmtree('thout%d' % rl)


################################################################################
# construct_bloomfilter
#
# Input
#  read_len:    Trimmed read length used to find filenames
#
# Output
#  multimap_bf: Bloom filter storing multimapping read headers (or None)
################################################################################
def construct_bloomfilter(read_len):
    # get multi-mapping headers
    subprocess.call('samtools view thout%d/accepted_hits.bam | grep -v -w "NH:i:1" | cut -f1 | sort -u -T tmp_sort > multimap.txt' % read_len, shell=True)

    # count multimappers
    multimap_count = int(subprocess.check_output('wc -l multimap.txt', shell=True).split()[0])

    if multimap_count == 0:
        multimap_bf = None
    else:
        # add reads
        multimap_bf = pybloomfilter.BloomFilter(multimap_count, 1e-3, 'multimap.bf')
        for line in open('multimap.txt'):
            multimap_bf.add(line.rstrip())

    return multimap_bf


################################################################################
# fastq_read_length
#
# Input
#  fastq_file: Fastq file name
#
# Output
#  read length
################################################################################
def fastq_read_length(fastq_file):
    if fastq_file[-2:] == 'gz':
        fastq_open = gzip.open(fastq_file)
    else:
        fastq_open = open(fastq_file)
    header = fastq_open.readline()
    seq = fastq_open.readline().rstrip()
    return len(seq)


################################################################################
# initial_fastq
#
# Input
#  fastq_files:    List of fastq file names.
#  read_finalized: Bit array describing whether the read alignment is done.
#  read_len:       Length to trim the reads to.
#
# Output
#  iter.fq:        New fastq file containing trimmed reads.
################################################################################
def initial_fastq(fastq_files, read_len, read_finalized):
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

            print >> out_fq, header.rstrip()
            print >> out_fq, seq[:read_len].rstrip()
            print >> out_fq, '+'
            print >> out_fq, qual[:read_len].rstrip()

            read_finalized.append(False)

            header = fq_open.readline()
        fq_open.close()

    out_fq.close()


################################################################################
# split_lost_multi
#
# Input
#  read_len:       Trimmed read length used to find filenames. We're working in
#                   the dir defined by read_len, by will use unmapped reads
#                   from read_len+1
#  write_all:      Write all multi-mappers regardless of unmapped_set
#
# Output
#  lost_multi.bam: BAM file filtered for multimapping reads lost in the next iter
################################################################################
def split_lost_multi(read_len, write_all=False):
    # store unmapped headers    
    unmapped_set = set()
    if write_all == False and os.path.isfile('thout%d/unmapped.bam' % (read_len+1)):
        for aligned_read in pysam.Samfile('thout%d/unmapped.bam' % (read_len+1), 'rb'):
            unmapped_set.add(aligned_read.qname)

    # open multimapping bam
    multi_bam = pysam.Samfile('thout%d/accepted_hits.bam' % read_len, 'rb')

    # initialize lost multi mapped read BAM file
    lost_multi_bam = pysam.Samfile('thout%d/lost_multi.bam' % read_len, 'wb', template=multi_bam)

    # print lost multis
    for aligned_read in multi_bam:
        if aligned_read.opt('NH') > 1 and write_all or aligned_read.qname in unmapped_set:
            lost_multi_bam.write(aligned_read)

    lost_multi_bam.close()


################################################################################
# split_iter_bam_bf
#
# Input
#  read_len:     Trimmed read length used to find filenames
#  multimap_bf:  Bloom filter storing multimapping read headers
#
# Output
#  unique.bam:   BAM file filtered for uniquely mapping reads
#  multimap.bam: BAM file filtered for multimapping reads
################################################################################
def split_iter_bam_bf(read_len, multimap_bf):
    # original BAM for header
    original_bam = pysam.Samfile('thout%d/accepted_hits.bam' % read_len, 'rb')

    # initialize split BAM files
    unique_bam = pysam.Samfile('thout%d/unique.bam' % read_len, 'wb', template=original_bam)
    multimap_bam = pysam.Samfile('thout%d/multimap.bam' % read_len, 'wb', template=original_bam)

    # parse and split
    for aligned_read in original_bam:
        if aligned_read.qname in multimap_bf:
            multimap_bam.write(aligned_read)
        else:
            unique_bam.write(aligned_read)

    unique_bam.close()
    multimap_bam.close()


################################################################################
# update_fastq
#
# Input
#  fastq_files:    List of fastq file names.
#  read_len:       Length to trim the reads to (and find prior multimaps).
#  read_finalized: Bit array describing whether the read alignment is done.
#  multimap_bf:    Bloom filter storing multimapping read headers.
#
# Output
#  iter.fq:        New fastq file containing the trimmed reads we want.
################################################################################
def update_fastq(fastq_files, read_len, read_finalized, multimap_bf):
    # initialize
    out_fq = open('iter.fq', 'w')

    # for each fastq file
    i = 0
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

            if not read_finalized[i]:
                if header[1:].split()[0] in multimap_bf:
                    print >> out_fq, header.rstrip()
                    print >> out_fq, seq[:read_len].rstrip()
                    print >> out_fq, mid.rstrip()
                    print >> out_fq, qual[:read_len].rstrip()
                else:
                    read_finalized[i] = True

            i += 1

            header = fq_open.readline()
        fq_open.close()

    out_fq.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
