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
#
# --restart mode is semi-tested to not crash, but not fully tested to be the
# same as the initial run.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bowtie index> <reads1> [<reads2> ...]'
    parser = OptionParser(usage)

    # shrinking options
    parser.add_option('-i', '--initial_seed', dest='initial_seed', type='int', default=18, help='Seed length to initially consider for aligning the reads [Default: %default]')
    parser.add_option('-m', '--max_length', dest='max_length', type='int', default=60, help='Maximum seed length to consider for aligning the reads [Default: %default]')
    parser.add_option('-r', '--restart', dest='restart', action='store_true', default=False, help='Use the existing directory structure to restart a failed run [Default: %default]')
    parser.add_option('-s', '--salvage', dest='salvage', action='store_true', default=False, help='Use the existing directory structure to salvage a failed run; i.e. do not align more. [Default: %default]')
    parser.add_option('-t', '--stop', dest='stop', type='float', default=0.0005, help='Stop when an iteration only adds this percentage of new unique reads [Default: %default]')

    # tophat options
    parser.add_option('-p', '--num_threads', dest='num_threads', type='int', default=2, help='# of TopHat threads to launch [Default: %default]')
    parser.add_option('-G','--GTF', dest='gtf_file', help='Reference GTF file')
    parser.add_option('--transcriptome-index', dest='tx_index', default='txome', help='Transcriptome bowtie2 index [Default: %default]')

    # output options
    parser.add_option('-o', dest='output_dir', default='.', help='Output directory [Default %default]')
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
    if not os.path.isdir(options.output_dir):        
        os.mkdir(options.output_dir)
    os.chdir(options.output_dir)

    # find read length
    full_read_length = fastq_read_length(fastq_files[0])

    # max it out
    full_read_length = min(options.max_length, full_read_length)

    # fastq bit array
    read_finalized = bitarray()

    if options.salvage or options.restart:
        # find latest thout
        read_len = options.initial_seed        
        while os.path.isfile('thout%d/accepted_hits.bam' % (read_len+1)):
            read_len += 1

        if options.restart:
            # compute unique reads
            total_unique = 0
            for i in range(options.initial_seed, read_len):
                total_unique += count_unique('thout%d/unique.bam' % i)

            # take the fastq through the same process
            restart_finalized(read_finalized, fastq_files, read_len, options.initial_seed)

    else:
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
        total_unique = 0

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
    while read_len <= full_read_length and not options.salvage:
        # construct bloom filter for last iteration multimaps
        multimap_bf = construct_bloomfilter(read_len-1)

        # if no multimappers, quit
        if multimap_bf == None:
            break

        # split last iteration unique and multimappers
        iter_unique = split_iter_bam_bf(read_len-1, multimap_bf)
        if total_unique == 0:
            unique_pct = 100.0
        else:
            unique_pct = iter_unique / float(total_unique)
        print >> sys.stderr, 'Iteration %d mapped %d additional reads uniquely, an increase of %.3f\%' % (read_len-1, iter_unique, 100.0*unique_pct)

        # consider stopping if the benefit is small
        if unique_pct < options.stop:
            print >> sys.stderr, 'Stopping early due to reduced benefit'
            break
        total_unique += iter_unique

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
    merge_cmd = 'samtools merge -f accepted_hits.bam %s' % ' '.join(bam_files)
    print >> sys.stderr, merge_cmd
    subprocess.call(merge_cmd, shell=True)

    # clean up
    if os.path.isfile('multimap.txt'):
        os.remove('multimap.txt')
    if os.path.isfile('multimap.bf'):
        os.remove('multimap.bf')
    if os.path.isdir('tmp_sort'):
        shutil.rmtree('tmp_sort')

    if not options.keep_tmp:
        if os.path.isfile('iter.fq'):
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
# count_unique
#
# Count the number of aligned reads in a BAM file, assuming no multimaps.
################################################################################
def count_unique(bam_file):
    count = 0

    bam_in = pysam.Samfile(bam_file, 'rb')
    for aligned_read in bam_in:
        count += 1
    bam_in.close()

    return count


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

            if read_len > 0:
                print >> out_fq, header.rstrip()
                print >> out_fq, seq[:read_len].rstrip()
                print >> out_fq, '+'
                print >> out_fq, qual[:read_len].rstrip()

            read_finalized.append(False)

            header = fq_open.readline()
        fq_open.close()

    out_fq.close()


################################################################################
# restart_finalized
#
# If we're restarting, we need to take read_finalized through the same
# process as before.
#
# Input
#  read_finalized: Bit array describing whether the read alignment is done.
#  fastq_files:    List of fastq file names.
#  last_len:       Last read length with a successful iteration.
#  initial_seed:   Initial read length seed.
#
# Output:
#  read_finalized: Updated and ready to progress.
################################################################################
def restart_finalized(read_finalized, fastq_files, last_len, initial_seed):
    # initialize
    initial_fastq(fastq_files, 0, read_finalized)

    # update read_finalized for each iteration
    for rl in range(initial_seed+1, last_len+1):
        # construct bloom filter for last iteration multimaps
        multimap_bf = construct_bloomfilter(rl-1)

        # update fastq for multimappers and grow
        update_fastq(fastq_files, 0, read_finalized, multimap_bf)

        print >> sys.stderr, 'Updated read_finalized for %d' % rl


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

    # NOTE:
    #
    #  The code below here has a bug. I'm pretty sure that I want the condition to be:
    #    if aligned_read.opt('NH') > 1 and (write_all or aligned_read.qname in unmapped_set):
    #
    #  As is, in the middle iterations, it's going to print any read in unmapped_set because
    #  the first condition is never met because write_all=False. So unique reads are printed,
    #  too. But unique reads did not continue on, so they will never be in unmapped.bam.
    #
    #  In the last iteration, it's going to print every read because write_all=True, so the
    #  first condition is always met, so we're not checking for being in the unmapped_set.
    #  That's what we want anyway.

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
#  unique_count: The number of uniquely mapped reads in this run
################################################################################
def split_iter_bam_bf(read_len, multimap_bf):
    # original BAM for header
    original_bam = pysam.Samfile('thout%d/accepted_hits.bam' % read_len, 'rb')

    # initialize split BAM files
    unique_bam = pysam.Samfile('thout%d/unique.bam' % read_len, 'wb', template=original_bam)
    multimap_bam = pysam.Samfile('thout%d/multimap.bam' % read_len, 'wb', template=original_bam)

    # parse and split
    unique_count = 0
    for aligned_read in original_bam:
        if aligned_read.qname in multimap_bf:
            multimap_bam.write(aligned_read)
        else:
            unique_bam.write(aligned_read)
            unique_count +=1

    unique_bam.close()
    multimap_bam.close()

    return unique_count


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
                    if read_len > 0:
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
