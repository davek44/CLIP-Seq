#!/usr/bin/env python
from optparse import OptionParser
import copy, os, subprocess, tempfile
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

    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error(usage)
    else:
        ref_gtf = args[0]
        clip_bam = args[1]

    # make a new gtf w/ unspliced RNAs
    pre_ref_gtf = prerna_gtf(ref_gtf)

    # run Cufflinks on new gtf file and control BAM
    p = subprocess.Popen('cufflinks --no-js-test -G %s %s' % (pre_ref_gtf,), shell=True)
    os.waitpid(p.pid,0)

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
# __main__
################################################################################
if __name__ == '__main__':
    main()
