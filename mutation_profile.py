#!/usr/bin/env python
from optparse import OptionParser
import os, pdb, string, subprocess, tempfile
import pysam

################################################################################
# mutation_profile.py
#
# Count the # of occurrences of each mutation type in a BAM file, possible
# limited to those reads overlapping a gff or bed file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fasta file> <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-f', dest='feature_file', help='Limit to reads overlapping these features in GFF or BED')
    (options,args) = parser.parse_args()
    
    if len(args) != 2:
        parser.error(usage)
    else:
        fasta_file = args[0]
        bam_file = args[1]

    # filter BAM using features
    if options.feature_file:
        bam_feat_fd, bam_feat_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        subprocess.call('intersectBed -abam %s -b %s > %s' % (bam_file, options.feature_file, bam_feat_file), shell=True)
        bam_file = bam_feat_file

    # load reference fasta
    fasta = pysam.Fastafile(fasta_file)

    # count mutations and nt's
    mutation_profile = {}
    acgt_content = {'A':0, 'C':0, 'G':0, 'T':0}

    bam_in = pysam.Samfile(bam_file, 'rb')
    for aligned_read in bam_in:
        if aligned_read.mapq > 0:
            # get read sequence
            align_seq = aligned_read.seq.upper()

            # correct for deleted nts
            deleted_ref_nts = sum([nt_count for code, nt_count in aligned_read.cigar if code in [2,3]])

            # get reference sequence
            ref_seq = fasta.fetch(reference=bam_in.references[aligned_read.tid], start=aligned_read.pos, end=aligned_read.aend+deleted_ref_nts).upper()

            # correct for multimap
            nh_tag = float(aligned_read.opt('NH'))

            ref_i = 0
            align_i = 0
            for code, nt_count in aligned_read.cigar:
                # match
                if code == 0:
                    for j in range(nt_count):
                        # acgt content
                        if aligned_read.is_reverse:
                            acgt_content[rc(ref_seq[ref_i+j])] += 1/nh_tag
                        else:
                            acgt_content[ref_seq[ref_i+j]] += 1/nh_tag

                        # mutations
                        if ref_seq[ref_i+j] != align_seq[align_i+j]:
                            if aligned_read.is_reverse:
                                mut_key = (rc(ref_seq[ref_i+j]), rc(align_seq[align_i+j]))
                            else:
                                mut_key = (ref_seq[ref_i+j], align_seq[align_i+j])
                            mutation_profile[mut_key] = mutation_profile.get(mut_key,0) + 1/nh_tag

                    # update indexes
                    ref_i += nt_count
                    align_i += nt_count

                # insertion
                elif code == 1:
                    # mutations
                    for j in range(nt_count):
                        if aligned_read.is_reverse:
                            mut_key = ('_', rc(align_seq[align_i+j]))
                        else:
                            mut_key = ('_', align_seq[align_i+j])
                        mutation_profile[mut_key] = mutation_profile.get(mut_key,0) + 1/nh_tag

                    # update indexes
                    align_i += nt_count

                # deletion
                elif code == 2:
                    # mutations
                    for j in range(nt_count):
                        if aligned_read.is_reverse:
                            mut_key = (rc(ref_seq[ref_i+j]), '_')
                        else:
                            mut_key = (ref_seq[ref_i+j], '_')

                    # update indexes
                    ref_i += nt_count

                # intron
                elif code == 3:
                    # update indexes
                    ref_i += nt_count

    # clean
    bam_in.close()
    if options.feature_file:
        os.close(bam_feat_fd)
        os.remove(bam_feat_file)

    # print acgt content
    print 'ACGT content:'
    for nt in ['A','C','G','T']:
        print '%s %9d' % (nt,acgt_content[nt])
    print ''

    # print raw stats
    print 'Raw mutations:'
    print_table(mutation_profile)

    # print normalized stats
    print 'Normalized mutations:'
    norm_mutation_profile = normalize_profile(mutation_profile, acgt_content)
    print_table(norm_mutation_profile)


################################################################################
# normalize_profile
#
# Normalize the mutation counts using the nt content of the sequencing reads.
# Leave insertions alone.
################################################################################
def normalize_profile(mutation_profile, acgt_content):
    nt_total = sum(acgt_content.values())

    norm_mutation_profile = {}

    nts = ['_','A','C','G','T']
    for nt1 in nts:
        # determine factor
        if nt1 == '_':
            factor = 1
        else:
            factor = 4.0*acgt_content[nt1]/nt_total

        # multiply
        for nt2 in nts:
            norm_mutation_profile[(nt1,nt2)] = factor*mutation_profile.get((nt1,nt2),0)

    return norm_mutation_profile


################################################################################
# print_table
################################################################################
def print_table(mutation_profile):
    nts = ['_','A','C','G','T']
    print '  %7s %7s %7s %7s %7s' % tuple(nts)
    for nt1 in nts:
        row_counts = [mutation_profile.get((nt1,nt2),0) for nt2 in nts]
        print '%1s %7d %7d %7d %7d %7d %7d' % tuple([nt1]+row_counts+[sum(row_counts)])

    col_sums = []
    for nt2 in nts:
        col_sums.append(sum([mutation_profile.get((nt1,nt2),0) for nt1 in nts]))
    print '  %7d %7d %7d %7d %7d' % tuple(col_sums)
    print ''


################################################################################
# rc
#
# Reverse complement sequence
################################################################################
def rc(seq):
    return seq.translate(string.maketrans("ATCGatcg","TAGCtagc"))[::-1]


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
