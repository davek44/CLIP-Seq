#!/usr/bin/env python
from optparse import OptionParser
import pdb
import pybedtools, pysam

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
    bam_bt = pybedtools.BedTool(bam_file)
    if options.feature_file:
        feature_bt = pybedtools.BedTool(options.feature_file)
        bam_filt_bt = bam_bt.intersect(feature_bt, u=True, s=True, f=0.5)
    else:
        bam_filt_bt = bam_bt

    # load reference fasta
    fasta = pysam.Fastafile(fasta_file)

    # count mutations and nt's
    mutation_profile = {}
    acgt_content = {'A':0, 'C':0, 'G':0, 'T':0}
    for align in bam_filt_bt.features():
        align_a = align.__str__().split('\t')

        # get read sequence
        align_seq = align_a[9].upper()        

        # count nt's
        for nt in align_seq:
            if nt != 'N':
                acgt_content[nt] += 1

        # parse cigar string
        cigar_str = align_a[5]
        cigar_a = parse_cigar(cigar_str)

        # load reference sequence
        deleted_ref_nts = sum([nt_count for nt_count,code in cigar_a if code in ['D','N']])
        ref_seq = fasta.fetch(reference=align.chrom, start=align.start, end=align.end+deleted_ref_nts).upper()

        #strand = get_strand(int(align_a[1]))
        #if strand == '+':
        #    continue

        # process cigar string
        ref_i = 0
        align_i = 0
        for nt_count, code in cigar_a:
            # match
            if code == 'M':
                for j in range(nt_count):
                    if ref_seq[ref_i+j] != align_seq[align_i+j]:
                        mut_key = (ref_seq[ref_i+j], align_seq[align_i+j])
                        mutation_profile[mut_key] = mutation_profile.get(mut_key,0) + 1
                ref_i += nt_count
                align_i += nt_count

            # insertion
            elif code == 'I':
                for j in range(nt_count):
                    mut_key = ('_',align_seq[align_i+j])
                    mutation_profile[mut_key] = mutation_profile.get(mut_key,0) + 1
                align_i += nt_count

            # deletion
            elif code == 'D':
                for j in range(nt_count):
                    mut_key = (ref_seq[ref_i+j],'_')
                    mutation_profile[mut_key] = mutation_profile.get(mut_key,0) + 1
                ref_i += nt_count

            # intron
            elif code == 'N':
                ref_i += nt_count

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
# get_strand
################################################################################
def get_strand(flag):
    strand = "+"
    if flag & (0x10):	# minus strand if true.
        strand = "-"		
    return strand

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
# parse_cigar
#
# Return a list of (nt count, cigar code) tuples from a cigar string.
################################################################################
def parse_cigar(cigar_str):
    cigar_a = []
    nt_count = 0
    for i in range(len(cigar_str)):
        if cigar_str[i].isdigit():
            nt_count = 10*nt_count + int(cigar_str[i])
        else:
            cigar_a.append((nt_count,cigar_str[i]))
            nt_count = 0
    return cigar_a

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
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
