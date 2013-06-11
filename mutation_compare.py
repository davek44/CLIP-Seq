#!/usr/bin/env python
from optparse import OptionParser
import math, pdb
import pybedtools, pysam

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2
grdevices = importr('grDevices')

################################################################################
# mutation_compare.py
#
# Print a table and plot a heatmap of enrichments/depletions of mutation types
# in a pair of datasets.
#
# Two normalization options:
#  1. # sequenced bp
#  2. # mutations
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <mut1 file> <mut2 file>'
    parser = OptionParser(usage)
    parser.add_option('-m', dest='mut_norm', action='store_true', default=False, help='Normalize by # mutations (as opposed to sequenced bp) [Default: %default]')
    parser.add_option('-o', dest='output_pdf', default='mut_cmp.pdf', help='Output pdf file for heatmap [Default: %default]')
    parser.add_option('-r', dest='raw', action='store_true', default=False, help='Use raw mutation counts (as opposed to normalized for ACGT content) [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error(usage)
    else:
        mut1_file = args[0]
        mut2_file = args[1]

    mutation_profile1, seq_bp1 = parse_mutations(mut1_file, options.raw)
    mutation_profile2, seq_bp2 = parse_mutations(mut2_file, options.raw)

    relative_mutation_profile = compute_relative_profile(mutation_profile1, seq_bp1, mutation_profile2, seq_bp2)

    print_table(relative_mutation_profile)

    # make plotting data structures
    nts = ['_','A','C','G','T']
    nts1 = []
    nts2 = []
    rel = []
    for nt1 in nts:
        for nt2 in nts:
            nts1.append(nt1)
            nts2.append(nt2)
            rel.append(relative_mutation_profile[(nt1,nt2)])

    nts1_r = ro.StrVector(nts1)
    nts2_r = ro.StrVector(nts2)
    rel_r = ro.FloatVector(rel)

    df = ro.DataFrame({'nt1':nts1_r, 'nt2':nts2_r, 'rel':rel_r})

    # plot
    '''
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='nt2', y='nt1', fill='rel') + \
        ggplot2.geom_tile() + \
        ggplot2.scale_x_discrete(mut2_file, limits=nts) + \
        ggplot2.scale_y_discrete(mut1_file, limits=nts) + \
        ggplot2.scale_fill_gradient('Enrichment 1/2')
    '''

    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='nt2', y='nt1', fill='rel') + \
        ggplot2.geom_tile() + \
        ggplot2.scale_x_discrete('Read') + \
        ggplot2.scale_y_discrete('Reference') + \
        ggplot2.scale_fill_gradient2('log2 enrichment', low='darkblue', mid='white', high='darkred')

    # save to file
    grdevices.pdf(file=options.output_pdf)
    gp.plot()
    grdevices.dev_off()


################################################################################
# compute_relative_profile
#
# Compute the relative enrichment/depletion between two mutation profiles
# after further normalizing each.
################################################################################
def compute_relative_profile(mut_prof1, norm1, mut_prof2, norm2):
    rel_mut_prof = {}
    for mut_key in mut_prof1.keys():
        if mut_prof1[mut_key] == 0 or mut_prof2[mut_key] == 0:
            rel_mut_prof[mut_key] = 0
        else:
            rel_mut_prof[mut_key] = math.log((mut_prof1[mut_key]/float(norm1)) / (mut_prof2[mut_key]/float(norm2))) / math.log(2)
    return rel_mut_prof


################################################################################
# parse_mutations
#
# Parse my mutation_profile.py output.
################################################################################
def parse_mutations(mut_file, raw):
    nts = ['_','A','C','G','T']

    mut_open = open(mut_file)

    # acgt content
    acgt_content = {}
    line = mut_open.readline()
    line = mut_open.readline()
    while line.rstrip():
        a = line.split()
        acgt_content[a[0]] = int(a[1])
        line = mut_open.readline()
    
    # raw mutation counts
    raw_mutation_profile = {}
    line = mut_open.readline()
    line = mut_open.readline()
    line = mut_open.readline()
    for i in range(len(nts)):        
        a = line.split()
        for j in range(len(nts)):
            raw_mutation_profile[(nts[i],nts[j])] = int(a[j+1])
        line = mut_open.readline()
    line = mut_open.readline()

    # norm mutation counts
    norm_mutation_profile = {}
    line = mut_open.readline()
    line = mut_open.readline()
    line = mut_open.readline()
    for i in range(len(nts)):
        a = line.split()
        for j in range(len(nts)):
            norm_mutation_profile[(nts[i],nts[j])] = int(a[j+1])
        line = mut_open.readline()
    line = mut_open.readline()

    mut_open.close()

    if raw:
        mutation_profile = raw_mutation_profile
    else:
        mutation_profile = norm_mutation_profile

    return mutation_profile, sum(acgt_content.values())


################################################################################
# print_table
################################################################################
def print_table(mutation_profile):
    nts = ['_','A','C','G','T']
    print '  %6s %6s %6s %6s %6s' % tuple(nts)
    for nt1 in nts:
        row_counts = [mutation_profile.get((nt1,nt2),0) for nt2 in nts]
        print '%1s %6.3f %6.3f %6.3f %6.3f %6.3f' % tuple([nt1]+row_counts)
    print ''


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
