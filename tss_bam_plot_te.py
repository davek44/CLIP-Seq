#!/usr/bin/env python
from optparse import OptionParser
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2
import math, os, pdb, shutil, subprocess, tempfile
import gff

grdevices = importr('grDevices')

################################################################################
# tss_bam_plot_te.py
#
# Plot read coverage in a BAM file surrounding the TSS's defined in a gtf file.
# Split plots by overlap with a transposable element family.
#
# The geometric mean idea is nice in theory to smooth outliers, but the IP
# sequencing will have more variance than the Control which throws off the
# the alignment # normalization.
################################################################################

hg19_reps_gff = '%s/research/common/data/genomes/hg19/annotation/repeatmasker/hg19.fa.out.tpf.gff' % os.environ['HOME']

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file> <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='control_bam_file', default=None, help='Control BAM file')
    parser.add_option('-d', dest='downstream', default=2000, type='int', help='TSS downstream [Default: %default]')
    parser.add_option('-g', dest='geo_mean', default=False, action='store_true', help='Plot coverage geometric means [Default: %default]')
    parser.add_option('-o', dest='out_prefix', default='tss', help='Output prefix [Default: %default]')
    parser.add_option('-u', dest='upstream', default=5000, type='int', help='TSS upstream [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) < 2:
        parser.error('Must provide gtf file and BAM file')
    else:
        gtf_file = args[0]
        bam_files = args[1:]

    # map genes to TEs
    gene_te = intersect_gene_te(gtf_file, options.upstream, options.downstream)

    # count tss
    tss_count = {}
    for gene in gene_te:
        for te in gene_te[gene]:
            tss_count[te] = tss_count.get(te,0) + 1

    # collect intervals
    tss_intervals = get_tss(gtf_file, options.upstream, options.downstream)    

    # process bam
    te_tss_cov = process_bams(bam_files, tss_intervals, gene_te, options.out_prefix, options.upstream, options.downstream, options.geo_mean)
    if options.control_bam_file:
        control_te_tss_cov = process_bams([options.control_bam_file], tss_intervals, gene_te, options.out_prefix, options.upstream, options.downstream, options.geo_mean)

    ############################################
    # output
    ############################################
    make_output(te_tss_cov, options.out_prefix, options.upstream, options.downstream)

    if options.control_bam_file:
        # normalize
        main_aligns = 0
        for bam_file in bam_files:
            main_aligns += float(subprocess.check_output('samtools view -c %s' % bam_file, shell=True))
        control_aligns = float(subprocess.check_output('samtools view -c %s' % options.control_bam_file, shell=True))

        if options.geo_mean:
            for te in te_tss_cov:
                te_tss_cov[te] = [1000000.0*math.exp(float(tc)/tss_count[te])/main_aligns for tc in te_tss_cov[te]]
            for te in control_te_tss_cov:
                control_te_tss_cov[te] = [1000000.0*math.exp(float(tc)/tss_count[te])/control_aligns for tc in control_te_tss_cov[te]]

        else:
            for te in te_tss_cov:
                te_tss_cov[te] = [1000000.0*tc/tss_count[te]/main_aligns for tc in te_tss_cov[te]]
            for te in control_te_tss_cov:
                control_te_tss_cov[te] = [1000000.0*tc/tss_count[te]/control_aligns for tc in control_te_tss_cov[te]]

        # plot subtraction
        sub_te_cov = {}
        for te in te_tss_cov:
            sub_te_cov[te] = [te_tss_cov[te][i]-control_te_tss_cov[te][i] for i in range(len(te_tss_cov[te]))]
        make_output(sub_te_cov, options.out_prefix+'_sub', options.upstream, options.downstream)

        # plot and
        make_output_and(te_tss_cov, control_te_tss_cov, options.out_prefix+'_and', options.upstream, options.downstream)
    

################################################################################
# get_tss
#
# Make a hash keyed by chromosome and valued with lists of (start,end,strand)
# TSS tuples. Sort by chromosome.
################################################################################
def get_tss(gtf_file, upstream, downstream):
    tss_intervals = {}

    genes = gff.read_genes(gtf_file)
    for gid in genes:
        g = genes[gid]
        if g.strand == '+':
            istart = g.exons[0].start-upstream
            iend = g.exons[0].start+downstream
        else:
            istart = g.exons[-1].end-downstream
            iend = g.exons[-1].end+upstream

        tss_intervals.setdefault(g.chrom,[]).append((istart,iend,g.strand,gid))

    for chrom in tss_intervals:
        tss_intervals[chrom].sort()

    return tss_intervals


################################################################################
# intersect_gene_te
#
# Map genes to the transposable elements they contain.
################################################################################
def intersect_gene_te(gtf_file, upstream, downstream):
    # focus on promoter
    tmp_fd, tmp_file = tempfile.mkstemp()
    gff.promoters(gtf_file, upstream, downstream, tmp_file)
    
    # intersect genes w/ repeats
    # hash transposon nt by gene
    gene_trans = {}
    p = subprocess.Popen('intersectBed -wo -a %s -b %s' % (tmp_file,hg19_reps_gff), shell=True, stdout=subprocess.PIPE)
    line = p.stdout.readline()
    while line:
        a = line.split('\t')

        # get names
        gene = gff.gtf_kv(a[8])['transcript_id']
        rep_kv = gff.gtf_kv(a[17])
        rep = rep_kv['repeat']
        fam = rep_kv['family']

        # add nt
        if gene not in gene_trans:
            gene_trans[gene] = {}
        gene_trans[gene][(rep,fam)] = gene_trans[gene].get((rep,fam),0) + int(a[18])
        gene_trans[gene][('*',fam)] = gene_trans[gene].get(('*',fam),0) + int(a[18])
        gene_trans[gene][('*','*')] = gene_trans[gene].get(('*','*'),0) + int(a[18])

        line = p.stdout.readline()
    p.communicate()

    # create a fake family for dTE-lncRNAs
    for line in open(gtf_file):
        a = line.split('\t')
        tid = gff.gtf_kv(a[8])['transcript_id']
        if tid not in gene_trans:
            gene_trans[tid] = {('n','n'):1}

    return gene_trans


################################################################################
# make_bed
#
# Make a bed file of the TSS intervals to give to samtools
################################################################################
def make_bed(tss_intervals, bed_file):
    bed_out = open(bed_file,'w')
    for chrom in tss_intervals:
        for (start,end,strand,gid) in tss_intervals[chrom]:
            print >> bed_out, '%s\t%d\t%d' % (chrom,start-1,end)
    bed_out.close()


################################################################################
# make_output
################################################################################
def make_output(te_tss_cov, out_prefix, upstream, downstream):
    # clean raw counts dir
    if os.path.isdir('%s_raw' % out_prefix):
        shutil.rmtree('%s_raw' % out_prefix)
    os.mkdir('%s_raw' % out_prefix)

    # dump raw counts to file
    for te in te_tss_cov:
        if te[0] in ['n','*','HERVH-int','L2a','AluSx','AluJb','MIRb','LTR7'] and te[1] in ['n','*','LINE/L1','SINE/Alu','LTR/ERV1','LTR/ERVL-MaLR','LINE/L2','LTR/ERVL','SINE/MIR','DNA/hAT-Charlie','LTR/ERVK','DNA/TcMar-Tigger']:
            raw_out = open('%s_raw/%s_%s.txt' % (out_prefix,te[0].replace('/','_'),te[1].replace('/','_')),'w')
            for i in range(-upstream,downstream+1):
                print >> raw_out, '%e\t%e' % (i, te_tss_cov[te][upstream+i])
            raw_out.close()

    # clean plot dirs
    if os.path.isdir('%s_plot' % out_prefix):
        shutil.rmtree('%s_plot' % out_prefix)
    os.mkdir('%s_plot' % out_prefix)

    # make data structures
    tss_i = ro.IntVector(range(-upstream,downstream+1))
    for te in te_tss_cov:
        if te[0] in ['n','*','HERVH-int','L2a','AluSx','AluJb','MIRb','LTR7'] and te[1] in ['n','*','LINE/L1','SINE/Alu','LTR/ERV1','LTR/ERVL-MaLR','LINE/L2','LTR/ERVL','SINE/MIR','DNA/hAT-Charlie','LTR/ERVK','DNA/TcMar-Tigger']:
            cov = ro.FloatVector(te_tss_cov[te])
            df = ro.DataFrame({'tss_i':tss_i, 'cov':cov})

            # construct full plot
            gp = ggplot2.ggplot(df) + \
                ggplot2.aes_string(x='tss_i', y='cov') + \
                ggplot2.geom_point() + \
                ggplot2.scale_x_continuous('TSS index') + \
                ggplot2.scale_y_continuous('Coverage')

            # plot to file
            grdevices.pdf(file='%s_plot/%s_%s.pdf' % (out_prefix,te[0].replace('/','_'),te[1].replace('/','_')))
            gp.plot()
            grdevices.dev_off()


################################################################################
# make_output_and
################################################################################
def make_output_and(te_tss_cov, control_te_tss_cov, out_prefix, upstream, downstream):
    # clean raw counts dir
    if os.path.isdir('%s_raw' % out_prefix):
        shutil.rmtree('%s_raw' % out_prefix)
    os.mkdir('%s_raw' % out_prefix)

    # dump raw counts to file
    for te in te_tss_cov:
        if te[0] in ['n','*','HERVH-int','L2a','AluSx','AluJb','MIRb','LTR7'] and te[1] in ['n','*','LINE/L1','SINE/Alu','LTR/ERV1','LTR/ERVL-MaLR','LINE/L2','LTR/ERVL','SINE/MIR','DNA/hAT-Charlie','LTR/ERVK','DNA/TcMar-Tigger']:
            raw_out = open('%s_raw/%s_%s.txt' % (out_prefix,te[0].replace('/','_'),te[1].replace('/','_')),'w')
            for i in range(-upstream,downstream+1):
                print >> raw_out, '%d\t%e\t%e' % (i, te_tss_cov[te][upstream+i], control_te_tss_cov[te][upstream+i])
            raw_out.close()

    # clean plot dirs
    if os.path.isdir('%s_plot' % out_prefix):
        shutil.rmtree('%s_plot' % out_prefix)
    os.mkdir('%s_plot' % out_prefix)

    # make data structures
    tss_i = ro.IntVector(2*range(-upstream,downstream+1))
    labels = ro.StrVector(['Main']*(upstream+downstream+1)+['Control']*(upstream+downstream+1))
    for te in te_tss_cov:
        if te[0] in ['n','*','HERVH-int','L2a','AluSx','AluJb','MIRb','LTR7'] and te[1] in ['n','*','LINE/L1','SINE/Alu','LTR/ERV1','LTR/ERVL-MaLR','LINE/L2','LTR/ERVL','SINE/MIR','DNA/hAT-Charlie','LTR/ERVK','DNA/TcMar-Tigger']:
            cov = ro.FloatVector(te_tss_cov[te] + control_te_tss_cov[te])
            df = ro.DataFrame({'tss_i':tss_i, 'cov':cov, 'label':labels})

            # construct full plot
            gp = ggplot2.ggplot(df) + \
                ggplot2.aes_string(x='tss_i', y='cov', colour='label') + \
                ggplot2.geom_point() + \
                ggplot2.scale_x_continuous('TSS index') + \
                ggplot2.scale_y_continuous('Coverage') + \
                ggplot2.scale_colour_discrete('')

            # plot to file
            grdevices.pdf(file='%s_plot/%s_%s.pdf' % (out_prefix,te[0].replace('/','_'),te[1].replace('/','_')))
            gp.plot()
            grdevices.dev_off()


################################################################################
# process_bams
#
# Count read coverage in a BAM file around the tss_intervals given. Hash
# counts by TE family.
################################################################################
def process_bams(bam_files, tss_intervals, gene_te, out_prefix, upstream, downstream, geo_mean):
    # initialize data structures
    ti_chrom = None
    ti_i = -1
    active_tss = []
    te_tss_cov = {}
    for gene in gene_te:
        for te in gene_te[gene]:
            if not te in te_tss_cov:
                te_tss_cov[te] = [0]*(upstream+downstream+1)

    make_bed(tss_intervals,'%s_tss.bed' % out_prefix)

    for bam_file in bam_files:
        # pileup
        proc = subprocess.Popen('samtools mpileup -l %s_tss.bed %s' % (out_prefix,bam_file), shell=True, stdout=subprocess.PIPE)
        line = proc.stdout.readline()
        while line:
            a = line.split()

            chrom = a[0]
            pos = int(a[1])
            cov = int(a[3])

            if chrom in tss_intervals:
                # update active TSS
                active_tss, ti_chrom, ti_i = update_active_tss(tss_intervals, ti_chrom, ti_i, active_tss, chrom, pos)

                # append to tss_cov
                for tss in active_tss:
                    (tchrom,tstart,tend,tstrand,tid) = tss
                    if tstrand == '+':
                        for te in gene_te[tid]:
                            if geo_mean:
                                te_tss_cov[te][pos-tstart] += math.log(cov)
                            else:
                                te_tss_cov[te][pos-tstart] += cov
                    else:
                        for te in gene_te[tid]:
                            if geo_mean:
                                te_tss_cov[te][tend-pos] += math.log(cov)
                            else:
                                te_tss_cov[te][tend-pos] += cov

            # next line
            line = proc.stdout.readline()

        # finish
        line = proc.communicate()[0]

    return te_tss_cov


################################################################################
# update_active_tss
#
# Maintain a set of TSS intervals which chrom and pos fall inside.
################################################################################
def update_active_tss(tss_intervals, ti_chrom, ti_i, active_tss, chrom, pos):
    # drop finished tss
    while active_tss and (active_tss[0][0] != chrom or active_tss[0][2] < pos):
        active_tss = active_tss[1:]

    # initialize tss_intervals indexes
    if ti_chrom != chrom:
        ti_chrom = chrom
        ti_i = 0

    # add opened tss
    while ti_i < len(tss_intervals[ti_chrom]):
        (tstart,tend,tstrand,tid) = tss_intervals[ti_chrom][ti_i]
        if pos < tstart:
            # before TSS: no new active
            break
        elif tstart <= pos <= tend:
            # inside TSS: add and move to next
            active_tss.append((ti_chrom,tstart,tend,tstrand,tid))
            ti_i += 1
        else:
            # past TSS: move to next
            ti_i += 1

    return active_tss, ti_chrom, ti_i


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
