__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import pybedtools
from pybedtools import BedTool
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os
import math
import load


#The purpose of this script is to take true positive MACS2 peak calls (as defined in get_true_pos.py script) and do
#fold change analysis on corresponding GRO-Seq data (+/- Nutlin in this case)


def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def interval_overlap(interval1,interval2):
    boolean = False
    chrom1,start1,stop1 = interval1
    chrom2,start2,stop2 = interval2
    if chrom1 == chrom2 and (stop2 >= start1 >= start2 or start2 <= stop1 <= stop2):
        boolean = True

    return boolean

def format_boxplot(bp):
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#1b9e77' )

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)

def run(Nutlin1,GRODMSO,GRONutlin,figuredir,filedir):
    x = BedTool(Nutlin1)
    y = BedTool(GRODMSO)
    z = BedTool(GRONutlin)

    a = x.map(y, c='4', o='sum')
    b = x.map(z, c='4', o='sum')

    F = plt.figure()
    ax = F.add_subplot(111)
    ax.set_title('Nutlin1hr vs. DMSO')
    ax.set_ylabel('Count')
    ax.set_xlabel('Log2 Fold Change (Nutlin/DMSO)')
    ax.hist([math.log(float(n[3])/float(m[3]),2) for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.'], bins =100)
    ax.set_xlim([-20,20])
    plt.axvline(0, color='red',linestyle='dashed')
    # ax.set_ylabel('Log2 Fold Change (Nutlin/DMSO)')
    # ax.set_xticklabels(['Nutlin/DMSO'])
    # bp = ax.boxplot([math.log(float(n[3])/float(m[3]),2) for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.'],patch_artist=True)
    # format_boxplot(bp)
    plt.savefig(figuredir + 'GRO_Analysis_Fold_Change_hist.png', dpi=1200)

    outfile = open(filedir + 'false_positives_GRO-Seq_fold_change.bed','w')
    for interval in [m[:3] for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.' and float(n[3])/float(m[3]) < 1]:
        outfile.write('\t'.join(interval) + '\n')
    outfile.close()

    outfile = open(filedir + 'p53_txn_fold_change.bed','w')
    for interval in sorted([(m[0],m[1],m[2],float(n[3])/float(m[3])) for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.'], key=lambda x: x[-1], reverse=True):
        outfile.write('\t'.join(interval[:3]) + '\t' + str(interval[-1]) + '\n')
    outfile.close()

def run2(Nutlin1,DMSO,GRODMSO,GRONutlin,figuredir,filedir):
    N1 = BedTool(Nutlin1)
    D = BedTool(DMSO)
    y = BedTool(GRODMSO)
    z = BedTool(GRONutlin)

    a = N1.map(y, c='4', o='sum')
    b = N1.map(z, c='4', o='sum')

    c = (N1-D).map(y, c='4', o='sum')
    d = (N1-D).map(z, c='4', o='sum')

    F = plt.figure()
    ax = F.add_subplot(121)
    ax.set_title('Nutlin 1Hr All Peaks')
    ax.set_ylabel('Count')
    ax.set_xlabel('Log2 Fold Change (Nutlin/DMSO)')
    ax.hist([math.log(float(n[3])/float(m[3]),2) for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.'], bins =50)
    ax.set_xlim([-20,20])
    plt.axvline(0, color='red',linestyle='dashed')
    ax = F.add_subplot(122)
    ax.set_title('Nutlin 1Hr - DMSO Peaks')
    ax.set_ylabel('Count')
    ax.set_xlabel('Log2 Fold Change (Nutlin/DMSO)')
    ax.hist([math.log(float(n[3])/float(m[3]),2) for m,n in zip(c,d) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.'], bins =50)
    ax.set_xlim([-20,20])
    plt.axvline(0, color='red',linestyle='dashed')
    # ax.set_ylabel('Log2 Fold Change (Nutlin/DMSO)')
    # ax.set_xticklabels(['Nutlin/DMSO'])
    # bp = ax.boxplot([math.log(float(n[3])/float(m[3]),2) for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.'],patch_artist=True)
    # format_boxplot(bp)
    plt.savefig(figuredir + 'GRO_Analysis_Fold_Change_hist.png', dpi=1200)

    outfile = open(filedir + 'false_positives_GRO-Seq_fold_change.bed','w')
    for interval in [m[:3] for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.' and float(n[3])/float(m[3]) < 1]:
        outfile.write('\t'.join(interval) + '\n')
    outfile.close()

    outfile = open(filedir + 'p53_txn_fold_change.bed','w')
    for interval in sorted([(m[0],m[1],m[2],float(n[3])/float(m[3])) for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.'], key=lambda x: x[-1], reverse=True):
        outfile.write('\t'.join(interval[:3]) + '\t' + str(interval[-1]) + '\n')
    outfile.close()


def chip_gro_scatter(Nutlin1,DMSO,DMSOCHIP,Nutlin1CHIP,GRODMSO,GRONutlin,figuredir,file2dir,genes,promoters):
    N1 = BedTool(Nutlin1)
    D = BedTool(DMSO)
    DC = BedTool(DMSOCHIP)
    NC = BedTool(Nutlin1CHIP)
    GD = BedTool(GRODMSO)
    GN1 = BedTool(GRONutlin)
    g = BedTool(genes).sort().merge()
    p = BedTool(promoters)
    # p = BedTool([(gene[0],gene[1],gene[1]) for gene in g])

    print len(g), len(p)

    a = (N1+p).map(DMSOCHIP, c=4, o="sum", null="0")
    b = (N1+p).map(Nutlin1CHIP, c=4, o="sum", null="0")
    c = (N1+p).intersect(g,wb=True).sort().map(GRODMSO, c=4, o="sum", null="0")
    d = (N1+p).intersect(g,wb=True).sort().map(GRONutlin, c=4, o="sum", null="0")

    print len(a), len(c)



    F = plt.figure()
    ax = F.add_subplot(121)
    ax.set_title('Gene Target')
    ax.set_ylabel('ChIP-Seq signal at promoter Log2(Nutlin1hr/DMSO)')
    ax.set_xlabel('GRO-Seq signal at gene target Log2(Nutlin1hr/DMSO)')
    ax.scatter([float(n[3])/float(m[3]) for m,n in zip(a,b) if n[3] != "0" and m[3] != "0"], [float(l[-1])/float(k[-1]) for l,k in zip(c,d) if l[-1] != "0" and k[-1] != "0"])

    a = (N1).map(DMSOCHIP, c=4, o="sum", null="0")
    b = (N1).map(Nutlin1CHIP, c=4, o="sum", null="0")
    c = (N1).map(GRODMSO, c=4, o="sum", null="0")
    d = (N1).map(GRONutlin, c=4, o="sum", null="0")

    ax2 = F.add_subplot(122)
    ax2.set_title('ChIP-Seq Peaks (non-gene associated)')
    ax2.set_ylabel('ChIP-Seq Log2(Nutlin1hr/DMSO)')
    ax2.set_xlabel('GRO-Seq Log2(Nutlin1hr/DMSO)')
    ax2.scatter([float(n[3])/float(m[3]) for m,n in zip(a,b) if n[3] != "0" and m[3] != "0"], [float(l[3])/float(k[3]) for l,k in zip(c,d) if l[3] != "0" and k[3] != "0"])

    plt.savefig(figuredir + 'ChIP_GRO_scatter.png', dpi=1200)



if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    file2dir = parent_dir(homedir) + '/files2/'
    figuredir = parent_dir(homedir) + '/figures/'

    Nutlin1 = filedir + 'Nutlin1Hr_peaks.merge.200.bed.true_positive.bed'
    GRODMSO = '/scratch/Users/joru1876/Allen2014_NutlinGRO/DMSO1Hr.mp.reflected.sorted.merge.BedGraph'
    GRONutlin = '/scratch/Users/joru1876/Allen2014_NutlinGRO/Nutlin1Hr.mp.reflected.sorted.merge.BedGraph'
    # run(Nutlin1,GRODMSO,GRONutlin,figuredir,filedir)

    Nutlin1 = file2dir + 'Nutlin1Hr_peaks.merge_200.bed'
    DMSO = file2dir + 'DMSO1Hr_peaks.merge_200.bed'
    # run2(Nutlin1,DMSO,GRODMSO,GRONutlin,figuredir,file2dir)

    genes = file2dir + 'refGene.bed'
    promoters = file2dir + 'refGene_hg19_TSS.bed'
    DMSOCHIP = '/projects/dowellLab/Taatjes/170420_K00262_0091_AHJTMHBBXX/trimmed/cat/bowtie2/sortedbam/genomecoveragebed/fortdf/DMSO1Hr_trimmed.fastq.bowtie2.sorted.reflected.BedGraph'
    Nutlin1CHIP = '/projects/dowellLab/Taatjes/170420_K00262_0091_AHJTMHBBXX/trimmed/cat/bowtie2/sortedbam/genomecoveragebed/fortdf/Nutlin1Hr_trimmed.fastq.bowtie2.sorted.reflected.BedGraph'
    chip_gro_scatter(Nutlin1,DMSO,DMSOCHIP,Nutlin1CHIP,GRODMSO,GRONutlin,figuredir,file2dir,genes,promoters)



