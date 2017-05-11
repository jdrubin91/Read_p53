__author__ = 'Jonathan Rubin'

import os
import sys
import pybedtools
from pybedtools import BedTool
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

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

def run(DMSO,Nutlin1,Nutlin3,DMSObedgraph,Nutlin1bedgraph,Nutlin3bedgraph,figuredir):
    D = BedTool(DMSO)
    N1 = BedTool(Nutlin1)
    N3 = BedTool(Nutlin3)
    Db = BedTool(DMSObedgraph)
    N1b = BedTool(Nutlin1bedgraph)
    N3b = BedTool(Nutlin3bedgraph)

    w1 = D
    w2 = N1-D
    w3 = N3-N1-D

    a = w1.map(Db, c='4', o='sum')
    b = w1.map(N1b, c='4', o='sum')
    c = w1.map(N3b, c='4', o='sum')

    F = plt.figure()
    ax = F.add_subplot(131)
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
    plt.savefig(figuredir + 'ChIP_Depth_Waves_bp.png', dpi=1200)




if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File/Figure directory
    filedir = parent_dir(homedir) + '/files/'
    file2dir = parent_dir(homedir) + '/files2/'
    figuredir = parent_dir(homedir) + '/figures/'

    DMSO = file2dir + 'DMSO1Hr_peaks.merge_200.bed'
    Nutlin1 = file2dir + 'Nutlin1Hr_peaks.merge_200.bed'
    Nutlin3 = file2dir + 'Nutlin3Hr_peaks.merge_200.bed'
    DMSObedgraph = '/projects/dowellLab/Taatjes/170420_K00262_0091_AHJTMHBBXX/trimmed/cat/bowtie2/sortedbam/genomecoveragebed/fortdf/DMSO1Hr_trimmed.fastq.bowtie2.sorted.reflected.BedGraph'
    Nutlin1bedgraph = '/projects/dowellLab/Taatjes/170420_K00262_0091_AHJTMHBBXX/trimmed/cat/bowtie2/sortedbam/genomecoveragebed/fortdf/Nutlin1Hr_trimmed.fastq.bowtie2.sorted.reflected.BedGraph'
    Nutlin3bedgraph = '/projects/dowellLab/Taatjes/170420_K00262_0091_AHJTMHBBXX/trimmed/cat/bowtie2/sortedbam/genomecoveragebed/fortdf/Nutlin3Hr_trimmed.fastq.bowtie2.sorted.reflected.BedGraph'
    run(DMSO,Nutlin1,Nutlin3,DMSObedgraph,Nutlin1bedgraph,Nutlin3bedgraph,figuredir)