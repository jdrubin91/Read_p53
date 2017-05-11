__author__ = 'Jonathan Rubin'

import os
import sys
import pybedtools
from pybedtools import BedTool
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

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

    a = w1.map(Db, c='4', o='sum',null=0)
    b = w1.map(N1b, c='4', o='sum',null=0)
    c = w1.map(N3b, c='4', o='sum',null=0)

    F = plt.figure()
    ax = F.add_subplot(131)
    ax.set_title('Wave1')
    ax.set_ylabel('Reads/Millions Mapped')
    ax.set_xticklabels(['DMSO','Nutlin 1hr','Nutlin 3hr'])
    bp = ax.boxplot([[m[3] for m in a],[n[3] for n in b],[m[3] for l in c]],patch_artist=True)
    format_boxplot(bp)

    a = w2.map(Db, c='4', o='sum',null=0)
    b = w2.map(N1b, c='4', o='sum',null=0)
    c = w2.map(N3b, c='4', o='sum',null=0)
    ax2 = F.add_subplot(131)
    ax2.set_title('Wave2')
    ax2.set_ylabel('Reads/Millions Mapped')
    ax2.set_xticklabels(['DMSO','Nutlin 1hr','Nutlin 3hr'])
    bp2 = ax2.boxplot([[m[3] for m in a],[n[3] for n in b],[m[3] for l in c]],patch_artist=True)
    format_boxplot(bp2)

    a = w3.map(Db, c='4', o='sum',null=0)
    b = w3.map(N1b, c='4', o='sum',null=0)
    c = w3.map(N3b, c='4', o='sum',null=0)
    ax3 = F.add_subplot(131)
    ax3.set_title('Wave3')
    ax3.set_ylabel('Reads/Millions Mapped')
    ax3.set_xticklabels(['DMSO','Nutlin 1hr','Nutlin 3hr'])
    bp3 = ax3.boxplot([[m[3] for m in a],[n[3] for n in b],[m[3] for l in c]],patch_artist=True)
    format_boxplot(bp3)


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