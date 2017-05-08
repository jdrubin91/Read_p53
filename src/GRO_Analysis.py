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
    bp = ax.boxplot([math.log(float(n[3])/float(m[3]),2) for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.'], bins =100)
    format_boxplot(bp)
    plt.savefig(figuredir + 'GRO_Analysis_Fold_Change.png', dpi=1200)

    outfile = open(filedir + 'false_positives_GRO-Seq_fold_change_bp.bed','w')
    for interval in [m[:3] for m,n in zip(a,b) if m[3] != 0 and n[3] != 0 and m[3] != '.' and n[3] != '.' and float(n[3])/float(m[3]) < 1]:
        outfile.write('\t'.join(interval) + '\n')




if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'

    Nutlin1 = filedir + 'Nutlin1Hr_peaks.merge.200.bed.true_positive.bed'
    GRODMSO = '/scratch/Users/joru1876/Allen2014_NutlinGRO/DMSO1Hr.mp.reflected.sorted.merge.BedGraph'
    GRONutlin = '/scratch/Users/joru1876/Allen2014_NutlinGRO/Nutlin1Hr.mp.reflected.sorted.merge.BedGraph'
    run(Nutlin1,GRODMSO,GRONutlin,figuredir,filedir)


