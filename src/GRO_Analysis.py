__author__ = 'Jonathan Rubin'

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os
import math
import load

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

def run(true_neg,Nutlin1):
    x = load.load_bed_full_intervals(true_neg)
    y = load.load_bed_full_intervals(Nutlin1)

    print len(y)
    for interval in x:
        for i in range(len(y)-1):
            # print x[i]
            chrom1,start1,stop1 = interval
            chrom2,start2,stop2 = y[i][:3]
            if chrom1 == chrom2 and start1 == start2 and stop1 == stop2:
                y.pop(i)

    print len(y)



if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'

    true_neg = filedir + 'true_negatives.txt'
    Nutlin1 = filedir + 'Nutlin1Hr_peaks.merge.200.bed'
    run(true_neg,Nutlin1)


