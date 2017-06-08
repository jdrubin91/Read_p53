__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import os
import sys
import pybedtools
from pybedtools import BedTool
from matplotlib import pyplot as plt
# from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})
import numpy as np

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def run(promoters, DMSO, Nutlin1, Nutlin3, figuredir):
    D = BedTool(DMSO)
    N1 = BedTool(Nutlin1)
    N3 = BedTool(Nutlin3)
    p = BedTool(promoters)

    N = 3
    ind = np.arange(N)
    width = 0.35
    p1 = plt.bar(ind,(len(D)-len(D+p),len(N1)-len(N1+p),len(N3)-len(N3+p)),width,color='blue')
    p2 = plt.bar(ind,(len(D+p),len(N1+p),len(N3+p)),width,color='red')
    plt.ylabel('Peaks')
    plt.title('Peak Promoter Association')
    plt.xticks(ind+(width/2), ('DMSO', 'Nutlin 1Hr', 'Nutlin 3Hr'))
    plt.yticks(np.arange(0, 10000, 1000))
    plt.legend((p1[0], p2[0]), ('Non-promoter', 'Promoter associated'))
    plt.savefig(figuredir + 'peak_promoter_overlap.png',dpi=1200)

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
    promoters = file2dir + 'refGene_hg19_TSS.bed'

    run(promoters, DMSO, Nutlin1, Nutlin3, figuredir)