__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import pybedtools
from pybedtools import BedTool
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os
import numpy as np
from numpy import random as rn
import load
from scipy import stats
import math
import time

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def find_nearest(array):
    l = len(array)
    result = list()
    for i in range(l):
        if i == 0:
            result.append(math.fabs(array[1][0] - array[i][1]))
        elif i == l-1:
            result.append(math.fabs(array[i][0] - array[i-1][0]))
        else:
            result.append(min(math.fabs(array[i][1] - array[i+1][0]),math.fabs(array[i][0] - array[i-1][1])))

    return result

def run(DMSO,Nutlin1,Nutlin3,P53,figuredir,file2dir):
    D = BedTool(DMSO)
    N1 = BedTool(Nutlin1)
    N3 = BedTool(Nutlin3)
    P = BedTool(P53).cut([0,1,2])
    
    start = time.time()

    w1 = (D+P).saveas(file2dir + 'Wave1.bed')
    w1rand = BedTool([P[i] for i in rn.randint(0,len(P),len(w1))]).sort()
    w2 = (N1+P-D).saveas(file2dir + 'Wave2.bed')
    w2rand = BedTool([P[i] for i in rn.randint(0,len(P),len(w2))]).sort()
    w3 = (N3+P-N1-D).saveas(file2dir + 'Wave3.bed')
    w3rand = BedTool([P[i] for i in rn.randint(0,len(P),len(w3))]).sort()

    # w1 = D+P
    # print rn.shuffle(list(P))[:len(w1)]
    # w1rand = BedTool(rn.shuffle(P)[:len(w1)]).sort()
    # w2 = N1+P-D
    # w2rand = BedTool(rn.shuffle(P)[:len(w2)]).sort()
    # w3 = N3+P-N1-D
    # w3rand = BedTool(rn.shuffle(P)[:len(w3)]).sort()

    end = time.time()

    print(end - start)

    a = w2.closest(w1, d=True)
    b = w2rand.closest(w1rand, d=True)
    c = w3.closest(w2, d=True)
    d = w3rand.closest(w2rand, d=True)


    w21 = list()
    w2r1r = list()
    w32 = list()
    w3r2r = list()

    for x in a:
        try:
            w21.append(math.log(float(x[-1]),10))
        except:
            w21.append(0)

    for x in b:
        try:
            w2r1r.append(math.log(float(x[-1]),10))
        except:
            w2r1r.append(0)

    for x in c:
        try:
            w32.append(math.log(float(x[-1]),10))
        except:
            w32.append(0)

    for x in d:
        try:
            w3r2r.append(math.log(float(x[-1]),10))
        except:
            w3r2r.append(0)

    print len(w21),len(w2r1r),len(w32),len(w3r2r)

    # w21 = [math.log(float(x[-1])) for x in a if float(x[-1]) != 0]
    # w2r1r = [math.log(float(x[-1])) for x in b if float(x[-1]) != 0]
    # w32 = [math.log(float(x[-1])) for x in c if float(x[-1]) != 0]
    # w3r2r = [math.log(float(x[-1])) for x in d if float(x[-1]) != 0]

    # print stats.ks_2samp(w21, w2r1r)
    # print stats.ks_2samp(w32, w3r2r)

    F = plt.figure()
    ax1 = F.add_subplot(221)
    ax1.set_title('Wave2 to Wave1 (pval: ' + str(stats.ks_2samp(w21, w2r1r)[1]) + ')')
    ax1.set_ylabel('Count')
    ax1.set_xlabel('Log 10 Distance (bp)')
    ax1.hist(w21,bins=np.arange(0, 8 + 0.2, 0.2),alpha=0.5,color='green')
    # ax1.set_xlim([0,500000])
    # ax1.set_ylim([0,600])
    # ax1.hist(w21,bins=np.arange(0, 18 + 0.2, 0.2))
    # ax1.set_xscale('log')

    # ax2.F.add_subplot(222)
    # ax2.set_title('Wave2rand to Wave1rand')
    # ax2.set_ylabel('Count')
    # ax2.set_xlabel('Distance (bp)')
    ax1.hist(w2r1r,bins=np.arange(0, 8 + 0.2, 0.2),alpha=0.5,color='red')
    ax1.legend(['Observed','Expected'],loc='upper left')
    # ax2.set_xlim([0,500000])
    # ax2.set_ylim([0,600])
    # ax1.hist(w2r1r,bins=np.arange(0, 18 + 0.2, 0.2))
    # ax2.set_xscale('log')

    ax2 = F.add_subplot(222)
    ax2.set_title('Cumulative distribution function')
    ax2.set_ylabel('CDF')
    ax2.set_xlabel('Log 10 Distance (bp)')
    # Use the histogram function to bin the data
    counts, bin_edges = np.histogram(w21, bins=np.arange(0, 18 + 0.2, 0.2), normed=True)
    counts_r, bin_edges_r = np.histogram(w2r1r, bins=np.arange(0, 18 + 0.2, 0.2), normed=True)
    # Now find the cdf
    cdf = np.cumsum(counts)
    cdf_r = np.cumsum(counts_r)
    # And finally plot the cdf
    plt.plot(bin_edges[1:], cdf, color='green')
    plt.plot(bin_edges_r[1:], cdf_r, color='red')
    ax2.legend(['Observed','Expected'],loc='upper left')

    ax3 = F.add_subplot(223)
    ax3.set_title('Wave3 to Wave2 (pval: ' + str(stats.ks_2samp(w32, w3r2r)[1]) + ')')
    ax3.set_ylabel('Count')
    ax3.set_xlabel('Log 10 Distance (bp)')
    ax3.hist(w32,bins=np.arange(0, 8 + 0.2, 0.2),alpha=0.5,color='green')
    # ax3.set_xlim([0,500000])
    # ax3.set_ylim([0,3500])
    # ax3.hist(w32,bins=np.arange(0, 18 + 0.2, 0.2))
    # ax3.set_xscale('log')

    ax4 = F.add_subplot(224)
    ax4.set_title('Cumulative distribution function')
    ax4.set_ylabel('CDF')
    ax4.set_xlabel('Log 10 Distance (bp)')
    # Use the histogram function to bin the data
    counts, bin_edges = np.histogram(w32, bins=np.arange(0, 18 + 0.2, 0.2), normed=True)
    counts_r, bin_edges_r = np.histogram(w3r2r, bins=np.arange(0, 18 + 0.2, 0.2), normed=True)
    # Now find the cdf
    cdf = np.cumsum(counts)
    cdf_r = np.cumsum(counts_r)
    # And finally plot the cdf
    plt.plot(bin_edges[1:], cdf, color='green')
    plt.plot(bin_edges_r[1:], cdf_r, color='red')
    ax4.legend(['Observed','Expected'],loc='upper left')

    # ax4 = F.add_subplot(224)
    # ax4.set_title('Wave3rand to Wave2rand')
    # ax4.set_ylabel('Count')
    # ax4.set_xlabel('Distance (bp)')
    ax3.hist(w3r2r,bins=np.arange(0, 8 + 0.2, 0.2),alpha=0.5,color='red')
    ax3.legend(['Observed','Expected'],loc='upper left')
    # ax4.set_xlim([0,500000])
    # ax4.set_ylim([0,3500])
    # ax3.hist(w3r2r,bins=np.arange(0, 18 + 0.2, 0.2))
    # ax4.set_xscale('log')



    plt.savefig(figuredir + 'Cluster_analysis.png', dpi=1200)

def nearest_neighbor(DMSO,Nutlin1,Nutlin3,P53):
    a = load.load_bed_points(DMSO)
    b = load.load_bed_points(Nutlin1)
    c = load.load_bed_points(Nutlin3)
    P = load.load_bed_points(P53)

    arand = set()
    while len(arand) < len(a):
        arand.add(rn.randint(0,len(P)))
    arand = list(arand)

    brand = set()
    while len(arand) < len(b):
        brand.add(rn.randint(0,len(P)))
    brand = list(brand)

    crand = set()
    while len(arand) < len(c):
        crand.add(rn.randint(0,len(P)))
    crand = list(crand)

    print len(arand),len(brand),len(crand)

    # F = plt.figure()
    # ax1 = F.add_subplot(221)
    # plt.hist([l for l in find_nearest(a)],bins=100)
    # # ax1.set_xlim([0,500])
    # ax1.set_title('DMSO')
    # ax1.set_ylabel('Count')
    # ax1.set_xlabel('Peak Distance')
    # ax2 = F.add_subplot(222)
    # plt.hist([m for m in find_nearest(b)],bins=100)
    # # ax2.set_xlim([0,500])
    # ax2.set_title('Nutlin 1hr')
    # ax2.set_ylabel('Count')
    # ax2.set_xlabel('Peak Distance')
    # ax3 = F.add_subplot(223)
    # plt.hist([n for n in find_nearest(c)],bins=100)
    # # ax3.set_xlim([0,500]) 
    # ax3.set_title('Nutlin 3hr')
    # ax3.set_ylabel('Count')
    # ax3.set_xlabel('Peak Distance')
    # # plt.savefig(figuredir + 'nearest_neighbor.png', dpi=1200)
    # plt.show()
    


if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    file2dir = parent_dir(homedir) + '/files2/'
    figuredir = parent_dir(homedir) + '/figures/'

    DMSO = file2dir + 'DMSO1Hr_peaks.merge_200.bed'
    Nutlin1 = file2dir + 'Nutlin1Hr_peaks.merge_200.bed'
    Nutlin3 = file2dir + 'Nutlin3Hr_peaks.merge_200.bed'
    P53 = file2dir + 'HO_P53_HUMAN.H10MO.B.pval-6.bed'
    P53 = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/P53_fimo.bed'
    P53 = '/scratch/Users/joru1876/P53_fimo.bed'
    run(DMSO,Nutlin1,Nutlin3,P53,figuredir,file2dir)

    P53 = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/P53_fimo.no_header.bed'
    P53 = '/scratch/Users/joru1876/P53_fimo.no_header.bed'
    # nearest_neighbor(DMSO,Nutlin1,Nutlin3,P53)
