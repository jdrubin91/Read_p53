__author__ = 'Jonathan Rubin'

# from pybedtools import BedTool
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os
import math
from scipy import spatial
import load

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

def peak_distance(bed1,bed2,bed3,figuredir):

    a = load.load_bed_points(bed1)
    b = load.load_bed_points(bed2)
    c = load.load_bed_points(bed3)

    F = plt.figure()
    ax1 = F.add_subplot(221)
    plt.hist([l for l in find_nearest(a) if l < 500],bins=100)
    ax1.set_xlim([0,500])
    ax1.set_title('DMSO')
    ax1.set_ylabel('Count')
    ax1.set_xlabel('Peak Distance')
    ax2 = F.add_subplot(222)
    plt.hist([m for m in find_nearest(b) if m < 500],bins=100)
    ax2.set_xlim([0,500])
    ax2.set_title('Nutlin 1hr')
    ax2.set_ylabel('Count')
    ax2.set_xlabel('Peak Distance')
    ax3 = F.add_subplot(223)
    plt.hist([n for n in find_nearest(c) if n < 500],bins=100)
    ax3.set_xlim([0,500])
    ax3.set_title('Nutlin 3hr')
    ax3.set_ylabel('Count')
    ax3.set_xlabel('Peak Distance')
    # plt.savefig(figuredir + 'peak_distance.png', dpi=1200)
    plt.show()

def width_depth_analysis(bed1,bed2,bed3,figuredir):
    a = load.load_bed_points(bed1)
    b = load.load_bed_points(bed2)
    c = load.load_bed_points(bed3)

    F = plt.figure()
    ax1 = F.add_subplot(221)
    plt.hist([l[1] - l[0] for l in a],bins=100)
    # ax1.set_xlim([0,1000])
    ax1.set_title('DMSO')
    ax1.set_ylabel('Count')
    ax1.set_xlabel('Peak Width')
    ax2 = F.add_subplot(222)
    plt.hist([l[1] - l[0] for l in b],bins=100)
    # ax2.set_xlim([0,1000])
    ax2.set_title('Nutlin 1hr')
    ax2.set_ylabel('Count')
    ax2.set_xlabel('Peak Width')
    ax3 = F.add_subplot(223)
    plt.hist([l[1] - l[0] for l in c],bins=100)
    # ax3.set_xlim([0,2000])
    ax3.set_title('Nutlin 3hr')
    ax3.set_ylabel('Count')
    ax3.set_xlabel('Peak Width')
    plt.savefig(figuredir + 'peak_width.png', dpi=1200)
    # plt.show()

def fold_change_analysis(DMSO,Nutlin1,Nutlin3):
    x = load.load_counts_file(DMSO)
    y = load.load_counts_file(Nutlin1)
    z = load.load_counts_file(Nutlin3)


    print len([b/a for a,b in zip(x,y) if a != 0 and b/a < 1])
    print len([b/a for a,b in zip(x,z) if a != 0 and b/a < 2])
    # F = plt.figure()
    # ax1 = F.add_subplot(121)
    # ax1.set_title('DMSO vs. Nutlin1hr')
    # ax1.set_ylabel('Count')
    # ax1.set_xlabel('Fold Change')
    # ax1.hist([b/a for a,b in zip(x,y) if a != 0 and b/a < 12],bins=100)
    # ax2 = F.add_subplot(122)
    # ax2.set_title('DMSO vs. Nutlin3hr')
    # ax2.set_ylabel('Count')
    # ax2.set_xlabel('Fold Change')
    # ax2.hist([b/a for a,b in zip(x,z) if a != 0 and b/a < 12],bins=100)
    # # format_boxplot(bp1)
    # # plt.show()
    # plt.savefig(figuredir + 'fold_change_hist.png', dpi=1200)


    # ax1 = F.add_subplot(111)
    # ax1.set_title('p53 Peaks Fold Changes')
    # ax1.set_ylabel('Log2(Nutlin/DMSO)')
    # ax1.set_xticklabels(['1hr','3hr'])
    # bp1 = ax1.boxplot([[math.log(b/a,2) for a,b in zip(x,y) if b != 0 and a != 0],[math.log(d/c,2) for c,d in zip(x,z) if d != 0 and c!= 0]],patch_artist=True)
    # # print bp1
    # format_boxplot(bp1)
    # plt.savefig(figuredir + 'fold_change_bp.png', dpi=1200)
    # plt.show()

def false_positive_overlap(DMSO,Nutlin1,Nutlin3,filedir):
    x = load.load_counts_file_full_intervals(DMSO)
    y = load.load_counts_file_full_intervals(Nutlin1)
    z = load.load_counts_file_full_intervals(Nutlin3)

    # print len([b[3]/a[3] for a,b in zip(x,y) if a[3] != 0 and b[3]/a[3] < 1])
    # print len([b[3]/a[3] for a,b in zip(x,z) if a[3] != 0 and b[3]/a[3] < 2])

    false1 = list()
    for a,b in zip(x,y):
        D = a[3] 
        N = b[3]
        try:
            fc = N/D
            if fc < 1:
                false1.append(a[:3])
        except:
            print "Error in following interval:", a, b

    false2 = list()
    for a,b in zip(x,z):
        D = a[3] 
        N = b[3]
        try:
            fc = N/D
            if fc < 2:
                false2.append(a[:3])
        except:
            print "Error in following interval:", a, b

    
    true_neg = list()
    for interval in false1:
        if interval in false2:
            true_neg.append(interval)

    print len(true_neg)

    # outfile = open(filedir + 'true_negatives.txt','w')
    # for interval in true_neg:
    #     outfile.write('\t'.join(interval) + '\n')

    for interval in true_neg:
        for i in range(len(x)-1):
            # print x[i]
            chrom1,start1,stop1 = interval
            chrom2,start2,stop2 = x[i][:3]
            if chrom1 == chrom2 and start1 == start2 and stop1 == stop2:
                x.pop(i)
                y.pop(i)
                z.pop(i)

    outfilex = open(filedir + DMSO.split('/')[-1].split('_')[0] + '_true_pos_counts.bed','w')
    outfiley = open(filedir + Nutlin1.split('/')[-1].split('_')[0] + '_true_pos_counts.bed','w')
    outfilez = open(filedir + Nutlin3.split('/')[-1].split('_')[0] + '_true_pos_counts.bed','w')
    for interval in x:
        outfilex.write('\t'.join(interval[:3]) + '\t' + str(interval[3]) + '\n')
    for interval in y:
        outfiley.write('\t'.join(interval[:3]) + '\t' + str(interval[3]) + '\n')
    for interval in z:
        outfilez.write('\t'.join(interval[:3]) + '\t' + str(interval[3]) + '\n')








if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'

    bed1 = filedir + 'DMSO_peaks.merge.200.bed'
    bed2 = filedir + 'Nutlin1Hr_peaks.merge.200.bed'
    bed3 = filedir + 'Nutlin3Hr_peaks.merge.200.bed'
    # peak_distance(bed1,bed2,bed3,figuredir)
    # width_depth_analysis(bed1,bed2,bed3,figuredir)
    DMSO = filedir + 'DMSO1Hr_AllPeaks_counts.bed'
    Nutlin1 = filedir + 'Nutlin1Hr_AllPeaks_counts.bed'
    Nutlin3 = filedir + 'Nutlin3Hr_AllPeaks_counts.bed'
    # fold_change_analysis(DMSO,Nutlin1,Nutlin3)
    false_positive_overlap(DMSO,Nutlin1,Nutlin3,filedir)

