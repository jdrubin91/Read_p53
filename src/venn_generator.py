__author__ = 'Jonathan Rubin'

import os
import sys
import pybedtools
from pybedtools import BedTool
from matplotlib import pyplot as plt
# from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})
import numpy as np
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def run(DMSO,Nutlin1,Nutlin3,figuredir):
    D = BedTool(DMSO)
    N1 = BedTool(Nutlin1)
    N3 = BedTool(Nutlin3)

    a = len(D) - len(D+N1-N3) - len(D+N3-N1) - len(N3+N1+D)
    b = len(N1) - len(N1+D-N3) - len(N1+N3-D) - len(N3+N1+D)
    c = len(D+N1-N3)
    d = len(N3) - len(N3+D-N1) - len(N3+N1-D) - len(N3+N1+D)
    e = len(N3+D-N1)
    f = len(N3+N1-D)
    g = len(N3+N1+D)
    

    plt.figure()
    v = venn3(subsets=(a, b, c, d, e, f, g), set_labels = ('DMSO', 'Nutlin 1hr', 'Nutlin 3hr'))
    # v.get_patch_by_id('100').set_alpha(1.0)
    # v.get_patch_by_id('100').set_color('white')
    # v.get_label_by_id('100').set_text('Unknown')
    # xy = v.get_label_by_id('100').get_position()
    # v.get_label_by_id('DMSO').set_text('')
    # c[0].set_lw(1.0)
    # c[0].set_ls('dotted')
    plt.title("Peak Overlaps")
    plt.text('DMSO', xy=xy)
    plt.savefig(figuredir + 'venn_diagram.png',dpi=1200)
    # plt.show()

def promoter_overlap(promoters, DMSO, Nutlin1, Nutlin3, figuredir):
    D = BedTool(DMSO)
    N1 = BedTool(Nutlin1)
    N3 = BedTool(Nutlin3)
    p = BedTool(promoters)

    AB = len(D+p)
    aB = 0
    Ab = len(D) - AB
    plt.figure()
    v = venn2(subsets(1,2,3), set_labels = ('DMSO','Promoters'))


def example():
    plt.figure(figsize=(4,4))
    v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'))
    print v
    v.get_patch_by_id('100').set_alpha(1.0)
    v.get_patch_by_id('100').set_color('white')
    v.get_label_by_id('100').set_text('Unknown')
    v.get_label_by_id('A').set_text('Set "A"')
    c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
    c[0].set_lw(1.0)
    c[0].set_ls('dotted')
    plt.title("Sample Venn diagram")
    plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
                 ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    plt.show()

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
    # run(DMSO,Nutlin1,Nutlin3,figuredir)

    promoters = file2dir + 'refGene_hg19_TSS.bed'

    promoter_overlap(promoters, DMSO, Nutlin1, Nutlin3, figuredir)

