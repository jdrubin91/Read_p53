__author__ = 'Jonathan Rubin'

import os
import pybedtools
from pybedtools import BedTool

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def run(DMSO,Nutlin1,Nutlin3,genes,file2dir):
    D = BedTool(DMSO)
    N1 = BedTool(Nutlin1)
    N3 = BedTool(Nutlin3)
    g = BedTool(genes).sort().merge(c=5,o="distinct")

    w1 = D
    w2 = N1-D
    w3 = N3-N1-D

    outfile1 = open(file2dir + 'wave1_gene_targets.txt','w')
    outfile2 = open(file2dir + 'wave2_gene_targets.txt','w')
    outfile3 = open(file2dir + 'wave3_gene_targets.txt','w')

    for item in [str(m[3]) for m in g.intersect(w1,u=True)]:
        outfile1.write(item + '\n')
    for item in [str(m[3]) for m in g.intersect(w2,u=True)]:
        outfile2.write(item + '\n')
    for item in [str(m[3]) for m in g.intersect(w3,u=True)]:
        outfile3.write(item + '\n')


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
    genes = file2dir + 'refGene_hg19_TSS.bed'

    run(DMSO,Nutlin1,Nutlin3,genes,file2dir)
