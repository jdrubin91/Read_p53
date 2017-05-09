__author__ = 'Jonathan Rubin'

import os
import sys
import pybedtools
from pybedtools import BedTool

#The purpose of this script is to perform MEME on different bed file intervals. For the purposes of
#p53 ChIP with Nutlin, this script will perform MEME on: True negatives, Global true positives, exclusive DMSO (truepos) sites,
#exclusive 1hr Nutlin (truepos) sites, and exclusive 3hr Nutlin (truepos) sites

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

def run_MEME(fastafile,outdir,scriptdir):
    os.system("qsub -v fastafile=" + fastafile + " output=" + outdir + " " + scriptdir + "MEME_run.sh")

def run(true_neg,true_pos,DMSO,Nutlin1,Nutlin3,hg19fasta,outdir,filedir,scriptdir):
    neg = BedTool(true_neg)
    negfasta = filedir + true_neg + '.fasta'
    pos = BedTool(true_pos)
    posfasta = filedir + true_pos + '.fasta'
    D = BedTool(DMSO)
    D_unique = filedir + DMSO + '_unique.fasta'
    N1 = BedTool(Nutlin1)
    N1_unique = filedir + Nutlin1 + '_unique.fasta'
    N3 = BedTool(Nutlin3)
    N3_unique = filedir + Nutlin3 + '_unique.fasta'
    hg19 = BedTool(hg19fasta)

    neg.sequence(fi=hg19fasta).saveas(negfasta)
    pos.sequence(fi=hg19fasta).saveas(posfasta)

    run_MEME(negfasta,outdir,scriptdir)
    run_MEME(posfasta,outdir,scriptdir)

    D-neg-N1-N3.sequence(fi=hg19fasta).saveas(D_unique)
    run_MEME(D_unique,outdir,scriptdir)
    
    N1-neg-D-N3.sequence(fi=hg19fasta).saveas(N1_unique)
    run_MEME(N1_unique,outdir,scriptdir)

    N3-neg-D-N1.sequence(fi=hg19fasta).saveas(N3_unique)
    run_MEME(N3_unique,outdir,scriptdir)



if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File/Figure directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'
    scriptdir = parent_dir(homedir) + '/scripts/'


    true_neg = filedir + 'true_negatives.txt'
    true_pos = filedir + 'All_Peaks.sorted.merge.bed.true_positive.bed'
    DMSO = filedir + 'DMSO_peaks.merge.200.bed'
    Nutlin1 = filedir + 'Nutlin1Hr_peaks.merge.200.bed'
    Nutlin3 = filedir + 'Nutlin3Hr_peaks.merge.200.bed'
    hg19fasta = '/scratch/Users/joru1876/hg19_reference_files/hg19_all.fa'
    outdir = parent_dir(homedir) + '/MEME/'
    run(true_neg,true_pos,DMSO,Nutlin1,Nutlin3,hg19fasta,outdir,filedir,scriptdir)