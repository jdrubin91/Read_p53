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
    os.system("qsub -v fastafile=" + fastafile + ",output=" + outdir + " " + scriptdir + "MEME_run.sh")

def run(DMSO,Nutlin1,Nutlin3,hg19fasta,outdir,filedir,scriptdir):
    # neg = BedTool(true_neg)
    # negfasta = true_neg + '.fasta'
    # pos = BedTool(true_pos)
    # posfasta = true_pos + '.fasta'
    D = BedTool(DMSO)
    D_unique = DMSO + '_unique.fasta'
    N1 = BedTool(Nutlin1)
    N1_unique = Nutlin1 + '_unique.fasta'
    N3 = BedTool(Nutlin3)
    N3_unique = Nutlin3 + '_unique.fasta'
    hg19 = BedTool(hg19fasta)

    print len(D), len(N1-D), len(N3-D-N1)

    # neg.sequence(fi=hg19fasta).save_seqs(negfasta)
    # pos.sequence(fi=hg19fasta).save_seqs(posfasta)

    # run_MEME(negfasta,outdir+'neg/',scriptdir)
    # run_MEME(posfasta,outdir+'pos/',scriptdir)

    # (D-neg).sequence(fi=hg19fasta).save_seqs(D_unique)
    # run_MEME(D_unique,outdir+'DMSO/',scriptdir)
    
    # (N1-neg-D).sequence(fi=hg19fasta).save_seqs(N1_unique)
    # run_MEME(N1_unique,outdir+'Nutlin1Hr/',scriptdir)

    # (N3-neg-D-N1).sequence(fi=hg19fasta).save_seqs(N3_unique)
    # run_MEME(N3_unique,outdir+'Nutlin3Hr/',scriptdir)

    # (D).sequence(fi=hg19fasta).save_seqs(D_unique)
    # run_MEME(D_unique,outdir+'DMSO/',scriptdir)
    
    # (N1-D).sequence(fi=hg19fasta).save_seqs(N1_unique)
    # run_MEME(N1_unique,outdir+'Nutlin1Hr/',scriptdir)

    # (N3-D-N1).sequence(fi=hg19fasta).save_seqs(N3_unique)
    # run_MEME(N3_unique,outdir+'Nutlin3Hr/',scriptdir)



if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File/Figure directory
    filedir = parent_dir(homedir) + '/files/'
    file2dir = parent_dir(homedir) + '/files2/'
    figuredir = parent_dir(homedir) + '/figures/'
    scriptdir = parent_dir(homedir) + '/scripts/'


    # true_neg = filedir + 'true_negatives.txt.wcut_8.bed'
    # true_pos = filedir + 'All_Peaks.sorted.merge.bed.true_positive.bed.wcut_8.bed'
    # DMSO = filedir + 'DMSO_peaks.merge.200.bed.wcut_8.bed'
    # Nutlin1 = filedir + 'Nutlin1Hr_peaks.merge.200.bed.wcut_8.bed'
    # Nutlin3 = filedir + 'Nutlin3Hr_peaks.merge.200.bed.wcut_8.bed'
    DMSO = file2dir + 'DMSO1Hr_peaks.merge_200.bed'
    Nutlin1 = file2dir + 'Nutlin1Hr_peaks.merge_200.bed'
    Nutlin3 = file2dir + 'Nutlin3Hr_peaks.merge_200.bed'
    hg19fasta = '/scratch/Users/joru1876/hg19_reference_files/hg19_all.fa'
    outdir = parent_dir(homedir) + '/MEME2/'
    run(DMSO,Nutlin1,Nutlin3,hg19fasta,outdir,filedir,scriptdir)