__author__ = 'Jonathan Rubin'

import os
import math
import pybedtools
from pybedtools import BedTool

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def remove_repeats(repeats,intervalfile,filedir):
    a = BedTool(intervalfile)
    b = BedTool(repeats)

    (a-b).saveas(filedir + 'Non_repeat_peaks.bed')
    (a+b).saveas(filedir + 'Repeat_peaks.bed')

if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'

    repeats = filedir + 'simple_repeat_regions.bed'
    intervalfile = filedir + 'All_Peaks.sorted.merge.bed'
    remove_repeats(repeats,intervalfile,filedir)
