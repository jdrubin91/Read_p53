__author__ = 'Jonathan Rubin'

from pybedtools import BedTool


#Takes as input two bed files and outputs to the specified location (and name). Options may be specified as a list ex. (options = [u=True,v=True,etc.])
def intersect(bed1,bed2,output,u=True):
    BedTool(bed1).intersect(BedTool(bed2),u=u).saveas(output)

#Takes as input a bed file and a bam file and outputs bedtools map (o='sum') to specified location (and name).
def count(bed,bam,output):
    Bedtool(bed).map(BedTool(bam),o='sum').saveas(output)

#Takes as input a bed file and outputs bedtools merge (with specified options) to output location (and name).
def merge(bed,output,d=0):
    BedTool(bed).merge(d=d).saveas(output)


