__author__ = 'Jonathan Rubin'

def load_bed_points(bed):
    x = list()
    with open(bed) as F:
        for line in F:
            line = line.strip().split()
            val1 = float(line[1])
            val2 = float(line[2])
            x.append((val1,val2))
    return x

def load_counts_file(file1):
    x = list()
    with open(file1) as F:
        for line in F:
            try:
                val = float(line.strip().split()[-1])
            except ValueError:
                val = 0
            x.append(val)
    return x

def load_counts_file_full_intervals(file1):
    x = list()
    with open(file1) as F:
        for line in F:
            line = line.strip().split()
            chrom,start,stop = line[:-1]
            try:
                val = float(line[-1])
            except ValueError:
                val = 0
            x.append((chrom,start,stop,val))
    return x

def load_bed_full_intervals(file1):
    x = list()
    with open(file1) as F:
            for line in F:
                chrom,start,stop = line.strip().split()
                x.append((chrom,start,stop))
    return x