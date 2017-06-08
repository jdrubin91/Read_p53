'''
	Questions: 
		are p53 ChIP-seq peaks (over motif sites)
		close to themselves then what we would 
		expect from the random motifs

'''

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np,math,time
'''
	bed file class
'''	

class peak:
	def __init__(self, chrom, x, ID):
		self.ID 		= ID #either 0 if it is a conserved peak or 1 if it is a unique peak	
		self.chrom 	= chrom
		self.x 		= x


class bed_file:
	def __init__(self, FILE=None,DICT=None,TUPLE=None):
		if FILE is not None:
			self.F 		= FILE
			self._load()
		if DICT is not None:
			self.G 		= DICT
		if TUPLE is not None:
			self.tuple 	= TUPLE
			self._load_tuple()
		self._stats()
	def _load_tuple(self):
		self.G 	= dict()
		for chrom,val in self.tuple:
			if chrom not in self.G:
				self.G[chrom]=list()
			self.G[chrom].append(float(val))
		for chrom in self.G:
			self.G[chrom].sort()
	def _load(self):
		self.G 	= dict()
		with open(self.F)as FH:
			for line in FH:
				chrom,start,stop = line.strip("\n").split("\t")[:3]
				x 	= 0.5*(float(start) + float(stop))
				if chrom not in self.G:
					self.G[chrom]=list()
				self.G[chrom].append(x)
		for chrom in self.G:
			self.G[chrom].sort()
	def _stats(self,SHOW=False):
		self.d 		= np.log10([ abs(l[i+1]-l[i]+1) for l in self.G.values() for i in range(0,len(l)-1) ])
		self.mean 	= np.mean(self.d)

		self.std 	= np.std(self.d)
		self.median = np.median(self.d)
		self.N 		= len(self.d)+len(self.G.keys())
		if SHOW:
			plt.hist(self.d,bins=int(math.sqrt(len(self.d))))
			plt.show()
		self.tPeaks = np.array([(l,x) for l in self.G.keys() for x in self.G[l]])
	def get_random_peaks(self,N):
		return self.tPeaks[np.random.choice(self.N, N)]
def intersect(BedFileA,BedFileB, window=200):
	chroms 	= set(BedFileA.G.keys()).intersection(set(BedFileB.G.keys()))
	I 			= dict()
	for chrom in chroms:
		j,N 		= 0,len(BedFileB.G[chrom])
		I[chrom] = list()
		for x in BedFileA.G[chrom]:
			while j < N and BedFileB.G[chrom][j] + window < x - window:
				j+=1
			if j < N and BedFileB.G[chrom][j] - window < x + window:
				'''Intersection'''
				I[chrom].append(x)
	return bed_file(DICT=I)

def simulate(base_case,test_case,N=pow(10,4),ax=None,title="",R=None):
	'''run N simulates'''
	BDS 			= list()
	for i in range(N):

		tuples 	= base_case.get_random_peaks(test_case.N)
		bd 	 	= bed_file(TUPLE=tuples)
		BDS.append(bd)
	'''viz'''
	a 		= "mean"
	BINS 	= 60
	if ax is None:
		F 	 	= plt.figure(tight_layout=True,figsize=(15,10))
		ax 	= plt.gca()
	ax.set_title(title)
	counts,edges 	= np.histogram([getattr(b,a) for b in BDS],bins=BINS)

	edges = 0.5*(edges[1:] + edges[:-1])
	width = (max(edges)-min(edges)) / BINS
	ax.bar(edges, counts,width=width,label="Simulated Data",edgecolor="white")
	ax.plot([getattr(test_case,a ),getattr(test_case,a )],[0,max(counts)],
					color="red", label="Observed Data")
	ax.set_ylabel("Frequency")
	if R is None:
		ax.legend(loc="best")
		R 	= getattr(test_case,a)-0.2,max(edges)
	ax.set_xlim(R)
	sns.despine()
	if ax is None:

		plt.show()
	return R
def viz_all_three():
	FILE_DMSO 	= "../files2/DMSO1Hr_peaks.merge_200.bed"
	FILE_Nut1h 	= "../files2/Nutlin1Hr_peaks.merge_200.bed"
	FILE_Nut3h 	= "../files2/Nutlin3Hr_peaks.merge_200.bed"
	FILE_motif 	= "../files2/pwmscan_hg19_25506_4823.bed"

	DMSO 			= bed_file(FILE=FILE_DMSO)
	Nut1h 		= bed_file(FILE=FILE_Nut1h)
	Nut3h 		= bed_file(FILE=FILE_Nut3h)
	motif 		= bed_file(FILE=FILE_motif)

	I_DMSO 		= intersect(DMSO,motif)
	I_Nut1h 		= intersect(DMSO,Nut1h)
	I_Nut3h 		= intersect(DMSO,Nut3h)

	'''viz all of them'''
	sns.set(font_scale=2.0),sns.set_style("ticks")
	F 				= plt.figure(figsize=(10,10),tight_layout=True)
	ax1,ax2,ax3 = F.add_subplot(3,1,1),F.add_subplot(3,1,2),F.add_subplot(3,1,3)
	'''run simulations'''

	R 				= simulate(motif,I_DMSO,ax=ax1,title="DMSO")
	simulate(motif,I_Nut1h,ax=ax2,title="Nutlin 1 Hour",R=R)
	simulate(motif,I_Nut3h,ax=ax3,title="Nutlin 3 Hour",R=R)

	ax3.set_xlabel("Average Distance to Nearest")
	plt.savefig("../figures/NearestNeighborSimulations.png")
	plt.show()




def main():
	viz_all_three()


if __name__ == "__main__":
	main()