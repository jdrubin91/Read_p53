'''
	Question:
		If you are a peak that is unique to either Nutlin
		1 hour relative to DMS0 or a unique peak in Nutlin 
		3 hour relative to 1 hour ARE YOU CLOSER than 
		you would expect by chance to a peak that is conserved
		between the two populations.
'''
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np,math,time
from collections import defaultdict
'''
	had to add this peak class to keep track of the type of peak
'''

class peak:
	def __init__(self, chrom, x, ID):
		self.ID 		= ID #either 0 if it is a conserved peak or 1 if it is a unique peak	or None
		self.chrom 	= chrom
		self.x 		= x


class bed_file:
	def __init__(self, FILE=None,DICT=None,LST=None):
		if FILE is not None:
			self.F 		= FILE
			self._load()
		if DICT is not None:
			self.G 		= DICT
		if LST is not None:
			self.LST 	= LST
			self._load_LST()
	def _load_LST(self):
		self.G 	= defaultdict(list)
		for pk in self.LST:
			self.G[pk.chrom].append((pk.x, pk))
		for chrom in self.G:
			self.G[chrom].sort()
			self.G[chrom]=[j for i,j in self.G[chrom]]
	def _load(self):
		self.G 	= dict()
		with open(self.F)as FH:
			for line in FH:
				chrom,start,stop = line.strip("\n").split("\t")[:3]
				x 	= 0.5*(float(start) + float(stop))
				if chrom not in self.G:
					self.G[chrom]=list()
				self.G[chrom].append((x, peak(chrom, x, None)))
		for chrom in self.G:
			self.G[chrom].sort()
			self.G[chrom]=[j for i,j in self.G[chrom]]
	def _stats(self,SHOW=False,peak_lbl=False):
		if not peak_lbl:
			self.d 		= np.log10([ abs(l[i+1].x-l[i].x+1) for l in self.G.values() for i in range(0,len(l)-1) ])
		else:
			self.d 		= list()
			for chrom in self.G:
				pks 	= self.G[chrom]
				i,N 	= 0, len(self.G[chrom])
				while i < N:
					j 	= i+1
					while j < N and pks[i].ID!=pks[j].ID:
						j+=1
					k 	= i-1
					while k > 0 and pks[i].ID!=pks[k].ID:
						k-=1
					if j < N and k > 0 and i!=k and i!=j:
						d1,d2 	= abs(pks[i].x-pks[j].x),abs(pks[i].x-pks[k].x)
						self.d.append(min(d1,d2)+1)
					i+=1
			self.d 		= np.array(self.d)
			self.d 		= np.log10(self.d[~np.isnan(self.d) & ~np.isinf(self.d)] )

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
def intersect(BedFileA,BedFileB, window=500,BOTH=False):
	chroms 	= set(BedFileA.G.keys()).intersection(set(BedFileB.G.keys()))
	I 			= dict()
	for chrom in chroms:
		j,N 		= 0,len(BedFileB.G[chrom])
		I[chrom] = list()
		for A in BedFileA.G[chrom]:
			while j < N and BedFileB.G[chrom][j].x + window < A.x - window:
				j+=1
			if j < N and BedFileB.G[chrom][j].x - window < A.x + window:
				'''Intersection'''
				I[chrom].append(peak(chrom,A.x,0))					
			elif j < N and BOTH:
				'''Unique to A; these are presumably the newly created p53 ChIP-sites'''
				I[chrom].append(peak(chrom,A.x,1))
	return bed_file(DICT=I)

def simulate(base_case, observed_case,simN=pow(10,4),ax=None,title="",R=None ):
	stat 		= "mean"
	'''compute the nearest neighbor statistics for the test case'''
	observed_case._stats(peak_lbl=True)

	BDS 		= list()
	pks 		= np.array([pk.ID for chrom in observed_case.G for pk in observed_case.G[chrom] ])
	motifs 	= np.array([(chrom, pk.x) for chrom in base_case.G for pk in base_case.G[chrom] ])
	S,N 		= sum(pks),len(pks)
	mN 		= len(motifs)
	'''run simN simulates'''
	for i in range(simN):
		'''
			1. randomly select N motifs
		'''
		idx 		= np.random.choice(mN,S)
		'''
			2. randomly label S of these motifs as conserved 
				and N-S of these motifs as new
		'''
		nucleates= np.random.permutation(np.append(np.zeros((N-S,)),np.ones((S,))))
		'''
			then compute stats
		'''
		nPeaks 	= [peak(chrom, float(x), nucleates[j]) for j,(chrom, x) in enumerate(motifs[idx]) ]
		cMotifs 	= bed_file(LST=nPeaks)
		cMotifs._stats(peak_lbl=True)
		'''
			remember for the histogram
		'''
		BDS.append(getattr(cMotifs, stat))
	BINS 	= 60
	if ax is None:
		F 	 	= plt.figure(tight_layout=True,figsize=(15,10))
		ax 	= plt.gca()

	ax.set_title(title)
	counts,edges 	= np.histogram(BDS,bins=BINS)

	edges = 0.5*(edges[1:] + edges[:-1])
	width = (max(edges)-min(edges)) / BINS
	ax.bar(edges, counts,width=width,label="Simulated Data",edgecolor="white")
	ax.plot([getattr(observed_case,stat ),getattr(observed_case,stat )],[0,max(counts)],
					color="red", label="Observed Data",lw=2)
	ax.set_ylabel("Frequency")
	if R is None:
		ax.legend(loc="best")
		R 	= getattr(observed_case,stat)-0.4,max(edges)
	ax.set_xlim(R)
	sns.despine()
	return R



def main():
	FILE_DMSO 	= "../files2/DMSO1Hr_peaks.merge_200.bed"
	FILE_Nut1h 	= "../files2/Nutlin1Hr_peaks.merge_200.bed"
	FILE_Nut3h 	= "../files2/Nutlin3Hr_peaks.merge_200.bed"
	FILE_motif 	= "../files2/pwmscan_hg19_25506_4823.bed"

	DMSO 			= bed_file(FILE=FILE_DMSO)
	Nut1h 		= bed_file(FILE=FILE_Nut1h)
	Nut3h 		= bed_file(FILE=FILE_Nut3h)
	motif 		= bed_file(FILE=FILE_motif)


	'''grab only peaks over p53 sites
	'''
	I_DMSO 		= intersect(DMSO,motif,BOTH=False)
	I_Nut1h 		= intersect(Nut1h,motif,BOTH=False)
	I_Nut3h 		= intersect(Nut3h,motif,BOTH=False)

	'''label the peaks as either conserved (0) or new (1)
	'''
	WAVE2 		= intersect(I_Nut1h,I_DMSO, BOTH=True)
	WAVE3 		= intersect(I_Nut3h,I_Nut1h, BOTH=True)



	sns.set(font_scale=2.0);sns.set_style("ticks")
	F 				= plt.figure(figsize=(15,10),tight_layout=True)
	ax1 			= F.add_subplot(2,1,1)
	ax2 			= F.add_subplot(2,1,2)
	R 				= simulate(motif,WAVE2,ax=ax1, title="Nutlin 1 hour unique Peaks\n(relative to DMSO)")
	R 				= simulate(motif,WAVE3,ax=ax2,R=R, title="Nutlin 3 hour unique Peaks\n(relative to Nutlin 1 hour)")
	ax2.set_xlabel("Average Distance between Conserved and New Peak")

	plt.savefig("../figures/NearestNeighborSimulationsConservedVsNewPeaks.png")
	plt.show()



if __name__ == "__main__":
	main()



