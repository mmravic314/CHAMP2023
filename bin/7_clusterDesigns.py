# clean up redundant fasta files
# cluster sequences with appropriate BLOSOM/PAM matrix, average linkage
# python ~/CHAMP/bin/clust_fasta.py ~/CHAMP/EGFR/bbDesign_CTERMrnd2/match_1/redund_m1.fasta ~/CHAMP/EGFR/bbDesign_CTERMrnd2/match_1/redund_m1-dMat.pkl

# 

import sys, os, numpy as np, Levenshtein

from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from itertools import combinations
from collections import defaultdict



#################### Global definitions, classes, functions  ####

class Model:

	def __init__(self,seq, model, ps, score, hbonds):
		self.seq 	= seq
		self.model 	= model
		self.ps 	= float( ps )
		self.score	= float( score )
		self.hbonds	= int( hbonds ) 
		self.arrayIn= 10000
		self.clusLab= 10000

	def __repr__(self):
		return '#%s %.3f %.1f#' % ( self.model, self.ps, self.score )

dmatrix_HASH  = {}
dmatrix_HASH[30] = matlist.blosum30; dmatrix_HASH[35] = matlist.blosum35; dmatrix_HASH[40] = matlist.blosum40
dmatrix_HASH[45] = matlist.blosum45; dmatrix_HASH[50] = matlist.blosum50; dmatrix_HASH[55] = matlist.blosum55
dmatrix_HASH[60] = matlist.blosum60; dmatrix_HASH[65] = matlist.blosum65; dmatrix_HASH[70] = matlist.blosum70
dmatrix_HASH[75] = matlist.blosum75; dmatrix_HASH[80] = matlist.blosum80; dmatrix_HASH[85] = matlist.blosum85
dmatrix_HASH[90] = matlist.blosum90; dmatrix_HASH[95] = matlist.blosum95; dmatrix_HASH[100] = matlist.blosum100


def score_Matrix(matrix, AA_tuple):
	if AA_tuple not in matrix:
		return matrix[(tuple(reversed(AA_tuple)))]
	else:
		return matrix[pair]

#################### Global defs end  ###########################


#################  MAIN  ########################################

## Read in file
models		 = {}
sequences 	 = defaultdict(list)

with open( sys.argv[1] ) as fin:
	
	for i in fin:
		#print(i)
		if i[0] == '>':

			model, ps, score, sc, hbonds, unsatPolars = tuple( i.split()[1:] )
			continue
		if len(i.rstrip()) <1: continue
		if len(i) > 20: 
			seq 	= i.rstrip()
			coreSeq = seq[3:-4]			# disclude variations in the last 3 polar AA at N & C-termini

			design 			= Model( seq, model, ps, score, hbonds )
			models[model] 	= design

			sequences[coreSeq].append( design )



# figure out the approximate mean sequence identity of the input fasta, comparing maximum 1000 design combos
seq_IDs = []
for a,b in combinations( sequences.keys(), 2 ):
	#print(a,b, Levenshtein.ratio(a,b) )
	seq_IDs.append( Levenshtein.ratio(a,b) )

	if len(seq_IDs) > 1000: 
		break

mean_seqID 		= np.mean(seq_IDs)* 100
global_seq_id 	= 5.0 * round(float(mean_seqID) / 5.0)

print('\n\n\t\tGlobal sequence ID approximated to be %.f, so use BLOSOM%d for clustering \n' % ( round(mean_seqID), round(global_seq_id) ) ) 
	
blosom = dmatrix_HASH[ global_seq_id  ]



ssd_blosum_similarity = []

for a,b in combinations( sequences.keys(), 2 ):

	score = 0
	for i, j in zip(a, b):
		try:
			score += blosom[(i, j)]
		except:
			score += blosom[(j, i)]
	ssd_blosum_similarity.append(score)


ssd_blosum_distance = np.max(ssd_blosum_similarity) - ssd_blosum_similarity 


'''
#### Write the distance matrix for clustering to pickle file once only
import cPickle as pic
pklPath = sys.argv[2]
if not os.path.exists( pklPath ):
	dMat = fillMatrix( unique )
	pic.dump( dMat, open( pklPath, 'wb' ) )
else:
	dMat = pic.load( open( pklPath, 'rb' ) )
'''

###### CLUSTERING #####

from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, cophenet
import time
ts = time.time()



linkage_mat = linkage( ssd_blosum_distance, method='ward')
clust 		= fcluster( linkage_mat, 1.15 )		# this clustering metric/threshold is going to be example specific

clusters 	= defaultdict(list) 
for n, s in enumerate(sequences.keys()):
	#print(s, clust[n]) 
	clusters[ clust[n] ].append( s )


print( '\n\tUnique clusters:', len(set(clust)) )
print ('\t\t# model-ID packstat dG_complex' )
	
for c, seq_list in sorted( clusters.items() ):

	num_designs = np.sum( [ len(sequences[s]) for s in seq_list ] )
	print ('>>> cluster %d with %s unique core TM sequences; %d design models/trajectories' % ( c, len(seq_list), np.sum(num_designs)   ) )
	
	#merge all designs, allow sorting by energy
	step = 0
	for s in sorted( seq_list, key=lambda x: np.min( [y.score for y in sequences[x] ]) ):

		print ('\t', s)
		for model in sorted( sequences[s], key=lambda x: x.score ):

			if step ==0:
				print('\t\t',model, '\t << best design in cluster')
			else:
				print('\t\t',model)
			step +=1

	'''
	all_models = []
	for s in sorted( seq_list ):
		print ('\t', s)
		for model in sorted( sequences[s], key=lambda x: x.score ):
			print('\t\t',model)
	'''

	print()
	#print('\n')

'''


linkage_mat = linkage( ssd_blosum_distance, method='complete')
for k in np.arange(1.0, 1.2, 0.002):

	clust 		= fcluster( linkage_mat, k )

	#print np.max(clust), k
	#print clust
	# organize clusters
	clusters 	= defaultdict(list)
	stp = 0
	for s in sequences.keys():

		clusters[ clust[stp] ].append( unique[s] )
		unique[s].clusLab = clust[stp]

		stp +=1
	for k, v in clusters.items():
		repSeq = sorted( v, key= lambda x: x.ps , reverse=True)[0]
		print (k, 'n=%d' % (len(v)), repSeq.ps, repSeq.seq)

'''




