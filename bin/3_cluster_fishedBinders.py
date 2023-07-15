from numpy import linalg as LA
import sys, os, numpy as np, time, pickle as pkl, math, shutil
from prody import *
from collections import defaultdict, Counter
from itertools import combinations
from scipy.cluster.hierarchy import linkage, fcluster


### 
# input 1:  directory of single helix PDBs to cluster
# input 2:  parallel or anti-parallel relative to target frag 1 'chain A', as single charater p or a (only lowercase is accepted)
# input 3: (optional):  cluster linkage critera (Å cut-off)

# Output:  clusters will be written in same directory as path given, input 1

## >> python ../bin/ClusterTM_domains.py fished_mid15_frag/fishedTM_Binds_interface/ a

############# < I/O ###########

path2PDBs  		= sys.argv[1]
par_v_aparKey 	= sys.argv[2]
par_v_aparDict 	= {'a':'ap', 'p':'pa'}

if par_v_aparKey not in par_v_aparDict.keys():
	print('\nOnly a or p are accepted as options for parallel vs antiparallel flag, input argument 2. Abort...\n')
	sys.exit()
else:
	par_v_aparFlg = par_v_aparDict[par_v_aparKey]
 
if len(sys.argv) > 3:
	clust_rmsd = float( sys.argv[3] )
else:
	clust_rmsd = 2.5	

############# I/O > ###########

# inputs must be numpy arrays
def rmsd_calc(A, B):
    d = A - B
    l = len(A)
    return np.sqrt((d*d).sum()/l)

dirname 			= os.path.dirname( path2PDBs )
eligible_windows	= {}

# Split PDB into all possible 12 residue windows. 
for pdb in os.listdir( path2PDBs ):
	if not 'pdb' in pdb: 
		continue

	flag = pdb[0:1]

	print (pdb, flag, par_v_aparFlg[0])

	if flag != par_v_aparFlg[0]: 
		print ('hey')
		continue


	else:

		pdb_path 	= os.path.join( path2PDBs , pdb )
		tmp_pdb 	= parsePDB(pdb_path)
		tmp_pdb_ca  = tmp_pdb.select('ca')
		eligible_windows[pdb] = []

		for i in np.arange( len( tmp_pdb_ca) - 14 ):
			resi_range 	= tmp_pdb_ca.getResnums()[i:i+13] 
			bb_window 	= tmp_pdb.select('bb resnum %s' % ' '.join([str(x) for x in resi_range]))
			resi_ID		= '%d_%d' % (resi_range[0],resi_range[-1])

			eligible_windows[pdb].append( bb_window  )

num_pdbs = len( eligible_windows.keys() )
print( 'clustering', num_pdbs, 'fished PDBs' )
TM_list_sorted = sorted( eligible_windows.keys() )

### pairwise alignment - compare all windows ( maybe some branch & bound in here )
# store as a pair of lists, one list to ID the two PDB windows compared & list of all RMSD - turn to Pandas later
## do pairwise comparison (incomplete non-square matrix)

# skip all this work if an existing RMSD calculation file present
path2_pkl_RMSDdata = os.path.join(path2PDBs, 'rmsd_data_%a.pkl' % par_v_aparKey)
path2_pkl_RMSDdata2= os.path.join(path2PDBs, 'rmsd_datav2_%a.pkl'% par_v_aparKey)


if os.path.exists( path2_pkl_RMSDdata ) and os.path.exists( path2_pkl_RMSDdata2 ):
	rmsd_data 					= pkl.load( open( path2_pkl_RMSDdata , 'rb') )
	rmsd_data_v2, comparisons 	= tuple(pkl.load( open( path2_pkl_RMSDdata2 , 'rb') ) )
	print('... retrieved RMSD matrix from previous run', path2_pkl_RMSDdata)

else:
	step = 1
	#rmsd_data, rmsd_index 	= [], []
	rmsd_data 			= {}
	time_init 			= time.time()
	num_calcs			= math.comb(num_pdbs, 2) #/ 2*(  math.factorial(num_pdbs - 2) )
	rmsd_data_v2 	= []
	comparisons		= []

	for ab in combinations( TM_list_sorted, 2):
		print(ab)

		a,b = ab
		comparisons.append(ab)
		best_matchRMSD		= 100

		if step in np.arange(0, num_calcs, 1000):
			print ('-- step %d of %d -- RMSD matrix calculation...  %.3f sec taken' % (step, num_calcs, time.time() - time_init) )
		#print(a,b, '-- step %d of %d -- RMSD matrix calculation...  %.3f sec taken' % (step, num_pdbs*(num_pdbs-1), time.time() - time_init) )
		for win_a in eligible_windows[a]:
			for win_b in eligible_windows[b]:

				resi_a 		= list( set( win_a.getResnums() ) )
				resi_b 		= list( set( win_b.getResnums() ) )
				resi_ID_a	= '%s-%d_%d' % (a[:-4], resi_a[0],resi_a[-1])
				resi_ID_b	= '%s-%d_%d' % (b[:-4], resi_b[0],resi_b[-1])

				#if sorted( [ resi_ID_a, resi_ID_b ] ) in rmsd_index: continue

				if len(win_a.getCoords()) != len(win_b.getCoords()): 	# hack skip fix for some coordinates missing backbone atoms?
					rmsd  = 10
				else:
					rmsd = rmsd_calc( win_a.getCoords(), win_b.getCoords() )
					
				if rmsd < best_matchRMSD:
					best_matchRMSD = rmsd
				#rmsd_data.append( rmsd )
				#rmsd_index.append( sorted([resi_ID_a, resi_ID_b] ) )
				#print ( resi_ID_a, resi_ID_b, rmsd )
				rmsd_data[ '_|_'.join( sorted( [resi_ID_a, resi_ID_b] ) ) ] = rmsd
		rmsd_data_v2.append(best_matchRMSD)

		step += 1
	print ('time elapsed: %.2f' % (time.time() - time_init) )

	## sort then store rmsdata
	rmsd_data = dict(sorted(rmsd_data.items(), key=lambda item: item[1]))
	#pkl.dump( rmsd_data , open( path2_pkl_RMSDdata , 'wb') )
	#pkl.dump( (rmsd_data_v2, comparisons), open( path2_pkl_RMSDdata2 , 'wb') )



################################
####   Attempt with heirachical clusterings
#### this results in redundency - because input PDB list has redundencies (for now)
################################
ssd_rmsd = np.array(rmsd_data_v2)		## pseudo-square distance matrix in reduced array format 
init_T = time.time()

linkMat 	= linkage( ssd_rmsd , method='complete', metric='euclidean')
h_clust 	= fcluster( linkMat, clust_rmsd, criterion='distance')
numClust 	= len( set(h_clust) )
print (round( time.time() - init_T, 3), 's to complete clustering')

print ('RMS cutoff at %.2f, Unique clusters found:' % clust_rmsd, numClust, '\n')

# reorganize & label clusters
TM_clusters 		= defaultdict(list)
TM_array_index_key 	= defaultdict(list)

for n in np.arange(len(TM_list_sorted)):
	TM_clusters[ h_clust[n] ].append( TM_list_sorted[n] )
	TM_array_index_key[ h_clust[n] ].append( TM_list_sorted[n] )
	#print ( n, h_clust[n], TM_list_sorted[n] )



singles_dir 	= os.path.join(  path2PDBs, '%s_singletons' % (par_v_aparFlg) )
if not os.path.exists( singles_dir  ):
			os.mkdir(singles_dir )

#within clusters, print RMSD matrix
# and fasta files with a local sequence alignment to arbitraty 'first' cluster member
step = 1	
for k,v in sorted( TM_clusters.items(), key=lambda x: len(x[1]), reverse=True ):
	clust_size = len( v )
	if clust_size > 1:
		clust_dir = os.path.join( path2PDBs, '%s_cluster_%d_n-%d' % (par_v_aparFlg, step, clust_size ) )
		print('\n> %s_node %d: \t members: %d\n'  % (par_v_aparFlg, step, clust_size ) )
	
		print (clust_dir, '____')
		if not os.path.exists( clust_dir ):
			os.mkdir(clust_dir)
		pdb_path  	=  os.path.join( path2PDBs , TM_list_sorted[k] )
		shutil.copy( pdb_path , clust_dir  )

		for p in v:
			print('\t', p)
			pdb_path  	=  os.path.join( path2PDBs , p )
			shutil.copy( pdb_path , clust_dir  )

	else:
		print('> singleton: %s' % TM_list_sorted[k])
		pdb_path  	=  os.path.join( path2PDBs , TM_list_sorted[k] )
		shutil.copy( pdb_path , singles_dir  )
	step +=1

# 56 + 33 + 19 + 7 + 7 + 6 + 30