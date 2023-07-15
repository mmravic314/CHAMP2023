
# analyze a directory full of PDB files output from rosetta. 
# extract the sequences, (single sequence can have multiple molecular models aka trajectories)
# cluster by BLOSOM matrix

import sys, os, numpy as np, subprocess as sp, seaborn as sns, shutil, pickle as pkl
from prody import *
from collections import defaultdict, Counter
#from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt



# print the compiled sequences to a fasta file, output to the same directory path as supplied (input 1)
# grab a series of Rosetta calculated molecular features (energies & filters) from model files
# compile features in header of sequences in fasta file 
# & also plot distribution of features for all sequences
# cluster the sequences, print out clusters - fasta file has assigned clusters in header

# input 1) 	Path to directory full of rosetta output PDB files, 
# input 2)  chain (single letter) or chains (comma separated list) within the model that are being designed which will be analyzed
# input 3)  Path to the 'top' designs passing all filters we calculated

## usage
# python ../../../bin/analyze_cluster_designsV1.py setup3_output/ X
# python ../../../bin/analyze_cluster_designsV1.py setup4_output/ X


#################

UnNatAA={}
UnNatAA["ALA"] = 'A'; UnNatAA["CYS"] = 'C'; UnNatAA["ASP"] = 'D'; UnNatAA["GLU"] = 'E'; UnNatAA["PHE"] = 'F';
UnNatAA["GLY"] = 'G'; UnNatAA["HIS"] = 'H'; UnNatAA["ILE"] = 'I'; UnNatAA["LYS"] = 'K';
UnNatAA["LEU"] = 'L'; UnNatAA["MET"] = 'M'; UnNatAA["ASN"] = 'N'; UnNatAA["PRO"] = 'P'; UnNatAA["GLN"] = 'Q';
UnNatAA["ARG"] = 'R'; UnNatAA["SER"] = 'S'; UnNatAA["THR"] = 'T'; UnNatAA["VAL"] = 'V'; UnNatAA["TRP"] = 'W'; UnNatAA["TYR"] = 'Y';
UnNatAA['ABA'] = 'A'; UnNatAA['CSO'] = 'C'; UnNatAA['CSD'] = 'C'; UnNatAA['CME'] = 'C';
UnNatAA['OCS'] = 'C'; UnNatAA["HSD"] = 'H'; UnNatAA['KCX'] = 'K'; UnNatAA['LLP'] = 'K';
UnNatAA['MLY'] = 'K'; UnNatAA['M3L'] = 'K'; UnNatAA['MSE'] = 'M'; UnNatAA['PCA'] = 'P'; UnNatAA['HYP'] = 'P';
UnNatAA['SEP'] = 'S'; UnNatAA['TPO'] = 'T'; UnNatAA['PTR'] = 'Y'


class Model:
		def __init__(self, name, chains, filepath):
			self.name		= name
			self.chains		= chains
			self.seqS		= ''
			self.filePath	= filepath
			self.cluster 	= 1000
			self.passed_filters	= 0
			self.score_feat = {'dG_separated':0 , 'delta_unsatHbonds':10, 'hbonds_int':0, 'holes-ABX_interface':10 , 'sc_value':0}
			self.perRes_scores 	= {} 

		def __repr__( self ):
			return self.name

		'''
		def Calc_num_lowProbability_rotamers(self, verbose=1):
			residues, fa_dun_rotatmer = [], []
			record = '\n'
			for r, score_list in sorted(self.perRes_scores.items(), key=lambda x: int(x[0].split('_')[1]) ):
				residues.append(r)
				rotamer_score = float( score_list['fa_dun'] )
				## mods for MET & Tyr   (always high)
				if r[:3] in ['MET', 'TYR' ]:
					rotamer_score = round( rotamer_score/1.5, 4) 

				fa_dun_rotatmer.append( rotamer_score )
				if rotamer_score > 2.25:
					record += '\t>>%s %.3f\n' % ( r, rotamer_score )
			self.num_lowProbability_rotamers = len([ x for x in fa_dun_rotatmer if x > 2.25])
			return record
		'''
#################

################# Main #################

input_dir	 = sys.argv[1]
dir_files	 = os.listdir(input_dir)
output_dir	 = sys.argv[3]


pdb_files = [x for x in dir_files if x[-3:] == 'pdb']

## chain(s) to analyze
chains2check = sys.argv[2]
if len(chains2check) > 1:
	if ',' not in chain2check:
		print('\n\t ERROR:  Format for arg2 (specifying chains) is one chain (single letter) or chains (comma separated list)\n')
		sys.exit()
	chains2check = chains2check.split(',')



#### Parse PDB files, get sequence of binder and rosetta model stats/scores
# save in pickle so parseing not repeated
chainSeqDict, all_designs, good_design_seqs = {}, [], defaultdict(list)
score_features = {'dG_separated':[] , 'delta_unsatHbonds':[], 'hbonds_int':[], 'packstat':[] , 'sc_value':[], 'complex_normalized': [], 'total_energy':[] }
res_indices, chainsChecked, residues_in_design = np.array([]), [], np.array([])


for design in pdb_files:
	path 		= os.path.join(input_dir, design)
	design_id	= '-'.join( design[:-4].split('_')[-2:] )
	pdb 		= parsePDB( path )
	desi 		= Model( design_id,  chains2check , path)

	for ch in chains2check:
		selection 		= pdb.select('ca chain %s' % ch )
		desi.seqS 	 	= ''.join(  selection.getSequence()  )
		desi.design_id  = design_id
		print( design_id, selection.getSequence() )

	all_designs.append(desi)

	inF 		= open(path, 'r')
	pdb_text 	= inF.readlines()
	perRes_score_parse_flg = 0
	for l in pdb_text[900:]:
		if len(l) < 1: continue
		line = l.rsplit()
		if len(line) < 2: continue

		if line[0] in score_features.keys():
			score_features[line[0]].append( float( line[1] ) )
			desi.score_feat[line[0]] = float( line[1] )
			continue

		if 'dG_separated' == line[0]:
			dG_complex 		= float( line[1] )
			desi.dG_complex = float( line[1] )
#			print('\n\n*******',l)

		if 'delta_unsatHbonds' == line[0]:
			unsatHbonds 		= int( line[1] )
			desi.unsatHbonds 	= int( line[1] )
#			print('\n\n*******',l)

		if 'hbonds_int' == line[0]:
			num_HB 		= int( line[1] )
			desi.num_HB = int( line[1] )
#			print('\n\n*******',l)

		if 'packstat' == line[0]:
			holes_interface 		= float( line[1] )
			desi.holes_interface 	= float( line[1] )
#			print('\n\n*******',l)

		if 'sc_value' == line[0]:
			ShapeComplement_interface 		= float( line[1] )
			desi.ShapeComplement_interface 	= float( line[1] )
#			print('\n\n*******',l)

		if 'pose' == line[0]:
			total_energy	= float( line[-1] )
			desi.total_energy	= float( line[-1] )
			score_features[ 'total_energy' ].append( float( line[-1] ) )
			desi.score_feat[ 'total_energy' ] = float( line[-1] )

		
		if '#END_POSE_ENERGIES_TABLE' == line[0]:
			break

		'''
		elif '#BEGIN_POSE_ENERGIES_TABLE' == line[0]:
			perRes_score_parse_flg = 1
		elif 'label' == line[0]:
			score_term_keys = line[1:]
		else:
			pass


		## gather all per-residue scores from energy table of chain of interest
		if perRes_score_parse_flg:
			aa_check 	= line[0].split('_')
			aa_basename	= aa_check[0][:3]
			if aa_basename in UnNatAA.keys():
				
				aa_num = int( aa_check[-1] )
				if aa_num in res_indices:

					if ':' in aa_check[0]:
						rosetta_termini_long = aa_check[0].split(':')[1] 
						mod = rename_Rosetta_termini[ rosetta_termini_long ]
						if 'Nterm' in rosetta_termini_long:
							aa_basename = mod + aa_basename
						else:
							aa_basename += mod 

					residue_ID  	= aa_basename + '_' +  str(aa_num)
					perRes_scores 	= line[1:]
					desi.perRes_scores[residue_ID] = {}
					for term, score in zip( score_term_keys, perRes_scores ):
						desi.perRes_scores[residue_ID][term] = score

				else:
					continue
		'''

	## Knobs to holes analysis



	# assign to design object
#	score_features = [ dG_complex, unsatHbonds, num_HB, holes_interface, ShapeComplement_interface ]
#	print (design_id, score_features, seq)

	inF.close()






## plot distribution of score features/values:
score_Axis 		= {'dG_separated':'Association Energy (REU)' , 'delta_unsatHbonds':"Unsaturated Buried Polar Groups", 'hbonds_int':"# Interface H-bonds", 'packstat':"RosettaHoles1" ,'total_energy':"Total Energy (REU)",  'complex_normalized':"Interface Energy normalized (REU)", 'sc_value':"Interface Shape Complementarity"}
score_limits 	= {'dG_separated':0 , 'holes-ABX_interface':10 , 'sc_value':0}
score_limit_flg = {'dG_separated':0 , 'holes-ABX_interface':0 , 'sc_value':2}

plt.rcParams["font.family"] = "Arial"
print('\n\n\t Summarizing scores for %d design trajectories' % len(pdb_files) )

for k,v in score_features.items():
	mean, stdev = round( np.mean(v), 2), round( np.std(v), 2 )

	print ('\n>>%s\t\taverage score: %.2f,\tstdev %.2f,\tmax: %.2f,\tmin: %.2f' %  (k, mean, stdev, np.max(v), np.min(v) ),  )

	if k == 'dG_separated':
		#limit = mean - 1*stdev 			# limit for dG interface score is mean score 
		score_limits[k] = mean
		print('limit: %0.2f' % mean)

	elif k == 'sc_value':				# limit for surface complementarity is 0.75
		score_limits[k] = mean
		print('limit: %0.2f' % mean )

	elif k == 'packstat':
		top10_packstat	= round( np.percentile(v, 90), 2)
		print ('packstat top 10 percent cutoff',  top10_packstat )
		score_limits[k] = top10_packstat

	else:
		print()


	## plot
	bin_number = int( round(  len(v)**0.5 , 0) )

	sns.histplot(data=v, bins=bin_number, stat='density', alpha= 0.6, kde=True, color='y', 
    edgecolor='black', linewidth=0.5, line_kws=dict(color='black', alpha=0.5,linewidth=2, label='KDE'))
	plt.gca().get_lines()[0].set_color('black') # manually edit line color due to bug in sns v 0.11.0
	plt.legend(frameon=False)

	plt.xlabel( score_Axis[k], fontsize=20)
	plt.xticks(  fontsize=16)

	output_image_path = os.path.join(  sys.argv[1],  '%s.png' % k )
	if not os.path.exists( output_image_path ):
		plt.savefig( output_image_path , dpi=300, bbox_inches='tight' )
	#plt.show()
	plt.clf()


##### plot packstat versus total energy, interface energy, and interface energy normalized by residue

plt.plot( score_features['total_energy'], score_features['packstat'], 'bo', markersize=3 )
plt.xlabel('Total Energy (REU)', fontsize=22)
plt.ylabel('RosettaHoles1 (Packstat)', fontsize=22)
plt.xticks( fontsize=16)
plt.yticks( fontsize=12)
plt.savefig( os.path.join(  sys.argv[1],  'TotalREU_vs_packstat.png' ) , dpi=300, bbox_inches='tight' )
#plt.show()
plt.clf()



plt.plot( score_features['dG_separated'], score_features['packstat'], 'bo', markersize=3 )
plt.xlabel('Association Energy (REU)', fontsize=22)
plt.ylabel('RosettaHoles1 (Packstat)', fontsize=22)
plt.xticks( fontsize=16)
plt.yticks( fontsize=12)
plt.savefig( os.path.join(  sys.argv[1],  'InterfaceREU_vs_packstat.png' ) , dpi=300, bbox_inches='tight' )
#plt.show()
plt.clf()




## gather the design models which pass the sequence filte, move them to new file
## Number of H-bonds must be >2 unsaturated buried groups
## SC fraction > 0.65
## dG separated at least > median
## 
score_keys = [ 'packstat', 'dG_separated', 'sc_value', 'hbonds_int', 'delta_unsatHbonds' ]

good_design_list 	= []
all_designs_fasta 	= ''
good_designs_fasta 	= ''
for d in all_designs:

	print ( d.seqS, desi.design_id, [x for x in d.score_feat.values() ] )
	all_designs_fasta += '> %s %.3f %.1f %.2f %d %d\n%s\n\n' % (d.design_id, d.score_feat['packstat'], d.score_feat['dG_separated'], d.score_feat['sc_value'], d.score_feat['hbonds_int'], d.score_feat['delta_unsatHbonds'],  d.seqS )
	for k, v in d.score_feat.items():
		#print (d.seqS, v)

		
		if k == 'packstat':
				if v > score_limits[k]:
					d.passed_filters += 1
					#print('pass packstat', v)

		elif k == 'dG_separated':
				if v < score_limits[k]:
					d.passed_filters += 1
					#print('pass dG', v)
		elif k == 'sc_value':
				if v < 0.65:
					d.passed_filters -=10
					#print('fail SC', v)
		else:
				pass
	### this setion for trying to score H-bonds & polar groups packing in protein interface but not Hbonding: "unsatisfied"

	# any designs with >3 polar buried unsatisfied groups rejected
	if d.score_feat['delta_unsatHbonds'] > 3: 			
			d.passed_filters -=10
			#print('fail unsatHbonds', k)

	# any designs with 0 polar buried unsatisfied groups are okay or with number of H-bonds at least double
	if d.score_feat['delta_unsatHbonds'] > (d.score_feat['hbonds_int'] / 2.0):
			d.passed_filters -=10
			#print('fail unsatHbonds/Hbond ratio', k)


	if d.passed_filters >= 2:
		# primary sequence filters:
		seq 		= d.seqS

		print ('\n\n*** good design\t\t>> %s \t|\t%s ' %  (d.seqS , d.filePath) ) 
		print ('#\t', ' '.join(score_keys))
		print ('\n', d, d.passed_filters, sorted( d.score_feat.items() ) )
		print (os.path.basename(d.filePath ) )
		score_summary = [ d.score_feat[x] for x in score_keys ]

		'''
		## any rotamer outliers?
		rotamer_details = d.Calc_num_lowProbability_rotamers()
		score_summary.append(d.num_lowProbability_rotamers)
		print ( '   # Questionable rotamers:', d.num_lowProbability_rotamers, rotamer_details )


		## DSSP calculation to potentially filter extra helix defects
		pout = sp.Popen('~/bin/mkdssp -i %s' % d.filePath, shell=1, stdout=sp.PIPE, stderr=sp.PIPE)
		SSE_helix_assignments, doubleHbonded_bb_atoms = [], []
		for line in pout.stdout.readlines():
			line 	= line.decode('UTF-8') 
			chain 	= line[11:12].rstrip()
	
			if chain not in chains2check:
				continue 
			resid 	= line[7:10].rstrip()
			resname = line[13:14].rstrip()
			resi_id = ( chain, resid )

			if resi_id in residues_in_design[4:-4]:
				H_bond_code = line[14:20].strip().split()
				
				if H_bond_code[0] == 'H':
					SSE_helix_assignments.append(1)
				else:
					pass

				if 'X' in H_bond_code[1]:
					doubleHbonded_bb_atoms.append(1)
				elif '>' in H_bond_code[1] and '<' in H_bond_code[1]:
					doubleHbonded_bb_atoms.append(1)
				else:
					pass

				#print(resi_id, resname,  H_bond_code, '--- %d %d' % (len(SSE_helix_assignments), len(doubleHbonded_bb_atoms) ), '\n')

		#d.non_helix_SSE_calls = len( residues_in_design[4:-4] ) - len( SSE_helix_assignments )
		#d.broken_bb_hbonds 	= len( residues_in_design[4:-4] ) - len( doubleHbonded_bb_atoms )
		#print ("   # residues\t non-helix assigments: %d\t backbone H-bonding disrupted %d" % (d.non_helix_SSE_calls, d.broken_bb_hbonds )	 )
				#print(line[:90], '\n')
		#score_summary.append( d.non_helix_SSE_calls )
		#score_summary.append( d.broken_bb_hbonds )
		'''
		good_design_list.append(d)
		good_designs_fasta += '> %s %.3f %.1f %.2f %d %d\n%s\n\n' % (d.design_id, d.score_feat['packstat'], d.score_feat['dG_separated'], d.score_feat['sc_value'], d.score_feat['hbonds_int'], d.score_feat['delta_unsatHbonds'],  d.seqS )
	
		#newname
		d.new_design_label = 'design_%s.pdb' % ( d.design_id )
		d.score_summary = score_summary
		#score_summary.append( d.new_design_label )

		# holes, dG, SC, Hbond_int, del_UnSat_HBs, lowProb_Rots, helixBreak_SSE, helixBreak_HB, newPath
		#good_design_seqs[d.seqS].append( score_summary )

		outfile = os.path.join(output_dir, d.new_design_label  )
		if not os.path.exists(outfile):
			shutil.copy( d.filePath,  outfile )


print ( '\n\n--------------\nfound %d designs passing filters, out of %d total' % ( len(good_design_list), len(all_designs) ) )

# write fasta files for all designs and good designs

all_fasta = open( os.path.join( sys.argv[1], 'all_designs.fasta' ), 'w')
good_fasta = open( os.path.join( sys.argv[3], 'good_designs.fasta' ), 'w')

all_fasta.write(all_designs_fasta)
good_fasta.write(good_designs_fasta)





