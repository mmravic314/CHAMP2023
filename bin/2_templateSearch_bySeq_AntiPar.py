# Marco Mravic DeGrado Lab UCSF Biophysics Nov 2015

## Search for user input seq pattern from database of membrane protein TM spans
## Then use membrane protein PDB to find other interacting TM span (tertiary/qarternary) contacting that sequence pattern
## Extract the structure, aligned to an input TM pdb containing the initial input sequence pattern

### **** NOTE ANTIPARALLEL VERSION OF CODE!!  Some changes to awareness/oreintation of membrane at lines XXX-XXX

# input 1: Directory of membrane protein PDB database list.txt, with one directory down used to find selectedfasta and selectedOPM
# input 2: target TM helix for alignment of modeled binder helices, must include INPUT SEQUENCE pattern
# input 3: reg ex to search sequences for those near the targeted motif on both input target & within TM-PDB
# input 4: output directory (write folder if doesn't already exist)

## Example command line
# python ./bin/2_templateSearch_bySeq_AntiPar.py ~/tmPDB_mar2020/40%_nonRedundant_finalPdbs.txt mEpoR_TM_threaded_0001.pdb [ASGC]\w\w\w\w\w\w[ASGC]\w\w\w\w\w\w[ASGC] smallX6small_x3_matches

### NOTE*****!!!! 
#
# At around line 140, there are two places where the assumption of parallel insertion into the membrane of both chains is made
# YOne should change this if the topology is not n2c for both chains... N-term to outer leaflet, etc ...
#
###
UnNatAA={}
UnNatAA["ALA"] = 'A'; UnNatAA["CYS"] = 'C'; UnNatAA["ASP"] = 'D'; UnNatAA["GLU"] = 'E'; UnNatAA["PHE"] = 'F';
UnNatAA["GLY"] = 'G'; UnNatAA["HIS"] = 'H'; UnNatAA["ILE"] = 'I'; UnNatAA["LYS"] = 'K';
UnNatAA["LEU"] = 'L'; UnNatAA["MET"] = 'M'; UnNatAA["ASN"] = 'N'; UnNatAA["PRO"] = 'P'; UnNatAA["GLN"] = 'Q';
UnNatAA["ARG"] = 'R'; UnNatAA["SER"] = 'S'; UnNatAA["THR"] = 'T'; UnNatAA["VAL"] = 'V'; UnNatAA["TRP"] = 'W'; UnNatAA["TYR"] = 'Y';
UnNatAA['ABA'] = 'A'; UnNatAA['CSO'] = 'C'; UnNatAA['CSD'] = 'C'; UnNatAA['CME'] = 'C';
UnNatAA['OCS'] = 'C'; UnNatAA["HSD"] = 'H'; UnNatAA['KCX'] = 'K'; UnNatAA['LLP'] = 'K';
UnNatAA['MLY'] = 'K'; UnNatAA['M3L'] = 'K'; UnNatAA['MSE'] = 'M'; UnNatAA['PCA'] = 'P'; UnNatAA['HYP'] = 'P';
UnNatAA['SEP'] = 'S'; UnNatAA['TPO'] = 'T'; UnNatAA['PTR'] = 'Y'

import sys, os, re, numpy as np
from collections import defaultdict
from prody import *

TMpdb_list 	= sys.argv[1]
targetTM 	= sys.argv[2]
regex 		= sys.argv[3]
fasta_dir 	= sys.argv[4]
TMpdb_dir 	= sys.argv[5]
outdir		= sys.argv[6]

inPDB 		= parsePDB( targetTM )
target_seq 	= ''.join([UnNatAA[x] for x in inPDB.select('ca').getResnames()] )

# find the input sequence pattern within the input TM domain
#	 locate the subset of amino acids / residue range of this region
target_match  = re.search(regex, target_seq)			## alpha-5
if target_match:
	target_pattern_indicies = ( target_match.start(), target_match.end() )
	print('input sequence pattern match:', target_match.start(), target_match.end(), target_seq[target_match.start() : target_match.end()]) 
else:
	print ('\n\n\t regex %s was not found within input target TM template %s \n' % (regex, target_seq) )
	sys.exit()

# locate the sequence pattern subset atoms in the input TM domain
target_matchAAs = inPDB.select('ca').copy()[target_match.start():target_match.end() ]
target_matchbb	= inPDB.select('bb').copy()[target_match.start():target_match.end() ]
# calculate the TM helical axis normal vector
dataMean        = target_matchbb.getCoords().mean(axis=0)
uu, dd, vv      = np.linalg.svd( target_matchbb.getCoords() - dataMean)
helix_axis_tar  = vv[0]





# locate filepaths of TM PDBs and fasta files to search the sequence pattern
#TMpdb_fa	= os.path.join( TMpdb_dir , 'selectedfasta')
#TMpdb_pdbs	= os.path.join( TMpdb_dir , 'selectedOPM')
if not os.path.exists(fasta_dir):
	print ('\n\n\t filepath %s does not exist\n' % fasta_dir)
if not os.path.exists(TMpdb_dir):
	print ('\n\n\t filepath %s does not exist\n' % TMpdb_dir)

print ('\n\n')
# look through fasta files of  non-redundent chains in membrane protein database to find TM segments with regex matches
loadedPDBs 		= {}
PDBs_w_match	= []
fasta_data_all	= {}
TM_seg_matches  = {}
for pdbfile in [x for x in os.listdir(TMpdb_list) if 'pdb' in x]:
	fastaPath = os.path.join( fasta_dir, '%s.fasta' % pdbfile[:-4] )
	#print(pdbfile, fastaPath)
	TM_segs 	= []
	match_Flg 	= 0
	chain_id 	= pdbfile[:-4]

	#enter fasta
	for line in open(fastaPath):
		if line[:2] == '>>': 
			TM_id = line[2:].rstrip()
		else:
			target_seq 		= line.rstrip()
			target_match  	= re.search(regex, target_seq)	
			TM_segs.append( ( TM_id, target_seq ) )
			if target_match:

				#if match found in file, start saving fasta data to file
				match_Flg 	= 1
				#print( '--->', TM_id, line.rstrip(), target_match.start(), target_match.end(), target_seq[target_match.start() : target_match.end()]) 
				match_index = (target_match.start(), target_match.end())

				TM_seg_matches[TM_id] = ( line.rstrip(), target_match.start(), target_match.end() )

				#TM_id2 		= TM_id.split()[0].split('_')[]
				pdb_id 		= TM_id.split()[0].split('_')[0]
				if pdb_id not in PDBs_w_match:
					PDBs_w_match.append( pdb_id.upper() )
				#res_range 	= [ int(x) for x in TM_id.split()[1].split('-') ]
				#ch 			= TM_id2[1]
				#chain_id 	= 'tm_%s_%s_%d' % ( TM_id2[0], ch , ord(ch) ) 

	if match_Flg:
		#print 'match', pdbfile[:-4], chain_id
		#print TM_segs
		fasta_data_all[ chain_id ] = TM_segs


# parse & store all fasta files of different chains for PDB's with a TM domain sequence pattern match
for f in os.listdir(fasta_dir):
	pdb, ch, = f[3:7], f[8]
	chain_id = f[:-6]
		
	if chain_id in fasta_data_all.keys():
		continue

	else:
		if pdb in PDBs_w_match:
			fastaPath = os.path.join( fasta_dir, '%s.fasta' % chain_id )
			TM_segs   = []

			for line in open(fastaPath):
				if line[:2] == '>>': 
					TM_id = line[2:].rstrip()
				else:
					target_seq 		= line.rstrip()
					TM_segs.append( ( TM_id, target_seq ) )
			fasta_data_all[ chain_id ] = TM_segs


fished_TMseqs = {}
if not os.path.exists(outdir):
	os.mkdir(outdir)
		
### look through matching TM domain, load full PDB,
#  further fish for additional TM domains within PDB that interactin in 3D with the query TM domain w/ matching sequence pattern
#  critera of fished interacting TM domain: Calpha atom w. 7-8 angstrom distance to each residue in sequence pattern
for TM_id, seqData in sorted(TM_seg_matches.items()):
	#print(TM_id, seqData)
	match_index = ( seqData[1], seqData[2] )

	TM_id2 		= TM_id.split()[0].split('_')
	res_range 	= [ int(x) for x in TM_id.split()[1].split('-') ]
	pdb_id, ch 	= TM_id2[0], TM_id2[1]
	#chain_id 	= 'tm_%s_%s_%d' % ( TM_id2[0], ch , ord(ch) ) 
		
	pdbPath 	= os.path.join( TMpdb_dir, 'tm_%s.pdb' % ( pdb_id ) )
	try:
		pdb 	= loadedPDBs[ pdb_id ]
	except:
		pdb 	= parsePDB(pdbPath)
		loadedPDBs[ pdb_id ] = pdb
				
			
	#print( "entering PDB", os.path.basename(pdbPath) )
	# find and select the 'target' TM matching the input sequence pattern from within the larger PDB
	match_TM_all 	= pdb.select('chain %s resnum %d to %d' % (ch, res_range[0],res_range[1] ) )
	matched_TMaas 	= pdb.select('chain %s resnum %d to %d' % (ch, res_range[0]+match_index[0], res_range[0] + match_index[1] -1 ) )
	seq_targ		= ''.join( [ UnNatAA[x] for x in matched_TMaas.select('ca').getResnames()] )

	#print ( np.arange(res_range[0]+match_index[0], res_range[0] + match_index[1] ) )
	#print ('chain %s resnum %d to %d' % (ch, res_range[0]+match_index[0], res_range[0] + match_index[1] -1 ))
	#print (seq_targ)


	#### This section should be changed to reflect sequence pattern in the query 
	#      could be coded in a more fancy way to be 100% general to resi index for the key residues within the input regex

	resi_range_in_match = np.arange(res_range[0]+match_index[0], res_range[0] + match_index[1] )
	pattern_aa_nums		= ( resi_range_in_match[0], resi_range_in_match[7], resi_range_in_match[14] )
	motif_TMaas			= pdb.select('ca chain %s resnum %d %s %d' % (ch, resi_range_in_match[0], resi_range_in_match[7], resi_range_in_match[14] ) )

	# look for... 
	ca_protein_notTM 	= pdb.select('ca and not Match', Match = match_TM_all )
	#print ( len(ca_protein_notTM), len(motif_TMaas), len(matched_TMaas) )

	close_CAs_ 			= [ ca_protein_notTM.select('within 7.7 of motif', motif=x.getCoords() ) for x in motif_TMaas]
	dMat 				= buildDistMatrix( motif_TMaas	, ca_protein_notTM)
	minDist = [np.min(x) for x in dMat]

	# Look for TM segments with contacts to all three Calpha atoms in the sequence pattern
	if all(i <= 7 for i in minDist):

		print "\n\n--> potential close interacting TM segment fished nearby atoms from matching query TM:", TM_id
		print(seq_targ, pattern_aa_nums, seq_targ[0], seq_targ[7], seq_targ[14] )
		
		# all chains involved in contacting the query/pattern residues.  Then, find any/all chain(s) in contact with all three
		all_chains	= [ set( x.getChids() ) for x in close_CAs_ ] 
		allclose_ch	= set.intersection(*all_chains)

		# for chains with residues contacting all pattern residues, find the full TM segment by residue index & select atoms
		for fished_ch in list(allclose_ch):

			#print 'matching TM sub-segment in chain', list(allclose_ch)[0], 'contacting all sequence pattern residues'
			
			chain_id_fished = 'tm_%s_%s_%d' % (pdb_id.upper(), fished_ch, ord(fished_ch)) 
			resi_from_ch	= [  x.select('chain %s' % fished_ch).getResnums() for x in close_CAs_ ] 
			# check if there's an atom remaining for list of fished interacting atoms to each motif aa, when selected by chain
			for check in resi_from_ch:
				if check is None:
					continue

			# find if a set of fished residues contacting the seq pattern residues come all from the same TM domain/span
			# versus contacts to the seq pattern residues coming from multiple spans 
			for TMs in fasta_data_all[chain_id_fished]:
				
				res_limits 	= [ int(x) for x in TMs[0].split()[1].split('-') ]
				resi_rangeX = np.arange( res_limits[0], res_limits[1] + 1)

				aa_in_TMseg = 0
				for a_set in resi_from_ch: 
					if len( list( set(a_set) & set(resi_rangeX) ) ) > 0:
						aa_in_TMseg += 1

				if aa_in_TMseg == 3:
					fished_TM 	= pdb.select('chain %s resnum %d to %d' % (fished_ch, resi_rangeX[0],resi_rangeX[-1] ) )
					segi 		= fished_TM.getSegnames()[0]

					moved_matchTM 	= matched_TMaas.select('ca').copy()
					if len(moved_matchTM) != len(target_matchAAs):
						continue
					transform 		= superpose( moved_matchTM, target_matchAAs)
					fished_TM_moved = applyTransformation( transform[1], fished_TM.copy() )


					# if the TM segment is <18 aa the reject/skip it
					fished_TMca		= fished_TM_moved.select('ca').copy()
					len_Fished		= len(fished_TMca)
					if len_Fished < 17: 
						print "TM skipped", "%s_%s%s" % (pdb_id.upper(), fished_ch, segi), "too short segment", len_Fished
						continue

					## find the "core" closest interacting segment of fished TM vs the target PDB; all possible 17aa possible windows
					dMat2 = buildDistMatrix( fished_TMca, target_matchAAs )

					index_range = np.arange(len_Fished)
					
					top_window, minVal = 0, 1000
					for n in np.arange(0, len_Fished - 16):
						#print n, n+16, fished_TMca[n].getResnum(), fished_TMca[n+16].getResnum(),  np.min( dMat2[n:n+16])
						tmp_minVal = np.min( dMat2[n:n+16])
						if tmp_minVal <minVal:
							minVal 		=  tmp_minVal
							top_window 	= (n, n+16)

					mini_resRange 		= fished_TMca[ top_window[0]: top_window[-1] ].getResnums()
					mini_fished_TM		= fished_TM_moved.select('resnum %d to %d' % (mini_resRange[0], mini_resRange[-1] ) )
					mini_fished_TMbb	= mini_fished_TM.select('bb')

					#calculate crossing angle, determine if parallel or antiparallel to input targeted helix
					# determine if parallel/antiparallel: Fit helix axis by SVD, compare to reference vector for input fragment 
					dataMean        = mini_fished_TMbb.getCoords().mean(axis=0)
					uu, dd, vv      = np.linalg.svd( mini_fished_TMbb.getCoords() - dataMean)
					helix_axis      = vv[0]
					cosine_angle = np.dot(helix_axis, helix_axis_tar) / ( np.linalg.norm(helix_axis) * np.linalg.norm(helix_axis_tar) )
					Crossing_angle = np.arccos(cosine_angle) * 180.0 / np.pi
					par_or_apar_flag = 'pa'
					if Crossing_angle > 90:
						par_or_apar_flag = 'ap'

					## write PDB of matched TM segment
					fishedTM_path 	= os.path.join(outdir, "%s_fished_%s_%s%s_%d-%d.pdb" % (par_or_apar_flag, pdb_id.upper(), fished_ch, segi, resi_rangeX[0], resi_rangeX[-1]) )
					miniTM_path 	= os.path.join(outdir, "%s_core-fished_%s_%s%s_%d-%d.pdb" % (par_or_apar_flag, pdb_id.upper(), fished_ch, segi, resi_rangeX[0], resi_rangeX[-1]) )
					pairTM_path 	= os.path.join(outdir, "%s_pairTM_%s_%s%s_%d-%d.pdb" % (par_or_apar_flag, pdb_id.upper(), fished_ch, segi, resi_rangeX[0], resi_rangeX[-1]) )

					writePDB(  fishedTM_path, fished_TM_moved )
					writePDB(  miniTM_path, mini_fished_TM )
					writePDB(  pairTM_path, moved_matchTM + fished_TM_moved )
					print "     ------> found bonafide interacting TM, crossing angle:", round(Crossing_angle, 1) , ", fished & printed:", par_or_apar_flag , fishedTM_path, '\n'
			
	else:
		continue


	

	









