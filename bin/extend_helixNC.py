import numpy as np, os, sys
from prody import *

# input 1: PDB fine of helical protein chains to modify
# input 2: input PDB ideal helix to use as source of atoms for extending input 1
# input 3: chains to motify, no spaces, e.g.:  A;  AB; ADE; XYZ;
# input 4: amino acids at N term to add; dumb for now, applied to all chains (maximum of 9aa)
# input 5: amino acids at C term to add 

# Note:  deletes the N & C term residue and replaces with helical alinine

#usage python ~/CHAMP/2023_design-champ-mEpoR/bin/extend_helixNC.py  model_1.pdb ~/bin/20_allAlaIdealHelix.pdb AB 5 5

pdb 	= parsePDB(sys.argv[1])
outfile = 'ext_' + os.path.basename(sys.argv[1])
print(outfile)
helix 	= parsePDB(sys.argv[2])
chains 	= sys.argv[3]

naa 	=  int(sys.argv[4] ) + 5 
caa 	=  int(sys.argv[5] ) + 5
#print ('resnums 1 to %d' % 4+int(sys.argv[4]) , 'resnums %d to 20' % 20-4-int(sys.argv[5]) )

# grab fragments of ideal helix: 4 aa + number of extension residues
#nExt 	= helix.select('resnum 1 to %d' % naa   )
#cExt 	= helix.select('resnum %d to 20' % caa  )

newChains = []

for ch in chains:
	
	chain 		= pdb.select('chain %s' % ch)
	resnums		= chain.select('ca').getResnums()

	# transform the n-terminal extension aligned to the n-term of helix, skipping target pdb's first aa
	bb_chain_n 	= chain.select( 'bb resnum %s to %s' % (resnums[1], resnums[4]) )
	mobile_Nbb 	= helix.select('bb resnum %d to %d' % (naa-3, naa) ).copy()
	#print( 'bb resnum %d to %d' % (naa-3, naa),  'bb resnum %s to %s' % (resnums[1], resnums[4]) )
	bb, trans	= superpose(mobile_Nbb,bb_chain_n)
	rmsd 		= calcRMSD(bb, bb_chain_n)
	if rmsd > 1.5: 
		print('\n\n\tunacceptable helical conformation of the target/input "helix" N termini... is it really a helical conformation? \n\t\t%.3f Angstrom RMSD to ideal helix\n' % rmsd)
		sys.exit()

	#print ( round( rmsd, 3) )

	new_nExt 	= applyTransformation( trans, helix.select('resnum 1 to %d' % (naa-4) ).copy() )


	# transform the c-terminal extension aligned to the c-term of helix, skipping target pdb's last aa
	bb_chain_c 	= chain.select( 'bb resnum %s to %s' % (resnums[-5], resnums[-2]) )
	mobile_cbb 	= helix.select('bb resnum 1 to 4' ).copy()
	#print( 'bb resnum %s to %s' % (resnums[-5], resnums[-2]) )

	bb, trans	= superpose(mobile_cbb, bb_chain_c)
	rmsd 		= calcRMSD(bb, bb_chain_c)
	if rmsd > 1.5: 
		print('\n\n\tunacceptable helical conformation of the target/input "helix" C termini... is it really a helical conformation? \n\t\t%.3f Angstrom RMSD to ideal helix\n' % rmsd)
		sys.exit()
	#print ( round( rmsd, 3) )
	new_cExt 	= applyTransformation( trans, helix.select('resnum 5 to %d' % caa).copy() ) 
	#print( 'resnum 5 to %d' % caa )


	target_trimmed 	= 	chain.select( 'resnum %s to %s' % (resnums[1], resnums[-2]) )
	#print( 'bb resnum %s to %s' % (resnums[1], resnums[-2]) )


	newChainPDB 	= new_nExt.copy() + target_trimmed.copy() 
	# renumber the new PDBs
	r = 1
	for resi in newChainPDB.iterResidues():
		#print(resi.getResnum())
		resi.setResnum(r)
		r+=1
	#	print(resi.getResnum())

	r = resnums[-1] + int(sys.argv[4]) + 1
	for resi in new_cExt.iterResidues():
		#print(resi.getResnum())
		resi.setResnum(r)
		r+=1
	#	print(resi.getResnum())

		# renumber the new PDB
	newChainPDB 	= newChainPDB + new_cExt.copy()

	
	# correct chain and segment IDs
	newChainPDB.setChids( [ch for a in newChainPDB.getChids()] )
	newChainPDB.setSegnames( [ch for a in newChainPDB.getSegnames()] )


	writePDB('tmpChain_%s.pdb' % ch, newChainPDB)
	if len(newChains) == 0:
		newChains = newChainPDB
	else:
		newChains = newChains + newChainPDB


writePDB(outfile, newChains)




