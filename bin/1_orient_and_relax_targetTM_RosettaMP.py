# Marco Mravic, Scripps 2021

## Input 1: pdb of ideal helix threaded with target sequence, e.g. mEpoR TM region, with first guess insertion by OPM PPM 
## Input 2: path to rosetta main
## input 3: path to Rosetta scripts XML with protocol 

## Example command line
# python ./bin/orient_and_relax_target_RosettaMP.py opm_eVgL_pent_clean.pdb ~/rosetta/ ../bin/tmBundle_Relax.xml 


import sys, os, subprocess as sp
from prody import *

inPDB		= sys.argv[1]
rosiBase 	= sys.argv[2]
protocolPth = sys.argv[3]
span_file   = sys.argv[4]
rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.macosclangrelease' )
rosiDB 		= os.path.join( rosiBase, 'database/')
#rosiOctPrl 	= os.path.join( rosiBase, 'source/src/apps/public/membrane_abinitio/octopus2span.pl')
#rosiSpanGen = os.path.join( rosiBase, 'source/bin/spanfile_from_pdb.macosclangrelease')
#rosiHoles 	= os.path.join( rosiBase, 'source/external/DAlpahBall/DAlphaBall.gcc')

pdb 	= parsePDB(inPDB)
#### Renumber input target to start at 1
step 	= 1
for i in pdb.iterResidues():
	i.setResnum(step)
	step +=1
writePDB(inPDB, pdb)



#spanF 		= inPDB[:-4] + '.span'
## Make span file from input PDB

#cmdSpan = [ rosiSpanGen, 
#'-database', rosiDB, 
#'-in:file:s', inPDB, '-ignore_zero_occupancy', 'false'
#]

#if not os.path.exists( spanF ):
#	sp.call( cmdSpan )
#sys.exit()



### Relax
cmd = [  rosiScrps, 
'-parser:protocol', protocolPth, 			# Path to Rosetta script (see above)
'-in:file:s', inPDB,							# Input PDB structure
'-nstruct', '1', 							# Generate 1  output model
'-mp:setup:spanfiles', span_file,	
#'-mp::setup::spans_from_structure', '1',	# Generate RosettaMP .span file of 'embedded' residues from input model coords
'-relax:jump_move', 'true', 				# Allow jumps (membrane residue) to move during relax
'-out:overwrite',
'-packing:pack_missing_sidechains', '1',	
'-ignore_zero_occupancy', 'false',
'-mp:lipids:composition', 'DPPC',			# default mammalian membrane-like thickness
'-keep_input_protonation_state', '1',
 ]


sp.call( cmd )

sys.exit()







##

'''  OPTIONS FROM 2015 PLOS paper
-parser:protocol membrane_relax.xml # Path to Rosetta script (see above)
-in:file:s 3PXO_tr_native.pdb # Input PDB structure
-nstruct 1000 # Generate 1000 models
-mp:setup:spanfiles 3PX0.span # Input spanfile
-mp:scoring:hbond # Turn on membrane hydrogen bonding
-relax:jump_move true # Allow jumps to move during relax
-packing:pack_missing_sidechains 0 # Wait to pack sidechains 
'''