# Marco Mravic, Scripps 2021

import sys, os, subprocess as sp, numpy as np
from prody import *

## Input 1: target structure	path of this will be used to find span file w/ same name
## Input 2: path to rosetta main
## input 3: path to Rosetta scripts XML with protocol 
## input 4: path to resfile
## input 5: path to constraints

## Example command line
## single CPU runs
# ./bin/5_runRosettaMP_Design.py design_input_interface.pdb ~/rosetta/ ./bin/5_design_TM_protocol1.xml ./resfile.txt constraint.txt run1_outputs/
# 
# multiple CPUs on same machine
# seq 1 3 | xargs -n 1 -P 3 ./bin/5_runRosettaMP_Design.py design_input_interface.pdb ~/rosetta/ ./bin/5_design_TM_protocol1.xml ./resfile.txt constraint.txt run1_outputs/
##
#
#python3 ./bin/5_runRosettaMP_Design.py design_input_interface.pdb ~/rosetta/ ./bin/5_design_TM_protocol2.xml ./resfile.txt constraint.txt run2_outputs/



inPDB		= sys.argv[1]
rosiBase 	= sys.argv[2]
protocolPth = sys.argv[3]
resfile_path= sys.argv[4]
cstfile_path= sys.argv[5]
outputDir 	= sys.argv[6]
if len(sys.argv) == 8:
	offset  = sys.argv[7]
else:
	offset  = 0

rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.macosclangrelease' )
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiScoreF 	= os.path.join( rosiBase, 'database/scoring/score_functions/aa_composition/')
rosiHoles 	= os.path.join( rosiBase, 'source/external/DAlpahBall/DAlphaBall.gcc')

'''spanF 		= inPDB[:-4] + '.span'
if not os.path.exists( spanF ):
	print('ERROR no span file at: %s \n\n' % spanF )
	sys.exit()
'''

#check the last file's suffix numbering.  Next files should be increasing index
if len( os.listdir(outputDir) ) == 0:
	pdbs= []
else:
	pdbs = [ x for x in os.listdir(outputDir) if 'pdb' in x ]

if len(pdbs) < 1:
	suffix = "_%s" % str(1).zfill(3)
else:
	suffices = [0]
	for x in pdbs:
		print (x, x.split('_') )
		suf = int(x.split('_')[-2])
		suffices.append(suf)
	suffix = "_%s" % str( np.max(suffices) + 1 + int(offset) ).zfill(3)



cmd = [  rosiScrps, 
'-parser:protocol', protocolPth, 			# Path to Rosetta script (see above)
'-in:file:s', inPDB,							# Input PDB structure
'-nstruct', '10', 							# Generate 1 model
#'-mp:setup:spanfiles', spanF,				# Input spanfile
'-mp::setup::spans_from_structure', '1',
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:overwrite',
'-packing:resfile', resfile_path , 
'-packing:pack_missing_sidechains', '1',
'-ignore_zero_occupancy', 'false',
'-mp:lipids:composition', 'POPC',
'-out:path:all', outputDir,
'-parser:script_vars', 'cst_file=%s' % cstfile_path, #'AAcomp_file=./TMab_AA_composition.comp'
'-aa_composition_setup_file', os.path.join( rosiScoreF, 'TMab_AA_comp_less2-MWY.comp' ),
'-holes:dalphaball', rosiHoles, 
'-out:suffix', suffix
 ]

sp.call( cmd )

