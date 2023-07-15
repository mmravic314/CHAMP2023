import sys, os, numpy as np, random, re
from collections import Counter
from prody import *
# Marco Mravic DeGrado Lab UCSF Biophysics Nov 2015

## Some manual alignment or eyeball used to chose residue for sequence threading
## Yes, this step could be automated, but lots of lines of code to do something obvious

# input 1: target helix with MEM record from earlier rosetta minimization
# input 2: Two helix pair idealized model
# input 3: residue number of input 2 chain A to begin for threading/aligning the target input 1 helix



## Example command line
# python ./bin/


################## global definitions  #################

# Three to One res ID converter includes unnatural/modified amino acids
UnNatAA={}
UnNatAA["ALA"] = 'A'; UnNatAA["CYS"] = 'C'; UnNatAA["ASP"] = 'D'; UnNatAA["GLU"] = 'E'; UnNatAA["PHE"] = 'F';
UnNatAA["GLY"] = 'G'; UnNatAA["HIS"] = 'H'; UnNatAA["ILE"] = 'I'; UnNatAA["LYS"] = 'K';
UnNatAA["LEU"] = 'L'; UnNatAA["MET"] = 'M'; UnNatAA["ASN"] = 'N'; UnNatAA["PRO"] = 'P'; UnNatAA["GLN"] = 'Q';
UnNatAA["ARG"] = 'R'; UnNatAA["SER"] = 'S'; UnNatAA["THR"] = 'T'; UnNatAA["VAL"] = 'V'; UnNatAA["TRP"] = 'W'; UnNatAA["TYR"] = 'Y';
UnNatAA['ABA'] = 'A'; UnNatAA['CSO'] = 'C'; UnNatAA['CSD'] = 'C'; UnNatAA['CME'] = 'C';
UnNatAA['OCS'] = 'C'; UnNatAA["HSD"] = 'H'; UnNatAA['KCX'] = 'K'; UnNatAA['LLP'] = 'K';
UnNatAA['MLY'] = 'K'; UnNatAA['M3L'] = 'K'; UnNatAA['MSE'] = 'M'; UnNatAA['PCA'] = 'P'; UnNatAA['HYP'] = 'P';
UnNatAA['SEP'] = 'S'; UnNatAA['TPO'] = 'T'; UnNatAA['PTR'] = 'Y'


### Ez_potential
#Implementation from 
#Ez, a Depth-dependent Potential for Assessing the Energies of Insertion of Amino Acid Side-chains into Membranes: Derivation and Applications to Determining the Orientation of Transmembrane and Interfacial Helices
#
#http://www.sciencedirect.com/science/article/pii/S0022283606012095
#
Ez_potential = {}
Ez_potential['ALA'] = [0, -0.29, 10.22, 4.67]
Ez_potential['ASP'] = [0, 1.19, 14.25, 8.98]
Ez_potential['GLU'] = [0, 1.30, 14.66, 4.16]
Ez_potential['PHE'] = [0, -0.65, 19.67, 7.12]   # mm mod, since this value looks to high and I only get TRP 
Ez_potential['GLY'] = [0, -0.01, 13.86, 6.00]
Ez_potential['HIS'] = [0, 0.75, 12.26, 2.77]
Ez_potential['ILE'] = [0, -0.56, 14.34, 10.69]
Ez_potential['LYS'] = [0, 1.66, 11.11, 2.09]
Ez_potential['LEU'] = [0, -0.64, 17.34, 8.61]
Ez_potential['MET'] = [0, -0.28, 18.04, 7.13]
Ez_potential['ASN'] = [0, 0.89, 12.78, 6.28]
Ez_potential['PRO'] = [0, 0.83, 18.09, 3.53]
Ez_potential['GLN'] = [0, 1.21, 10.46, 2.59]
Ez_potential['ARG'] = [0, 1.55, 9.34, 4.68]
Ez_potential['SER'] = [0, 0.10, 13.86, 6.00]
Ez_potential['THR'] = [0, 0.01, 13.86, 6.00]
Ez_potential['VAL'] = [0, -0.47, 11.35, 4.97]
#equation 3
Ez_potential['TRP'] = [1, -0.65, 11.65, 7.20] # mm mod to -.65, since this value looks to high and I only get TRP 
Ez_potential['TYR'] = [1, -0.42, 13.04, 6.20]


def compute_Ez(r, z):
    if r not in Ez_potential:
        return 0
    equation,E0,zm,n = Ez_potential[r]

    if equation == 0:
        dEz =  E0/(1+(z/zm)**n)
    elif equation == 1:
        dEz = E0*np.exp(-(z-zm)**2/2/n**2)

    return dEz

# mock distribution of lipid-friendly amino acids in database of TM domains, found in Mravic et al 2019 Science. 
# suggested by Elazar et al (2022, eLife) that native like distributions give more functional designs in practice
TM_AAs		= [ 'A', 'G', 'I', 'L', 'F', 'S', 'T', 'V' ]
weightAA	= [0.16, 0.07, 0.14, 0.24, 0.11, 0.07, 0.07, 0.14]
weightAA	= [0.15, 0.05, 0.14, 0.35, 0.10, 0.05, 0.05, 0.12]

# original champ lipid-facing distrubtion
TM_AA_old	= [ 'A', 'I', 'F', 'V', 'L' ]
weightAAold	= [ 0.1, 0.1, 0.1, 0.1, 0.6 ]

################## global definitions  #################





################## main ops  #################



aligned_targetTM 	= parsePDB(sys.argv[1])
poly_ala_interface 	= parsePDB(sys.argv[2])
pattern_start	 	= int(sys.argv[3])

ca_interfaceA, ca_targetTM	= poly_ala_interface.select('ca chain A'), aligned_targetTM.select('ca')
#bb_interfaceA, bb_targetTM	= poly_ala_interface.select('bb chain A'), aligned_targetTM.select('bb')

length_target 		= len(ca_targetTM)
if length_target > 30 or length_target + pattern_start > len(ca_interfaceA):
	print ("input target and input residue offset error - too long for input interface dimer template")
	sys.exit()

if aligned_targetTM.getChids()[0] != 'A':
	print ("error detected... input target pdb not chain A")
	print ("to make it simple, the input target TM helix needs to be chain 'A'\n")


# align the ideal helix-helix interface model to the input aligned target TM domain
# transform the input expected complex chain A with target, with alignment excluding N & c-term less ideal helix geometries
# use user input for aligning the TM domain register between the target & template dimer/interface/complex model

firstResi, lastResi 	= ca_targetTM.getResnums()[0], ca_targetTM.getResnums()[-1]
bb_interfaceA_mobile 	= poly_ala_interface.select( 'ca chain A resnum %d to %d' % (pattern_start + 1 , pattern_start + length_target -2) )
bb_target 			 	= aligned_targetTM.select( 'ca resnum %d to %d' % ( firstResi+1, lastResi-1  ) )


transform 				= superpose(bb_interfaceA_mobile, bb_target  )
if calcRMSD(bb_interfaceA_mobile, bb_target  ) > 1:
	print ("poor alignment of input helix with ideal helix in dimer/interface template")
	sys.exit()



# create new merged molecular model w/ input target TM domain and chain B binder polyalanine 'CHAMP'
# trim the chain B of the docked / positioned binder (chain B) so it's within the expected membrane region. 

maxZ   			= ca_targetTM.copy()[0].getCoords()[2], ca_targetTM.copy()[-1].getCoords()[2]

good_resi_range = []
for r in poly_ala_interface.select('ca chain B').copy():
#	print r.getCoords(), r.getResnum(), 
	if r.getCoords()[2] < np.min(maxZ) - 1 or r.getCoords()[2] < np.max(maxZ)+1:
		good_resi_range.append( r.getResnum() )
#		print  '+'
#	else:
#		print 

chainB = poly_ala_interface.select('chain B resnum %d to %d ' % ( good_resi_range[0], good_resi_range[-1] )  ).copy()

# renumber chainB champ starting at residue 1
for n, r in enumerate( chainB.iterResidues(), 1 ):
	r.setResnum(n)

#combine input target & positioned CHAMP
newPDB = aligned_targetTM.select('protein').copy() + chainB
writePDB( 'design_input_interface.pdb',newPDB )


#### write constraints file for the target input PDB, chain A.  



#### write residue file 
#  using inter-helix distance to define lipid-facing vs interfacial / designable / repacking residues
newPDB 		= newPDB.copy()
champ_ca 	= newPDB.select('ca chain B').copy()
first_resnum_target = newPDB.select('ca chain A').copy()[0].getResnum()
nontermini_resnums	= newPDB.select('ca chain A').copy()[2:-3].getResnums()

target_BBsc 	= newPDB.select('chain A name N C CA O CB CG')
target_cen 		= calcCenter( newPDB.select('bb chain A') )[2]
design_flg_list = []
champ_Z_forRes	= []

index_champ = np.arange( len(champ_ca) )
first_six_r, last_six_r  = champ_ca[0:6].copy().getResnums(), champ_ca[-6:].copy().getResnums()
first_six, last_six  = index_champ[0:6], index_champ[-6:]
n4aa_EZ, c4aa_EZ = {}, {}


step = 0

cstStr = ''
txtA	= 'start\n'
seq 	= [ '-' for x in index_champ  ]

for r in newPDB.iterResidues():
	
	# write resfile based on distance for target chain 'A'; 
	# allow repacking for potential interacting residues (NATAA); fix roatmer for lipid-facing (NATRO) 
	if r.getChid() == 'A':
		if np.min( buildDistMatrix(r, chainB ) ) < 7:
			txtA += '%d %s NATAA\n' % ( r.getResnum(), r.getChid() )
		else:
			txtA += '%d %s NATRO\n' % ( r.getResnum(), r.getChid() )
	
		##  write constraint file for target
		#   skip if first AA
		ca 		= r.select('name CA')[0]
		coords 	= ca.getCoords()
		if ca.getResnum() in nontermini_resnums:
			cstStr += 'CoordinateConstraint CA %d%s CA %d%s %f %f %f HARMONIC 0.0 0.3\n' % ( ca.getResnum(),'A' ,first_resnum_target, 'A' , coords[0], coords[1], coords[2] )
		else:
			pass

		continue

	elif r.getChid() == 'B' and r.getResname() != 'MEM':
		design_flg 	= 0
		CA_CB 		= r.select( 'name CA CB ')

		notGly_flag	= len(CA_CB) - 1
		dMat = buildDistMatrix( CA_CB , target_BBsc )
		closest_SC_dist = np.min(dMat)

		# distance and sidechain orientation metrics to define champ residues as potentially interfacial for design
		if notGly_flag:
				CB_vs_CA_dmat 	= np.mean( dMat[1] - dMat[0] )		# is the CB of sidechain closer to interface than CA?
				design_flg 		= 0

				if closest_SC_dist < 6:
					design_flg = 1
				elif closest_SC_dist < 8 and CB_vs_CA_dmat < -0.2:
					design_flg = 1
				elif closest_SC_dist < 6 and CB_vs_CA_dmat < 0:
					design_flg = 1
				else:
					pass

		else:			# glycines
			if closest_SC_dist < 6:
				design_flg = 1
	else:
		continue


	# for champ lipid-facing residues, use Z depth to define residue identity 
	# find 'best' AA location for 1 Trp and 1 Tyr based on EZ potential at water-lipid interfacial region
	# all aa's in near bulk water, define as Lys
	# other lipid-facing residues in core of bilayer or interfacial region, define by semi-randomization
	# check to make sure no common motifs (small-x3-small, small-x6-small) are accidently placed

	design_flg_list.append( design_flg )
	seq[step] = '_'

	Z 		=  calcCenter(CA_CB)[2] - target_cen 
	Z_abs	= np.fabs( Z )
	champ_Z_forRes.append(Z)

	# assign non-interface non-designable residues based on depth dependent position in bilayer
	if not design_flg:
		# calc probability score of aromatic/cationic AA at that residue's Z depth with EZ potential
		probs 	= [ round( compute_Ez(x, Z_abs), 4 ) for x in ['TRP', 'TYR', 'ARG', 'LYS', 'PHE'] ]
		probs.append(Z)

		# in polar bulk water region of bilayer
		if Z_abs  > 16  :					# if outside hydrophobic core or headgroup region, default to Lys 
			if Z > 0: 
				seq[step] = 'K'
			else:
				seq[step] = 'E'

		# if well inside hydrophobic core, default to semi-randomized AA
		elif Z_abs < 12.5 :						
				#	print(r, step, '-', round(Z, 2) )
			pass

		# if at headgroup region, calculate EZ potential for F, Y, W
		elif step in first_six:		
				n4aa_EZ[step] =  probs
				#	print( r, step, '-',round(Z, 1), probs)

		elif step in last_six:			
				n4aa_EZ[step] =  probs
				#	print( r, step, '-',round(Z, 1), probs)
		else:
				#	print(r, step, '-', round(Z, 2) )
				pass
	else:
		
		#  if designable & in water-exposed region include polar amino acids in allowed residues
		if Z_abs > 16:
			seq[step] = 'Z'
		else:
			seq[step] = 'X'

	step += 1



# put exactly 1 Trp & 1 Tyr per TM domain, located at the 'optimal' depth in the bilayer interface region
# optimized using the EZ potential depth-dependent score for each potential residues in the CHAMP
W_top = sorted( n4aa_EZ.items(), key=lambda x: x[1][0] )[0][0]
seq[W_top] = 'W'
del n4aa_EZ[W_top]

Y_top = sorted( n4aa_EZ.items(), key=lambda x: x[1][1] )[1][0]
seq[Y_top] = 'Y'
del n4aa_EZ[Y_top]

#print seq


numRandom_AA = len( [ x for x in seq if x == '_'] )


## weighed randomization of remaining lipid-facing residues based on naturally occuring TM domain residue frequency
# perform while loop trials creation of sequences. quit only upon successful after sequence motif checks
# i.e. quit if well known sequence repeat motifs for homo-oligomerix self-assembly are found 
# no G/S/AxxxG/S/A motifs or 'G/A/SxxxxxxG/S/A', G/S/TxxxT/G/S,   
# LxxI(V/T/C/I)xxL is also rejected
# also, no more than 4 Leu in the lipid facing, no more than 1 of each small residue
known_motif_tegex = []
randomTMaa_success_flg = 0
#while randomTMaa_success_flg == 0:
while randomTMaa_success_flg == 0:
		step = 0
		random_AAs = random.choices(TM_AAs	,  weights=tuple(weightAA), k=numRandom_AA)

		seq2 = []
		for i in seq :
			if i=='_':
				seq2.append( random_AAs[step] )
				step +=1
			else:
				seq2.append( i )

		seq2_str = ''.join( seq2 ) 

		## sequence counter filters
		## no more than 1 of S,T,G; no more than 4 L
		counts = Counter( seq2 )
		if counts['S'] > 2 or counts['G'] > 1 or counts['T'] > 2 or counts['L'] > 4 or counts['I'] > 3 or counts['F'] < 2:
			#print (counts.most_common(), '--------------BAD SEQ, count filter broken!!! \n')
			continue

		## sequence pattern filters
		match = re.search(r"[GAS]\w\w\w[GAS]", seq2_str )
		if match:
			#print ( '--------------BAD SEQ, seq pattern 1 filter broken!!! \n')
			continue

		match = re.search(r"[GST]\w\w\w[GST]", seq2_str )
		if match:
			#print ( '--------------BAD SEQ, seq pattern 2 filter broken!!! \n')
			continue

		match = re.search(r"[GSA]\w\w\w\w\w\w[GSA]", seq2_str )
		if match:
			#print ( '--------------BAD SEQ, seq pattern 3 filter broken!!! \n')
			continue
		
		randomTMaa_success_flg = 1


print ('\t', seq2_str)
####  write resfile & constraint file

## you might want to reconsider the amino acid group for the designable set 
# e.g. to include polars if there's polars in the target; allow Gly if the input motif isn't already all small resi's

txtB = ''
for aa, r in zip(seq2_str, champ_ca.iterResidues() ):


	if aa == 'Z':
		print('%d %s ALLAAxc\n' % ( r.getResnum(), r.getChid()  ))
		txtB += '%d %s ALLAAxc\n' % ( r.getResnum(), r.getChid()  )
	
	elif aa == 'X':
		print('%d %s PIKAA ASTVIMLF\n' % ( r.getResnum(), r.getChid()  ))
		txtB += '%d %s PIKAA ASTVIMLF\n' % ( r.getResnum(), r.getChid() )

	else:
		print('%d %s PIKAA %s\n' % ( r.getResnum(), r.getChid(), aa ))
		txtB += '%d %s PIKAA %s\n' % ( r.getResnum(), r.getChid(), aa )

resfile_F = open('resfile.txt', 'w')
resfile_F.write( txtA + txtB )
resfile_F.close()


cstfile_F = open('constraint.txt', 'w')
cstfile_F.write( cstStr )
cstfile_F.close()
print(cstStr)







