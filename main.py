import networkx as nx
import math
import optparse
import os
from optparse import OptionParser
import itertools
from copy import deepcopy
import sys


'''
Takes as input a file name denoting a molecule .str file, and this then generates a graph 
of that molecule with the nodes as atom names, node attributes as atom types and partial charges,
and edges denoting bonds.
'''
def generateGraph(molecule):
	molecule_file = open(molecule, "rb")

	graph = nx.Graph()

	for line in molecule_file.readlines():

		if (line[0:4] == "ATOM"):
			attributes = line.split()
			atom_name = attributes[1]
			atom_type = attributes[2]
			atom_charge = attributes[3]
			graph.add_node(atom_name, atom_type=atom_type, charge=float(atom_charge))

		elif(line[0:5] == "BOND "):
			atom1 = line.split()[1]
			atom2 = line.split()[2]
			graph.add_edge(atom1, atom2)

	return graph

'''
For the given str file (molecule), this function generates a 4-tuple of dictionaries containing
all the parameters (bonds, anles, dihedrals, improper dihedrals) for that molecule. Each parameter type
has associated with it a dictionary, with keys being atom types, and values being actual parmeters.
For example, for angle parameters, there will be a dictionary with a key being something like
(CG2DC1, CG2DC1, CG2O2) and a value being something like (60.00, 120.00), in the exact same order
as that given in CHARMM str input files. 
'''
def getParameters(molecule):
	molecule_file = open(molecule, "rb")

	lines = molecule_file.readlines()
	i = 0
	line = lines[0]

	exist_bonds = True
	while (line[0:5] != "BONDS"):
		if (line == "END"):
			print "There are no bond length parameters"
			exist_bonds = False
			break
		i += 1
		line = lines[i]


	i += 1
	bond_line = deepcopy(i)

	exist_angles = True
	while (line[0:6] != "ANGLES"):
		if (line == "END"):
			print "There are no angle parameters"
			exist_angles = False
			break
		i += 1
		line = lines[i]

	i += 1
	angle_line = deepcopy(i)

	exist_dihedrals = True
	while (line[0:9] != "DIHEDRALS"):
		if (line == "END"):
			print "There are no dihedral parameters"
			exist_dihedrals = False
			break
		i += 1
		line = lines[i]

	i += 1
	dihedral_line = deepcopy(i)

	exist_impropers = True
	while (line[0:9] != "IMPROPERS"):
		if (line == "END"):
			print "There are no improper parameters"
			exist_angles = False
			break
		i += 1
		line = lines[i]

	i += 1
	improper_line = deepcopy(i)

	while (line[0:3] != "END" and line[0:9] != "NONBONDED"):
		i += 1
		line = lines[i]

	end_line = deepcopy(i)

	bonds = {}
	angles = {}
	dihedrals = {}
	impropers = {}

	if exist_bonds:
		for j in range(bond_line, angle_line-2):
			if (lines[j][0] == "!"): continue
			bond = lines[j].split()
			if len(bond) == 0: continue
			atoms = (bond[0], bond[1])
			parms = (bond[2], bond[3])
			bonds[atoms] = [float(k) for k in parms]

	if exist_angles:
		for j in range(angle_line, dihedral_line-2):
			if (lines[j][0] == "!"): continue
			angle = lines[j].split()
			if len(angle) == 0: continue
			atoms = (angle[0], angle[1], angle[2])
			parms = (angle[3], angle[4])
			angles[atoms] = [float(k) for k in parms]

	if exist_dihedrals:
		for j in range(dihedral_line, improper_line-2):
			if (lines[j][0] == "!"): continue
			dihedral = lines[j].split()
			if len(dihedral) == 0: continue
			atoms = (dihedral[0], dihedral[1], dihedral[2], dihedral[3])
			parms = (dihedral[4], dihedral[5], dihedral[6])
			dihedrals[atoms] = [float(k) for k in parms]

	if exist_impropers:
		for j in range(improper_line, end_line-2):
			if (lines[j][0] == "!"): continue
			improper = lines[j].split()
			if len(improper) == 0: continue
			atoms = (improper[0], improper[1], improper[2], improper[3])
			parms = (improper[4], improper[5], improper[6])
			impropers[atoms] = [float(k) for k in parms]
	return (bonds, angles, dihedrals, impropers)



#Calculates the total charge of the given molecule 
def totalCharge(molecule):
	total_charge = 0.0

	for atom in molecule.nodes():
		total_charge += molecule.node[atom]["charge"]

	return total_charge

#For a given atom in a given graph, returns the hydrogens bonded to that atom
def getHydrogenNeighbors(molecule, atom):
	atom_neighbors = set(molecule.neighbors(atom))
	atom_hydrogens = set()
	for neighbor in atom_neighbors:
		if(neighbor[0] == "H"):
			atom_hydrogens.add(neighbor)
	return atom_hydrogens

#For a given set of atoms in a molecule, this function subtracts the charges of the hydrogens
#bonded to those atoms from their partial charges. This is used  in the last stage of partial charge
#stitching.
def subtractHydrogens(molecule, atoms):
	for atom in atoms:
		atom_hydrogens = getHydrogenNeighbors(molecule, atom)
		for hydrogen in atom_hydrogens:
			molecule.node[atom]["charge"] -= molecule.node[hydrogen]["charge"]

	return molecule

'''
One of the primary functions of the code. updateCharges takes as input the "stitched" molecule, 
a set of fragments that will comrpise it, and non-hydrogen atoms that are common between fragments.
'''

def updateCharges(stitched, molecules, common_atoms):
	for atom in common_atoms:
		stitched.node[atom]["charge"] = 0.0
	
	for i in range(1, len(molecules)):
		molecule = molecules[i]
		common_hydrogens = set()

		for atom in (common_atoms.intersection(molecule.nodes())):
			atom_neighbors = set(molecule.neighbors(atom))
			for neighbor in atom_neighbors:
				if(neighbor[0] == "H"):
					common_hydrogens.add(neighbor)

		for atom in molecule.nodes():
			if(not (atom in common_atoms) and not (atom in common_hydrogens)):
				stitched.node[atom]["charge"] = molecule.node[atom]["charge"]
			elif(atom in common_atoms):
				atom_neighbors = set(molecule.neighbors(atom))
				atom_hydrogens = getHydrogenNeighbors(molecule, atom)
				stitched.node[atom]["charge"] += molecule.node[atom]["charge"]
				for hydrogen in atom_hydrogens:
					stitched.node[atom]["charge"] += molecule.node[hydrogen]["charge"]

	stitched = subtractHydrogens(stitched, common_atoms)
	return(stitched)

'''
The following function, updateParameters, takes as input the individual fragments, the original molecule,
a list of 4-tuples of dictionaries containing parameters. 

Summary of "parameters" data structure (a list of 4-tuples of dictionaries mapping tuples to tuples):

Parameters = [original_molecule_params, frag_1_params, frag_2_params,...,frag_n_params]
where:
e.g.: frag_1_params = (bonds, angles, dihedrals, impropers)
where:
e.g.: angles = {(atomtype1, atomtype2, atomtype3) --> (angle_constant, equilibribium angle)}

This function goes through the parameters listed for the original molecule, maps the atom types originally assigned
to atom types given to the equivalent atom in the individual fragments, then assigns the equivalent parameter terms
from the fragments to the original molecule. This gives the original molecule a new set of parameters consisting mostly
of parameters from the fragments (except near the atoms at which bonds were formed between fragments), which are returned
as stitched_parameters, a 4-tuple of dictionaries.

'''

def updateParameters(molecules, parameters, stitched):
	for i in range(1, len(molecules)):
		fragment = molecules[i]
		for frag_atom in fragment.nodes():
			matching = filter(lambda x: frag_atom[0:3] in x, stitched.nodes())
			for stitched_atom in matching:
				stitched.node[stitched_atom]["atom_type"] = fragment.node[frag_atom]["atom_type"]
	
	stitched_parameters = parameters[0]
	for j in range(0, len(stitched_parameters)):
		for parameter_key in stitched_parameters[j].keys():
			for i in range(1, len(molecules)):
				fragment = molecules[i]
				fragment_parameters = parameters[i]

				if parameter_key in fragment_parameters[j].keys():
					stitched_parameters[j][parameter_key] = fragment_parameters[j][parameter_key]

				elif reversed(parameter_key) in fragment_parameters[j].keys():
					stitched_parameters[j].pop(parameter_key, None)
					stitched_parameters[j][reversed(parameter_key)] = fragment_parameters[j][reversed(parameter_key)]

	return (stitched, stitched_parameters)

'''
Writes all new charges, topology, and parameters to a new .str file with user-provided name, 
and provides option to output all params including those in CGenFF, or to skip those. 
'''

def writeToFile(original_file, stitched, parameters, cgenff, output_all_params, new_title):
	if(output_all_params.lower() == "true"):
		print "Outputting all bond/angle/dihedral params, even those already in CGenFF"
	else:
		print "Outputting only those bond/angle/dihedral params not already in CGenFF"


	molecule_file = open(original_file, "rb")
	atom_names = []
	bond_lines = []
	improper_lines = []
	for line in molecule_file.readlines():
		if (line[0:4] == "ATOM"):
			atom_names.append(line.split()[1])
		if (line[0:5] == "BOND "):
			bond_lines.append(line)
		if (line[0:5] == "IMPR "):
			improper_lines.append(line)

			
	new_file = open("%s" %new_title, 'wb')

	for atom in atom_names:
		atom_type = stitched.node[atom]["atom_type"]
		charge = stitched.node[atom]["charge"]
		new_file.write("ATOM %s" %atom)

		nspaces = 7 - len(atom)
		for i in range(0,nspaces):
			new_file.write(" ")
		new_file.write(atom_type)

		if charge > 0: nspaces = 8 - len(atom_type)
		if charge < 0: nspaces = 7 - len(atom_type)
		for i in range(0,nspaces):
			new_file.write(" ")
		new_file.write("%.3f\n" %charge)

	new_file.write("\n")

	for bond in bond_lines:
		new_file.write(bond)
	for improper in improper_lines:
		new_file.write(improper)
	new_file.write("\n")

	new_file.write("END \n \n ! Begin Parameters: \n \n ")

	for j in range(0, len(parameters)):
		if j == 0:
			new_file.write("BONDS\n")
			for parm_key in parameters[j].keys():
				if output_all_params.lower() == "true":
					if parm_key in cgenff[j].keys(): continue
				atom = ""
				for k in range(0,len(parm_key)):
					atom = parm_key[k]
					new_file.write(atom)
					for l in range(0,7-len(atom)):
						new_file.write(" ")
				new_file.write(" ")
				new_file.write("%.2f" %parameters[j][parm_key][0])
				for l in range(0,5):
					new_file.write(" ")
				new_file.write("%.4f\n" %parameters[j][parm_key][1])
			new_file.write("\n")

		if j == 1:
			new_file.write("ANGLES\n")
			for parm_key in parameters[j].keys():
				if output_all_params.lower() == "true":
					if parm_key in cgenff[j].keys(): continue
				atom = ""
				for k in range(0,len(parm_key)):
					atom = parm_key[k]
					new_file.write(atom)
					for l in range(0,7-len(atom)):
						new_file.write(" ")
				new_file.write("  ")
				new_file.write("%.2f" %parameters[j][parm_key][0])
				for l in range(0,4):
					new_file.write(" ")
				new_file.write("%.2f\n" %parameters[j][parm_key][1])
			new_file.write("\n")

		if j == 2:
			new_file.write("DIHEDRALS\n")
			for parm_key in parameters[j].keys():
				if output_all_params.lower() == "true":
					if parm_key in cgenff[j].keys(): continue
				atom = ""
				for k in range(0,len(parm_key)-1):
					atom = parm_key[k]
					new_file.write(atom)
					for l in range(0,7-len(atom)):
						new_file.write(" ")
				new_file.write(atom)
				for l in range(0,10-len(atom)):
					new_file.write(" ")
				new_file.write("%.4f" %parameters[j][parm_key][0])
				for l in range(0,2):
					new_file.write(" ")
				new_file.write("%d" %parameters[j][parm_key][1])
				phase = parameters[j][parm_key][2]

				for l in range(0,8-len(str(int(parameters[j][parm_key][2])))):
					new_file.write(" ")
				new_file.write("%.2f\n" %parameters[j][parm_key][2])
			new_file.write("\n")

		if j == 3:
			new_file.write("IMPROPERS\n")
			for parm_key in parameters[j].keys():
				if output_all_params.lower() == "true":
					if parm_key in cgenff[j].keys(): continue
				atom = ""
				for k in range(0,len(parm_key)-1):
					atom = parm_key[k]
					new_file.write(atom)
					for l in range(0,7-len(atom)):
						new_file.write(" ")
				new_file.write(atom)
				for l in range(0,10-len(atom)):
					new_file.write(" ")
				new_file.write("%.4f" %parameters[j][parm_key][0])
				for l in range(0,2):
					new_file.write(" ")
				new_file.write("%d" %parameters[j][parm_key][1])
				phase = parameters[j][parm_key][2]

				for l in range(0,8-len(str(int(parameters[j][parm_key][2])))):
					new_file.write(" ")
				new_file.write("%.2f\n" %parameters[j][parm_key][2])
			new_file.write("\n")
	new_file.write("END\n")
	new_file.write("RETURN")				

				

PROG = os.path.basename(os.path.splitext(__file__)[0])

parser = OptionParser()
parser.add_option("-o", "--original_ligand", type="string", help="enter the str file for the entire ligand")
parser.add_option("-f", "--fragment", action="append", type="string", help="enter the str file for fragment. each fragment has its own -f flag")
parser.add_option("-s", "--stitched_name", type="string", help="enter the name of the file to which your new stitched molecule will be saved")
parser.add_option("-p", "--use_frag_params", type="string", help="enter TRUE if you want to use bonded parameters from fragments instead of original str",default="false")
parser.add_option("-a", "--output_all_params", type="string", help="enter TRUE if you want to output all parameters, including those already in CGenFF", default="false")

opts, args = parser.parse_args()

molecule_files = []

if(opts.original_ligand):
	molecule_files.append(opts.original_ligand)
else:
	print("No original molecule provided")
	sys.exit()

if(opts.fragment):
	for fragment in opts.fragment:
		molecule_files.append(fragment)
else:
	print("No provided fragments")
	sys.exit()

molecules = []

#print(molecule_files)

for molecule in molecule_files:
	molecules.append(generateGraph(molecule))

for i in range(0, len(molecules)):
	print("Total charge in %s is %.3f" %(molecule_files[i], totalCharge(molecules[i])))

#Finds all atoms common to those in the different fragments

common_atoms = set()

for comb in itertools.combinations(range(1,len(molecules)), 2):
	molecule_i = molecules[comb[0]]
	molecule_j = molecules[comb[1]]
	intersection = set.intersection(set(molecule_i.nodes()), set(molecule_j.nodes()))
	common_atoms = set.union(common_atoms, intersection)

hydrogens = set()
for atom in common_atoms:
	if(atom[0] == 'H'):
		hydrogens.add(atom)

#common_atoms now contains only non-hydrogen atoms 

common_atoms -= hydrogens

#stitched denotes the new version of the parent molecule that is being parameterized 

stitched = deepcopy(molecules[0])

stitched = updateCharges(stitched, molecules, common_atoms)


parameters = list()

#reads in bonded parameters from each of the user-listed files

for i in range(0, len(molecule_files)):
	parameters.append(getParameters(molecule_files[i]))

stitched_parameters = parameters[0]

if(opts.use_frag_params.lower() == "true"):
	print "Using bond/angle/dihedral params from the individual fragments when possible"
	(stitched, stitched_parameters) = updateParameters(molecules, parameters, stitched)
else:
	print "Using bond/angle/dihedral params from original .str file"

cgenff_parameters = getParameters("../toppar/par_all36_cgenff.prm")

writeToFile(opts.original_ligand, stitched, stitched_parameters, cgenff_parameters, opts.output_all_params, opts.stitched_name)

print("Total charge in final stitched molecule is %.3f" %(totalCharge(stitched)))
print("DONE.")


