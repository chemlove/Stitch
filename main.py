import networkx as nx
import math
import optparse
import os
from optparse import OptionParser
import itertools
from copy import deepcopy


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
		print(stitched.node["CAM"]["charge"])

	stitched = subtractHydrogens(stitched, common_atoms)
	return(stitched)

def writeToFile(original_file, stitched, new_title):
	molecule_file = open(original_file, "rb")
	atom_names = []
	for line in molecule_file.readlines():
		if (line[0:4] == "ATOM"):
			atom_names.append(line.split()[1])
			
	new_file = open("%s" %new_title, 'wb')

	for atom in atom_names:
		atom_type = stitched.node[atom]["atom_type"]
		charge = stitched.node[atom]["charge"]
		new_file.write("ATOM %s" %atom)
		for i in range(0,4-(len(atom)-3)):
			new_file.write(" ")
		new_file.write(atom_type)

		for i in range(0,4-(len(atom_type)-4)):
			new_file.write(" ")
		new_file.write("%.3f\n" %charge)
				

PROG = os.path.basename(os.path.splitext(__file__)[0])

parser = OptionParser()
parser.add_option("-o", "--original_ligand", type="string", help="enter the str file for the entire ligand")
parser.add_option("-f", "--fragment", action="append", type="string", help="enter the str file for fragment. each fragment has its own -f flag")
parser.add_option("-s", "--stitched_name", type="string", help="enter the name of the file to which your new stitched molecule will be saved")

opts, args = parser.parse_args()
print "arguments:", args
print "options:", opts

molecule_files = []

molecule_files.append(opts.original_ligand)

if(opts.fragment):
	for fragment in opts.fragment:
		molecule_files.append(fragment)

molecules = []

#print(molecule_files)

for molecule in molecule_files:
	molecules.append(generateGraph(molecule))

for molecule in molecules:
	print(totalCharge(molecule))

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

common_atoms -= hydrogens

stitched = deepcopy(molecules[0])

stitched = updateCharges(stitched, molecules, common_atoms)

writeToFile(opts.original_ligand, stitched, opts.stitched_name)

print("\n \n")
print(stitched.node["CAM"])
print(stitched.node["CAH"])
print(stitched.node["CBD"])
print(totalCharge(stitched))



