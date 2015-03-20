import networkx as nx
import math
import optparse
import os
from optparse import OptionParser

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

def totalCharge(molecule):
	total_charge = 0.0

	for atom in molecule.nodes():
		total_charge += molecule.node[atom]["charge"]

	return total_charge




PROG = os.path.basename(os.path.splitext(__file__)[0])

parser = OptionParser()
parser.add_option("-o", "--original_ligand", action="append", type="string", help="enter the str file for the entire ligand")
parser.add_option("-f", "--fragment", action="append", type="string", help="enter the str file for fragment. each fragment has its own -f flag")

opts, args = parser.parse_args()
print "arguments:", args
print "options:", opts

molecule_files = []

molecule_files.append(opts.original_ligand[0])

if(opts.fragment):
	for fragment in opts.fragment:
		molecule_files.append(fragment)

molecules = []

#print(molecule_files)

for molecule in molecule_files:
	molecules.append(generateGraph(molecule))

for molecule in molecules:
	print(totalCharge(molecule))


	