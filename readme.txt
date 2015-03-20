Stitch (formerly named Lilo) is a program for stitching together fragments of a parent molecule into a larger whole. The workflow
for generating input for Stitch is documented in "lilo_spec.*", located in this directory. As an overview, first you might have a 
.mol2 file of a ligand you want to parameterize with CHARMM36 (e.g., beta-FOA). You would then submit this mol2 file to ParamChem,
but will almost certainly be hit with very high penalties for partial charges and for bonded parameters. This is largely due to 
suboptimal ligand fragmentation done internally by ParamChem. To improve parameters, a sensible approach is to use your human knowledge
and intuition to make your own ligand fragmentation, and then submit each individual fragment to ParamChem. Now you will have several
.str (stream) files containing parameters both for your original, parent ligand as well as for the individual fragments of that ligand. 
Stitch automates the process of re-assembling those fragments with significantly lower penalties into a new assembled system.

Currently (this will be improved over time), Stitch takes as input:

flag -o "original_file_name" denoting the .str file containing ParamChem's original output for the parent molecule
a series of flag -f "fragment_i_file_name" where "i" is a number ranging from 1 to n, where n is the number of fragments.
a flag -s "title", the title of the file you would like to output

and Stitch currently outputs:

A text file containing rows in format:

ATOM ATOM_NAME ATOM_TYPE NEW_PARTIAL_CHARGE

which you can then copy paste into the original .str file for the parent ligand to update the partial charges.

SAMPLE USAGE (from test_files directory)

python ../main.py -o bfoa_paramchem.str -f bf0_lys_frag1.str -f bfna_nonbonded_methyl_frag2.str -f bfna_frag_3_alkene.str -s test_output.top 

********
In the next version, we will add: 

Automatic generation of a whole new .str file so you won't have to copy-paste anything

Assignation of atom types and bonded parameters from the fragments to the newly stitched ligand.

Please email enf@stanford.edu if you have questions.

Happy simulating :-)
Evan
