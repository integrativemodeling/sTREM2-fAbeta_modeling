## \example rmf/pdb.py
# Write a PDB to an RMF file.
#

from __future__ import print_function
import IMP.atom
import IMP.rmf
import RMF
import sys


m = IMP.Model()
pdbname = sys.argv[1]
# Create a new IMP.atom.Hierarchy from the contents of the pdb file
h = IMP.atom.read_pdb(pdbname, m,  IMP.atom.CAlphaPDBSelector())
h.set_name('System')
tfn = pdbname.split('.')[0] + '.rmf3'

print("File name is", tfn)

# open the file, clearing any existing contents
rh = RMF.create_rmf_file(tfn)

# add the hierarchy to the file
IMP.rmf.add_hierarchies(rh, [h])

# add the current configuration to the file as frame 0
IMP.rmf.save_frame(rh)
