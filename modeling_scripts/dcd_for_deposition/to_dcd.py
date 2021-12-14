#!/usr/bin/env python

# Adapted for the deposition of the 

"""This is a simple script to convert an IMP cluster (multiple single-frame
   RMF files) into DCD (CHARMM/NAMD) trajectories, with each bead or residue
   represented as a single 'atom', in residue number order (the same order
   as in the mmCIF file).

   We use the updated version of MDTools that's bundled with Chimera, so you
   may need to change sys.path accordingly, below.
"""

from __future__ import print_function
import os
import ast
import sys
sys.path.append('/opt/chimera-1.14-1.fc31/share/Trajectory/DCD/MDToolsMarch97/')
import md
import IMP.pmi.output
import IMP.rmf
import RMF

class IMPCluster(object):
    def __init__(self, directory, list_file):
        self.directory = directory
        with open(list_file) as fh:
            self.all_rmfs = [os.path.basename(x.rstrip(' \r\n'))
                             for x in fh.readlines()]

    def get_size(self):
        """Get number of models in the cluster"""
        return len(self.all_rmfs)

    def get_rmf_file(self, model_num):
        """Get path to RMF file for the given model in the cluster"""
        return os.path.join(self.directory, self.all_rmfs[model_num])


class DCDOutput(object):
    """Dump a series of IMP Hierarchies to a DCD file."""

    def __init__(self, fname):
        self.fname = fname
        self._d = None

    def dump(self, mhs):
        coords = self._get_coords(mhs)
        if self._d is None:
            self._init_dcd(len(coords))
        assert(len(coords) == len(self._ag.atoms))
        for i, coord in enumerate(coords):
            # Everything should be non-atomic for now
            assert(coord[1] is None)
            a = self._ag.atoms[i]
            a.x, a.y, a.z = coord[0]
        self._d.append()

    def _init_dcd(self, n_coords):
        self._ag = md.AtomGroup()
        for i in range(n_coords):
            self._ag.atoms.append(md.Atom())
        self._d = md.DCDWrite(self.fname, self._ag)

    def _get_coords(self, mhs):
        assert(len(mhs) == 1)
        o = IMP.pmi.output.Output(atomistic=True)
        name = 'dcd-output'
        o.dictionary_pdbs[name] = mhs[0]
        o._init_dictchain(name, mhs[0])
        coords, geometric_center = o.get_particle_infos_for_pdb_writing(name)
        return coords


if len(sys.argv) != 5:
    print("Usage: %s [all models directory] [cluster list file] "
          "[output DCD file] [percentage of the Cluster]" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

d = DCDOutput(sys.argv[3])
cluster = IMPCluster(sys.argv[1], sys.argv[2])

cluster_size = cluster.get_size()
if sys.argv[4] < 1:
    print("fifth argument is percentage of the cluster to extract, should be a number between 1-100")
    sys.exit()
else:
    percentage_to_extract = float(sys.argv[4])
for i in range(int(cluster_size*percentage_to_extract/100)):
    m = IMP.Model()
    r = RMF.open_rmf_file_read_only(cluster.get_rmf_file(i))
    mhs = IMP.rmf.create_hierarchies(r, m)
    m.update()
    print("Adding coordinates for model %d of %d" % (i + 1, cluster_size))
    d.dump(mhs)
