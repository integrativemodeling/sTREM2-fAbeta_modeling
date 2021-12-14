from multiprocessing import Process
import time
import resource
import numpy as np
import scipy
import random
import math
import sys, os, glob
import time

import scipy as sp
from scipy import spatial

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi
import IMP.pmi.analysis
import os
import sys
from math import sqrt
import glob


rmf_files = glob.glob("*.rmf3")
#rmf_files = ['0.rmf3', '2584.rmf3',  '2728.rmf3',  '2872.rmf3',  '3015.rmf3',  '315.rmf3',   '3303.rmf3',  '441.rmf3',   '586.rmf3',  '72.rmf3',   '874.rmf3']
for rmf_ref in rmf_files:
    print (rmf_ref)
    m = IMP.Model()
    h_ref = IMP.pmi.analysis.get_hiers_from_rmf(m,0,rmf_ref)[0]
    print(h_ref)
    s0 = IMP.atom.Selection(h_ref, resolution=1)
    print(s0)
    pdb_name = rmf_ref.strip(".rmf3")+".cif"
    o = IMP.pmi.output.Output()
    o.init_pdb(pdb_name,h_ref, mmcif=True)
    o.write_pdb(pdb_name)
    del o
