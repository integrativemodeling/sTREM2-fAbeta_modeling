import IMP
import IMP.atom
import IMP.isd
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.io
import IMP.pmi.io.crosslink
import sys
import os

# Imports needed to use ProtocolOutput
import IMP.pmi.mmcif
import ihm
import ihm.location
import ihm.model
import ihm.cross_linkers
import ihm.dumper




"""
Author: Dibyendu Mondal
email: dibyendu@salilab.org
"""

#########################################################################################
#-------------------------- Some useful functions --------------------------------------
#########################################################################################


def create_fibril_transforms(z_displacement, len_protofilament):
    """
    A function to compute transformation matrix to elongate the Abeta fibril.
    Assuming z-axis as the fibril growth axis
    :param z_displacement: gap length between two fibril chain layers along z axis (in Angs) [default: 4.8] 
    :param len_protofilament: number of fibril chain layers [default: 15]
    :return: transformation matrix
    """
    transforms = []
    half_len = int(len_protofilament//2) 
    if len_protofilament % 2 ==0:
        for t in range(-half_len, half_len, 1):
            if t != 0: #t = 0 will be the original copy no need to store the vector 
                trans_even = IMP.algebra.Transformation3D(IMP.algebra.Vector3D(0, 0, ((t * z_displacement))))
                transforms.append(trans_even)
                print("Transforms:", t, " | ", trans_even)
    else:
        for t in range(-1 * half_len, half_len + 1, 1):
            if t !=0:
                trans_odd = IMP.algebra.Transformation3D(IMP.algebra.Vector3D(0, 0, ((t * z_displacement))))
                transforms.append(trans_odd)
                print("Transforms:", t, " | ", trans_odd)
       

    return transforms

#########################################################################################
# ------------------- Define all data path and other variables -------------------------
#########################################################################################

data_dir = "../data"
pdb_dir = data_dir + "/" + "pdb/"
fasta_dir = data_dir + "/" + "fasta/"
xl_data_NHSF = data_dir + "/" + "xl/bs3_nshsf_mixed.csv"  # all xl-ms data in one file (BS3 and NSHSF)
outdir = "./only_xl/run"
output = '/output'

# restraints variables
xl_weight = 1.0
xl_length = 22.5
xl_slope = 0.02
spherical_barrier_radius = 50
excld_vol_trem2_abeta_weight = 10

# sampling parameters
max_translation = 1.0
max_rotation = 0.5
max_shuff_translation = 150.0
replica_temp_min = 1.5
replica_temp_max = 2.5
SA_temp_min = 1.0
SA_temp_max = 2.5
mc_steps = 10
num_frames = 50000
optimize_flex_beads_steps = 500

#########################################################################################
# -----------------------------CREATE SYSTEM AND STATE ---------------------------------
#########################################################################################

mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
bs = IMP.pmi.macros.BuildSystem(mdl)


if '--mmcif' in sys.argv:
# Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput()
    po.system.title = ('Integrative structure of Abeta-sTREM2 complex')
#    po.system.citations.append(ihm.Citation.from_pubmed_id(29531062))
    s.add_protocol_output(po)

st = s.create_state()


# Store the FASTA sequences in a dictionary
sequences = IMP.pmi.topology.Sequences(fasta_dir + "TREM2_ABfibril.fasta")

# Generate an ACTIVE Abeta Trimer
# Create a molecule from 1 trimer (chain A) of the Amyloid fibril (2LMP) & copy it to build the ACTIVE trimer.

Abeta = st.create_molecule("Abeta", sequence=sequences["Abeta40"], chain_id='A')
AbetaB = Abeta.create_copy('B')
AbetaC = Abeta.create_copy('C')
AbA = Abeta.add_structure(pdb_dir + "2LMP_D_J_P.pdb", chain_id='A', offset=0)
AbB = AbetaB.add_structure(pdb_dir + "2LMP_D_J_P.pdb", chain_id='B', offset=0)
AbC = AbetaC.add_structure(pdb_dir + "2LMP_D_J_P.pdb", chain_id='C', offset=0)

# Add representation to the monomers in the ACTIVE Abeta trimer

# Chain A monomer
Abeta.add_representation(AbA, resolutions=[1], color=0.1)
Abeta.add_representation(Abeta.get_non_atomic_residues(), resolutions=[1], color=0.2)

# Chain B monomer
AbetaB.add_representation(AbB, resolutions=[1], color=0.1)
AbetaB.add_representation(AbetaB.get_non_atomic_residues(), resolutions=[1], color=0.3)

# Chain C monomer
AbetaC.add_representation(AbC, resolutions=[1], color=0.1)
AbetaC.add_representation(AbetaC.get_non_atomic_residues(), resolutions=[1], color=0.4)

# Fibril extension
# Create the transforms needed to build up the fibril from symmetry z-displacement translation

fibril_transforms = create_fibril_transforms(z_displacement=4.8, len_protofilament=15)  
print('Transforms  || ', len(fibril_transforms), " || ", fibril_transforms)

# Create a clone of the trimer fibril unit for every transform
# Collect all trimer fibril units created in a list

fibril_mols = [Abeta, AbetaB, AbetaC]
protofil_mols = [[Abeta], [AbetaB], [AbetaC]]

abeta_chains = ['A', 'B', 'C']
additional_chains = [['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14', 'A15'],
                       ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12', 'B13', 'B14', 'B15'],
                       ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15']]

# somehow multiletter chain can not be handled, switching back to one-letter code
#additional_chains = [['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'],
#                       ['B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'],
#                       ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']]

for i in range(0, len(abeta_chains)):
    for j in range(0, len(fibril_transforms)):
        # print("Creating a clone with this chain id:", additional_chains[i][j])
        clone = fibril_mols[i].create_clone(chain_id=additional_chains[i][j])
        protofil_mols[i].append(clone)
        fibril_mols.append(clone)

# -----------------------TREM2-----------------------
# Create a molecule for 1 TREM2 Extra Cellular Domain (based on 5UD7, chain A)
trem2 = st.create_molecule('TREM2', sequence=sequences["TREM2"])

T = trem2.add_structure(pdb_dir + "TREM2_A.pdb",
                        chain_id='X',
                        offset=0)  # offset -18
trem2.add_representation(trem2.get_atomic_residues(), resolutions=[1], color=0.5)

# _____ Build the system & degrees of freedom _______________________________________________
root_hier = s.build()
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

# ----------- DEFINE THE SYSTEM RIGID BODIES & Set up connectivity restraints -----------------------------
# Lists for collecting molecules/particles
TREM2_mol = []
active_tri_mol = []
fib_ext_mol = []
shuffle_exclude_rbs = []

# Lists useful for collecting restraints
output_objects = []
sample_objects = []
crs = []

# get a dictionary of the molecules in the system w/ KEY = molecule name
moldict = st.get_molecules()
print(moldict)

for molname in moldict:
    for mol in moldict[molname]:
        if 'TREM2' in molname:
            atomic = mol.get_atomic_residues()
            dof.create_rigid_body(atomic,
                                  max_trans=max_translation,
                                  max_rot=max_rotation,
                                  resolution='all')
            TREM2_mol.append(mol)

        elif 'Abeta' in molname:
            crA = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
            crA.add_to_model()
            output_objects.append(crA)
            sample_objects.append(crA)
            crs.append(crA)

            atomic = mol.get_atomic_residues()
            dof.create_rigid_body(atomic,
                                  max_trans=0,
                                  max_rot=0,
                                  resolution='all')
            dof.create_flexible_beads(mol.get_non_atomic_residues(),
                                      max_trans=max_translation,
                                      resolution=1)
            active_tri_mol.append(mol)


#########################################################################################
# ------------------------Composite Constrains for Fibril-------------------------------
#########################################################################################

# Constrain the fibril trimer unit copies with z-trans. symmetry ### This is very important
for i in range(0, len(abeta_chains)):
    for t in range(0, len(fibril_transforms)):
        threefold_trans = fibril_transforms[t]
        dof.constrain_symmetry(protofil_mols[i][0], protofil_mols[i][t + 1], threefold_trans)
mdl.update()

#########################################################################################
# -----------------------------ADD RESTRAINTS -----------------------------------------
#########################################################################################

# 1. External Barrier Restraint

# get the center of mass of the fibril then
# set up an an external barrier sphere with radius equal to 
# the length of the fibril from the center of mass to
# exclude docking from fibril ends

Fibril_sel = IMP.atom.Selection(root_hier, molecule="Abeta")
fib_particles = Fibril_sel.get_selected_particles()
if len(fib_particles) == 0:
    print("COM not set up. Cannot select protein %s)" & Abeta)

fibril_COM = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(mdl), fib_particles)
fibril_COM_coor = IMP.core.XYZ(fibril_COM).get_coordinates()

eb = IMP.pmi.restraints.basic.ExternalBarrier(hierarchies=root_hier, radius=spherical_barrier_radius, center=fibril_COM_coor, label='barrier')
eb.add_to_model()
output_objects.append(eb)

# 2. Excluded Volume Restraint

ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=TREM2_mol, resolution=1)
ev1.add_to_model()
ev1.set_label('TREM2')
output_objects.append(ev1)


ev2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=TREM2_mol,
                                                              other_objects=[fibril_mols],
                                                              resolution=1)
ev2.add_to_model()
ev2.rs.set_weight(excld_vol_trem2_abeta_weight)
ev2.set_label('all')
output_objects.append(ev2)

# 3. Cross-linking Restraint

xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()

# Set up XL restraint for all XLs
crosslink_restraints = []
xldb_NHSF = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_NHSF.create_set_from_file(file_name=xl_data_NHSF, converter=xldbkc)
xlrN = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=root_hier,            # Must pass the root hierarchy to the system
    database=xldb_NHSF,    # The cross-link database.
    length=xl_length,               # The cross-linker plus side chain length
    resolution=1,                   # The resolution at which to evaluate the cross-link
    slope=xl_slope,
    label='NHSF',                   # This adds a linear term to the scoring function to bias cross-links towards each other
    weight=xl_weight)               # Scaling factor for the restraint score.
xlrN.add_to_model()
output_objects.append(xlrN)
crosslink_restraints.append(xlrN)

#########################################################################################
# ----------------------------- SAMPLING -----------------------------------------------
#########################################################################################
shuffle_exclude_rbs = dof.get_rigid_bodies()[:-1]



if '--mmcif' in sys.argv:
    num_frames=20
    start = 0
    end = 1
    optimize_flex_beads_steps = 1

run_num = 1
global_output_directory = outdir + str(run_num) + output
print(global_output_directory)

IMP.pmi.tools.shuffle_configuration(root_hier, excluded_rigid_bodies=shuffle_exclude_rbs,
        max_translation=max_shuff_translation)
dof.optimize_flexible_beads(optimize_flex_beads_steps)
out = IMP.pmi.output.Output()

rex = IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        output_objects=output_objects,
        monte_carlo_temperature=1.0,
        simulated_annealing=True,
        simulated_annealing_minimum_temperature=SA_temp_min,
        simulated_annealing_maximum_temperature=SA_temp_max,
        simulated_annealing_minimum_temperature_nframes=200,
        simulated_annealing_maximum_temperature_nframes=20,
        replica_exchange_minimum_temperature=replica_temp_min,
        replica_exchange_maximum_temperature=replica_temp_max,
        number_of_best_scoring_models=0,
        monte_carlo_steps=mc_steps,
        number_of_frames=num_frames,
        global_output_directory=global_output_directory,
        test_mode=True)
rex.execute_macro()

po.finalize()

sp = po.system

print(dir(s))
#print(s.citations)
with open('initial.cif', 'w') as fh:
    ihm.dumper.write(fh, [sp])

# Datasets for XL-MS restraint
for r in sp.restraints:
    if isinstance(r, ihm.restraint.CrossLinkRestraint):
        print("XL-MS dataset at:", r.dataset.location.path)
        print("Details:", r.dataset.location.details)

# Get last step of last protocol (protocol is an 'orphan' because
# we haven't used it for a model yet)
last_step = sp.orphan_protocols[-1].steps[-1]


# Correct number of output models to account for multiple runs
last_step.num_models_end = 3000000

# Get last protocol in the file
protocol = po.system.orphan_protocols[-1]
# State that we filtered the 3000000 frames
analysis = ihm.analysis.Analysis()
protocol.analyses.append(analysis)
analysis.steps.append(ihm.analysis.ClusterStep(
                      feature='RMSD', num_models_begin=3000000,
                      num_models_end=100))
