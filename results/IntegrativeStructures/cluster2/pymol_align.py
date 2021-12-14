import pymol
import glob
import time
pymol.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
# pymol.pymol_argv = ['pymol']  #with GUI visualization
import pickle

from itertools import combinations
from itertools import permutations
import collections
import sys
import subprocess
""""
# The Algorithma
1. Load in all of the complex molecules in mmCIF format.
2. The TREM2 is chainless, assign chain 'X'. [x]
3. Decompose the Abeta fibril to alpha, beta, gamma components and get chain assignments.
4. delete the flexible beads/residues, so as to do superfitting. []
5. Find the combination of chains, which when superimposed, get the lowest RMSD of the TREM2 and the chosen TREM2 chains.

"""

def super_fit(mol1, mol2, mol3):

    # Currently assigned chains to Abeta molecules.
    # abeta_chains = ['A', 'G', 'M']
    # modified abeta_chains = ['A', 'B', 'G'] as alpha, beta, gamma

    obj_list = pymol.cmd.get_object_list(selection='(all)')

    if mol1 in obj_list:
        pass
    else:
        pymol.cmd.load(mol1, object=mol1.strip(".pdb"))
        mol1 = mol1.strip(".pdb")
        assigned_chains = pymol.cmd.get_chains(mol1)

        # [1]TREM2 doesn't chain id, assign chain "X" to it.
        # alter (chain A),chain='B'
        # chain " ", chain='X'
        pymol.cmd.alter('chain " "', 'chain = "X"')

        # [2]Delete flexible beads within the Abeta
        # the flexible residues range from
        for chain in assigned_chains[1:]:  # Exclude TREM2 chain
            selection_string = "resi 1-8 in chain " + chain
            # print(chain, selection_string)
            pymol.cmd.remove(selection_string)

        # [3] Rename individual faces of the fibril to A0,A1.., B0,B1.., G0,G1...

    alpha_face_mdl1 = ['A', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q']
    beta_face_mdl1  = ['B', 'S', 'T', 'U', 'V', 'W', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
    gamma_face_mdl1 = ['C', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x']


    # [4] Find the optimal RMSD superposition between pairs ([ [[A0, B0, G0, X], [A1, B1, G1, X]],[ ...).
    # combinations of three chains at a time.
    # 15 * 15 = 225 different positions to calculate RMSD.
    # First get the basic rmsd calculation between chain calculation
    # sele D, obj1 and chain A8 + chain B8 + chain G8 + chain X
    # pair_fit A, A

    pymol.cmd.load(mol2, object=mol2.strip(".cif"))
    mol2 = mol2.strip(".cif")
    assigned_chains = pymol.cmd.get_chains(mol2)

    # [1]TREM2 doesn't chain id, assign chain "X" to it.
    # alter (chain A),chain='B'
    # chain " ", chain='X'
    pymol.cmd.alter('chain " "', 'chain = "X"')

    # [2]Delete flexible beads within the Abeta
    # the flexible residues range from
    for chain in assigned_chains[1:]:  # Exclude TREM2 chain
        selection_string = "resi 1-8 in chain " + chain
        # print(chain, selection_string)
        pymol.cmd.remove(selection_string)

    # [3] Rename individual faces of the fibril to A0,A1.., B0,B1.., G0,G1...

    alpha_face_mdl2 = ['A', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14']
    beta_face_mdl2  = ['B', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12', 'B13', 'B14']
    gamma_face_mdl2 = ['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14']

    pymol.cmd.load(mol3, object=mol3.strip(".pdb"))
    mol3 = mol3.strip(".pdb")
    assigned_chains = pymol.cmd.get_chains(mol3)
    pymol.cmd.alter('chain " "', 'chain = "X"')





    # [4] Find the optimal RMSD superposition between pairs ([ [[A0, B0, G0, X], [A1, B1, G1, X]],[ ...).
    # combinations of three chains at a time.
    # 15 * 15 = 225 different positions to calculate RMSD.
    # First get the basic rmsd calculation between chain calculation
    # sele D, obj1 and chain A8 + chain B8 + chain G8 + chain X
    # pair_fit A, A
    chain_modl1 = [alpha_face_mdl1, beta_face_mdl1, gamma_face_mdl1]
    chain_modl2 = [alpha_face_mdl2, beta_face_mdl2, gamma_face_mdl2]

    unique_chain = ['A', 'B', 'C']
    i = 0
#    for i in range(0, len(alpha_face_mdl2)):
            # Translate Alpha, beta and gamma chains
    pymol.cmd.pair_fit(mol1 + '//' + unique_chain[0] + '/9-40/CA', mol2 + '//'+unique_chain[0] + '/9-40/CA',
                                     mol1 + '//'+ unique_chain[1]+ '/9-40/CA', mol2 + '//'+unique_chain[1] + '/9-40/CA',
                                     mol1 + '//'+ unique_chain[2] + '/9-40/CA', mol2 + '//'+unique_chain[2] + '/9-40/CA')
    rmsd = pymol.cmd.rms_cur(mol1 + '//' + unique_chain[0] + '/9-40/CA', mol2 + '//'+unique_chain[0] + '/9-40/CA')
#            if rmsd <=0.0:
#                break
                

    pymol.cmd.pair_fit(mol3 + '//X/21-129/CA', mol2 + '//X/21-129/CA')
    Abeta_name = mol1 + "_"+mol2+"_.pdb"
    TREM2_name = mol3 + "_"+mol2+"_.pdb"
    pymol.cmd.save(Abeta_name, mol1)
    pymol.cmd.save(TREM2_name, mol3)
#    os.system('cat Abeta_name TREM2_name >> final_name')
    return Abeta_name, TREM2_name


if __name__ == '__main__':
    file_name_pairs = []
    pymol.finish_launching()
    pymol.cmd.feedback("disable", "all", "everything")
    # time.sleep(10)
#    all_cifs = glob.glob("*.cif")
    mol1 = 'Abeta_fibril.pdb'
    mol3 = 'TREM2_A.pdb'
    mol2 = sys.argv[1]
    print(mol2)
    final_name = mol2.split('.cif')[0]+"_final.pdb"
    print(final_name)
    Abeta_name, TREM2_name = super_fit(mol1, mol2, mol3)
    file_name_pairs.append([Abeta_name, TREM2_name])

    # Python program to
    # demonstrate merging
    # of two files

    data = data2 = ""

    # Reading data from file1
    with open(Abeta_name) as fp:
        data = fp.read()

    # Reading data from file2
    with open(TREM2_name) as fp:
        data2 = fp.read()

    # Merging 2 files
    # To add the data of file2
    # from next line
    data += "\n"
    data += data2

    with open (final_name, 'w') as fp:
        fp.write(data)
    subprocess.run(['rm', Abeta_name, TREM2_name])
    pymol.cmd.delete(mol2)


