import IMP
import IMP.atom
import IMP.rmf
import RMF
import glob
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import multiprocessing as mp
from matplotlib import rc
import sys

# Set the global font to be DejaVu Sans, size 10 (or any other sans-serif font of your choice!)
rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':20, 'weight': 'bold'})

# Set the font used for MathJax - more on this later
rc('mathtext',**{'default':'regular'})

def get_rmfs_coordinates(str_file, chain_names, subunit_name):

#    conform = []
    num = 1
    masses = []
    radii = []
    ps_names = []

    models_name = []
   
#    print(str_file)
    models_name.append(str_file)

    m = IMP.Model()
    inf = RMF.open_rmf_file_read_only(str_file)
    h = IMP.rmf.create_hierarchies(inf, m)[0]
    IMP.rmf.load_frame(inf, 0)

    pts = []

    s0 =[]
    if subunit_name == 'Abeta':
        for indx_1 in range(len(chain_names)):
            for indx_2 in range(len(chain_names[0])):
                s0 += IMP.atom.Selection(h, resolution=1,molecule=subunit_name, chain_id = chain_names[indx_1][indx_2]).get_selected_particles()[8:]
    else:
        s0 = IMP.atom.Selection(h, resolution=1, molecule=subunit_name).get_selected_particles()
#    print(s0)
    for leaf in s0:

        p=IMP.core.XYZR(leaf)
        pts.append([p.get_coordinates()[i] for i in range(3)])

    if num == 0:
        masses.append(IMP.atom.Mass(leaf).get_mass())
        radii.append(p.get_radius())
        mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(leaf))
        # traverse up the Hierarchy to get copy number of the molecule that the bead belongs to
        copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(leaf))

        if IMP.atom.Fragment.get_is_setup(leaf): #TODO not tested on non-fragment systems
            residues_in_bead = IMP.atom.Fragment(leaf).get_residue_indexes()

            ps_names.append(mol_name+"_"+str(min(residues_in_bead))+"_"+str(max(residues_in_bead))+"_"+str(copy_number))

        else:
            residue_in_bead = str(IMP.atom.Residue(leaf).get_index())

            ps_names.append(mol_name+"_"+residue_in_bead+"_"+residue_in_bead+"_"+str(copy_number))

    conform = pts
    pts = []
    num = num + 1

    return ps_names, masses, radii, np.array(conform), models_name


def slice_con_map(contact_map):
    con_map_list = []
    len_chain = int(contact_map.shape[0]/3)
    con_map_A = contact_map[:len_chain,]
    con_map_list.append([con_map_A, con_map_A.ravel().sum()])
    con_map_B = contact_map[len_chain:2*(len_chain),]
    con_map_list.append([con_map_B,con_map_B.ravel().sum()])
    con_map_C = contact_map[2*(len_chain):,]
    con_map_list.append([con_map_C, con_map_C.ravel().sum()])
    con_map_list.sort(key=lambda x:x[1])
    return(con_map_list[-1][0])


def main_worker_contact_map(file_name): 
    ps_names_Abeta, masses_Abeta, radii_Abeta, conform_Abeta, models_name_Abeta = get_rmfs_coordinates(file_name, chain_names, subunit_name='Abeta')
    ps_names_TREM2, masses_TREM2, radii_TREM2, conform_TREM2, models_name_TREM2 = get_rmfs_coordinates(file_name, chain_names, subunit_name='TREM2')
    distance_mat = distance_matrix(conform_Abeta, conform_TREM2)
    # a contact exist if the distance is less than 12 A
    b = distance_mat < 12.0
    contact_map = b.astype(int)
    contact_map_to_check = slice_con_map(contact_map).ravel()
    return contact_map_to_check


def main_worker2_dist_mat(index, data):
    padding_width = index + 1
    data_1 = data[index:index+1,]
    data_rest =data[index+1:,]
    difference_data = np.abs((data_rest - data_1).sum(axis=1))
    padded_dif_data = np.pad(difference_data, (padding_width, 0), 'constant', constant_values = (0,0))
    return padded_dif_data


# My models contain a Fibril with C3 symmetry and translational 
# There are 15 3-fold chains

chain_names_starts = ["A", "B", "C"] # 0 degrees rotation along C3 axis
chain_names = [[], [], []]

#assuming there are 15 chains and the chains are numbered as !!!!! CHECK YOUR .rmf3 FILES !!!!!

#[['A1', 'B1', 'C1'], ['A2', 'B2', 'C2'], ['A3', 'B3', 'C3'], ['A4', 'B4', 'C4'], ['A5', 'B5', 'C5'], ['A6', 'B6', 'C6'], ['A7', 'B7', 'C7'], ['A', 'B', 'C'], ['A8', 'B8', 'C8'], 
#['A9', 'B9', 'C9'], ['A10', 'B10', 'C10'], ['A11', 'B11', 'C11'], ['A12', 'B12', 'C12'], ['A13', 'B13', 'C13'], ['A14', 'B14', 'C14']]
for e in range(len(chain_names_starts)):
    for f in range(0, 15):
        if f < 7:
                chain_names[e].append(chain_names_starts[e]+str(f+1))

        elif f >= 8:
                chain_names[e].append(chain_names_starts[e]+str(f))

        else:
                chain_names[e].append(chain_names_starts[e])


#print(chain_names)

list_rmf = sorted(glob.glob("*.rmf3"))
pool = mp.Pool(int(mp.cpu_count()/2))
#print(pool)

num_structures = len(list_rmf)
a_list = pool.map(main_worker_contact_map, [file_name for file_name in list_rmf])
pool.close()

array_shape = a_list[0].shape
sum_ = np.zeros(array_shape)
for f in a_list:
    sum_ += f
avg = sum_/num_structures
col_label = []

# just ploting residue from 21 to 130
for i in range(21, 130):
    col_label.append(str(i))

df_contact = pd.DataFrame(avg.reshape(480, 109), columns=col_label)
img = plt.pcolor(df_contact)
plt.xticks(np.arange(0.5, len(df_contact.columns), 20), ['21', '41', '61', '81', '101', '121'])
plt.yticks(np.arange(15, len(df_contact.index), 32), ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14'])
cbar = plt.colorbar(img)
cbar.minorticks_on()
plt.tight_layout()
file_name = "clusterXX_contact_map_avg" + ".pdf"
plt.savefig(file_name)
print('done')

