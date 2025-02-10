import numpy as np
from scipy.spatial.distance import cdist
from pymatgen.core.structure import IStructure, Structure, Molecule
from pymatgen.core.surface import SlabGenerator, miller_index_from_sites
from pymatgen.core.composition import Composition
from pymatgen.io.cif import CifWriter
import sys, glob, os

# set binary elements
binary = ['Pt', 'Ni']

def get_num_atom_cfg():
    ele_1_num, ele_2_num, name = [], [], []
    for file in os.listdir():
        array, size = [], []
        if file.endswith('.cfg'):
            idx = file.find('.cfg')
            name.append(file[:idx])
            with open(file, 'r') as f:
                lines=f.readlines()
                ele_1, ele_2 = 0, 0
                for line in lines:
                    array.append(line.split())
                for i in range(len(array)):
                    if len(array[i]) != 0:
                        if array[i][0] == 'Size':
                            size.append(int(array[i+1][0]))
                        if len(array[i]) == 5:
                            if array[i][1] == '0':
                                ele_1 += 1
                            if array[i][1] == '1':
                                ele_2 += 1
            ele_1_num.append(ele_1)
            ele_2_num.append(ele_2)
    # write to csv file
    with open('struct.csv', 'w') as f:
        f.write('Configuration'+ ',' + 'Pt' + ',' + 'Ni' + '\n')
        for i, j, k in zip(name, ele_1_num, ele_2_num):
            f.write(str(i) + ',' + str(j) + ',' + str(k))
            f.write('\n')

def get_num_atom_cif():
    ele_1_num, ele_2_num, name = [], [], []
    for file in os.listdir():
        array, size = [], []
        if file.endswith('.cif'):
            print('file',file)
            name.append(file)
            structure = Structure.from_file(file)
            formula = structure.composition.formula
            composition = Composition(formula)
            atom_counts = composition.get_el_amt_dict()
            ele_1_num.append(atom_counts[binary[0]])
            ele_2_num.append(atom_counts[binary[1]])
            
    # write to csv file
    with open('struct.csv', 'w') as f:
        f.write('Configuration'+ ',' + f'{binary[0]}' + ',' + f'{binary[1]}' + '\n')
        for i, j, k in zip(name, ele_1_num, ele_2_num):
            f.write(str(i) + ',' + str(j) + ',' + str(k))
            f.write('\n')        
            

def coords_latt(SS_cif):
    structure = Structure.from_file(SS_cif)
    lattice = structure.lattice.matrix
    coords = structure.cart_coords
    return coords, lattice

def calculate_cluster_center(cluster):
    center = np.mean(cluster, axis=0)
    return center

def count_atoms_within_radius(cluster, radius):
    # Calculate the center of the cluster
    center = np.mean(cluster, axis=0)
    
    # Calculate the distance from each atom to the center
    distances = cdist(cluster, center.reshape(1, -1)).flatten()
    
    # Count the number of atoms within the specified radius
    num_atoms_within_radius = len(cluster[distances <= radius + np.finfo(float).eps])
    
    return num_atoms_within_radius

def get_atoms_and_distances(cluster, num_core_atom):
    # Calculate the center of the cluster
    center = np.mean(cluster, axis=0)
    
    # Calculate the distance from each atom to the center
    distances = cdist(cluster, center.reshape(1, -1)).flatten()
    
    # Get the indices that would sort the distances in ascending order
    indices_sorted = np.argsort(distances)
    
    # Get the sorted distances and coordinates of the atoms
    sorted_distances = distances[indices_sorted]
    sorted_coordinates = cluster[indices_sorted]
    core_coords = sorted_coordinates[0:num_core_atom]
    shell_coords = sorted_coordinates[num_core_atom:]
    return shell_coords, core_coords

get_num_atom_cif()
with open('struct.csv', 'r') as f:
    lines = f.readlines()
    nn = lines[0].split(',')
    ele_1, ele_2 = lines[0].split(',')[1], lines[0].split(',')[2]
    ele_2 = ele_2.replace("\n", "")
    for line in lines[1:]:
        line = line.split(',')
        name, num_ele_1, num_ele_2 = line[0], int(float(line[1])), int(float(line[2]))
        idx = name.find('.cif')
        newname = name[:idx]
        coords, lattice = coords_latt(name)
        # assign ele_1_core and ele_2_shell
        shell_coords_1, core_coords_1= get_atoms_and_distances(coords, num_ele_1)
        core_1 = [ele_1]*len(core_coords_1)
        shell_1 = [ele_2]*len(shell_coords_1)
        coreshell_1 = core_1 + shell_1
        coords_1=np.concatenate((core_coords_1, shell_coords_1))
        # Create the Structure object
        structure_1 = Structure(lattice, species=coreshell_1, coords=coords_1, coords_are_cartesian=True)
        # Write the CIF file
        cif_writer_1 = CifWriter(structure_1)
        cif_writer_1.write_file(f'{newname}_{ele_1}_core_{ele_2}_shell.cif')
    
        # assign ele_1_core and ele_2_shell
        shell_coords_2, core_coords_2= get_atoms_and_distances(coords, num_ele_2)
        core_2 = [ele_2]*len(core_coords_2)
        shell_2 = [ele_1]*len(shell_coords_2)
        coreshell_2 = core_2 + shell_2
        coords_2=np.concatenate((core_coords_2, shell_coords_2))

        # Create the Structure object
        structure_2 = Structure(lattice, species=coreshell_2, coords=coords_2, coords_are_cartesian=True)
        # Write the CIF file
        cif_writer_2 = CifWriter(structure_2)
        cif_writer_2.write_file(f'{newname}_{ele_2}_core_{ele_1}_shell.cif')

os.system(f"mkdir {ele_1}_core_{ele_2}_shell")
os.system(f'mv *{ele_1}_core_{ele_2}_shell.cif {ele_1}_core_{ele_2}_shell')

os.system(f"mkdir {ele_2}_core_{ele_1}_shell")
os.system(f'mv *{ele_2}_core_{ele_1}_shell.cif {ele_2}_core_{ele_1}_shell')
