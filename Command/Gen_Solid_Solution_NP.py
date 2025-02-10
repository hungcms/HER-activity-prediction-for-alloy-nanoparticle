import numpy as np
from scipy.spatial.distance import cdist
from pymatgen.core.structure import IStructure, Structure, Molecule
from pymatgen.core.surface import SlabGenerator, miller_index_from_sites
from pymatgen.io.cif import CifWriter
import sys, glob, os
from pymatgen.analysis.phase_diagram import ConvexHull

def coords_latt(SS_cif):
    structure = Structure.from_file(SS_cif)
    structure.make_supercell([15,15,15], to_unit_cell=True)
    lattice = structure.lattice.matrix
    coords = structure.cart_coords
    # lattice = structure_cluster.lattice.matrix
    return coords, lattice, structure

def get_input():
    if (len(sys.argv) == 3):
        ss_cif = sys.argv[1]
        radius = sys.argv[2]
    else:
        print("Missing the parameter")
        print("Usage: python [command] [ss_cif] [cluster radius]")
        sys.exit()
    return ss_cif, radius

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

def get_atoms_and_distances(supercell, radius, structure):
    # Calculate the center of the cluster
    center = np.mean(supercell, axis=0)
    
    # Calculate the distance from each atom to the center
    distances = cdist(supercell, center.reshape(1, -1)).flatten()
    
    # Get the indices that would sort the distances in ascending order
    indices_sorted = np.argsort(distances)
    sorted_distances = distances[indices_sorted]
    sorted_coordinates = supercell[indices_sorted]

    # ====new center======
    center_new = sorted_coordinates[0]
    # Calculate the distance from each atom to the center
    distances_new = cdist(supercell, center_new.reshape(1, -1)).flatten()
    # Get the indices that would sort the distances in ascending order
    indices_sorted_new = np.argsort(distances_new)
    sorted_distances_new = distances_new[indices_sorted_new]
    sorted_coordinates_new = supercell[indices_sorted_new]

    symbol_list = [structure.species[i].symbol for i in indices_sorted_new]
    final_symbol = []
    final_coords = []
    for i, j, g in zip(sorted_distances_new, sorted_coordinates_new, symbol_list):
        if float(i) <= float(radius):
            final_symbol.append(g)
            final_coords.append(j)
        else:
            pass   
    return final_symbol, final_coords

ss_cif, radius = get_input()
if ss_cif != 'all':
    coords, lattice, structure = coords_latt(ss_cif)

    symbol_list, final_coords = get_atoms_and_distances(coords, radius, structure)
    # Create the Structure object
    as_latt = [[60,0,0], [0,60,0], [0,0,60]] # define the lattice parameters for SNPs
    molecule = Molecule(species=symbol_list, coords=final_coords)
    structure = molecule.get_boxed_structure(a=60, b=60, c=60)
    # Write the CIF file
    cif_writer = CifWriter(structure)
    cif_writer.write_file(f'SS_radius_{radius}.cif') # custom filename of the CIF file
elif ss_cif == 'all':
    for file in os.listdir():
        if file.endswith('.cif'):
                coords, lattice, structure = coords_latt(file)
                symbol_list, final_coords = get_atoms_and_distances(coords, radius, structure)
                # Create the Structure object
                as_latt = [[60,0,0], [0,60,0], [0,0,60]]
                molecule = Molecule(species=symbol_list, coords=final_coords)
                structure = molecule.get_boxed_structure(a=60, b=60, c=60)
                # Write the CIF file
                cif_writer = CifWriter(structure)
                idx = file.find('.cif')
                new_name = file[:idx] + f'_{radius}nm'
                cif_writer.write_file(f'{new_name}.cif') # custom filename of the CIF file

os.system('mkdir file')
os.system(f'mv *_{radius}nm.cif file')

directory = os.getcwd()
cif_files = [f for f in os.listdir(directory) if f.endswith('.cif')]

new_names = []
# Iterate through each CIF file
for index, cif_file in enumerate(cif_files):
    # Parse the CIF file to get the chemical formula
    parser = os.path.join(directory, cif_file)
    structure = Structure.from_file(parser)
    chemical_formula = structure.composition

    # Generate the new file name
    new_name = f'{chemical_formula}_{index}.cif'
    new_names.append(new_name)
    # Rename the CIF file by composition
    os.rename(os.path.join(directory, cif_file), os.path.join(directory, new_name))

    print(f'Renamed {cif_file} to {new_name}')