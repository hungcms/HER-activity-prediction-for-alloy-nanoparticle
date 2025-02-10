from ase.build import molecule
import ase.io
from ase.build import surface 
from pymatgen.core.structure import Structure
from pymatgen.core.structure import Molecule
# from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.io.cif import CifWriter
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
import os
import itertools

def read_cif_file(file_path):
    atoms = ase.io.read(file_path, format='cif')
    return atoms

def get_surface_layer(atoms):    
    z_coords = [site[2] for site in atoms.get_positions()]
    sorted_z_coords = sorted(set(z_coords), reverse=True)
    for atom in atoms:
        if round(atom.position[2],2) >= round(sorted_z_coords[0],2):
            # print(atom.position[2])
            atom.tag = 1
        else:
            pass
    return atoms 

def write_to_defective_slab(atoms, name):
    index_to_del = []
    for atom in atoms:
        # print(atom.tag)
        if atom.tag == 1:
            index_to_del.append(atom.index)
        else:
            pass 

    idx = name.find('slab.cif')                                                                   
    new_name = file[0:idx] + f'_defective_'

    combinations = []
    for L in range(len(index_to_del) + 1):
        for subset in itertools.combinations(index_to_del, L):
            combinations.append(subset)
    combinations = [s for s in combinations if s and len(s)!= len(index_to_del)]
    count = 0 

    for index in combinations:
        atoms_tmp = atoms.copy()
        del atoms_tmp[list(index)]
        atoms_tmp.write(f'{new_name}{count}.cif', format='cif')
        count += 1
        
    return None

path_to_cif = os.getcwd()
for file in os.listdir():
    if file.endswith('.cif'):  
        file_path = f"{path_to_cif}/{file}"                                   
        atoms = read_cif_file(file_path)
        atoms = get_surface_layer(atoms)
        write_to_defective_slab(atoms, file)
        
