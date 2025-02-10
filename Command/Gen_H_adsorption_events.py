from pymatgen.core.structure import Structure
from pymatgen.core.structure import Molecule
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.io.cif import CifWriter
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
import os

# Create a Structure object for your surface slab

def sort_2_layer(structure):
    z_coords = [site.coords[2] for site in structure]
    sorted_z_coords = sorted(set(z_coords), reverse=True)
    top_two_layers = []
    for site in structure:
        if site.coords[2] >= sorted_z_coords[0]:
            top_two_layers.append(site.coords)
        elif site.coords[2] >= sorted_z_coords[1]:
            top_two_layers.append(site.coords)

    return top_two_layers

def sort_surface_layer(structure):
    z_coords = [site.coords[2] for site in structure]
    sorted_z_coords = sorted(set(z_coords), reverse=True)
    surface_layers = []

    for site in structure:

        if round(site.coords[2],1) >= round(sorted_z_coords[0],1):
            surface_layers.append(site.coords)
        else:
            pass

    return surface_layers

path_to_cif = os.getcwd()
os.chdir(path_to_cif)
#Read file function
def read_cif_file(file_path):
    structure = Structure.from_file(file_path)
    return structure


def calculate_average(coords):
    """Calculate the average of a list of 3D coordinates."""
    # Convert the list of coordinates to a NumPy array
    coords_array = np.array(coords)
    # Calculate the mean of each column (i.e. x, y, z)
    average = np.mean(coords_array, axis=0)
    return average.tolist()
    
#Read file
for file in os.listdir():
    if file.endswith('.cif'):
        file_path = f"{path_to_cif}/{file}"
        slab_structure = read_cif_file(file_path)
        surface_layer = sort_surface_layer(slab_structure)
        molecule = Molecule(['H'], [[0, 0, 0]])

        for k in range(len(surface_layer)):
            surface_layer[k][2] = surface_layer[k][2] + 1.5

        for count, site in enumerate(surface_layer):
            slab_structure_ = AdsorbateSiteFinder(slab_structure).add_adsorbate(molecule, site, translate=True)
            idx = file.find('.cif')
            new_name = file[0:idx] + f'{count}_ontop'
            CifWriter(slab_structure_).write_file(f'{new_name}.cif')


os.system('mkdir file')
os.system('mv *ontop.cif file')


