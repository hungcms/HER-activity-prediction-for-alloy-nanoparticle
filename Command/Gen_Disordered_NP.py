from ase.io import read, write
import random
import numpy as np
import sys, os, glob

def get_input():
    if (len(sys.argv) != 2):
        print("Missing the parameter, Please give me more")
        print("Usage: python [command] [filename1]")
        sys.exit()
    else:
        if sys.argv[1] == 'all':
            path = os.getcwd()
            cif_name = glob.glob(path + '/*.cif')
        else:
            cif_name = sys.argv[1]
    return cif_name


def get_symbols(input_name, metal_1, metal_2):
    # Read input CIF file
    input_structure = read(input_name)
    # Extract atomic symbols and positions
    symbols = input_structure.get_chemical_symbols()
    positions = input_structure.get_positions()

    # Count number of elemental atoms
    n_m1 = 0
    n_m2 = 0
    for symbol in symbols:
        if symbol == metal_1:
            n_m1 += 1
        elif symbol == metal_2:
            n_m2 += 1
    # Generate random permutation of Au and Pt atoms
    symbols = [metal_1] * n_m1 + [metal_2] * n_m2
    return symbols, input_structure

metal_list = ['Au', 'Pt'] # elements in binary alloy nanoparticles
name_temp = get_input()
if type(name_temp) is list:
    for name in name_temp:
        symbols, input_structure = get_symbols(name, metal_list[0], metal_list[1])
        for n,m  in enumerate(range(5)):
            random.shuffle(symbols)
            new_structure = input_structure.copy()
            new_structure.set_chemical_symbols(symbols)
            # Write new structure to CIF file
            write(name[0:name.find('.cif')]+ f'_{n}_ran.cif', new_structure)
            # write(f'output.cif', new_structure)
else:
    symbols, input_structure = get_symbols(name_temp, metal_list[0], metal_list[1])
    for n,m  in enumerate(range(5)):
        random.shuffle(symbols)
        new_structure = input_structure.copy()
        new_structure.set_chemical_symbols(symbols)
        # Write new structure to CIF file
        write(name_temp[0:name_temp.find('.cif')]+ f'_{n}_ran.cif', new_structure)
        # write(f'output.cif', new_structure)
        
os.system('mkdir random')                                                                                 
os.system('mv *ran.cif random')  
