import numpy as np
import sys, os, glob
from pymatgen.core.structure import IStructure, Structure
from pymatgen.core.surface import SlabGenerator, miller_index_from_sites
from pymatgen.io.cif import CifWriter


def gen_slab(filename, miller, min_slab=10, min_vacuum=10):
    structure = Structure.from_file(filename)
    slab_generator = SlabGenerator(initial_structure=structure, miller_index=miller, min_slab_size=min_slab, min_vacuum_size=min_vacuum, 
                                   center_slab=True, reorient_lattice=True, primitive=True, in_unit_planes=False, lll_reduce=True)
    slabs =  slab_generator.get_slabs()
    return slabs

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

miller_list = [(1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1), (1, 1, 1)]
name_temp = get_input()
if type(name_temp) is list:
    for i in name_temp:
        name = i
        for miller in miller_list:
            slabs = gen_slab(name, miller=miller)
            for k, slab in enumerate(slabs):
                CifWriter(struct=slab).write_file(name[0:name.find('.cif')]+ f'_{k}_{miller[0]}{miller[1]}{miller[2]}_slab.cif') 
else:
    for miller in miller_list:
        slabs = gen_slab(name_temp, miller=miller)
        for slab in slabs:
            CifWriter(struct=slab).write_file(name_temp[0:name_temp.find('.cif')]+'_slab.cif')
