import os, sys
from tracemalloc import start
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element

if len(sys.argv) == 1:
    print("""
[1]: cif files to cfg
[2]: vasprun to cfg
[3]: cfg to cifs
[4]: cfg to cifs -n=select number,....
""")
    quit()

elt_index_map = {"Pt": 0, "Ni": 1, 'H': 2}
index_elt_map = {0: "Pt", 1: "Ni", 2: 'H'}

def pmg_structure_to_cfg_block_1(structure):
    sites = structure.sites
    lat = structure.lattice.matrix
    string_block = ""
    string_block +=' Size\n'
    string_block +='    %d\n' % len(sites)
    string_block +=' Supercell\n'
    string_block +='    %13.6f %13.6f %13.6f\n' % (lat[0][0], lat[0][1], lat[0][2])
    string_block +='    %13.6f %13.6f %13.6f\n' % (lat[1][0], lat[1][1], lat[1][2])
    string_block +='    %13.6f %13.6f %13.6f\n' % (lat[2][0], lat[2][1], lat[2][2])
    string_block +=' AtomData:  id type       cartes_x      cartes_y      cartes_z\n'
    for i, s in enumerate(sites):
        coords = s.coords
        specie_type = elt_index_map[str(s.specie)]
        string_block +='           %3d  %3d  %13.6f %13.6f %13.6f\n' % (i+1, specie_type, coords[0], coords[1], coords[2])
        
    return string_block

def pmg_structure_to_cfg_block_2(structure, forces):
    sites = structure.sites
    lat = structure.lattice.matrix
    string_block = ""
    string_block +=' Size\n'
    string_block +='    %d\n' % len(sites)
    string_block +=' Supercell\n'
    string_block +='    %13.6f %13.6f %13.6f\n' % (lat[0][0], lat[0][1], lat[0][2])
    string_block +='    %13.6f %13.6f %13.6f\n' % (lat[1][0], lat[1][1], lat[1][2])
    string_block +='    %13.6f %13.6f %13.6f\n' % (lat[2][0], lat[2][1], lat[2][2])
    string_block +=' AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n'
    for i, s in enumerate(sites):
        coords = s.coords
        specie_type = elt_index_map[str(s.specie)]
        string_block +='           %3d  %3d  %13.6f %13.6f %13.6f  %11.6f %11.6f %11.6f\n' % (i+1, specie_type, coords[0], coords[1], coords[2], forces[i][0], forces[i][1], forces[i][2])
        
    return string_block

cfg_filename = "gen.cfg"
if sys.argv[1] == "1":
    cifs = [f for f in sys.argv[2:] if '.cif' in f]
    cifs.sort()
    print(cifs)
    f = open(cfg_filename, 'w')
    for cif in cifs:
        name=cif.replace("_poscar.cif","")
        f.write('BEGIN_CFG\n')
        structure = Structure.from_file(cif)
        structure_block_string = pmg_structure_to_cfg_block_1(structure)
        f.write(structure_block_string)
        f.write('Feature\tStructure %s\n' % name)
        f.write('END_CFG\n\n')
    f.close()
    
elif sys.argv[1] == "2":
    from pymatgen.io.vasp.outputs import Vasprun
    vaspruns = [f for f in sys.argv[2:] if 'vasprun.xml' in f]
    print(vaspruns)
    f = open(cfg_filename, 'w')
    for filepath in vaspruns:
        run = Vasprun(filepath, parse_dos=False, parse_eigen=False, parse_potcar_file=False)
        structures = run.structures
        energies = []
        forces = []
        for istep in run.ionic_steps:
            if istep["e_wo_entrp"] != istep["electronic_steps"][-1]["e_0_energy"]:
                print("CAUTION TO PARSING VASPRUN!!!!")
                energies.append(istep["e_wo_entrp"])
            else:
                energies.append(istep["electronic_steps"][-1]["e_0_energy"])
                forces.append(istep["forces"])
        for i in range(len(structures)):
            f.write('BEGIN_CFG\n')
            structure_block_string = pmg_structure_to_cfg_block_2(structures[i], forces[i])
            f.write(structure_block_string)
            f.write(' Energy\n')
            f.write(' %23.12f\n' % energies[i])
            f.write('END_CFG\n\n')
        
elif sys.argv[1] == "3":
    cfg_filename = [f for f in sys.argv[2:] if '.cfg' in f][0]
    print(cfg_filename)

    file_lines = open(cfg_filename, 'r').readlines()
    beginblocks = []
    endblocks = []
    featureblocks = []
    for i, line in enumerate(file_lines):
        if 'BEGIN_CFG' in line:
            beginblocks.append(i)
        elif 'END_CFG' in line:
            endblocks.append(i)
        elif 'Feature   Structure' in line:
            featureblocks.append(i)

    for i in range(len(beginblocks)):
        lattice = []
        species = []
        coords = []
        species_2 = []
        coords_2 = []
        for j, line in enumerate(file_lines[beginblocks[i]:endblocks[i]]):
            if j == 4:
                line = line.replace('\n', '') 
                spl = line.split()
                v_x = [float(spl[0]), float(spl[1]), float(spl[2])]
                lattice.append(v_x)
            if j == 5:
                line = line.replace('\n', '') 
                spl = line.split()
                v_y = [float(spl[0]), float(spl[1]), float(spl[2])]
                lattice.append(v_y)
            if j == 6:
                line = line.replace('\n', '') 
                spl = line.split()
                v_z = [float(spl[0]), float(spl[1]), float(spl[2])]
                lattice.append(v_z)
            if j >= 8 and len(line.split()) == 5:
                line = line.replace('\n', '') 
                spl = line.split()
                elt = index_elt_map[int(spl[1])]
                species.append(str(elt))
                x, y, z = float(spl[2]), float(spl[3]), float(spl[4])
                coords.append([x, y, z]) 
            if j >= 8 and len(line.split()) == 8:
                line = line.replace('\n', '') 
                spl = line.split()
                elt = index_elt_map[int(spl[1])]
                species_2.append(str(elt))
                x, y, z = float(spl[2]), float(spl[3]), float(spl[4])
                coords_2.append([x, y, z]) 
            
        for n, k in enumerate(featureblocks):
            if n == i:
                Feature = file_lines[k]
                Feature = Feature.split()[-1].replace('\n','')

        structure = Structure(lattice, species, coords, coords_are_cartesian=True)
        cfg_filename = cfg_filename.split("/")[-1]
        filename = "%s_%02d.cif" % (Feature, i+1)
        print(filename)
        try:
            structure.to('cif', filename)
        except:
            structure = Structure(lattice, species_2, coords_2, coords_are_cartesian=True)
            structure.to('cif', filename)
        
elif sys.argv[1] == "4":
    cfg_filename = [f for f in sys.argv[2:] if '.cfg' in f][0]
    print(cfg_filename)
    for arg in sys.argv[2:]:
        if '-n=' in arg:
            start_list = []
            startnumber = arg.split('=')[-1]
            if ',' in startnumber:
                startnumber = startnumber.split(',')
                for sn in startnumber:
                    start_list.append(int(sn))
            elif '-' in startnumber:
                startnumber = startnumber.split('-')
                for i in range(int(startnumber[0]), int(startnumber[1])+1):
                    start_list.append(i)
        else:
            pass

    file_lines = open(cfg_filename, 'r').readlines()
    beginblocks = []
    endblocks = []
    featureblocks = []
    for i, line in enumerate(file_lines):
        if 'BEGIN_CFG' in line:
            beginblocks.append(i)
        elif 'END_CFG' in line:
            endblocks.append(i)
        elif 'Feature   Structure' in line:
            featureblocks.append(i)
    
    for i in range(len(beginblocks)):
        if i in start_list:
            lattice = []
            species = []
            coords = []
            species_2 = []
            coords_2 = []
            for j, line in enumerate(file_lines[beginblocks[i]:endblocks[i]]):
                if j == 4:
                    line = line.replace('\n', '') 
                    spl = line.split()
                    v_x = [float(spl[0]), float(spl[1]), float(spl[2])]
                    lattice.append(v_x)
                if j == 5:
                    line = line.replace('\n', '') 
                    spl = line.split()
                    v_y = [float(spl[0]), float(spl[1]), float(spl[2])]
                    lattice.append(v_y)
                if j == 6:
                    line = line.replace('\n', '') 
                    spl = line.split()
                    v_z = [float(spl[0]), float(spl[1]), float(spl[2])]
                    lattice.append(v_z)
                if j >= 8 and len(line.split()) == 5:
                    line = line.replace('\n', '') 
                    spl = line.split()
                    elt = index_elt_map[int(spl[1])]
                    species.append(str(elt))
                    x, y, z = float(spl[2]), float(spl[3]), float(spl[4])
                    coords.append([x, y, z]) 
                if j >= 8 and len(line.split()) == 8:
                    line = line.replace('\n', '') 
                    spl = line.split()
                    elt = index_elt_map[int(spl[1])]
                    species_2.append(str(elt))
                    x, y, z = float(spl[2]), float(spl[3]), float(spl[4])
                    coords_2.append([x, y, z]) 
             
            for n, k in enumerate(featureblocks):
                if n == i:
                    Feature = file_lines[k]
                    Feature = Feature.split()[-1].replace('\n','')

            structure = Structure(lattice, species, coords, coords_are_cartesian=True)
            cfg_filename = cfg_filename.split("/")[-1]
            filename = "%s_%02d.cif" % (Feature, i+1)
            print(filename)
            try:
                structure.to('cif', filename)
            except:
                structure = Structure(lattice, species_2, coords_2, coords_are_cartesian=True)
                structure.to('cif', filename)
