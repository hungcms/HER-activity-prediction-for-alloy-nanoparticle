import os
import sys
import math
from pymatgen.core.structure import Structure
import matplotlib.pyplot as plt

def get_inputs():
    if not (len(sys.argv) == 10 or len(sys.argv) == 8):
        print("Wrong input")
        print("usage: command_name cfg_file_ss cfg_file_core_shell_1 cfg_file_name_core_shell_2 \
              cfg_file_name_random cfg_file_name_all first_element second_element y_min y_max")
        sys.exit()
    if len(sys.argv) == 10:
        name_1 = sys.argv[1]
        name_2 = sys.argv[2]
        name_3 = sys.argv[3]
        name_4 = sys.argv[4]
        name_5 = sys.argv[5]
        element_1 = sys.argv[6]
        element_2 = sys.argv[7]
        y_min = float(sys.argv[8])
        y_max = float(sys.argv[9])
    if len(sys.argv) == 8:
        name_1 = sys.argv[1]
        name_2 = sys.argv[2]
        name_3 = sys.argv[3]
        name_4 = sys.argv[4]
        name_5 = sys.argv[5]
        element_1 = sys.argv[6]
        element_2 = sys.argv[7]
        y_min = -0.5
        y_max = 0.5
    return name_1, name_2, name_3, name_4, name_5, element_1, element_2, y_min, y_max

def get_file_string_array(file_name):
    try:
        file = open(file_name, "r")
    except IOError:
        print('Error: file (%s) not found!\n' % (file_name))
        sys.exit()
    lines = file.readlines()
    file.close()
    array = []
    for line in lines:
        array.append(line.split())
    return array

def selection_sort_energy(arr):
    for i in range(len(arr)):
        minimum = i

        for j in range(i + 1, len(arr)):
            if arr[j][1] < arr[minimum][1]:
                minimum = j

        arr[minimum], arr[i] = arr[i], arr[minimum]

    return arr

def selection_sort_ratio(arr):
    for i in range(len(arr)):
        minimum = i

        for j in range(i + 1, len(arr)):
            if arr[j][0] < arr[minimum][0]:
                minimum = j

        arr[minimum], arr[i] = arr[i], arr[minimum]

    return arr

def for_energy(ratio, Energy, size, super_cell, ratio_raw):
    E_pure_atom_1 = []
    E_pure_atom_2 = []
    for i, k, l in zip(ratio, Energy, size):
        if i == 1:
            E_pure_atom_1.append(k/l)
        if i == 0:
            E_pure_atom_2.append(k/l)
    atom_1_value = sum(E_pure_atom_1) / len(E_pure_atom_1)
    atom_2_value = sum(E_pure_atom_2) / len(E_pure_atom_2)
    
    formation = []
    for i, k, l in zip(size, Energy, ratio_raw):
        formation.append((float(k) - l[0] *
                        atom_1_value - l[1]*atom_2_value)/float(i))

    return formation
   
   
def get_info(file_name):
    array = get_file_string_array(file_name)
    
    Energy, size, index, struc, conf, begin_line, end_line = [], [], [], [], [], [], []
    
    for i in range(len(array)):
        if len(array[i]) != 0:
            if array[i][0] == 'Size':
                size.append(int(array[i+1][0]))
            if array[i][0] == "AtomData:":
                index.append(i)
            if array[i][0] == "Energy":
                Energy.append(float(array[i+1][0]))
            if array[i][0] == "Feature" and array[i][1] == "Structure":
                struc.append([array[i][2]])
            if array[i][0] == "BEGIN_CFG":
                begin_line.append(i)
            if array[i][0] == "END_CFG":
                end_line.append(i)
        else:
            pass
    new_size = list(dict.fromkeys(size))
    super_cell = [size[i]/min(new_size) for i in range(len(size))]
    
    for i, k in zip(begin_line, end_line):
        conf.append([i, k])
        
    ratio = []
    ratio_raw = []
    for g, v in zip(size, index):
        atom_1 = 0
        atom_2 = 0
        for k in range(v+1, g+v+1):
            if array[k][1] == "0":
                atom_1 += 1
            else:
                atom_2 += 1

        if atom_1 == 0 and atom_2 == g:
            ratio.append(0)
            raw = [0, atom_2]
            ratio_raw.append(raw)
        elif atom_1 == g and atom_2 == 0:
            ratio.append(1)
            raw = [atom_1, 0]
            ratio_raw.append(raw)
        else:
            raw = [atom_1, atom_2]
            ratio_raw.append(raw)
            ratio_bimetallic = float(atom_1/(atom_1 + atom_2))
            ratio.append(ratio_bimetallic)
            
    info_conf = []
    for_ene = for_energy(ratio, Energy, size, super_cell, ratio_raw)
    for a, b, c, d, e in zip(ratio, for_ene, size, struc, conf):
        info_conf.append([a, b, c, d, e])

    return info_conf, ratio

def min_ene(file_name):
    info_conf, ratio = get_info(file_name)
    ratio_new = list(dict.fromkeys(ratio))
    list_info = [[] for _ in range(len(ratio_new))]
    
    for i in range(len(ratio_new)):
        for k in range(len(info_conf)):
            if info_conf[k][0] == ratio_new[i]:
                list_info[i].append(info_conf[k])
        selection_sort_energy(list_info[i])

    final_info = []

    for i in range(len(list_info)):
        final_info.append(list_info[i][0])

    return final_info


def det(p1, p2, p3):
	return (p2[0]-p1[0])*(p3[1]-p1[1]) \
            - (p2[1]-p1[1])*(p3[0]-p1[0])

def hull_point(file_name):
    connect_sorta_energy = selection_sort_energy(min_ene(file_name))
    lowest = connect_sorta_energy[0]
    connect = selection_sort_ratio(min_ene(file_name))
    for i in range(len(connect)):
        degree = math.degrees(math.atan2(connect[i][1]-lowest[1], connect[i][0]-lowest[0]))
        connect[i].append(degree)
    
    for i in range(len(connect)):
        min = i
        for j in range(i + 1, len(connect)):
            if connect[j][5] < connect[min][5]:
                min = j
        connect[min], connect[i] = connect[i], connect[min]
    
    hull = [lowest, connect[1]]

    for i in connect[2:]:
        while det(hull[-2], hull[-1], i) <= 0:
            del hull[-1]
        hull.append(i)

    hull_point = []
    for i in range(len(hull)):
        if hull[i][1] <= 0:
            hull_point.append(hull[i])
    
    return hull_point

def write_cif(file_name, atom_1, atom_2):
    array = get_file_string_array(file_name)
    final_info = hull_point(file_name)
    index_atom = {0: f"{atom_1}", 1: f"{atom_2}"}
    for i in range(len(final_info)):
        lattice, coords, species = [], [], []
        start, end = final_info[i][4][0], final_info[i][4][1]
        for v, j in enumerate(array[start:end]):
            if v == 4 or v == 5 or v == 6:
                lattice.append([float(j[0]), float(j[1]), float(j[2])])
            if 8 <= v < (8 + int(final_info[i][2])):
                coords.append([float(j[2]), float(j[3]), float(j[4])])
                if j[1] == "0":
                    species.append(index_atom[0])
                else:
                    species.append(index_atom[1])
        from pymatgen.io.cif import CifWriter
        structure = Structure(lattice, species, coords, coords_are_cartesian= True)
        writer = CifWriter(structure)
        writer.write_file(f"{final_info[i][3][0]}.cif")

def write_txt(file_name):
    info = hull_point(file_name)
    with open("ene.txt", "wt") as f:
        f.write("Configuration".ljust(30) + "Ratio".ljust(20) + "Number of Atom".ljust(20) + "Formation energy".ljust(20) + "\n")
        for i in range(len(info)):
            f.write(str(info[i][3][0]).ljust(30) + str(float("{:.3f}".format(info[i][0]))).ljust(20) + str(info[i][2]).ljust(20) + str(float("{:.3f}".format(info[i][1]))).ljust(20) + "\n")


def plot_convexhull(filename_1, filename_2, filename_3, filename_4, filename_5):
    info_1, ratio_1 = get_info(filename_1)
    info_2, ratio_2 = get_info(filename_2)
    info_3, ratio_3 = get_info(filename_3)
    info_4, ratio_4 = get_info(filename_4)
    x1, y1, x2, y2, x3, y3, x4, y4 = [], [], [], [], [], [], [], []
    for i in range(len(info_1)):
        x1.append(info_1[i][0])
        y1.append(info_1[i][1])
    for i in range(len(info_2)):
        x2.append(info_2[i][0])
        y2.append(info_2[i][1])
    for i in range(len(info_3)):
        x3.append(info_3[i][0])
        y3.append(info_3[i][1])
    for i in range(len(info_4)):
        x4.append(info_4[i][0])
        y4.append(info_4[i][1])
    hull_point_fn = selection_sort_ratio(hull_point(filename_5))
    a ,b = [], []
    for i in range(len(hull_point_fn)):
        a.append(hull_point_fn[i][0])
        b.append(hull_point_fn[i][1])


    plt.figure(figsize=(5.3,3.8))    
    
    plt.scatter(x1, y1, color = "green", s=20, marker = "s", label='Order')
    plt.scatter(x2, y2, color = "purple", s=20, marker = "^", label=f'{atom_1}@{atom_2}')
    plt.scatter(x3, y3, color = "orange",s=20, marker = "h", label=f'{atom_2}@{atom_1}')
    plt.scatter(x4, y4, color = "blue", s=20, marker = "o", label = 'Disorder')
    
    plt.hlines(y=0, xmin=0, xmax=1, linestyle="dashed", color='black')
    plt.plot(a, b, marker= "o", markersize=12, mfc="none", color='r')
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [0, 3, 1, 2]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=10)
    plt.xlabel(atom_1 + ' ' + 'concentration', fontsize=12)
    plt.ylabel("Formation Energy (eV/atom)",fontsize=12)
    plt.xlim(0, 1)
    plt.ylim(y_min, y_max)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    plt.savefig(f"convex_hull_{atom_1}_{atom_2}.jpg", dpi=600)
    plt.show()

file_name_1, file_name_2, file_name_3, file_name_4, file_name_5, atom_1, atom_2, y_min, y_max = get_inputs()

write_txt(file_name_5)
write_cif(file_name_5, atom_1, atom_2)
plot_convexhull(file_name_1, file_name_2,  file_name_3, file_name_4, file_name_5)
