from molSimplify.Classes.mol3D import *
import glob, os, time
import numpy as np
from molSimplify.Classes.globalvars import globalvars
from molSimplify.Classes.ligand import *
from molSimplify.Scripts.geometry import vecangle, distance
from openpyxl import load_workbook
from openpyxl import Workbook
import json

### Check whether the complex form an octahedral.
### Written by Chenru Duan
### upload: 3/7/2018
### update: 3/7/2018

dict_oct_check_loose = {'rmsd_max': 0.4, 'atom_dist_max': 0.7,
                        'num_coord_metal': 6, 'oct_angle_devi_max': 15,
                        'dist_del_eq': 0.45, 'dist_del_ax': 0.6,
                        'dist_del_eq_ax': 1.2}

dict_oct_check_st = {'rmsd_max': 0.3, 'atom_dist_max': 0.5,
                     'num_coord_metal': 6, 'oct_angle_devi_max': 10,
                     'dist_del_eq': 0.35, 'dist_del_ax': 0.4,
                     'dist_del_eq_ax': 0.8}  # default cutoff

dict_staus = {'good': 1, 'bad': 0}


## input: a xyz file
## output: a mol3D object.
def create_mol_with_xyz(_file_in):
    my_mol = mol3D()
    my_mol.readfromxyz(_file_in)
    return my_mol


## description: select the closest elements in input_arr compared
##              to the target array. Is used to screen atoms that
##              might construct an octehedral structure.
## input: two array of float. dim(input_arr) >= dim(target_arr)
## output: output_arr that is most similar to target_arr,
##         summation for the difference for output_arr and target_arr
def comp_angle_arr(input_arr, target_arr):
    output_arr = []
    del_arr = []
    for idx, ele in enumerate(input_arr):
        del_abs = []
        for _ele in target_arr:
            del_abs.append(abs(ele - _ele))
        min_del = min(del_abs)
        del_arr.append([min_del, idx])
    del_arr.sort()
    sum_del = 0
    for idx in range(len(target_arr)):
        posi = del_arr[idx][1]
        sum_del += del_arr[idx][0]
        output_arr.append(input_arr[posi])
    output_arr.sort()
    sum_del = sum_del / len(target_arr)
    return output_arr, sum_del


def sort_sec_ele(ele):
    return ele[1]


## standard 1: number of coordination for metal.
## For octehegral, if num_coord_metal < 6, it will be considered as
## problemetic. For num_coord_metal > 6, it might be caused by the
## ambiguity of C-N ring in the getBondedAtom function. We will use
## other metric (oct_comp) to see whether some of these catoms can
## construct an Oct.
## input: a xyz file
## output: coordination number for metal, and their indexes.
def get_num_coord_metal(file_in):
    my_mol = create_mol_with_xyz(_file_in=file_in)
    metal_ind = my_mol.findMetal()[0]
    metal_coord = my_mol.getAtomCoords(metal_ind)
    # print('metal index:', metal_ind)
    print('metal coordinate:', metal_coord)
    catoms = my_mol.getBondedAtoms(ind=metal_ind)
    print('coordination number:', catoms)
    ## standard 1: number of coordination for metal.
    num_coord_metal = len(catoms)
    return num_coord_metal, catoms


## standard 2: ligand changes (rmsd, max atom distance)
## for each ligand in the complex. Need input for
## original gemoetry.
## input: optimized and original xyz file.
## output: two scalar of maximum rmsd for ligands, and the
##         maximum distance change in ligands.
def ligand_comp_org(file_in, file_init_geo):
    liglist, liglist_init, flag_match = match_lig_list(file_in, file_init_geo)
    print('lig_list:', liglist)
    print('lig_list_init:', liglist_init)
    if flag_match:
        rmsd_arr, max_atom_dist_arr = [], []
        for idx, lig in enumerate(liglist):
            print('----This is %d th piece of ligand.' % (idx + 1))
            posi_shift = 2
            start_posi = posi_shift + lig[0]
            end_posi = posi_shift + lig[len(lig) - 1]
            lig_init = liglist_init[idx]
            start_posi_init = posi_shift + lig_init[0]
            end_posi_init = posi_shift + lig_init[len(lig_init) - 1]
            with open(file_in, 'r') as fo:
                lines = fo.readlines()[start_posi:end_posi + 1]
                with open('tmp.xyz', 'w') as foo:
                    foo.write(''.join(lines))
            with open(file_init_geo, 'r') as fo:
                lines = fo.readlines()[start_posi_init:end_posi_init + 1]
                with open('tmp_org.xyz', 'w') as foo:
                    foo.write(''.join(lines))
            tmp_mol = create_mol_with_xyz('tmp.xyz')
            tmp_org_mol = create_mol_with_xyz('tmp_org.xyz')
            print('# atoms: %d, init: %d' % (tmp_mol.natoms, tmp_org_mol.natoms))
            mol0, U, d0, d1 = kabsch(tmp_org_mol, tmp_mol)
            rmsd = tmp_mol.rmsd(tmp_org_mol)
            rmsd_arr.append(rmsd)
            print('rmsd:', rmsd)
            atom_dist_max = tmp_mol.maxatomdist(tmp_org_mol)
            max_atom_dist_arr.append(atom_dist_max)
            print('atom_dist_max', atom_dist_max)
        rmsd_max = max(rmsd_arr)
        atom_dist_max = max(max_atom_dist_arr)
    else:
        rmsd_max, atom_dist_max = -1, -1
    return rmsd_max, atom_dist_max


def match_lig_list(file_in, file_init_geo):
    flag_match = True
    my_mol = create_mol_with_xyz(_file_in=file_in)
    # print('natoms: ', my_mol.natoms)
    liglist, ligdents, ligcons = ligand_breakdown(my_mol)
    init_mol = create_mol_with_xyz(_file_in=file_init_geo)
    liglist_init, ligdents_init, ligcons_init = ligand_breakdown(init_mol)
    liglist_atom = [[my_mol.getAtom(x).symbol() for x in ele]
                    for ele in liglist]
    liglist_init_atom = [[init_mol.getAtom(x).symbol() for x in ele]
                         for ele in liglist_init]
    print(liglist_atom, liglist_init_atom)
    liglist_shifted = []
    ## --- match opt geo with init geo---
    # for ele in liglist_atom:
    #     try:
    #         posi = liglist_init_atom.index(ele)
    #         liglist_init_shifted.append(liglist_init[posi])
    #         liglist_init_atom.pop(posi)
    #         liglist_init.pop(posi)
    #     except ValueError:
    #         print('Ligands cannot match!')
    #         flag_match = False

    ## --- match init geo with opt geo---
    for ele in liglist_init_atom:
        try:
            posi = liglist_atom.index(ele)
            liglist_shifted.append(liglist[posi])
            liglist_atom.pop(posi)
            liglist.pop(posi)
        except ValueError:
            print('Ligands cannot match!')
            flag_match = False
    return liglist_shifted, liglist_init, flag_match


## standard 3: catoms structure compared to a perfect octahedral.
## Note: The same logic can be used to judge whether a complex is
##       Tetrahedron or panor square.
## input: a xyz file.
## output: the summation for the angle deviation for the closest
##         Oct structure. Candidates are from GetBondedAtom function
##         for the metal. A metric for the distance deviation from the
##         pefect Oct, including the diff in eq ligands, ax ligands and
##         eq-ax ligands.
def oct_comp(file_in):
    my_mol = create_mol_with_xyz(_file_in=file_in)
    num_coord_metal, catoms = get_num_coord_metal(file_in=file_in)
    # metal_ind = my_mol.findMetal()[0]
    metal_coord = my_mol.getAtomCoords(my_mol.findMetal()[0])
    catom_coord = []
    angle_ref = [90, 90, 90, 90, 180]
    theta_arr, oct_dist = [], []
    for atom in catoms:
        coord = my_mol.getAtomCoords(atom)
        catom_coord.append(coord)
    for idx1, coord1 in enumerate(catom_coord):
        delr1 = (np.array(coord1) - np.array(metal_coord)).tolist()
        theta_tmp = []
        for idx2, coord2 in enumerate(catom_coord):
            if idx2 != idx1:
                delr2 = (np.array(coord2) - np.array(metal_coord)).tolist()
                theta = vecangle(delr1, delr2)
                theta_tmp.append(theta)
        # print('For idx %d, theta array is:'%idx1, theta_tmp)
        out_theta, sum_del = comp_angle_arr(input_arr=theta_tmp, target_arr=angle_ref)
        # print('The adjusted array is: ', out_theta)
        theta_arr.append([catoms[idx1], sum_del, out_theta])

    theta_arr.sort(key=sort_sec_ele)
    theta_trunc_arr = theta_arr[0:6]
    # print('truncated theta array:', theta_trunc_arr)
    theta_trunc_arr_T = list(map(list, zip(*theta_trunc_arr)))
    oct_catoms = theta_trunc_arr_T[0]
    oct_angle_devi = theta_trunc_arr_T[1]
    oct_angle_all = theta_trunc_arr_T[2]
    print('Summation of deviation angle for catoms:', oct_angle_devi)
    print('Angle for catoms:', oct_angle_all)
    for atom in oct_catoms:
        coord = catom_coord[catoms.index(atom)]
        dist = distance(coord, metal_coord)
        oct_dist.append(dist)
    oct_dist.sort()
    if (oct_dist[2] - oct_dist[1]) > (oct_dist[4] - oct_dist[3]):
        dist_ax, dist_eq = oct_dist[:2], oct_dist[2:]  # ax dist is smaller
    else:
        dist_ax, dist_eq = oct_dist[4:], oct_dist[:4]  # eq dist is smaller
    print('dist:', dist_eq, dist_ax)
    dist_del_eq = max(dist_eq) - min(dist_eq)
    dist_del_ax = max(dist_ax) - min(dist_ax)
    dist_del_eq_ax = max(abs(max(dist_eq) - min(dist_ax)), abs(max(dist_ax) - min(dist_eq)))
    oct_dist_del = [dist_del_eq, dist_del_ax, dist_del_eq_ax]
    print('distance difference for catoms to metal (eq, max, eq_ax):', oct_dist_del)
    return oct_angle_devi, oct_dist_del


## See whether a complex is Oct or not.
## output: flag: 1 is good and 0 is bad
##         flag_list: if structure is bad, which test it fails
##         dict_oct_info: values for each metric we check.
def IsOct(file_in, file_init_geo=None, dict_check=dict_oct_check_st):
    num_coord_metal, catoms = get_num_coord_metal(file_in)

    if file_init_geo != None:
        rmsd_max, atom_dist_max = ligand_comp_org(file_in, file_init_geo)
    else:
        rmsd_max, atom_dist_max = -1, -1
    if num_coord_metal >= 6:
        oct_angle_devi, oct_dist_del = oct_comp(file_in)
    else:
        oct_angle_devi, oct_dist_del = [-1, -1], [-1, -1, -1]
    dict_oct_info = {}
    dict_oct_info['num_coord_metal'] = num_coord_metal
    dict_oct_info['rmsd_max'] = rmsd_max
    dict_oct_info['atom_dist_max'] = atom_dist_max
    dict_oct_info['oct_angle_devi_max'] = max(oct_angle_devi)
    dict_oct_info['dist_del_eq'] = oct_dist_del[0]
    dict_oct_info['dist_del_ax'] = oct_dist_del[1]
    dict_oct_info['dist_del_eq_ax'] = oct_dist_del[2]
    print('dict_oct_info', dict_oct_info)
    flag_list = []
    ## ---Adjust cutoff value---
    # if dict_oct_info['rmsd_max'] > 0.3:
    #     flag_list.append('rmsd_max')
    # if dict_oct_info['atom_dist_max'] > 0.5:
    #     flag_list.append('atom_dist_max')
    # if dict_oct_info['num_coord_metal'] < 6:
    #     flag_list.append('num_coord_metal')
    # if dict_oct_info['oct_angle_devi_max'] > 10:
    #     flag_list.append('oct_angle_devi_max')
    # if dict_oct_info['dist_del_eq'] > 0.35:
    #     flag_list.append('dist_del_eq')
    # if dict_oct_info['dist_del_ax'] > 0.4:
    #     flag_list.append('dist_del_ax')
    # if dict_oct_info['dist_del_eq_ax'] > 0.8:
    #     flag_list.append('dist_del_eq_ax')

    for key, values in dict_check.items():
        if dict_oct_info[key] > values:
            flag_list.append(key)

    if not len(flag_list):
        flag_oct = 1  # good structure
        flag_list = 'None'
    else:
        flag_oct = 0
        flag_list = ', '.join(flag_list)
    return flag_oct, flag_list, dict_oct_info


## input: _path: path for opt geo
##         path_init_geo: path for init geo
## generate and xlsx file of ['unique_num', 'oxstate', 'spinmult', 'flag_oct', 'flag_list', 'dict_oct_info']
def loop_structure(_path, path_init_geo):
    charac = ['unique_num', 'oxstate', 'spinmult', 'flag_oct', 'flag_list']
    info_tot = []
    for dirpath, dirs, files in sorted(os.walk(_path)):
        for name in sorted(files):
            if name.split('.')[1] == 'xyz':
                unique_num, oxstate, spinmult = name.split('_')[0], name.split('_')[4][1:], name.split('_')[5][2:]
                print(unique_num, oxstate, spinmult)
                file_in = '%s/%s' % (dirpath, name)
                file_init_geo = find_file_with_unique_num(path_init_geo, unique_num)
                # file_init_geo = gen_file_with_name(path_init_geo, name)
                print('!!!!file info:!!!!!', file_in, file_init_geo)
                if os.path.exists(file_init_geo):
                    flag_oct, flag_list, dict_oct_info = IsOct(file_in, file_init_geo)
                else:
                    print('No init_geo!', name)
                    flag_oct, flag_list, dict_oct_info = IsOct(file_in)
                # dict_oct_info = json.dumps(dict_oct_info) # dictionary to string
                _c, _dict = [], []
                for key, value in dict_oct_info.items():
                    _c.append(key)
                    _dict.append(value)
                info = [unique_num, oxstate, spinmult, flag_oct, flag_list] + _dict
                info_tot.append(info)
    charac += _c  # append dict_oct_info
    write_list_2_xlsx(charac, info_tot)
    return 0

## find the file name by matching unique_ID and spin mulplicity
def find_file_with_unique_num(_path, unique_num):
    filename = None
    for dirpath, dirs, files in sorted(os.walk(_path)):
        for name in sorted(files):
            if name.split('_')[0] == str(unique_num):
                filename = '%s/%s' % (dirpath, name)
                break
    return filename


def gen_file_with_name(path_init_geo, name_opt):
    name_opt = name_opt.split('_')
    name_opt = '_'.join(name_opt[:len(name_opt) - 1])
    name_init = '%s/%s_mols.xyz' % (path_init_geo, name_opt)
    return name_init


def write_list_2_xlsx(charac, txt_list, filename_save='new.xlsx', title='new'):
    wb_write = Workbook()
    ws_write = wb_write.active
    ws_write.title = title

    ws_write.append(charac)
    for info in txt_list:
        ws_write.append(info)
    wb_write.save(filename=filename_save)


def final_status(filename, sheetname):
    charac, txt_list = read_xlsx_2_list(filename, sheetname)
    for idx, ele in enumerate(txt_list):
        flag_oct = ele[charac.index('flag_oct')]
        run_status = ele[charac.index('run_status')]
        if flag_oct == '1' and run_status == 'Y':
            _status = 1
        elif flag_oct == '0':
            _status = 0
        else:
            _status = 'DN'
        txt_list[idx].append(_status)
    charac.append('_status')
    write_list_2_xlsx(charac, txt_list, filename_save='new_status.xlsx', title='new')
    return 0


##########################
## ---batch processing----
_path = './fe_optimized_geometry'
path_init_geo = './fe_init_geo'
# loop_structure(_path, path_init_geo)

## --For single complex-----
# file_in = '1694_co_eq_40_ax1_47_ax2_49_2_2.xyz'
# file_org = '1694_co_eq_40_ax1_47_ax2_49_2_2_origin.xyz'
file_in = '10016_fe_smi16_ligcharge-4_c2_sm5_optgeo.xyz'
file_org = '10016_fe_smi16_ligcharge-4_c2_sm5_mols.xyz'
flag_oct, flag_list, dict_oct_info = IsOct(file_in, file_org)
print('=====molecule: %s, flag_oct: %d' % (file_in, flag_oct))
if not flag_oct:
    print('=====flag_list:', flag_list)
    print('=====dict_oct_info:', dict_oct_info)
