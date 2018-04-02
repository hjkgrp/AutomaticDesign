from molSimplify.Classes.mol3D import *
import glob, os, time
import numpy as np
from molSimplify.Classes.globalvars import globalvars
from molSimplify.Classes.ligand import *
from molSimplify.Scripts.geometry import vecangle, distance

### Check whether the complex form an octahedral.
### Written by Chenru Duan
### upload: 3/7/2018
### update: 3/19/2018

# from openpyxl import load_workbook
# from openpyxl import Workbook
# import json

### Check whether the complex form an octahedral.
### Written by Chenru Duan
### upload: 3/7/2018
### update: 3/16/2018

dict_oct_check_loose = {'rmsd_max': 0.4, 'atom_dist_max': 0.6,
                        'num_coord_metal': 6, 'oct_angle_devi_max': 15,
                        'dist_del_eq': 0.45, 'max_del_sig_angle': 27,
                        'dist_del_all':1.2}

dict_oct_check_st = {'rmsd_max': 0.3, 'atom_dist_max': 0.45,
                     'num_coord_metal': 6, 'oct_angle_devi_max': 12,
                     'dist_del_eq': 0.35, 'dist_del_all': 1,
                     'max_del_sig_angle': 22.5}  # default cutoff

dict_staus = {'good': 1, 'bad': 0}

oct_angle_ref = [[90, 90, 90, 90, 180] for x in range(6)]


## input: a xyz file
## output: a mol3D object.
def create_mol_with_xyz(_file_in):
    my_mol = mol3D()
    my_mol.readfromxyz(_file_in)
    return my_mol


### Deprecated
# def comp_angle_arr(input_arr, target_arr):
#     output_arr = []
#     del_arr = []
#     print('!!!!!!!!', input_arr)
#     for idx, ele in enumerate(input_arr):
#         del_abs = []
#         for _ele in target_arr:
#             del_abs.append(abs(ele - _ele))
#         min_del = min(del_abs)
#         del_arr.append([min_del, idx])
#     del_arr.sort()
#     sum_del = 0
#     for idx in range(len(target_arr)):
#         posi = del_arr[idx][1]
#         sum_del += del_arr[idx][0]
#         output_arr.append(input_arr[posi])
#     output_arr.sort()
#     sum_del = sum_del / len(target_arr)
#     print('++++++', output_arr)
#     return output_arr, sum_del

## description: select the closest elements in input_arr compared
##              to the target array. Is used to screen atoms that
##              might construct an octehedral structure. Used to select
##              m targeting catoms out of n candidates (n > m).
## input: two array of float. dim(input_arr) >= dim(target_arr)
## output: output_arr that is most similar to target_arr after filtering,
##         summation for the difference for output_arr and target_arr
def comp_two_angle_array(input_angle, target_angle):
    _i = input_angle[1][:]
    # print('target_angle', _t)
    del_act = []
    output_angle = []
    for idx, ele in enumerate(target_angle):
        del_arr = []
        for _idx, _ele in enumerate(_i):
            del_arr.append([abs(ele - _ele), _idx, _ele])
        del_arr.sort()
        posi = del_arr[0][1]
        _i.pop(posi)
        # print('!!!input_a:', input_angle[1])
        del_act.append(del_arr[0][0])
        output_angle.append(del_arr[0][2])
    max_del_angle = max(del_act)
    sum_del = sum(del_act) / len(target_angle)
    #### older version
    # del_arr = []
    # for idx, ele in enumerate(input_angle[1]):
    #     del_abs = []
    #     # print('ele:', ele)
    #     _t =
    #     for _ele in target_angle:
    #         del_abs.append(abs(ele - _ele))
    #         # catoms.append(ele[0])
    #     min_del = min(del_abs)
    #     del_arr.append([min_del, idx])
    # del_arr.sort()
    # sum_del = 0
    # output_angle = []
    # max_del_angle = 0
    # for idx in range(len(target_angle)):
    #     posi = del_arr[idx][1]
    #     # print('posi:', posi)
    #     sum_del += del_arr[idx][0]
    #     if del_arr[idx][0] > max_del_angle:
    #         max_del_angle = del_arr[idx][0]
    #     output_angle.append(input_angle[1][posi])
    # output_angle.sort()
    # sum_del = sum_del / len(target_angle)
    return output_angle, sum_del, max_del_angle


## description: Given the target_angle, choose the input_angle that has
##              the smallest angle deviation in input_array. In this process,
##              the input_angle has already been filtered.
## input:       input_arr: array of input_angle, target_angle: array of angle you want.
## output:      output_angle: the input_angle that has the smallest deviation compared
##              to target_arr. del_angle: the deviation. catoms: the connecting antom for
##              the output_angle.
def comp_angle_pick_one_best(input_arr, target_angle):
    del_arr = []
    for ii, input_angle in enumerate(input_arr):
        out_angle, sum_del, max_del_angle = comp_two_angle_array(input_angle, target_angle)
        del_arr.append([sum_del, ii, max_del_angle, out_angle])
    del_arr.sort()
    posi = del_arr[0][1]
    del_angle = del_arr[0][0]
    output_angle = input_arr[posi][1]
    catoms = input_arr[posi][0]
    max_del_sig_angle = del_arr[0][2]
    # print('!!!!posi', input_arr[posi])
    # print('!!!!out:', del_angle, del_arr[0][3])
    input_arr.pop(posi)
    return output_angle, del_angle, catoms, max_del_sig_angle


## description: Loop the target_angle in target_arr.
## output: three lists corresponding to the outputs of comp_angle_pick_one_best.
def loop_target_angle_arr(input_arr, target_arr):
    output_arr = []
    sum_del = []
    catoms_arr = []
    max_del_sig_angle_arr = []
    for idx, ele in enumerate(target_arr):
        output_angle, del_angle, catoms, max_del_sig_angle = comp_angle_pick_one_best(input_arr, ele)
        output_arr.append(output_angle)
        sum_del.append(del_angle)
        catoms_arr.append(catoms)
        max_del_sig_angle_arr.append(max_del_sig_angle)
    return output_arr, sum_del, catoms_arr, max(max_del_sig_angle_arr)


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
    catoms = my_mol.getBondedAtomsOct(ind=metal_ind)
    print('coordinations: ', catoms, len(catoms))
    ## standard 1: number of coordination for metal.
    num_coord_metal = len(catoms)
    return num_coord_metal, catoms


## standard 2: ligand changes (rmsd, max atom distance)
## for each ligand in the complex. Need input for
## original gemoetry.
## input: optimized and original xyz file.
## output: two scalar of maximum rmsd for ligands, and the
##         maximum distance change in ligands.
def ligand_comp_org(file_in, file_init_geo, flag_deleteH=True, flag_loose=False,
                    flag_lbd=True):
    liglist, liglist_init, flag_match = match_lig_list(file_in, file_init_geo,
                                                       flag_loose, flag_lbd)
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
                    foo.write('number of atoms not necessary\n')
                    foo.write('tmp file for ligands comparison\n')
                    foo.write(''.join(lines))
            with open(file_init_geo, 'r') as fo:
                lines = fo.readlines()[start_posi_init:end_posi_init + 1]
                with open('tmp_org.xyz', 'w') as foo:
                    foo.write('number of atoms not necessary\n')
                    foo.write('tmp file for ligands comparison\n')
                    foo.write(''.join(lines))
            tmp_mol = create_mol_with_xyz('tmp.xyz')
            tmp_org_mol = create_mol_with_xyz('tmp_org.xyz')
            print('# atoms: %d, init: %d' % (tmp_mol.natoms, tmp_org_mol.natoms))
            if flag_deleteH:
                tmp_mol.deleteHs()
                tmp_org_mol.deleteHs()
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
        rmsd_max, atom_dist_max = 'lig_mismatch', 'lig_mismatch'
    return rmsd_max, atom_dist_max


## Match the ligend list generated by ligand_breakdown.
## useful for cases where the init geo and opt geo have the
## different ligands arrangement (when you only have the opt geo
## but still want init geo)
def match_lig_list(file_in, file_init_geo, flag_loose, flag_lbd=True):
    flag_match = True
    my_mol = create_mol_with_xyz(_file_in=file_in)
    # print('natoms: ', my_mol.natoms)
    print('In match, flag_loose', flag_loose)
    init_mol = create_mol_with_xyz(_file_in=file_init_geo)
    liglist_init, ligdents_init, ligcons_init = ligand_breakdown(init_mol)
    if flag_lbd: ## Also do ligand breakdown for opt geo
        liglist, ligdents, ligcons = ligand_breakdown(my_mol, flag_loose)
    else: ## ceate/use the liglist, ligdents, ligcons of initial geo as we just wanna track them down
        liglist, ligdents, ligcons = liglist_init[:], ligdents_init[:], ligcons_init[:]
    liglist_atom = [[my_mol.getAtom(x).symbol() for x in ele]
                    for ele in liglist]
    liglist_init_atom = [[init_mol.getAtom(x).symbol() for x in ele]
                         for ele in liglist_init]
    print('ligand_list opt in symbols:', liglist_atom)
    print('ligand_list init in symbols: ', liglist_init_atom)
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
def oct_comp(file_in, angle_ref=oct_angle_ref, catoms_arr=None):
    my_mol = create_mol_with_xyz(_file_in=file_in)
    num_coord_metal, catoms = get_num_coord_metal(file_in=file_in)
    # metal_ind = my_mol.findMetal()[0]
    metal_coord = my_mol.getAtomCoords(my_mol.findMetal()[0])
    catom_coord = []
    if not catoms_arr == None:
        catoms = catoms_arr
        num_coord_metal = len(catoms_arr)
    theta_arr, oct_dist = [], []
    for atom in catoms:
        coord = my_mol.getAtomCoords(atom)
        catom_coord.append(coord)
    th_input_arr = []
    for idx1, coord1 in enumerate(catom_coord):
        delr1 = (np.array(coord1) - np.array(metal_coord)).tolist()
        theta_tmp = []
        for idx2, coord2 in enumerate(catom_coord):
            if idx2 != idx1:
                delr2 = (np.array(coord2) - np.array(metal_coord)).tolist()
                theta = vecangle(delr1, delr2)
                theta_tmp.append(theta)
        th_input_arr.append([catoms[idx1], theta_tmp])
        # print('For idx %d, theta array is:'%idx1, theta_tmp)
        # out_theta, sum_del = comp_angle_arr(input_arr=theta_tmp, target_arr=angle_ref)
        # print('The adjusted array is: ', out_theta)
        # theta_arr.append([catoms[idx1], sum_del, out_theta])
    th_output_arr, sum_del_angle, catoms_arr, max_del_sig_angle = loop_target_angle_arr(th_input_arr, angle_ref)
    print('th:', th_output_arr)
    print('sum_del:', sum_del_angle)
    print('catoms_arr:', catoms_arr)
    print('catoms_type:', [my_mol.getAtom(x).symbol() for x in catoms_arr])
    for idx, ele in enumerate(th_output_arr):
        theta_arr.append([catoms_arr[idx], sum_del_angle[idx], ele])
    # theta_arr.sort(key=sort_sec_ele)
    # theta_trunc_arr = theta_arr[0:6]
    theta_trunc_arr = theta_arr
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
    # print('!!!oct_dist:', oct_dist)
    try:  ### For Oct
        dist_del_arr = np.array([oct_dist[3] - oct_dist[0], oct_dist[4] - oct_dist[1], oct_dist[5] - oct_dist[2]])
        min_posi = np.argmin(dist_del_arr)
        if min_posi == 0:
            dist_eq, dist_ax = oct_dist[:4], oct_dist[4:]
        elif min_posi == 1:
            dist_eq, dist_ax = oct_dist[1:5], [oct_dist[0], oct_dist[5]]
        else:
            dist_eq, dist_ax = oct_dist[2:], oct_dist[:2]
    except IndexError:  ## For one empty site
        if (oct_dist[3] - oct_dist[0]) > (oct_dist[4] - oct_dist[1]):
            dist_ax, dist_eq = oct_dist[:1], oct_dist[1:]  # ax dist is smaller
        else:
            dist_ax, dist_eq = oct_dist[4:], oct_dist[:4]  # eq dist is smaller
    dist_del_all = oct_dist[-1] - oct_dist[0]
    print('dist:', dist_eq, dist_ax)
    dist_del_eq = max(dist_eq) - min(dist_eq)
    dist_del_ax = max(dist_ax) - min(dist_ax)
    dist_del_eq_ax = max(abs(max(dist_eq) - min(dist_ax)), abs(max(dist_ax) - min(dist_eq)))
    oct_dist_del = [dist_del_eq, dist_del_ax, dist_del_eq_ax, dist_del_all]
    print('distance difference for catoms to metal (eq, ax, eq_ax):', oct_dist_del)
    return oct_angle_devi, oct_dist_del, max_del_sig_angle, catoms_arr


def Oct_inspection(file_in, file_init_geo=None, catoms_arr=None, dict_check=dict_oct_check_st,
                   std_not_use=[], angle_ref=oct_angle_ref, flag_loose=True, flag_lbd=False):
    if catoms_arr == None:
        print('Error, must have ctoms! If not, please use IsOct.')
        quit()
    elif len(catoms_arr) != 6:
        print('Error, must have 6 connecting atoms for octahedral.')
        quit()
    num_coord_metal = 6
    oct_angle_devi, oct_dist_del, max_del_sig_angle = [-1, -1], [-1, -1, -1, -1], -1
    rmsd_max, atom_dist_max = -1, -1
    if not file_init_geo == None:
        print('!!!Inspection,flag_loose:', flag_loose)
        rmsd_max, atom_dist_max = ligand_comp_org(file_in, file_init_geo, flag_loose=flag_loose,flag_lbd=flag_lbd)
    if not rmsd_max == 'lig_mismatch':
        oct_angle_devi, oct_dist_del, max_del_sig_angle, catoms_arr = oct_comp(file_in, angle_ref, catoms_arr)
    else:
        num_coord_metal = -1
        rmsd_max, atom_dist_max = -1, -1
        print('!!!!!Should always match. WRONG!!!!!')
        quit()
    dict_oct_info = {}
    dict_oct_info['num_coord_metal'] = num_coord_metal
    dict_oct_info['rmsd_max'] = rmsd_max
    dict_oct_info['atom_dist_max'] = atom_dist_max
    dict_oct_info['oct_angle_devi_max'] = max(oct_angle_devi)
    dict_oct_info['max_del_sig_angle'] = max_del_sig_angle
    dict_oct_info['dist_del_eq'] = oct_dist_del[0]
    dict_oct_info['dist_del_ax'] = oct_dist_del[1]
    dict_oct_info['dist_del_eq_ax'] = oct_dist_del[2]
    dict_oct_info['dist_del_all'] = oct_dist_del[3]
    print('dict_oct_info', dict_oct_info)
    for ele in std_not_use:
        dict_oct_info[ele] = 'banned_by_user'
    flag_list = []
    for key, values in dict_check.items():
        if not dict_oct_info[key] == 'banned_by_user':
            if dict_oct_info[key] > values:
                flag_list.append(key)
    if not len(flag_list):
        flag_oct = 1  # good structure
        flag_list = 'None'
    else:
        flag_oct = 0
        flag_list = ', '.join(flag_list)
        print('------bad structure!-----')
        print('flag_list:', flag_list)
    return flag_oct, flag_list, dict_oct_info


## See whether a complex is Oct or not.
## output: flag: 1 is good and 0 is bad
##         flag_list: if structure is bad, which test it fails
##         dict_oct_info: values for each metric we check.
def IsOct(file_in, file_init_geo=None, dict_check=dict_oct_check_st,
          std_not_use=[], angle_ref=oct_angle_ref, flag_catoms=False,
          catoms_arr=None):
    num_coord_metal, catoms = get_num_coord_metal(file_in)
    if not catoms_arr == None:
        catoms = catoms_arr
        num_coord_metal = len(catoms_arr)

    # if file_init_geo != None and num_coord_metal>=6:
    #     rmsd_max, atom_dist_max = ligand_comp_org(file_in, file_init_geo)
    # else:
    #     rmsd_max, atom_dist_max = -1, -1
    oct_angle_devi, oct_dist_del, max_del_sig_angle = [-1, -1], [-1, -1, -1, -1], -1
    rmsd_max, atom_dist_max = -1, -1
    catoms_arr = catoms
    if num_coord_metal >= 6:
        if not file_init_geo == None:
            rmsd_max, atom_dist_max = ligand_comp_org(file_in, file_init_geo)
        if not rmsd_max == 'lig_mismatch':
            num_coord_metal = 6
            oct_angle_devi, oct_dist_del, max_del_sig_angle, catoms_arr = oct_comp(file_in, angle_ref, catoms_arr)
        else:
            num_coord_metal = -1
            rmsd_max, atom_dist_max = -1, -1
    dict_oct_info = {}
    dict_oct_info['num_coord_metal'] = num_coord_metal
    dict_oct_info['rmsd_max'] = rmsd_max
    dict_oct_info['atom_dist_max'] = atom_dist_max
    dict_oct_info['oct_angle_devi_max'] = max(oct_angle_devi)
    dict_oct_info['max_del_sig_angle'] = max_del_sig_angle
    dict_oct_info['dist_del_eq'] = oct_dist_del[0]
    dict_oct_info['dist_del_all'] = oct_dist_del[3]
    print('dict_oct_info', dict_oct_info)
    for ele in std_not_use:
        dict_oct_info[ele] = 'banned_by_user'
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
        if not dict_oct_info[key] == 'banned_by_user':
            if dict_oct_info[key] > values:
                flag_list.append(key)
    if num_coord_metal < 6:
        flag_list.append('num_coord_metal')
    # if num_coord_metal == -1:
    #     flag_list.remove('rmsd_max')
    #     flag_list.remove('atom_dist_max')
    ## Case when the num_coord_metal > 6 but still forms a octahedral.
    # if ('num_coord_metal' in flag_list) and (not 'oct_angle_devi_max' in flag_list) and \
    #         (not 'dist_del_eq' in flag_list) and (not 'dist_del_ax' in flag_list) and \
    #         (not 'dist_del_eq_ax' in flag_list):
    #     num_coord_metal = 6
    #     flag_list.remove('num_coord_metal')

    if not len(flag_list):
        flag_oct = 1  # good structure
        flag_list = 'None'
    else:
        flag_oct = 0
        flag_list = ', '.join(flag_list)
        print('------bad structure!-----')
        print('flag_list:', flag_list)
    if not flag_catoms:
        return flag_oct, flag_list, dict_oct_info
    else:
        return flag_oct, flag_list, dict_oct_info, catoms_arr


def IsStructure(file_in, file_init_geo=None, dict_check=dict_oct_check_st,
                std_not_use=[], angle_ref=oct_angle_ref, num_coord=6):
    num_coord_metal, catoms = get_num_coord_metal(file_in)

    if file_init_geo != None:
        rmsd_max, atom_dist_max = ligand_comp_org(file_in, file_init_geo)
    else:
        rmsd_max, atom_dist_max = -1, -1
    if num_coord_metal >= num_coord:
        struct_angle_devi, struct_dist_del, max_del_sig_angle, catoms_arr = oct_comp(file_in, angle_ref)
    else:
        struct_angle_devi, struct_dist_del, max_del_sig_angle = [-1, -1], [-1, -1, -1, -1], -1
    dict_struct_info = {}
    dict_struct_info['num_coord_metal'] = num_coord_metal
    dict_struct_info['rmsd_max'] = rmsd_max
    dict_struct_info['atom_dist_max'] = atom_dist_max
    dict_struct_info['struct_angle_devi_max'] = max(struct_angle_devi)
    dict_struct_info['max_del_sig_angle'] = max_del_sig_angle
    dict_struct_info['dist_del_eq'] = struct_dist_del[0]
    dict_struct_info['dist_del_all'] = struct_dist_del[3]
    print('dict_struct_info', dict_struct_info)
    for ele in std_not_use:
        dict_struct_info[ele] = 'banned_by_user'
    flag_list = []

    for key, values in dict_check.items():
        if not dict_struct_info[key] == 'banned_by_user':
            if dict_struct_info[key] > values:
                flag_list.append(key)
    ## Case when the num_coord_metal > 6 but still forms a structahedral.
    if ('num_coord_metal' in flag_list) and (not 'struct_angle_devi_max' in flag_list) and \
            (not 'dist_del_eq' in flag_list) and (not 'dist_del_ax' in flag_list) and \
            (not 'dist_del_eq_ax' in flag_list):
        num_coord_metal = num_coord
        flag_list.remove('num_coord_metal')

    if not len(flag_list):
        flag_struct = 1  # good structure
        flag_list = 'None'
    else:
        flag_struct = 0
        flag_list = ', '.join(flag_list)
        print('------bad structure!-----')
        print('flag_list:', flag_list)
    return flag_struct, flag_list, dict_struct_info


## input: _path: path for opt geo
##         path_init_geo: path for init geo
## output: list of info: ['unique_num', 'oxstate', 'spinmult', 'flag_oct', 'flag_list', 'dict_oct_info']
def loop_structure(_path, path_init_geo):
    charac = ['unique_num', 'oxstate', 'spinmult', 'flag_oct', 'flag_list']
    info_tot = []
    for dirpath, dirs, files in sorted(os.walk(_path)):
        for name in sorted(files):
            if name.split('.')[1] == 'xyz':
                unique_num, oxstate, spinmult = name.split('_')[0], name.split('_')[4][1:], name.split('_')[5][2:]
                print(unique_num, oxstate, spinmult)
                file_in = '%s/%s' % (dirpath, name)
                ## ---You may add the function to find initial geo here.
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
    # write_list_2_xlsx(charac, info_tot)
    return charac, info_tot


## find the file name by matching unique_ID and spin mulplicity
def find_file_with_unique_num(_path, unique_num):
    filename = None
    for dirpath, dirs, files in sorted(os.walk(_path)):
        for name in sorted(files):
            if name.split('_')[0] == str(unique_num):
                filename = '%s/%s' % (dirpath, name)
                break
    return filename


## Only use this function when they have roughly have the same naming.
def gen_file_with_name(path_init_geo, name_opt):
    name_opt = name_opt.split('_')
    name_opt = '_'.join(name_opt[:len(name_opt) - 1])
    name_init = '%s/%s_mols.xyz' % (path_init_geo, name_opt)
