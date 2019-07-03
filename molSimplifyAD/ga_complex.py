import glob, math, numpy, subprocess, os, random, shutil
import sys, shlex, time, datetime
from molSimplifyAD.ga_tools import *
from molSimplify.Classes.mol3D import *
#####IMPORTS ADDED BY ADITYA
from molSimplify.Classes.atom3D import *
from molSimplify.Informatics.RACassemble import *
from molSimplify.Classes.globalvars import *
from molSimplify.Scripts.io import *
from molSimplify.Scripts.structgen import *
############################


class octahedral_complex:
    def __init__(self, ligands_list):
        self.free_sites = [1, 2, 3, 4, 5, 6]
        ### mark 3x bidentate
        self.three_bidentate = False
        self.ready_for_assembly = False
        self.has_vacant = False
        self.ligands_list = get_ligands()
        self.gene_template = get_gene_template()
        self.metals_list = get_metals()
        if isKeyword("oxocatalysis"):
            self.metals_list = get_metals(first_row=True)
        self.ox_list = get_ox_states()
        self.core = 'U'
        self.symclass = isKeyword("symclass")
        if self.gene_template['legacy']:
            self.ox = 'U'
            self.ax_dent = False
            self.eq_dent = False
            self.eq_ligands = list()
            self.ax_ligands = list()
            self.ax_inds = list()
            self.eq_inds = list()
            # GA_run = get_current_GA()
        else:
            if self.gene_template['ox']:
                self.ox = 'U'
            if self.gene_template['spin']:
                self.spin = 'U'
            self.dent_list = False
            self.ligands = list()
            self.inds = list()
            self.occs = list()
        self.ahf = int(isKeyword("exchange"))  # % HFX, B3LYP = 20

    def random_gen(self):
        if self.gene_template['legacy']:
            self._get_random_metal()
            self._get_random_ox()
            self._get_random_equitorial()
            self._get_random_axial()
            self._name_self()
        else:
            self._get_random_metal()
            if self.gene_template['ox']:
                self._get_random_ox()
            if self.gene_template['spin']:
                self._get_random_spin()
            self._get_random_ligands()
            self._name_self()

    def _name_self(self):
        ## this creates the gene for a given complex
        # GA_run = get_current_GA()
        if self.gene_template['legacy']:
            self.name = "_".join(
                [str(self.core), str(self.ox), str(self.eq_inds[0]), str(self.ax_inds[0]), str(self.ax_inds[1]),
                 str(self.ahf)])
        else:
            namelist = [str(self.core)]
            if self.gene_template['ox']:
                namelist.append(str(self.ox))
            if self.gene_template['spin']:
                namelist.append(str(self.spin))
            for i in self.inds:
                namelist.append(str(i))
            namelist.append(str(self.ahf))
            self.name = "_".join(namelist)

    def copy(self, partner):
        self.core = partner.core
        self.three_bidentate = partner.three_bidentate
        self.ahf = partner.ahf
        if self.gene_template['legacy']:
            self.ox = partner.ox
            self.ax_dent = partner.ax_dent
            self.ax_inds = partner.ax_inds
            self.ax_ligands = partner.ax_ligands
            self.eq_dent = partner.eq_dent
            self.eq_inds = partner.eq_inds
            self.eq_oc = partner.eq_oc
            self.eq_ligands = partner.eq_ligands
        else:
            if self.gene_template['ox']:
                self.ox = partner.ox
            if self.gene_template['spin']:
                self.spin = partner.spin
            self.dent_list = partner.dent_list
            self.ligands = partner.ligands
            self.inds = partner.inds
            self.occs = partner.occs
        self._name_self()

    def _get_random_metal(self):
        n = len(self.metals_list)
        metal_ind = numpy.random.randint(low=0, high=n)
        # if isOxocatalysis():
        #     while self.metals_list[metal_ind] in ['cr','co']:
        #         metal_ind = numpy.random.randint(low = 0,high = n)
        #         if self.metals_list[metal_ind] in ['fe','mn']:
        #             break
        self.core = metal_ind

    def _get_random_ox(self):
        possible_ox_states = get_ox_states()
        self.ox = numpy.random.choice(possible_ox_states)

    def _get_random_spin(self):
        local_spin_dictionary = spin_dictionary()
        metal_list = get_metals()
        # metal_key = metal_list[self.metal]
        metal_key = metal_list[self.core]
        these_states = local_spin_dictionary[metal_key][self.ox]
        self.spin = numpy.random.choice(these_states)
        print('This is the spin '+ str(self.spin)+' for metal '+str(metal_key)+' and ox '+str(self.ox))

    def _get_random_equitorial(self):
        ### choose the equitorial ligand
        n = len(self.ligands_list)
        eq_ind = numpy.random.randint(low=0, high=n)
        if isKeyword('oxocatalysis'):
            smiles = False
            if hasattr(self.ligands_list[eq_ind][0], '__iter__'):
                smiles = True
            while (self.ligands_list[eq_ind][0] in ['x', 'hydroxyl', 'oxo']) or (
                    smiles and self.ligands_list[eq_ind][0][0] in ['[O--]', '[OH-]']):
                eq_ind = numpy.random.randint(low=0, high=n)
                if (self.ligands_list[eq_ind][0] not in ['x', 'hydroxyl', 'oxo']) or (
                        smiles and self.ligands_list[eq_ind][0][0] in ['[O--]', '[OH-]']):
                    break
            # print('CHOSEN RANDOM EQLIG: ',self.ligands_list[eq_ind])
        # print('choosing '+str(eq_ind) + ' from '+str(n))
        # now we test if it is a SMILEs or molsimplify ligand    
        # print('name  is ' +str(self.ligands_list[eq_ind][0]))
        eq_ligand_properties = self.ligands_list[eq_ind][1]
        self.eq_dent = eq_ligand_properties[0]
        self.eq_oc = int(4 / self.eq_dent)
        self.eq_ligands = [self.ligands_list[eq_ind][0] for i in range(0, self.eq_oc)]
        self.eq_inds = [eq_ind]

    def _get_random_axial(self):
        ### choose axial ligands:
        n = len(self.ligands_list)
        self.ready_for_assembly = False
        self.ax_ligands = list()
        self.ax_inds = list()
        while not self.ready_for_assembly:
            ax_ind = numpy.random.randint(low=0, high=n)
            if isKeyword('oxocatalysis'):
                if hasattr(self.eq_ligands[0], '__iter__'):
                    oxo = find_ligand_idx('[O--]')
                    dent = self.ligands_list[ax_ind][0][1]
                else:
                    oxo = find_ligand_idx('oxo')
                    dent = self.ligands_list[ax_ind][1]
                while (ax_ind == oxo) or (dent > 1):  # No bidentate axials
                    ax_ind = numpy.random.randint(low=0, high=n)
                    if (ax_ind != int(oxo)):
                        break
                # print('CHOSEN RANDOM AXLIG: ',self.ligands_list[ax_ind])
            ax_ligand_properties = self.ligands_list[ax_ind][1]
            ax_dent = ax_ligand_properties[0]
            if ax_dent > 1 and not isKeyword('oxocatalysis'):
                if ((self.eq_dent == 2) and (ax_dent == 2) and (len(self.ax_ligands) == 0)):
                    three_bidentate = True
                    self.ax_ligands = [self.ligands_list[ax_ind][0], self.ligands_list[ax_ind][0]]
                    self.ax_inds = [ax_ind, ax_ind]
                    self.ready_for_assembly = True
                    self.ax_dent = 2
                    self.ax_oc = [2]
                else:
                    self.ready_for_assembly = False
            elif ax_dent == 1:
                # GA_run = get_current_GA()
                if isKeyword("symclass") == "strong" and not isKeyword('oxocatalysis'):
                    self.ax_ligands = [self.ligands_list[ax_ind][0], self.ligands_list[ax_ind][0]]
                    self.ax_inds = [ax_ind, ax_ind]
                    if (len(self.ax_ligands) == 2):
                        self.ax_dent = 1
                        self.ax_oc = [1, 1]
                        self.ready_for_assembly = True
                elif isKeyword('oxocatalysis'):
                    if hasattr(self.eq_ligands[0], '__iter__'):
                        oxo = find_ligand_idx('[O--]')
                    else:
                        oxo = find_ligand_idx('oxo')
                    self.ax_ligands = [self.ligands_list[ax_ind][0], self.ligands_list[oxo][0]]
                    self.ax_inds = [ax_ind, oxo]
                    if (len(self.ax_ligands) == 2):
                        self.ax_dent = 1
                        self.ax_oc = [1, 1]
                        self.ready_for_assembly = True
                    # print('AXLIGS ARE ', self.ax_ligands)
                else:
                    #                    This section is intended to allow for vacant axial sites
                    #                    It is not currently implemented
                    #                    if ((ax_ind == 0) and (not has_zero)):
                    #                        has_zero = True
                    #                        self.ax_ligands.append(self.ligands_dict[ax_ind][0])
                    #                        self.ax_inds.append(ax_ind)
                    #                        self.geo = "spy"
                    #
                    #                    elif ax_ind !=  0:
                    self.ax_ligands.append(self.ligands_list[ax_ind][0])
                    self.ax_inds.append(ax_ind)
                    if (len(self.ax_ligands) == 2):
                        self.ax_dent = 1
                        self.ax_oc = [1, 1]
                        # checks for consistent ordering
                        self.ax_inds.sort(reverse=True)
                        self.ax_ligands = [self.ligands_list[self.ax_inds[0]][0], self.ligands_list[self.ax_inds[1]][0]]
                        # print('after sorting, ax ligands are')
                        # print(self.ax_ligands,self.ax_inds)
                        # time.sleep(5)
                        self.ready_for_assembly = True
                    else:
                        self.ready_for_assembly = False

    def _get_random_ligands(self):
        n = len(self.ligands_list)
        self.ready_for_assembly = False
        self.ligands = list()
        self.inds = list()
        symclass = isKeyword("symclass")
        hyd = find_ligand_idx('hydroxyl')
        oxo = find_ligand_idx('oxo')
        x = find_ligand_idx('x')
        while not self.ready_for_assembly:
            ## get lig
            ind = numpy.random.randint(low=0, high=n)
            if isKeyword('oxocatalysis'):
                while ind in [x, hyd, oxo]:
                    ind = numpy.random.randint(low=0, high=n)

            ligand_properties = self.ligands_list[ind][1]
            lig_dent = ligand_properties[0]
            if symclass in ['weak', 'strong']:  # Switch this to be consistent with SMU
                if len(self.ligands) == 0:
                    if lig_dent in [1, 2, 4]:
                        # print('Found eqlig: ', self.ligands_list[ind])
                        self.ligands += 4 * [self.ligands_list[ind][0]]
                        self.inds += 4 * [ind]
                        eqdent = lig_dent
                elif lig_dent == 1:
                    if symclass == 'weak':
                        self.ligands.append(self.ligands_list[ind][0])
                        self.inds.append(ind)
                        if isKeyword('oxocatalysis'):
                            self.ligands.append(self.ligands_list[oxo][0])
                            self.inds.append(oxo)
                    elif symclass == 'strong':
                        self.ligands += 2 * [self.ligands_list[ind][0]]
                        self.inds += 2 * [ind]
                elif lig_dent == 2 and eqdent == 2 and len(self.ligands) == 4 and False and not isKeyword(
                        'oxocatalysis'):  # Triple bidentate handling, may change later
                    self.ligands += 2 * [self.ligands_list[ind][0]]
                    self.inds += 2 * [ind]
            if len(self.ligands) == 6:
                self.ready_for_assembly = True

        self.ligand_sort()

    def enforce_consistency(self, fixed_inds):
        # fixed inds is what we will hold constant, other spots will be randomized to enforce fixed inds
        used_dent = 0
        dent_valid = False
        hyd = find_ligand_idx('hydroxyl')
        oxo = find_ligand_idx('oxo')
        x = find_ligand_idx('x')
        if self.symclass in ['weak','strong']:
            if not len(set(self.inds[0:4])) == 1:
                print('Fails weak symclass check. ', self.inds[0:4],fixed_inds)
                # print(len(set(self.inds[0:4])))
                # print(self.inds)
                for f_ind in fixed_inds:  ## Loop over the positions that must be fixed
                    if f_ind < 4:
                        self.inds[0:4] = 4 * [self.inds[f_ind]]
                        self.ligands[0:4] = 4 * [self.ligands_list[self.inds[f_ind]][0]]
                        new_eq_dent = self.ligands_list[self.inds[f_ind]][1][0]

                        ## we want to check that the axial ligands are compatible with
                        ## the fixed equ index f_ind. Therefore, we will get the denticity
                        ## of the non-fixed axial ligands
                        axial_indices_that_can_be_changed = [i for i in range(4, 6) if i not in fixed_inds]
                        axial_ligands_that_can_be_changed = [self.inds[i] for i in axial_indices_that_can_be_changed]
                        all_ax_dent = [self.ligands_list[i][1][0] for i in axial_ligands_that_can_be_changed]
                        if (max(all_ax_dent) == 2 and not new_eq_dent == 2) or (isKeyword('oxocatalysis') and max(all_ax_dent)!=1):
                            while max(all_ax_dent) == 2 or (isKeyword('oxocatalysis') and max(all_ax_dent)!=1):
                                if self.symclass == 'strong':
                                    ax_inds = numpy.random.randint(low=0, high=n, size=1)
                                    ax_inds = [ax_inds, ax_inds]
                                else:
                                    ax_inds = numpy.random.randint(low=0, high=n, size=2)
                                all_ax_dent = [self.ligands_list[i][1][0] for i in ax_inds]
                                if isKeyword('oxocatalysis'):
                                    for ax in ax_inds:
                                        if ax in [x, hyd, oxo]:
                                            print('catalytic moiety, randomizing.')
                                            ax_inds = numpy.random.randint(low = 0,high = n,size=2)
                                all_ax_dent = [self.ligands_list[i][1][0] for  i in ax_inds]
                            self.inds[4:6] = ax_inds
                            self.ligands[4:6] = [self.ligands_list[i][0] for i in ax_inds]
                            if isKeyword('oxocatalysis'):
                                self.inds[5] = oxo
                                self.ligands[5] = self.ligands_list[oxo][0]
                    elif f_ind >= 4:
                        this_ax_dent = self.ligands_list[self.inds[f_ind]][1][0]
                        if this_ax_dent == 2:
                            self.inds[4:6] = 2 * [self.inds[f_ind]]
                            self.ligands[4:6] = 2 * [self.ligands_list[self.inds[f_ind]][0]]
                            eq_dent = self.ligands_list[self.inds[0]][1][0]
                            while eq_dent != 2:
                                eq_ind = numpy.random.randint(low=0, high=n)
                                eq_ligand_properties = self.ligands_list[eq_ind][1]
                                eq_dent = eq_ligand_properties[0]
                            self.inds[0:4] = 4 * [eq_ind]
                            self.ligands[0:4] = 4 * [self.ligands_list[eq_ind][0]]
                        this_ax_dent  = self.ligands_list[self.inds[f_ind]][1][0]
                        elif this_ax_dent > 1 and isKeyword('oxocatalysis'):
                            print('ILLEGAL BIDENTATE AXIAL!')
                            sardines
                        else:
                            # print(self.inds[[i for i in range(4,5) if i not in fixed_inds]])
                            other_ax_dent = \
                                self.ligands_list[self.inds[[i for i in range(4, 6) if i not in fixed_inds]][0]][1][0]
                            while other_ax_dent != 1:
                                ax_ind = numpy.random.randint(low=0, high=n)
                                ax_ligand_properties = self.ligands_list[ax_ind][1]
                                other_ax_dent = ax_ligand_properties[0]
                            self.inds[[i for i in range(4, 5) if i not in fixed_inds]] = ax_ind
                            self.ligands[[i for i in range(4, 5) if i not in fixed_inds]] = self.ligands_list[ax_ind][0]

            if self.symclass == 'strong':
                if (not len(set(self.inds[4:6])) == 1):
                    print('Fails strong symclass check.')
                    for f_ind in fixed_inds:  ## Loop over the positions that must be fixed
                        if f_ind >= 4:
                            self.inds[4:6] = 2 * [self.inds[f_ind]]
                            self.ligands[4:6] = 2 * [self.ligands_list[self.inds[f_ind]][0]]
        self.ligand_sort()

    def ligand_sort(self):
        self.inds[4:6].sort(reverse=True)
        templist = [self.inds[0], self.inds[2]]
        templist.sort(reverse=True)
        self.inds[0], self.inds[2] = templist[0], templist[1]
        templist = [self.inds[1], self.inds[3]]
        templist.sort(reverse=True)
        self.inds[1], self.inds[3] = templist[0], templist[1]
        self.ligands = [self.ligands_list[x][0] for x in self.inds]
        self.lig_occs = 6 * [0]
        ind = 0
        while ind < 6:
            ligand_properties = self.ligands_list[self.inds[ind]][1]
            lig_dent = ligand_properties[0]
            # print(ligand_properties)
            self.lig_occs[ind] = 1
            ind += lig_dent

    def examine(self):
        print("name is " + self.name)
        if self.gene_template['legacy']:
            print("eq", self.eq_ligands, self.eq_inds)
            print("axial", self.ax_ligands, self.ax_inds)
        else:
            print("ligands", self.ligands, self.inds)

    def encode(self, gene):
        self.random_gen()
        ll = gene.split("_")
        ll = [int(item) for item in ll]
        self.replace_metal(ll[0])
        if self.gene_template['legacy']:
            # print('ll is '+str(ll))
            self.replace_ox(ll[1])
            self.replace_equitorial([ll[2]])
            self.replace_axial(ll[3:5])
            self.ahf = ll[5]
        else:
            inds = []
            # print(ll)
            current_index = 1
            if self.gene_template['ox']:
                self.replace_ox(ll.pop(current_index))
            if self.gene_template['spin']:
                self.replace_spin(ll.pop(current_index))
            while len(ll) > 2:
                # print(ll)
                inds.append(ll.pop(current_index))
            self.ahf = ll.pop(current_index)
            # print('This is inds', inds)
            self.replace_ligands(inds)
        self._name_self()

    def replace_metal(self, new_metal_ind):
        self.core = new_metal_ind

    def replace_ox(self, new_ox):
        self.ox = new_ox

    def replace_spin(self, new_spin):
        self.spin = new_spin

    def replace_ahf(self, new_ahf):
        self.ahf = new_ahf

    def replace_ligands(self, inds):
        self.ready_for_assembly = False
        self.ligands = list()
        self.inds = list()
        for ind in inds:
            ## get lig
            self.ligands.append(self.ligands_list[ind][0])
            self.inds.append(ind)
        self.ready_for_assembly = True
        self.ligand_sort()
        if self.gene_template['legacy']:
            self.eq_inds = [self.inds[0]]
            self.ax_inds = self.inds[4:6]
            self.replace_equitorial(self.eq_inds)
            self.replace_axial(self.ax_inds)
        self._name_self()

    def replace_equitorial(self, new_eq_ind):
        # print('in repcoding, setting eq to ' + str(new_eq_ind))
        eq_ligand_properties = self.ligands_list[new_eq_ind[0]][1]
        self.eq_dent = eq_ligand_properties[0]
        self.eq_oc = int(4 / self.eq_dent)
        self.eq_ligands = [self.ligands_list[new_eq_ind[0]][0] for i in range(0, self.eq_oc)]
        self.eq_inds = new_eq_ind
        if (self.ax_dent == 1) or ((self.ax_dent == 2) and (self.eq_dent == 2) and not isKeyword(
                'oxocatalysis')):  # No triple bidentate in oxocat
            ## everything is ok!
            if (self.ax_dent == 2):
                self.three_bidentate = True
        else:  ## this complex cannot exist. keeping  equitorial,
            ## regenerating axial ligands
            # print("complex with" + str(self.eq_ligands) + " and " + str(self.ax_ligands) + " cannot exist, randomizing axial")
            self._get_random_axial()
        self._name_self()

    def replace_axial(self, new_ax_ind):
        self.ax_ligands = list()
        self.ax_inds = list()
        n = len(self.ligands_list)
        # print(new_ax_ind)
        self.three_bidentate = False
        for i, indices in enumerate(new_ax_ind):
            # print('trying to add ',indices )
            ax_ligand_properties = self.ligands_list[indices][1]
            ax_dent = ax_ligand_properties[0]
            # print('ax dent ' + str(ax_dent))
            if (ax_dent > 1):
                if (ax_dent == 2) and (i == 0):
                    three_bidentate = True
                    self.ax_ligands = [self.ligands_list[indices][0], self.ligands_list[indices][0]]
                    self.ax_inds = [indices, indices]
                    self.ax_dent = 2
                    self.three_bidentate = True
                    break
                else:
                    print('impossible, ax_dent  = ' + str(ax_dent))
            else:
                self.ax_ligands.append(self.ligands_list[indices][0])
                self.ax_inds.append(indices)
                self.ax_dent = 1

        if (self.three_bidentate) and not (self.eq_dent == 2):
            ## this complex cannot exist
            ## regenerating equitorial ligands
            # print("complex with" + str(self.eq_ligands) + " and " + str(self.ax_ligands) +
            #      " cannot exist, regenerating equitorial ligands")
            self.ready_for_assembly = False
            while not self.ready_for_assembly:
                eq_ind = numpy.random.randint(low=0, high=n)
                eq_ligand_properties = self.ligands_list[eq_ind][1]
                eq_dent = eq_ligand_properties[0]
                if (eq_dent == 2):
                    self.eq_dent = 2
                    self.eq_inds = [eq_ind]
                    self.eq_oc = int(4 / eq_dent)
                    self.eq_ligands = [self.ligands_list[eq_ind][0] for i in range(0, self.eq_oc)]
                    self.ready_for_assembly = True
        self._name_self()

    def exchange_ligands(self, partner, eq_swap):
        child = octahedral_complex(self.ligands_list)
        child.copy(self)  # copies this parent

        print("swapping from", partner.name, " to ", self.name)
        self.examine()
        if self.gene_template['legacy']:
            if eq_swap:
                print("swapping equitorial " + str(child.eq_inds) + ' -> ' + str(partner.eq_inds))
                child.replace_equitorial(partner.eq_inds)
            else:
                print("swapping axial" + str(child.ax_inds) + ' -> ' + str(partner.ax_inds))
                child.replace_axial(partner.ax_inds)
        else:
            if eq_swap:
                print("swapping equitorial " + str(child.inds[0:4]) + ' -> ' + str(partner.inds[0:4]))
                child.replace_ligands(partner.inds[0:4] + child.inds[4:6])
                child.enforce_consistency(range(0, 4))
            else:
                print("swapping axial" + str(child.inds[4:6]) + ' -> ' + str(partner.inds[4:6]))
                child.replace_ligands(child.inds[0:4] + partner.inds[4:6])
                child.enforce_consistency(range(4, 6))
                # range(4, 6) says that the axial should be swapped and held constant.

        child.examine()
        child._name_self()
        return child

    def exchange_metal(self, partner):
        child = octahedral_complex(self.ligands_list)
        child.copy(self)  # copies this parent
        print("swapping from", partner.name, " to ", self.name)
        self.examine()
        child.replace_metal(partner.core)
        child._name_self()
        return child

    def exchange_ox(self, partner):
        child = octahedral_complex(self.ligands_list)
        child.copy(self)  # copies this parent
        print("swapping from", partner.name, " to ", self.name)
        self.examine()
        child.replace_ox(partner.ox)
        child._name_self()
        return child

    def mutate(self):
        ## mutates either the axial
        ## or equitorial ligand a random
        # GA_run = get_current_GA()
        child = octahedral_complex(self.ligands_list)
        child.copy(self)  # copies this parent
        n = len(self.ligands_list)
        self.examine()
        fixed_inds = []
        if self.gene_template['legacy']:
            lig_to_mutate = random.randint(0, 4)  # Metal mutation allowed
            print('I think this is 3x bidentate: ', self.three_bidentate, self.ax_dent)
            print('in mutate')
            print('lig to mutate is ' + str(lig_to_mutate))
            if (lig_to_mutate == 0):
                print("mutating equitorial ligand")
                rand_ind = numpy.random.randint(low=0, high=n)
                smiles = False
                if hasattr(self.ligands_list[rand_ind][0], '__iter__'):
                    smiles = True
                if isKeyword('oxocatalysis'):
                    while (self.ligands_list[rand_ind][0] in ['x', 'hydroxyl', 'oxo']) or (
                            smiles and self.ligands_list[rand_ind][0][0] in ['[O--]', '[OH-]']):
                        rand_ind = numpy.random.randint(low=0, high=n)
                        if self.ligands_list[rand_ind][0] not in ['x', 'hydroxyl', 'oxo'] or (
                                smiles and self.ligands_list[rand_ind][0][0] not in ['[O--]', '[OH-]']):
                            break
                    print('MUTATED EQUATIORIAL IS ', self.ligands_list[rand_ind][0])
                child.replace_equitorial([rand_ind])
            elif (lig_to_mutate == 1) or (lig_to_mutate == 2):
                print('mutating axial ligand')
                if isKeyword('oxocatalysis'):
                    lig_to_mutate = 1  # Always keep lig_to_mutate = 1 since do not want to mutate axial moiety
                ready_flag = False
                while not ready_flag:
                    new_ax_list = list()
                    rand_ind = numpy.random.randint(low=0, high=n)
                    if isKeyword('oxocatalysis'):
                        smiles = False
                        if hasattr(self.ligands_list[rand_ind][0], '__iter__'):
                            smiles = True
                        while (self.ligands_list[rand_ind][0] in ['x', 'hydroxyl', 'oxo']) or (
                                smiles and self.ligands_list[rand_ind][0][0] in ['[O--]', '[OH-]']) or \
                                (self.ligands_list[rand_ind][1][0] > 1):  # No bidentate axials
                            rand_ind = numpy.random.randint(low=0, high=n)
                            if ((self.ligands_list[rand_ind][0] not in ['x', 'hydroxyl', 'oxo']) or (
                                    smiles and self.ligands_list[rand_ind][0][0] not in ['[O--]', '[OH-]'])) and \
                                    self.ligands_list[rand_ind][1][0] == 1:
                                break
                    ax_ligand_properties = self.ligands_list[rand_ind][1]
                    ax_dent = ax_ligand_properties[0]
                    if (ax_dent == self.ax_dent):
                        if (lig_to_mutate == 1):
                            print("mutating axial 1 ")
                            if isKeyword('symclass') == "strong" and not isKeyword('oxocatalysis'):
                                new_ax_list = [rand_ind, rand_ind]
                            else:
                                new_ax_list = [rand_ind, self.ax_inds[1]]
                                if not isKeyword('oxocatalysis'):
                                    new_ax_list.sort(reverse=True)
                        elif (lig_to_mutate == 2):
                            print("mutating axial 2 ")
                            if isKeyword('symclass') == "strong":
                                new_ax_list = [rand_ind, rand_ind]
                            else:
                                if not isKeyword('oxocatalysis'):
                                    new_ax_list = [self.ax_inds[0], rand_ind]
                                    new_ax_list.sort(reverse=True)
                                else:
                                    new_ax_list = [self.ax_inds[0], rand_ind]

                        child.ax_dent = 1
                        child.three_bidentate = False
                        ready_flag = True
                    elif (ax_dent == 2) and (self.ax_dent == 1):
                        ## here, we want to add a bidentate but
                        ## need to swap the second ligand too
                        print("swapping both axial ")
                        new_ax_list = [rand_ind, rand_ind]
                        child.ax_dent = 2
                        child.three_bidentate = True
                        ready_flag = True
                print("trying to add " + str(new_ax_list))
                child.replace_axial(new_ax_list)
            elif (lig_to_mutate == 3):  ## metal mutation
                print('mutating metal')
                child._get_random_metal()
                print('mutating ox')
                child._get_random_ox()
            child._name_self()
            child.examine()
        else:
            message = 'old ligands ' + '/'.join([str(i) for i in child.ligands])
            pos_to_mutate = np.random.choice([i for i in range(0, 7) if not i in fixed_inds])
            print(pos_to_mutate)
            if pos_to_mutate < 6:
                ready_for_assembly = False
                while not ready_for_assembly:
                    new_lig = np.random.randint(0, n)
                    new_dent = self.ligands_list[new_lig][1][0]
                    print(new_lig, new_dent)
                    if pos_to_mutate >= 4:
                        if new_dent >= 2:
                            ready_for_assembly = False
                        elif isKeyword('oxocatalysis') and new_dent != 1:
                            ready_for_assembly = False
                        else:
                            print('entered here!')
                            ready_for_assembly = True
                    else:
                        ready_for_assembly = True
                child.inds[pos_to_mutate] = new_lig
                fixed_inds += [pos_to_mutate]
                child.enforce_consistency(fixed_inds)
                newmessage = 'new ligands ' + '/'.join([str(i) for i in child.ligands])
                logger(setup_paths()['state_path'],
                       str(datetime.datetime.now()) + '   ' + str('ligand mutation performed at site ') + str(
                           pos_to_mutate))
            else:
                message = 'old metal/ox/spin ' + "/".join([str(child.core)])
                child._get_random_metal()
                newmessage = 'new metal/ox/spin ' + "/".join([str(child.core)])
                if self.gene_template['ox']:
                    message += '/' + str(child.ox)
                    child._get_random_ox()
                    newmessage += '/' + str(child.ox)
                if self.gene_template['spin']:
                    message += '/' + str(child.spin)
                    child._get_random_spin()
                    newmessage += '/' + str(child.spin)
                ####MUTATE OX???
                logger(setup_paths()['state_path'],
                       str(datetime.datetime.now()) + '   ' + str('entered metal section...'))
            logger(setup_paths()['state_path'], str(datetime.datetime.now()) + '   ' + message + ' --> ' + newmessage)
            child.ligand_sort()
            child._name_self()
            child.examine()
        return child

    def generate_geometry_legacy(self, prefix, spin, path_dictionary, rundirpath, gen):
        # get path info
        ligloc_cont = True
        # set metal properties:
        this_metal = self.metals_list[self.core]

        mol_name = prefix + self.name + "_" + str(spin)

        smicat = False  # holder for connection atoms calls for SMILES ligands
        if self.ax_dent == 1:
            # assemble SMILEs ligands
            liglist = ""  # empty string
            for eq_lig in self.eq_ligands:
                if not hasattr(eq_lig, '__iter__'):  # test if SMILES:
                    # print('molSimplify ligand in eq position  '+ str(eq_lig))
                    liglist += " " + str(eq_lig).strip("'[]'")
                elif hasattr(eq_lig, '__iter__'):  # this is the mark of SMILES strings:
                    # print('SMILEs ligand in eq position  '+ str(eq_lig))
                    liglist += " " + "'" + str(eq_lig[0]) + "'"
                    if not smicat:  # false on first hit
                        smicat = " [" + str(eq_lig[1]).replace("'", "")  # cat list
                    else:
                        smicat += ",  " + str(eq_lig[1]).replace("'", "")  # cat list
            # print(' after eq-shell, liglist is ' + liglist)
            # print(' after eq-shell, smicat is ' + str(smicat))
            for ax_lig in self.ax_ligands:
                # print('processing ligand ' + str(ax_lig))
                if not hasattr(ax_lig, '__iter__'):  # test if SMILES:
                    # print('molSimplify ligand in ax position  '+ str(ax_lig))
                    liglist += " " + str(ax_lig).strip("'[]'")
                elif hasattr(ax_lig, '__iter__'):  # this is the mark of SMILES strings:
                    # print('SMILEs ligand in ax position  '+ str(ax_lig))
                    liglist += " " + "'" + str(ax_lig[0]) + "'"
                    if not smicat:  # false on first hit
                        smicat = " [ " + str(ax_lig[1]).replace("'", "")  # cat list
                    else:
                        smicat += ",  " + str(ax_lig[1]).replace("'", "")  # cat list
            # liglist.replace("'","")
            if smicat:
                smicat += "]"
            # print(' after ax-shell, liglist is ' + liglist)
            # print(' after ax-shell, smicat is ' + str(smicat))
            ligloc = 1
            ligalign = 0
        elif self.ax_dent == 2:
            half = len(self.eq_ligands) / 2
            liglist = ""  # empty string
            # first eq ligand 
            eq_lig = self.eq_ligands[0]
            # print('processing ligand ' + str(eq_lig))
            if not hasattr(eq_lig, '__iter__'):  # test if SMILES:
                # print('molSimplify ligand in eq position  '+ str(eq_lig))
                liglist += " " + str(eq_lig).strip("'[]'")
            elif hasattr(eq_lig, '__iter__'):  # this is the mark of SMILES strings:
                # print('SMILEs ligand in eq position  '+ str(eq_lig))
                liglist += " " + str(eq_lig[0])
                if not smicat:  # false on first hit
                    smicat = " [" + str(eq_lig[1]).replace("'", "")  # cat list
                else:
                    smicat += ",  " + str(eq_lig[1]).replace("'", "")  # cat list
            # print(' after eq-shell 1, liglist is ' + liglist)
            # print(' after eq-shell 1, smicat is ' + str(smicat))

            # only one axial ligand allowed
            ax_lig = self.ax_ligands[0]
            # print('processing ligand ' + str(ax_lig))
            if not hasattr(ax_lig, '__iter__'):  # test if SMILES:
                # print('molSimplify ligand in ax position  '+ str(ax_lig))
                liglist += " " + str(ax_lig).strip("'[]'")
            elif hasattr(ax_lig, '__iter__'):  # this is the mark of SMILES strings:
                # print('SMILEs ligand in ax position  '+ str(ax_lig))
                liglist += " " + str(ax_lig[0])
                if not smicat:  # false on first hit
                    smicat = " [" + str(ax_lig[1]).replace("'", "")  # cat list
                else:
                    smicat += ",  " + str(ax_lig[1]).replace("'", "")  # cat list

            # second eq ligand 
            eq_lig = self.eq_ligands[1]
            # print('processing ligand ' + str(eq_lig))
            if not hasattr(eq_lig, '__iter__'):  # test if SMILES:
                # print('molSimplify ligand in eq position  '+ str(eq_lig))
                liglist += " " + str(eq_lig).strip("'[]'")
            elif hasattr(eq_lig, '__iter__'):  # this is the mark of SMILES strings:
                # print('SMILEs ligand in eq position  '+ str(eq_lig))
                liglist += " " + str(eq_lig[0])
                if not smicat:  # false on first hit
                    smicat = " [" + str(eq_lig[1]).replace("'", "")  # cat list
                else:
                    smicat += ",  " + str(eq_lig[1]).replace("'", "")  # cat list
            # print(' after eq-shell 2, liglist is ' + liglist)
            # print(' after eq-shell 2, smicat is ' + str(smicat))
            liglist.replace("'", "")
            if smicat:
                smicat += "]"
            #            liglist = (str([str(element).strip("'[]'") for element in (self.eq_ligands[:half])] +
            #             [str(self.ax_ligands[0]).strip("'[]'")] +
            #              [str(element).strip("'[]'") for element in (self.eq_ligands[half:])]).strip("[]")).replace("'", "")
            ligloc = 'false'
            ligalign = 0
        ## disable force field opt
        ## if not DFT
        if isKeyword('DFT'):
            ff_opt = 'A'
        else:  # optional denticity based FF control
            if max([self.ax_dent, self.eq_dent]) == 1:
                ff_opt = 'A'
            else:
                ff_opt = 'A'
        if smicat:
            ff_opt = 'BA'
        ## get custom exchange fraction
        # this_GA = get_current_GA()
        use_old_optimizer = isKeyword('old_optimizer')
        exchange = isKeyword('exchange')
        optimize = isKeyword('optimize')

        if optimize:
            # print(' setting up GEO optimization ')
            rty = 'minimize'
            scrpath = "scr/geo/gen_" + str(gen) + '/' + mol_name
        else:
            rty = 'energy'
            scrpath = "scr/sp/gen_" + str(gen) + '/' + mol_name

        geometry = "oct"

        ## set paths for generation
        ms_dump_path = path_dictionary["molsimplify_inps"] + 'ms_output.txt'
        ms_error_path = path_dictionary["molsimplify_inps"] + 'ms_errors.txt'
        jobpath = path_dictionary["job_path"] + mol_name + '.in'
        inpath = path_dictionary["infiles"] + mol_name + '.in'

        geometry_path = path_dictionary["initial_geo_path"] + '/' + mol_name + '.xyz'

        ## check if already exists:
        geo_exists = os.path.isfile(path_dictionary["initial_geo_path"] + mol_name + '.xyz')

        # Initialize ANN results dictionary
        ANN_results = {}
        if not (geo_exists):
            print('generating ' + str(mol_name) + ' with ligands ' + str(self.eq_ligands) + ' and' + str(
                self.ax_ligands))
            try:
                # if True:
                ### Add protection for brackets "()" in bash commands.
                liglist = protect_lig_bash(liglist)
                with open(ms_dump_path, 'a') as ms_pipe:
                    with open(ms_error_path, 'a') as ms_error_pipe:
                        call = " ".join(
                            ["molsimplify ", '-core ' + this_metal, '-lig ' + str(liglist), '-ligocc 1,1,1,1,1,1',
                             '-rundir ' + "'" + rundirpath.rstrip("/") + "'", '-keepHs yes,yes,yes,yes,yes,yes',
                             '-jobdir', 'temp',
                             '-coord 6', '-ligalign ' + str(ligalign), '-ligloc ' + str(ligloc), '-calccharge yes',
                             '-name ' + "'" + mol_name + "'",
                             '-geometry ' + geometry, '-spin ' + str(spin), '-oxstate ' + str(self.ox),
                             '-exchange ' + str(exchange),
                             '-qccode TeraChem', '-runtyp ' + rty, "-ffoption " + ff_opt, ' -ff UFF'])
                        if smicat:
                            call += ' -smicat ' + smicat

                        if isKeyword('oxocatalysis'):
                            call += ' -qoption dftd,d3 -qoption min_maxiter,1100'
                        print(call)
                        p2 = subprocess.Popen(call, stdout=ms_pipe, stderr=ms_error_pipe, shell=True)
                        p2.wait()

                assert (os.path.isfile(rundirpath + 'temp' + '/' + mol_name + '.molinp'))
                shutil.move(rundirpath + 'temp' + '/' + mol_name + '.molinp',
                            path_dictionary["molsimplify_inps"] + '/' + mol_name + '.molinp')
                shutil.move(rundirpath + 'temp' + '/' + mol_name + '.xyz', geometry_path)
            except:
                print('Error: molSimplify failure when calling ')
                print(call)
                sys.exit()

            # if this_GA.config['symclass']=="strong":
            properties = ['split', 'split_dist' 'homo', 'homo_dist', 'gap', 'gap_dist', 'oxo', 'oxo_dist',
                          'hat', 'hat_dist', 'oxo20', 'oxo20_dist', 'homo_empty', 'homo_empty_dist']
            with open(rundirpath + 'temp' + '/' + mol_name + '.report') as report_f:
                for line in report_f:
                    for prop in properties:
                        current_prop = line.split(',')[0]
                        if current_prop == prop:
                            print('****')
                            print(line)
                            prop_val = float(line.split(",")[1])
                            ANN_results.update({prop: float(line.split(",")[1])})
                            print('ANN_' + prop + ' is ' + "{0:.2f}".format(prop_val))
                if len(list(set(properties).difference(ANN_results.keys()))) > 0 and not isKeyword('DFT'):
                    for i in properties:
                        if i not in ANN_results.keys():
                            ANN_results.update(
                                {i: float(10000)})  # Chosen to be arbitrarily large to reduce the fitness value to 0.
                            print(str(i) + ' set to 10000 in ANN_results, chosen so that the fitness goes to 0. The key was not present.')
                        else:
                            print(str(i) + ' set to ' + str(ANN_results[i]) + ' since the key was present')

            if isKeyword('oxocatalysis') and 'oxo' in liglist and isKeyword('DFT'):  # Subbing in 1.65 as Oxo BL
                print('Modifying initial oxo geom file ' + mol_name + '.xyz to have oxo BL 1.65')
                geo_ref_file = open(path_dictionary["initial_geo_path"] + '/' + mol_name + '.xyz', 'r')
                lines = geo_ref_file.readlines()
                geo_ref_file.close()
                geo_replacement = open(path_dictionary["initial_geo_path"] + '/' + mol_name + '.xyz', 'w')
                adjusted_lines = lines[:-1] + [lines[-1][:-9] + '1.650000\n']
                geo_replacement.writelines(adjusted_lines)
                geo_replacement.close()

            shutil.move(rundirpath + 'temp' + '/' + mol_name + '.report',
                        path_dictionary["ms_reps"] + '/' + mol_name + '.report')

            ## write the job file
            if os.path.isfile(rundirpath + 'temp/' + mol_name + '.in'):
                with open(jobpath, 'w') as newf:
                    with open(rundirpath + 'temp/' + mol_name + '.in', 'r') as oldf:
                        for line in oldf:
                            if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line):
                                newf.writelines(line)
                    if isKeyword("molscontrol"):
                        newf.writelines("scrdir " + "./" + "\n")
                    else:
                        newf.writelines("scrdir " + scrpath + "\n")
                    os.remove(rundirpath + 'temp/' + mol_name + '.in')
            else:
                print(self.ox, self.spin)
                self.write_jobfile(jobpath, runtype=rty)

            ### check if ligands in old optimizer list
            old_optimizer_list = get_old_optimizer_ligand_list()
            # use_old_optimizer = False
            for ligs in (self.eq_ligands + self.ax_ligands):
                if ligs in old_optimizer_list:
                    use_old_optimizer = True
            ### make an infile!
            create_generic_infile(jobpath, restart=False, use_old_optimizer=use_old_optimizer)
            flag_oct, _, _ = self.inspect_initial_geo(geometry_path)
            if flag_oct == 0:
                print('Bad initial geometry. Setting all of the fitness values to 0 so it is not used.')
                print('All ANN dictkeys set to 10000, chosen so that the fitness goes to 0.')
                ANN_results = {k: float(10000) for k in ANN_results}
                for i in properties:
                    ANN_results.update({i: float(10000)})  # Chosen to be arbitrarily large to reduce the fitness value to 0.
                    print(str(i) + ' set to 10000 in ANN_results, chosen so that the fitness goes to 0.')
        should_geo_exist = not (isKeyword('no_geo')) #will be false if DFT is off.
        #Initialize ANN results dictionary
        ANN_results = {}
        properties = ['split','split_dist' 'homo','homo_dist','gap','gap_dist','oxo','oxo_dist',
                        'hat','hat_dist','oxo20','oxo20_dist','homo_empty','homo_empty_dist']
        if not (geo_exists) or not (should_geo_exist):
            print('generating '+ str(mol_name) + ' with ligands ' + str(self.eq_ligands) + ' and'  + str(self.ax_ligands))
            try:
            # if True:
                ### Add protection for brackets "()" in bash commands.
                liglist = protect_lig_bash(liglist)
                with open(ms_dump_path,'a') as ms_pipe:
                    with open(ms_error_path,'a') as ms_error_pipe:
                        call = " ".join(["molsimplify " ,'-core ' + this_metal,'-lig ' +str(liglist),'-ligocc 1,1,1,1,1,1',
                                 '-rundir ' +"'"+ rundirpath.rstrip("/")+"'",'-keepHs yes,yes,yes,yes,yes,yes','-jobdir','temp',
                                 '-coord 6','-ligalign '+str(ligalign),'-ligloc ' + str(ligloc),'-calccharge yes','-name '+"'"+mol_name+"'",
                                 '-geometry ' + geometry,'-spin ' + str(spin),'-oxstate '+ str(self.ox), '-exchange '+str(exchange),
                                 '-qccode TeraChem','-runtyp '+rty,"-ffoption "+ff_opt,' -ff UFF'])
                        if smicat:
                            call += ' -smicat ' + smicat

                        if isKeyword('oxocatalysis'):
                            call += ' -qoption dftd,d3 -qoption min_maxiter,1100'
                        if not should_geo_exist:
                            call += ' -reportonly True'
                        print(call)
                        p2 = subprocess.Popen(call,stdout = ms_pipe,stderr=ms_error_pipe, shell=True)
                        p2.wait()

                assert(os.path.isfile(rundirpath + 'temp'+'/' + mol_name + '.molinp'))
                shutil.move(rundirpath + 'temp'+'/' + mol_name + '.molinp', path_dictionary["molsimplify_inps"]+'/' + mol_name + '.molinp')
                if should_geo_exist:
                    shutil.move(rundirpath + 'temp'+'/' + mol_name + '.xyz', geometry_path)
            except:
                print('Error: molSimplify failure when calling ')
                print(call)
                sys.exit()

            with open(rundirpath + 'temp' +'/' + mol_name + '.report') as report_f:
                for line in report_f:
                    for prop in properties:
                        current_prop = line.split(',')[0]
                        if current_prop == prop:
                            print('****')
                            print(line)
                            prop_val = float(line.split(",")[1])
                            ANN_results.update({prop:float(line.split(",")[1])})
                            print('ANN_'+prop+' is ' +"{0:.2f}".format(prop_val))
                if len(list(set(properties).difference(ANN_results.keys())))>0 and not isKeyword('DFT'):
                    for i in properties:
                        if i not in ANN_results.keys():
                            ANN_results.update({i:float(10000)}) #Chosen to be arbitrarily large to reduce the fitness value to 0.
                            print(str(i)+ ' set to 10000 in ANN_results, chosen so that the fitness goes to 0. The key was not present.')
                        else:
                            print(str(i)+ ' set to '+str(ANN_results[i])+' since the key was present')

            if isKeyword('oxocatalysis') and 'oxo' in liglist and isKeyword('DFT'): #Subbing in 1.65 as Oxo BL
                print('Modifying initial oxo geom file '+ mol_name + '.xyz to have oxo BL 1.65')
                geo_ref_file = open(path_dictionary["initial_geo_path"] +'/'+ mol_name + '.xyz','r')
                lines = geo_ref_file.readlines()
                geo_ref_file.close()
                geo_replacement = open(path_dictionary["initial_geo_path"] +'/'+ mol_name + '.xyz','w')
                adjusted_lines = lines[:-1]+[lines[-1][:-9]+'1.650000\n']
                geo_replacement.writelines(adjusted_lines)
                geo_replacement.close()

            shutil.move(rundirpath + 'temp' +'/' + mol_name + '.report', path_dictionary["ms_reps"] +'/'+ mol_name + '.report')

            ## write the job file
            if should_geo_exist:
                with open(jobpath,'w') as newf:
                    with open(rundirpath + 'temp/' + mol_name + '.in','r') as oldf:
                        for line in oldf:
                            if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line):
                                newf.writelines(line)
                    newf.writelines("scrdir " +scrpath + "\n")
                os.remove(rundirpath + 'temp/' + mol_name + '.in')

                ### check if ligands in old optimizer list
                old_optimizer_list = get_old_optimizer_ligand_list()
                # use_old_optimizer = False
                for ligs in (self.eq_ligands+self.ax_ligands):
                    if ligs in old_optimizer_list:
                        use_old_optimizer = True
                ### make an infile!
                create_generic_infile(jobpath,restart=False,use_old_optimizer=use_old_optimizer)
                flag_oct, _, _ = self.inspect_initial_geo(geometry_path)
                if flag_oct == 0:
                    print('Bad initial geometry. Setting all of the fitness values to 0 so it is not used.')
                    print('All ANN dictkeys set to 10000, chosen so that the fitness goes to 0.')
                    ANN_results = {k:float(10000) for k in ANN_results}
                    for i in properties:
                        ANN_results.update({i:float(10000)}) #Chosen to be arbitrarily large to reduce the fitness value to 0.
                        print(str(i)+ ' set to 10000 in ANN_results, chosen so that the fitness goes to 0.')
            else:
                with open(jobpath,'w') as newf:
                    newf.writelines('ANN Prediction Directly from Molecular Graph')
                flag_oct = 1
        else:
            flag_oct = 1
        sorted_ANN_results = {}
        if len(ANN_results.keys()) > 0:
            for key in sorted(ANN_results.iterkeys()):
                sorted_ANN_results[key] = ANN_results[key]
        return jobpath, mol_name, sorted_ANN_results, flag_oct

    def generate_geometry(self, prefix, ox, spin, path_dictionary, rundirpath, gen, tmcdoc=False):
        # get path info
        ligloc_cont = True
        # set metal properties:
        this_metal = self.metals_list[self.core]

        mol_name = prefix + jobname_from_parts(metal=self.core, ox=ox, lig_inds=self.inds, ahf=self.ahf, spin=spin)

        smicat = False  # holder for connection atoms calls for SMILES ligands
        # print(self.lig_occs)
        purified_ligands = [self.ligands[i] for i in range(6) if self.lig_occs[i] > 0]
        # print(purified_ligands)

        liglist, smicat = SMILES_converter(purified_ligands)
        ligalign = 0
        # print(liglist,smicat)
        if self.lig_occs == [1, 0, 1, 0, 1, 0]:
            ligloc = 'false'
        else:
            ligloc = 1
        ligalign = 0

        ## disable force field opt
        ## if not DFT
        if isKeyword('DFT'):
            ff_opt = 'A'
        else:
            ff_opt = 'no'
        if smicat:
            ff_opt = 'no'
        ## get custom exchange fraction
        # this_GA = get_current_GA()
        use_old_optimizer = isKeyword('old_optimizer')
        exchange = isKeyword('exchange')
        optimize = isKeyword('optimize')

        if optimize:
            # print(' setting up GEO optimization ')
            rty = 'minimize'
            scrpath = "scr/geo/gen_" + str(gen) + '/' + mol_name
        else:
            rty = 'energy'
            scrpath = "scr/sp/gen_" + str(gen) + '/' + mol_name

        geometry = "oct"

        ## set paths for generation
        ms_dump_path = path_dictionary["molsimplify_inps"] + 'ms_output.txt'
        ms_error_path = path_dictionary["molsimplify_inps"] + 'ms_errors.txt'
        jobpath = path_dictionary["job_path"] + mol_name + '.in'
        inpath = path_dictionary["infiles"] + mol_name + '.in'

        geometry_path = path_dictionary["initial_geo_path"] + '/' + mol_name + '.xyz'

        ## check if already exists:
        geo_exists = os.path.isfile(path_dictionary["initial_geo_path"] + mol_name + '.xyz')

        # Initialize ANN results dictionary
        ANN_results = {}
        if not (geo_exists):
            print('generating ' + str(mol_name) + ' with ligands ' + str(liglist))
            try:
                # if True:
                with open(ms_dump_path, 'a') as ms_pipe:
                    with open(ms_error_path, 'a') as ms_error_pipe:
                        call = " ".join(["molsimplify ", '-core ' + this_metal, '-lig ' + str(liglist),
                                         '-ligocc ' + ','.join([str(i) for i in self.lig_occs if i > 0]),
                                         '-rundir ' + "'" + rundirpath.rstrip("/") + "'",
                                         '-keepHs yes,yes,yes,yes,yes,yes', '-jobdir', 'temp',
                                         '-coord 6', '-ligalign ' + str(ligalign), '-ligloc ' + str(ligloc),
                                         '-calccharge yes', '-name ' + "'" + mol_name + "'",
                                         '-geometry ' + geometry, '-spin ' + str(spin), '-oxstate ' + str(ox),
                                         '-exchange ' + str(exchange),
                                         '-qccode TeraChem', '-runtyp ' + rty, "-ffoption " + ff_opt, ' -ff UFF'])
                        if smicat:
                            call += ' -smicat ' + smicat
        properties = ['split','split_dist' 'homo','homo_dist','gap','gap_dist','oxo','oxo_dist',
                        'hat','hat_dist','oxo20','oxo20_dist','homo_empty','homo_empty_dist']
        if not (geo_exists):
            print('generating '+ str(mol_name) + ' with ligands ' + str(liglist))
            # try:
            if True:
                with open(ms_dump_path,'a') as ms_pipe:
                    with open(ms_error_path,'a') as ms_error_pipe:
                        call = " ".join(["molsimplify " ,'-core ' + this_metal,'-lig ' +str(liglist),'-ligocc '+','.join([str(i) for i in self.lig_occs if i>0]),
                                 '-rundir ' +"'"+ rundirpath.rstrip("/")+"'",'-keepHs yes,yes,yes,yes,yes,yes','-jobdir','temp',
                                 '-coord 6','-ligalign '+str(ligalign),'-ligloc ' + str(ligloc),'-calccharge yes','-name '+"'"+mol_name+"'",
                                 '-geometry ' + geometry,'-spin ' + str(spin),'-oxstate '+ str(ox), '-exchange '+str(exchange),
                                 '-qccode TeraChem','-runtyp '+rty,"-ffoption "+ff_opt,' -ff UFF'])
                        if smicat:
                            call += ' -smicat ' + smicat

                        if isKeyword('oxocatalysis'):
                            call += ' -qoption dftd,d3 -qoption min_maxiter,1100'
                        print(call)
#                            p2 = subprocess.call(call,stdout = ms_pipe,stderr=ms_error_pipe, shell=True)
                        p2 = subprocess.Popen(call,stdout = ms_pipe,stderr=ms_error_pipe, shell=True)
                        p2.wait()

                assert(os.path.isfile(rundirpath + 'temp'+'/' + mol_name + '.molinp'))
                shutil.move(rundirpath + 'temp'+'/' + mol_name + '.molinp', path_dictionary["molsimplify_inps"]+'/' + mol_name + '.molinp')
                shutil.move(rundirpath + 'temp'+'/' + mol_name + '.xyz', geometry_path)
            # except:
            #         print('Error: molSimplify failure when calling ')
            #         print(call)
            #         sys.exit()

            #if this_GA.config['symclass']=="strong":

            with open(rundirpath + 'temp' +'/' + mol_name + '.report','r') as report_f:
                for line in report_f:
                    print(line)
                    for prop in properties:
                        current_prop = line.split(',')[0]
                        dist_var = prop+'_dist'
                        if current_prop == prop:
                            print('****')
                            print(line)
                            prop_val = float(line.split(",")[1])
                            ANN_results.update({prop:float(line.split(",")[1])})
                            print('ANN_'+prop+' is ' +"{0:.2f}".format(prop_val))
                if len(list(set(properties).difference(ANN_results.keys())))>0 and not isKeyword('DFT'):
                    for i in properties:
                        if i not in ANN_results.keys():
                            ANN_results.update({i:float(10000)}) #Chosen to be arbitrarily large to reduce the fitness value to 0.
                            print(str(i)+ ' set to 10000 in ANN_results, chosen so that the fitness goes to 0. The key was not present.')
                        else:
                            print(str(i)+ ' set to '+str(ANN_results[i])+' since the key was present')
            if isKeyword('oxocatalysis') and 'oxo' in liglist and isKeyword('DFT'):  # Subbing in 1.65 as Oxo BL
                print('Modifying initial oxo geom file ' + mol_name + '.xyz to have oxo BL 1.65')
                geo_ref_file = open(path_dictionary["initial_geo_path"] + '/' + mol_name + '.xyz', 'r')
                lines = geo_ref_file.readlines()
                geo_ref_file.close()
                geo_replacement = open(path_dictionary["initial_geo_path"] + '/' + mol_name + '.xyz', 'w')
                adjusted_lines = lines[:-1] + [lines[-1][:-9] + '1.650000\n']
                geo_replacement.writelines(adjusted_lines)
                geo_replacement.close()

            shutil.move(rundirpath + 'temp' + '/' + mol_name + '.report',
                        path_dictionary["ms_reps"] + '/' + mol_name + '.report')

            flag_oct, _, _ = self.inspect_initial_geo(geometry_path)
            if flag_oct == 0:
                print('Bad initial geometry. Setting all of the fitness values to 0 so it is not used.')
                print('All ANN dictkeys set to 10000, chosen so that the fitness goes to 0.')
                ANN_results = {k: float(10000) for k in ANN_results}
                for i in properties:
                    ANN_results.update({i: float(10000)})  # Chosen to be arbitrarily large to reduce the fitness value to 0.
                    print(str(i) + ' set to 10000 in ANN_results, chosen so that the fitness goes to 0.')

            ## write the job file
            if os.path.isfile(rundirpath + 'temp/' + mol_name + '.in'):
                with open(jobpath, 'w') as newf:
                    with open(rundirpath + 'temp/' + mol_name + '.in', 'r') as oldf:
                        for line in oldf:
                            if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line):
                                newf.writelines(line)
                            if "charge" in line:
                                charge = int(line.split()[-1])
                            if "spinmult" in line:
                                spinmult = int(line.split()[-1])
                    if isKeyword("molscontrol"):
                        newf.writelines("scrdir " + "./" + "\n")
                    else:
                        newf.writelines("scrdir " + scrpath + "\n")
                os.remove(rundirpath + 'temp/' + mol_name + '.in')
            else:
                raise FileNotFoundError(rundirpath + 'temp/' + mol_name + '.in' + 'does not exist.')

            if tmcdoc:
                print("Using geometries from the database.")
                flag_oct = 1
                self.generate_from_db(tmcdoc, mol_name, path_dictionary, charge, spinmult)


            ### check if ligands in old optimizer list
            old_optimizer_list = get_old_optimizer_ligand_list()
            # use_old_optimizer = False
            for ligs in (liglist):
                if ligs in old_optimizer_list:
                    use_old_optimizer = True
            ### make an infile!
            create_generic_infile(jobpath, restart=False, use_old_optimizer=use_old_optimizer)
            flag_oct, _, _ = self.inspect_initial_geo(geometry_path)
            if flag_oct == 0:
                print('Bad initial geometry. Setting all of the fitness values to 0 so it is not used.')
                print('All ANN dictkeys set to 10000, chosen so that the fitness goes to 0.')
                ANN_results = {k:float(10000) for k in ANN_results}
                for i in properties:
                    ANN_results.update({i:float(10000)}) #Chosen to be arbitrarily large to reduce the fitness value to 0.
                    print(str(i)+ ' set to 10000 in ANN_results, chosen so that the fitness goes to 0.')
        else:
            flag_oct = 1
        sorted_ANN_results = {}
        if len(ANN_results.keys()) > 0:
            for key in sorted(ANN_results.iterkeys()):
                sorted_ANN_results[key] = ANN_results[key]
        return jobpath, mol_name, sorted_ANN_results, flag_oct

    def geo_free_RAC_generator(self,prefix, ox,spin,path_dictionary,rundirpath,gen):
        ##### The geo free RACs get generated on a choice of genes
        ligloc_cont = True
        licores = getlicores()
        if self.gene_template['legacy']:
            inds = 4*self.eq_inds+self.ax_inds
        else:
            inds = self.inds
        this_metal = self.metals_list[self.core]
        mol_name = prefix + jobname_from_parts(metal=self.core, ox=ox, lig_inds=inds, ahf=self.ahf, spin=spin)
        jobpath = path_dictionary["job_path"]  + mol_name + '.in'
        ligands = [str(self.ligands_list[i][0]) for i in inds]
        ligalign = 0
        ligs  = []
        cons  = []
        emsg, complex3D = False, []
        occs0 = []      # occurrences of each ligand
        toccs = 0       # total occurrence count (number of ligands)
        catsmi = []     # SMILES ligands connection atoms
        smilesligs = 0  # count how many smiles strings
        cats0 = []      # connection atoms for ligands
        dentl = []      # denticity of ligands
        connected = []  # indices in core3D of ligand atoms connected to metal
        frozenats = []  # atoms to be frozen in optimization
        freezeangles = False # custom angles imposed
        MLoptbds = []   # list of bond lengths
        rempi = False   # remove dummy pi orbital center of mass atom
        net_lig_charge = 0 
        backbatoms = []
        batslist = []
        bats = []
        metal_mol = mol3D()
        metal_mol.addAtom(atom3D(str(this_metal).capitalize()))
        ### CURRENTLY only works for molsimplify ligands... ###
        for i,name in enumerate(ligands):     
            this_mol, emsg = lig_load(name, licores)
            net_lig_charge += this_mol.charge
            this_mol.convert2mol3D()
            this_lig  = ligand(mol3D(), [],this_mol.denticity)
            this_lig.mol = this_mol
            this_con = this_mol.cat
            ligs.append(this_lig)
            cons.append(this_con)
            if name not in licores.keys():
                cats0.append([0])
                dent_i = len(cats0[-1])
                smilesligs += 1
            else:
                cats0.append(False)
            # otherwise get denticity from ligands dictionary
                if 'pi' in licores[name][2]:
                    dent_i = 1
                else:
                    if isinstance(licores[name][2], (str, unicode)):
                        dent_i = 1
                    else:
                        dent_i = int(len(licores[name][2]))
            # oc_i = int(ligoc[i]) if i < len(ligoc) else 1
            # occs0.append(0)         # initialize occurrences list
            dentl.append(dent_i)    # append denticity to list
            # loop over occurrence of ligand i to check for max coordination
            # for j in range(0,oc_i):
            #     occs0[i] += 1
            #     toccs += dent_i
        # sort by descending denticity (needed for adjacent connection atoms)
        ligandsU,occsU,dentsU = ligs,occs0,dentl # save unordered lists
        indcs = range(0,len(ligs)) #assuming no ligalign
        lig_instances = [ligs[i] for i in indcs]  # sort ligands list
        dents = [dentl[i] for i in indcs]   # sort denticities list
        #### CURRENTLY ASSUMES EQ, AX1, AX2 in that ORDER, CANNOT HANDLE DENTICITIES OUT OF 1, 2, 4 #####
        eq_number = int(4/dents[0])
        eq_cons = eq_number*[cons[0]]
        eq_ligs = eq_number*[lig_instances[0]]
        
        if dents[-2] == 2:
            ax_ligs = 1*[lig_instances[-2]]
            ax_cons = 1*[cons[-2]]
        else:
            ax_ligs = [lig_instances[-2],lig_instances[-1]]
            ax_cons = [cons[-2],cons[-1]]
        
        custom_ligand_dict = {"eq_ligand_list":eq_ligs,
                          "ax_ligand_list":ax_ligs,
                          "eq_con_int_list":eq_cons,
                          "ax_con_int_list":ax_cons}
        ox_modifier = {str(this_metal).capitalize(): -ox}
        core3D = assemble_connectivity_from_parts(metal_mol,custom_ligand_dict)
        descriptor_names, descriptors = get_descriptor_vector(core3D, custom_ligand_dict=custom_ligand_dict, ox_modifier=ox_modifier)
        descriptor_names += ["charge_lig"]
        descriptors += [net_lig_charge]
        try:
            descriptor_names, descriptors = create_OHE(descriptor_names, descriptors, str(this_metal).capitalize(), ox)
        except:
            print('Not adding one hot encoding descriptors...')
        if float(self.ahf) > 1:
            alpha = float(self.ahf)/float(100)
        else:
            alpha = self.ahf
        descriptor_names += ['alpha']
        descriptors += [alpha]
        descriptor_names += ['ox']
        descriptors += [ox]
        descriptor_names += ['spin']
        descriptors += [spin]
        flag_oct = True
        return descriptor_names, descriptors, flag_oct, mol_name,jobpath
        
    def inspect_initial_geo(self,geometry_path):
        ## this function contains the logic for inspecting a
        ## initial geo file and reporting if there are problems with it
        mol = mol3D()  # load blank mol3D()
        if os.path.isfile(geometry_path):
            mol.readfromxyz(geometry_path)
        flag_oct, flag_list, dict_oct_info = mol.IsOct()
        # flag_H = not mol.closest_H_2_metal()[0]
        # flag_oct = flag_oct and flag_H
        return flag_oct, flag_list, dict_oct_info

    def assemble_constraints(self, ox, spin):
        id_keys = ["metal", "ox", "spin", "lig1", "lig2", "lig3", "lig4", "lig5",
                   "lig6"]  ## Note: Here alpha is not considered.
        purified_ligands = [self.ligands[i] for i in range(6) if self.lig_occs[i] > 0]
        print(purified_ligands)
        _liglist, smicat = SMILES_converter(purified_ligands)
        liglist = [ele.strip("'") for ele in _liglist.split()]

        query_constraints = {}
        query_constraints.update({"metal": self.metals_list[self.core]})
        query_constraints.update({"ox": ox})
        query_constraints.update({"spin": spin})
        for ii in range(6):
            query_constraints.update({"lig%d" % (ii + 1): liglist[ii]})
        return query_constraints

    def generate_from_db(self, tmcdoc, mol_name, path_dictionary, charge, spinmult):
        numatoms = len(tmcdoc["init_geo"].split('\n')) - 1
        init_geopath = path_dictionary["initial_geo_path"] + mol_name + '.xyz'
        opt_geopath = path_dictionary["optimial_geo_path"] + mol_name + '.xyz'
        scrpath = path_dictionary["scr_path"] + mol_name
        outpath = path_dictionary["geo_out_path"] + mol_name + '.out'
        with open(init_geopath, "w") as fo:
            fo.write("%d\n" % numatoms)
            fo.write("====Geometry adopted from the database====\n")
            fo.write(tmcdoc["init_geo"])
        with open(opt_geopath, "w") as fo:
            fo.write("%d\n" % numatoms)
            fo.write("====Geometry adopted from the database====\n")
            fo.write(tmcdoc["opt_geo"])
        ensure_dir(scrpath)
        with open(scrpath + '/optim.xyz', "w") as fo:
            fo.write("%d\n" % numatoms)
            fo.write("====Geometry adopted from the database====\n")
            fo.write(tmcdoc["opt_geo"])
        with open(outpath, "w") as fo:
            fo.write("=======This outfile is FAKE and generated artificialy======\n")
            fo.write("Total charge:    %d\n" % int(charge))
            fo.write("Spin multiplicity: %d\n" % int(spinmult))
            fo.write("FINAL ENERGY: %f a.u.\n" % float(tmcdoc["energy"]))
            fo.write("-=#=-     Optimization Converged.     -=#=-\n")
            fo.write("Total processing time: 0.00 sec\n")


