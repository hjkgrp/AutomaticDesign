# Looks into the ANN_results.csv files
##Return file with list of genes and their ANN distances.
##Returns files with genes, their fitnesses, and their frequences per generation.
import csv
import os
from molSimplifyAD.ga_tools import *
#from molSimplifyAD.ga_tools import get_run_dir
from molSimplifyAD.ga_get_general import get_gen_npool


##Finds the ANN distances, writes them to a dictionary and then to a file
def find_distances():
    # print('in find distances')
    lastgen, npool = get_gen_npool(isKeyword('rundir'))
    gene_dist_dict = dict()
    gene_prop_dict = dict()
    gene_name_dict = dict()
    #GA_run = get_current_GA()
    runtype = isKeyword("runtype")
    # print('find dist part 1')
    # lastgen is the last generation that has been run
    for generation in range(lastgen + 1):
        ANN_dir = isKeyword('rundir') + "ANN_output/gen_" + str(generation) + "/ANN_results.csv"
        emsg, ANN_dict = read_ANN_results_dictionary(ANN_dir)
        for keys in list(ANN_dict.keys()):
            #gene, _, _, metal, ox, eqlig, axlig1, axlig2, _, _, _, spin, spin_cat, ahf, _, _ = translate_job_name(keys)
            translate_dict = translate_job_name(keys)
            gene = translate_dict['gene']
            gen = translate_dict['gen']
            slot = translate_dict['slot']
            metal = translate_dict['metal']
            ox = translate_dict['ox']
            liglist = translate_dict['liglist']
            gene_template = get_gene_template()
            if gene_template['legacy']:
                eqlig = liglist[0]
                axlig1 = liglist[1]
                axlig2 = liglist[2]
            gene = translate_dict['gene']
            indlist = translate_dict['indlist']
            spin = translate_dict['spin']
            spin_cat = translate_dict['spin_cat']
            ahf = translate_dict['ahf']
            base_name = translate_dict['basename']
            base_gene = translate_dict['basegene']
            chem_name = translate_dict['chem_name']
            split_energy = float(ANN_dict[keys]['split'])
            if runtype in ['homo','gap']:
                if (split_energy > 0 and int(spin)<=3) or (split_energy < 0 and int(spin)>3):
                    this_prop = float(ANN_dict[keys][runtype])
                    this_dist = float(ANN_dict[keys][runtype + '_dist'])
                    # geneName = "_".join(keys.split('_')[4:10])
                    metal = get_metals()[metal]
                    chem_name = translate_dict['chem_name']
                    if gene in list(gene_dist_dict.keys()):
                        pass
                    else:
                        gene_dist_dict.update({gene: this_dist})
                        gene_prop_dict.update({gene: this_prop})
                        gene_name_dict.update({gene: chem_name})
            elif runtype in ['oxo','hat']:
                # if (spin_cat == 'HS' or (get_metals()[metal] == 'cr' and int(spin) == 2)):
                if gene_template['legacy']:
                    if (spin_cat == isKeyword('spin_constraint')):
                        # print('Entered into HAT and OXO statement because HIGH SPIN')
                        print('Entered into HAT and OXO statement because LOW SPIN')
                        this_prop = float(ANN_dict[keys][runtype])
                        this_dist = float(ANN_dict[keys][runtype + '_dist'])
                        # keys_temp = keys.split('_')[4:10]
                        # geneName = "_".join(keys_temp)
                        metal = get_metals()[metal]
                        # chem_name = translate_dict['chem_name']
                        print((chem_name+' logged in dictionary'))
                        # print('genename is '+geneName)
                        if gene in list(gene_dist_dict.keys()):
                            print(('SKIPPING GENENAME '+str(gene)))
                            pass
                        else:
                            gene_dist_dict.update({gene: this_dist})
                            gene_prop_dict.update({gene: this_prop})
                            gene_name_dict.update({gene: chem_name})
                    elif (isKeyword('spin_constraint') == 'HS') and get_metals()[metal].lower() == 'cr' and ox == 5:
                        print('Cr(V) does not exist in HS')
                        metal = get_metals()[metal]
                        chem_name = translate_dict['chem_name']
                        gene_dist_dict.update({gene: 10000})
                        gene_prop_dict.update({gene: 10000})
                        gene_name_dict.update({gene: chem_name})
                else:
                    this_prop = float(ANN_dict[keys][runtype])
                    this_dist = float(ANN_dict[keys][runtype + '_dist'])
                    gene_dist_dict.update({gene: this_dist})
                    gene_prop_dict.update({gene: this_prop})
                    gene_name_dict.update({gene: chem_name})
            elif runtype == 'split':
                this_prop = float(ANN_dict[keys][runtype])
                this_dist = float(ANN_dict[keys][runtype + '_dist'])
                # geneName = "_".join(keys.split('_')[4:10])
                metal = get_metals()[metal]
                chem_name = translate_dict['chem_name']
                print(chem_name)
                if gene in list(gene_dist_dict.keys()):
                    pass
                else:
                    gene_dist_dict.update({gene: this_dist})
                    gene_prop_dict.update({gene: this_prop})
                    gene_name_dict.update({gene: chem_name})
            elif type(runtype) == list: #Currently only supports spin dependent properties with a spin constraint
                this_prop = []
                this_dist = []
                if gene_template['legacy']:
                    if spin_cat == isKeyword('spin_constraint'): #Constraining this to a single spin state.
                        for run in runtype:
                            this_prop.append(float(ANN_dict[keys][run]))
                            this_dist.append(float(ANN_dict[keys][run + '_dist']))
                        geneName = "_".join(keys.split('_')[4:10])
                        metal = get_metals()[metal]
                        chem_name = translate_dict['chem_name']
                        if geneName in list(gene_dist_dict.keys()):
                            pass
                        else:
                            gene_dist_dict.update({geneName: this_dist})
                            gene_prop_dict.update({geneName: this_prop})
                            gene_name_dict.update({geneName: chem_name})
                        print('Multiple factors in fitness (get distances)')
                    elif (isKeyword('spin_constraint') == 'HS') and get_metals()[metal].lower() == 'cr' and ox == 5:
                        print('Cr(V) does not exist in HS')
                        geneName = "_".join(keys.split('_')[4:10])
                        metal = get_metals()[metal]
                        chem_name = translate_dict['chem_name']
                        if geneName in list(gene_dist_dict.keys()):
                            pass
                        else:
                            gene_dist_dict.update({geneName: [10000,10000]})
                            gene_prop_dict.update({geneName: [10000,10000]})
                            gene_name_dict.update({geneName: chem_name})
                else:
                    # print('in here')
                    # print('gene is ',gene)
                    for run in runtype:
                        this_prop.append(float(ANN_dict[keys][run]))
                        this_dist.append(float(ANN_dict[keys][run + '_dist']))
                    # print('after lists')
                    chem_name = translate_dict['chem_name']
                    # print(chem_name)
                    # print(this_prop, this_dist)
                    if gene in list(gene_dist_dict.keys()):
                        pass
                    else:
                        gene_dist_dict.update({gene: this_dist})
                        gene_prop_dict.update({gene: this_prop})
                        gene_name_dict.update({gene: chem_name})

    ## Writes genes and distances to a .csv file
    write_path = isKeyword('rundir') + "statespace/all_distances.csv"
    if not os.path.isfile(write_path):
        open(write_path, 'w').close()
    emsg = write_dictionary(gene_dist_dict, write_path)
    if emsg:
        print(emsg)
    print('Now printing all dictionaries (dist, prop, name)')
    print(gene_dist_dict)
    print(gene_prop_dict)
    print(gene_name_dict)
    return gene_dist_dict, npool, gene_prop_dict, gene_name_dict


def mean_distances(gene_dist_dict):
    lastgen, npool = get_gen_npool(isKeyword('rundir'))
    mean_dist_dict = dict()
    dist_dict = gene_dist_dict
    print(dist_dict)
    dist_sum = 0
    curr_gen = 0
    read_path = isKeyword('rundir') + "statespace/all_results.csv"
    with open(read_path, 'r') as fi:
        list_of_lines = fi.readlines()
        for line in list_of_lines:
            gen, gene, fitness, freq = line.split(",")
            gen = int(gen)
            freq = int(freq)
            if gen != curr_gen:
                mean_dist = dist_sum / npool
                mean_dist_dict.update({curr_gen: mean_dist})
                curr_gen += 1
                dist_sum = 0
            if type(dist_dict[gene]) == list:
                if 10000 in dist_dict[gene]:
                    npool -= 1 #subtract one, the 10000 does not count
                    continue
                else:
                    dist_sum += np.mean(dist_dict[gene])*int(freq) #else average the two distance metrics
            else:
                dist_sum += float(dist_dict[gene]) * int(freq)
        mean_dist = dist_sum / npool
        mean_dist_dict.update({curr_gen: mean_dist})
    fi.close()

    write_path = isKeyword('rundir') + "statespace/mean_distances.csv"
    with open(write_path, 'w') as fi:
        emsg = write_dictionary(mean_dist_dict, write_path)
        if emsg:
            print(emsg)


# Uses the same directory as get_general, which is get_run_dir() from ga_tools
def format_distances():
    gene_dist_dict, npool, _, _ = find_distances()
    print('DONE with find distances')
    mean_distances(gene_dist_dict)
    print('DONE with mean distances')
