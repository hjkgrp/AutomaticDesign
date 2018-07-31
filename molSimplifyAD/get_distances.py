# Looks into the ANN_results.csv files
##Return file with list of genes and their ANN distances.
##Returns files with genes, their fitnesses, and their frequences per generation.
import csv
import os
from molSimplifyAD.ga_tools import *
from molSimplifyAD.ga_tools import get_run_dir
from molSimplifyAD.ga_get_general import _get_gen_npool


##Finds the ANN distances, writes them to a dictionary and then to a file
def _find_distances():
    lastgen, npool = _get_gen_npool(get_run_dir())
    gene_dist_dict = dict()
    gene_prop_dict = dict()
    gene_name_dict = dict()
    # lastgen is the last generation that has been run
    for generation in range(lastgen + 1):
        ANN_path = get_run_dir() + "ANN_ouput/gen_" + str(generation) + "/ANN_results.csv"

        if os.path.isfile(ANN_path):
            with open(ANN_path, 'r') as fi:
                # lol stands for list_of_lines
                lol = fi.readlines()
                GA_run = get_current_GA()
                runtype = GA_run.config["runtype"]
                for i, line in enumerate(lol):
                    if i == 0:
                        line_list = line.strip('\n').split(",")
                        prop_idx = line_list.index(runtype)
                        dist_idx = line_list.index(runtype + '_dist')
                        spin_idx = line_list.index('split')
                        continue
                    line_list = line.split(",")
                    job = line_list[0]
                    # print('job:', job)
                    gene, _, _, metal, ox, eqlig, axlig1, axlig2, _, _, _, spin, _, ahf, _, _ = translate_job_name(job)
                    split_energy = float(line_list[spin_idx])
                    spinmult = int(job.split('_')[-1])
                    if (split_energy > 0 and spinmult <= 3) or (split_energy < 0 and spinmult > 3):
                        prop = float(line_list[prop_idx])
                        dist = float(line_list[dist_idx])
                        # print(job, prop, dist, split_energy, spinmult)
                        ll = job.split("_")

                        geneName = "_".join(ll[4:10])
                        metal = get_metals()[metal]
                        chem_name = '_'.join(
                            [str(metal), str(ox), 'eq', str(eqlig), 'ax1', str(axlig1), 'ax2', str(axlig2), str(ahf),
                             str(spin)])
                        if geneName in gene_dist_dict.keys():
                            pass
                        else:
                            gene_dist_dict.update({geneName: dist})
                            gene_prop_dict.update({geneName: prop})
                            gene_name_dict.update({geneName: chem_name})
            fi.close()

    ## Writes genes and distances to a .csv file
    write_path = get_run_dir() + "statespace/all_distances.csv"
    if not os.path.isfile(write_path):
        open(write_path, 'w').close()
    emsg = write_dictionary(gene_dist_dict, write_path)
    if emsg:
        print(emsg)

    return gene_dist_dict, npool, gene_prop_dict, gene_name_dict


def _mean_distances(gene_dist_dict):
    lastgen, npool = _get_gen_npool(get_run_dir())
    mean_dist_dict = dict()
    dist_dict = gene_dist_dict
    dist_sum = 0
    curr_gen = 0
    read_path = get_run_dir() + "statespace/all_results.csv"
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
            dist_sum += float(dist_dict[gene]) * int(freq)
        mean_dist = dist_sum / npool
        mean_dist_dict.update({curr_gen: mean_dist})
    fi.close()

    write_path = get_run_dir() + "statespace/_mean_distances.csv"
    with open(write_path, 'w') as fi:
        emsg = write_dictionary(mean_dist_dict, write_path)
        if emsg:
            print(emsg)


# Uses the same directory as get_general, which is get_run_dir() from ga_tools
def format_distances():
    gene_dist_dict, npool, _, _ = _find_distances()
    _mean_distances(gene_dist_dict)
