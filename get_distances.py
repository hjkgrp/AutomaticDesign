#Looks into the ANN_results.csv files
##Return file with list of genes and their ANN distances.
##Returns files with genes, their fitnesses, and their frequences per generation.
import csv
import os
from ga_tools import *
from get_general import _get_gen_npool,data_dir

##Finds the ANN distances, writes them to a dictionary and then to a file
def _find_distances():
    lastgen,npool = _get_gen_npool(data_dir)
    gene_dist_dict = dict()

    #lastgen is the last generation that has been run
    for generation in xrange(lastgen+1):
        ANN_path = data_dir + "ANN_ouput/gen_" + str(generation) + "/ANN_results.csv"

        if os.path.isfile(ANN_path):
            with open(ANN_path, 'r') as fi:
                #lol stands for list_of_lines
                lol = fi.readlines()

                for line in lol:
                    gene, energy, distance = line.split(",")
                    ll = gene.split("_")
                    geneName = "_".join(ll[4:9])
                    distance = float(distance.strip('\n'))
                    if geneName in gene_dist_dict.keys():
                        pass
                    else:
                        gene_dist_dict.update({geneName:distance})
                        print geneName + " : " + str(distance)

            fi.close()

    ## Writes genes and distances to a .csv file
    write_path = data_dir + "statespace/all_distances.csv"
    if not os.path.isfile(write_path):
            open(write_path,'w').close()
    emsg = write_dictionary(gene_dist_dict,write_path)
    if emsg:
            print(emsg)

    return gene_dist_dict,npool

def _mean_distances(gene_dist_dict):
    lastgen,npool = _get_gen_npool(data_dir)
    mean_dist_dict = dict()
    dist_dict = gene_dist_dict
    dist_sum = 0
    curr_gen = 0

    read_path = data_dir + "statespace/all_results.csv"
    with open(read_path, 'r') as fi:
        list_of_lines = fi.readlines()
        for line in list_of_lines:
            gen, gene, fitness, freq = line.split(",")
            gen = int(gen)
            freq = int(freq)
            if gen != curr_gen:
                mean_dist = dist_sum / npool
                mean_dist_dict.update({curr_gen:mean_dist})
                curr_gen += 1
                dist_sum = 0
            dist_sum += float(dist_dict[gene])*int(freq)
        mean_dist = dist_sum / npool
        mean_dist_dict.update({curr_gen:mean_dist})
    fi.close()

    write_path = data_dir + "statespace/_mean_distances.csv"
    with open(write_path, 'w') as fi:
        emsg = write_dictionary(mean_dist_dict,write_path)
        if emsg:
                print(emsg)

# Uses the same directory as get_general, which is get_run_dir() from ga_tools
gene_dist_dict,npool = _find_distances()
_mean_distances(gene_dist_dict)
