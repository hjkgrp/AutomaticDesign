#Looks into the ANN_results.csv files
##Return list of genes and their ANN distances
import csv
from ga_tools import *

#Need number of generations (not sure about npool yet)
def _get_gen_npool():
    with open(get_run_dir()+"statespace/current_status.csv",'r') as fi:
        reader = csv.reader(fi)
        for row in reader:
            if row[0] == "gen":
                lastgen = int(row[1].strip('\n'))
            elif row[0] == "npool":
                npool = int(row[1].strip('\n'))
    return lastgen,npool

##Find the ANN distances, write to dictionary and then file
def _find_distances():
    lastgen,npool = _get_gen_npool()
    gene_dist_dict = dict()

    #lastgen is the last generation that has been run
    for generation in xrange(lastgen+1):
        ANN_path = get_run_dir() + "ANN_ouput/gen_" + str(generation) + "/ANN_results.csv"

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

    ## Write genes and distances to a .csv file
    write_path = get_run_dir() + "statespace/distances.csv"
    if not os.path.isfile(write_path):
            open(write_path,'w').close()
    emsg = write_dictionary(gene_dist_dict,write_path)
    if emsg:
            print(emsg)

_find_distances()
