import csv
from molSimplifyAD.ga_tools import *

class gene:
    def __init__(self,name,fitness,frequency,generation):
        self.name = name
        self.fitness = '{0:.12g}'.format(float(fitness))
        self.frequency = frequency
        self.generation = generation
    def _long(self):
        return 'The gene {0} has a fitness {1} and appears {2} time(s).'.format(self.name, self.fitness, str(self.frequency))

#Find a gene by name in a given list.
def _find_gene(geneName, list):
    for i in xrange(len(list)):
        if list[i].name == geneName:
            return i
    return -1

#Return generation and npool from current_status.csv
def _get_gen_npool(data_dir):
    with open(data_dir+"statespace/current_status.csv",'r') as fi:
        reader = csv.reader(fi)
        for row in reader:
            if row[0] == "gen":
                lastgen = int(row[1].strip('\n'))
            elif row[0] == "npool":
                npool = int(row[1].strip('\n'))
    return lastgen,npool

#Change mode of writing for outputs.
def _write_mode(generation):
    if generation == 0:
        return 'w'
    else:
        return 'a'

#Write gene, fitness, and frequency to text file.
def _write_all_txt(base_path, generation, end_results):
    txt_results_path = base_path + "all_results.txt"
    with open(txt_results_path,_write_mode(generation)) as fo:
        out = "\nGen: " + str(generation) + " [gene, fitness, frequency]\n\n"
        fo.write(out)
        for i in xrange(len(end_results)):
            fo.write(end_results[i]._long()+"\n")
    fo.close()

#Write gene, fitness, and frequency to .csv file.
def _write_all_csv(base_path, generation, end_results):
    csv_results_path = base_path + "all_results.csv"
    with open(csv_results_path,_write_mode(generation)) as fo:
        writer = csv.writer(fo)
        for i in range(len(end_results)):
            t = end_results[i]
            writer.writerow( (t.generation,t.name,t.fitness,t.frequency) )
    fo.close()
def _write_split_csv(base_path, generation, split_results):
    csv_results_path = base_path + "_final_ANN_split_results.csv"
    with open(csv_results_path,_write_mode(generation)) as fo:
        writer = csv.writer(fo)
        writer.writerow((generation,split_results))
    fo.close()

#Write generation, mean fitness, and number of unique genes to .csv file.
def _gen_gene_fitness_csv(base_path, generation, end_results,sumt):
    lastgen,npool = _get_gen_npool(get_run_dir())
    mean_fitness = sumt/float(npool)
    csv_results_path = base_path + "_generation_meanFitness_diversity.csv"
    with open(csv_results_path,_write_mode(generation)) as fo:
        writer = csv.writer(fo)
        writer.writerow( (generation,mean_fitness,len(end_results)))
    fo.close()

def _human_readable_csv(base_path, generation, end_results):
    from molSimplifyAD.get_distances import _find_distances
    gene_dist_dict, _, gene_prop_dict, gene_name_dict = _find_distances()
    csv_results_path = base_path + "human_readable_results.csv"
    with open(csv_results_path,_write_mode(generation)) as fo:
        writer = csv.writer(fo)
        if int(generation) == 0:
            writer.writerow( ('Generation','Gene','Chem Name','Fitness','Property','Distance','Frequency') )
        else:
            writer.writerow(('\n'))
        for i in range(len(end_results)):
            t = end_results[i]
            writer.writerow( (t.generation,t.name,gene_name_dict[t.name],t.fitness,gene_prop_dict[t.name], gene_dist_dict[t.name],t.frequency) )
    fo.close()
##########################################################################################
#Find unique genes and their frequencies by name in current_genes and their fitness from gene_fitness. Output to a text file named results.txt
def _get_freq_fitness(lastgen, npool):
    ## if not DFT, get all current ANN splitting energies
    if not isDFT():
            full_gene_info = dict()
            ANN_prop_dict = dict()
            GA_run = get_current_GA()
            runtype = GA_run.config["runtype"]
            for generation in xrange(lastgen+1):
                ANN_dir = get_run_dir() + "ANN_ouput/gen_"+str(generation)+"/ANN_results.csv"
                emsg, ANN_dict = read_ANN_results_dictionary(ANN_dir)
                for keys in ANN_dict.keys():
                    this_gene = "_".join(keys.split("_")[4:10])
                    if runtype == "redox":
                        this_gene = "_".join(keys.split("_")[4:9])                    
                    this_prop = float(ANN_dict[keys][runtype])
                    this_dist = float(ANN_dict[keys][runtype+'_dist'])
                    if not(this_gene in ANN_prop_dict.keys()):
                        ANN_prop_dict.update({this_gene:this_prop})
    
    
    for generation in range(lastgen+1):
        end_results = []
        current_gene_list = list()

        #First find all unique genes and add to list from current_genes.csv.
        ## If gene is already present, increase its frequency by 1.
        base_path = get_run_dir()+"statespace/"
        read_path = base_path + "gen_"+str(generation)+"/gene_fitness.csv"
        fi = open(read_path,'r')
        print("opened Gen: " + str(generation))

        for line in fi:
            #print(fi.readline())
            geneName = line.split(",")[0]
            geneName = geneName.strip('\n')
            # print('GET GENERAL GENENAME:',geneName)
            current_gene_list.append(geneName)
            index = _find_gene(geneName, end_results)
            if index >= 0:
                end_results[index].frequency += 1
            else:
                end_results.append(gene(geneName, 0, 1,generation))

        fi.close()

        #Second, find fitness values of genes, add to list, and calculate mean fitness.
        sumt = 0

        read_path = base_path + "gen_"+str(generation)+"/gene_fitness.csv"
        with open(read_path, 'r') as fi:
            list_of_lines = fi.readlines()
            for line in list_of_lines:
                geneName, fitness = line.split(",")

                fitness = fitness.strip('\n')
                index = _find_gene(geneName, end_results)
                if index >= 0:
                    temp = end_results[index]
                    sumt += temp.frequency * float(fitness)
                    temp.fitness = format(float(fitness), '.12f')
                    print(temp._long())

        fi.close()

        #Third, output the unique genes and their fitness values to .txt and .csv files.
        _write_all_txt(base_path, generation, end_results)
        _write_all_csv(base_path, generation, end_results)
        _gen_gene_fitness_csv(base_path, generation, end_results,sumt)
        _human_readable_csv(base_path, generation, end_results)
        
        # Fourth, recover actual splitting energies only in ANN case
        if not isDFT():
            mean_split = 0.0 
            count = 0
            print('ANN keys are ')
            print(ANN_prop_dict.keys())
            for geneName in current_gene_list:
                
                if geneName in ANN_prop_dict.keys():
                    mean_split += abs(ANN_prop_dict[geneName])
                    count += 1
                else:
                    print('Error: expected ' + geneName + ' to be in ANN results...')
                    pass
            mean_split = mean_split/float(count)
            ## write
            _write_split_csv(base_path, generation, mean_split)  



def format_freqeuncies():
    lastgen,npool = _get_gen_npool(get_run_dir())
    _get_freq_fitness(lastgen, npool,)
