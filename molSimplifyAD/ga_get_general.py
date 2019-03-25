import csv
from molSimplifyAD.ga_tools import *


class gene:
    def __init__(self, name, fitness, frequency, generation, fit_func, dist_param):
        self.name = name
        self.fitness = '{0:.12g}'.format(float(fitness))
        self.frequency = frequency
        self.generation = generation
        self.dist_param = str(dist_param)
        self.fitness_function = str(fit_func)

    def _long(self):
        return 'The gene {0} has a fitness {1} and appears {2} time(s).'.format(self.name, self.fitness,str(self.frequency))


# Find a gene by name in a given list.
def _find_gene(geneName, list):
    for i in xrange(len(list)):
        if list[i].name == geneName:
            return i
    return -1


# Return generation and npool from current_status.csv
def _get_gen_npool(data_dir):
    with open(data_dir + "statespace/current_status.csv", 'r') as fi:
        reader = csv.reader(fi)
        for row in reader:
            if row[0] == "gen":
                lastgen = int(row[1].strip('\n'))
            elif row[0] == "npool":
                npool = int(row[1].strip('\n'))
        print('lastgen is ', lastgen)
    return lastgen, npool


# Change mode of writing for outputs.
def _write_mode(generation):
    if generation == 0:
        return 'w'
    else:
        return 'a'


# Write gene, fitness, and frequency to text file.
def _write_all_txt(base_path, generation, end_results):
    txt_results_path = base_path + "all_results.txt"
    with open(txt_results_path, _write_mode(generation)) as fo:
        out = "\nGen: " + str(generation) + " [gene, fitness, frequency]\n\n"
        fo.write(out)
        for i in xrange(len(end_results)):
            fo.write(end_results[i]._long() + "\n")
    fo.close()


# Write gene, fitness, and frequency to .csv file.
def _write_all_csv(base_path, generation, end_results):
    csv_results_path = base_path + "all_results.csv"
    with open(csv_results_path, _write_mode(generation)) as fo:
        writer = csv.writer(fo)
        for i in range(len(end_results)):
            t = end_results[i]
            writer.writerow((t.generation, t.name, t.fitness, t.frequency))
    fo.close()


def _write_prop_csv(base_path, generation, prop_results):
    csv_results_path = base_path + "_final_ANN_prop_results.csv"
    with open(csv_results_path, _write_mode(generation)) as fo:
        writer = csv.writer(fo)
        writer.writerow((generation, prop_results))
    fo.close()


# Write generation, mean fitness, and number of unique genes to .csv file.
def _gen_gene_fitness_csv(base_path, generation, end_results, sumt):
    lastgen, npool = _get_gen_npool(isKeyword('rundir'))
    mean_fitness = sumt / float(npool)
    csv_results_path = base_path + "_generation_meanFitness_diversity.csv"
    with open(csv_results_path, _write_mode(generation)) as fo:
        writer = csv.writer(fo)
        writer.writerow((generation, mean_fitness, len(end_results)))
    fo.close()


def _human_readable_csv(base_path, generation, end_results):
    print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ HR CSV $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    from molSimplifyAD.get_distances import _find_distances
    from molSimplifyAD.ga_main import *
    gene_dist_dict, _, gene_prop_dict, gene_name_dict = _find_distances()
    csv_results_path = base_path + "human_readable_results.csv"
    print('AT HRCSV!!!!')
    with open(csv_results_path, _write_mode(generation)) as fo:
        writer = csv.writer(fo)
        # print('didnt fail here')
        if int(generation) == 0:
            writer.writerow(('Generation', 'Gene', 'Chem Name', 'Fitness', 'Property', 'Distance', 'Frequency','Fitness Function','Dist Parameter'))
        else:
            writer.writerow(('\n'))
            # print('didnt fail here2')
        for i, val in enumerate(end_results):
            # t = end_results[i]
            # print('didnt fail here3')
            # print(i, len(end_results))
            # print(gene_name_dict)
            # print([x.fitness for x in end_results])
            t = val
            # print(t.generation) #CR V NOT GETTING ASSIGNED TO THE GENE NAME DICT--> problematic!!!
            # print(t.name)
            # print(t.fitness_function)
            # print(t.name in gene_name_dict.keys())
            # print(gene_name_dict[t.name])
            # print(t.fitness)
            # print(gene_prop_dict[t.name])
            # print(gene_dist_dict[t.name])
            # print(t.frequency)
            # print(t.dist_param)
            # print('-!-!-!-!-!-!-!-')
            writer.writerow((t.generation, t.name, gene_name_dict[t.name], t.fitness, gene_prop_dict[t.name],
                             gene_dist_dict[t.name], t.frequency, t.fitness_function, t.dist_param))
    # print('hereeee')
    fo.close()


def _write_summary_csv(base_path, generation, sum_results):
    from molSimplifyAD.get_distances import _find_distances
    gene_dist_dict, _, gene_prop_dict, gene_name_dict = _find_distances()
    csv_results_path = base_path + "results_summary.csv"
    with open(csv_results_path, _write_mode(generation)) as fo:
        writer = csv.writer(fo)
        if int(generation) == 0:
            writer.writerow(('Generation', 'Gene', 'Chem Name', 'Fitness', 'Property', 'Distance'))
        else:
            writer.writerow(('\n'))
        for i in range(len(sum_results)):
            t = sum_results[i]
            writer.writerow((t.generation, t.name, gene_name_dict[t.name], t.fitness, gene_prop_dict[t.name],
                             gene_dist_dict[t.name]))
    fo.close()


##########################################################################################
# Find unique genes and their frequencies by name in current_genes and their fitness from gene_fitness. Output to a text file named results.txt
def _get_freq_fitness(lastgen, npool):
    ## if not DFT, get all current ANN splitting energies
    if not isKeyword('DFT'):
        full_gene_info = dict()
        ANN_prop_dict = dict()
        #GA_run = get_current_GA()
        runtype = isKeyword("runtype")
        for generation in xrange(lastgen + 1):
            print('------GENERATION IN XRANGE IS '+str(generation))
            ANN_dir = isKeyword('rundir') + "ANN_ouput/gen_" + str(generation) + "/ANN_results.csv"
            emsg, ANN_dict = read_ANN_results_dictionary(ANN_dir)
            for keys in ANN_dict.keys():
                # print('got here')
                print(ANN_prop_dict.keys())
                translate_dict = translate_job_name(keys)
                this_gene = translate_dict['gene']
                metal = translate_dict['metal']
                ox = translate_dict['ox']
                spin = translate_dict['spin']
                spin_cat = translate_dict['spin_cat']
                # this_gene = "_".join(keys.split("_")[4:10])
                # if runtype == "redox":
                #     this_gene = "_".join(keys.split("_")[4:9])
                # elif runtype in ['oxo','hat']:
                #     this_gene = "_".join(keys.split("_")[4:10])
                split_energy = ANN_dict[keys]['split']
                if runtype in ['split']:
                    this_prop = float(ANN_dict[keys][runtype])
                    this_dist = float(ANN_dict[keys][runtype + '_dist'])
                    ANN_prop_dict.update({this_gene: this_prop})
                elif runtype in ['oxo','hat']:
                    # print('got here 3')HR CSVHR CSV
                    if not (this_gene in ANN_prop_dict.keys()):
                        # print('got here 2')
                        this_prop = float(ANN_dict[keys][runtype])
                        this_dist = float(ANN_dict[keys][runtype + '_dist'])
                        # if spin_cat == 'HS' or (metal == 'cr' and int(spin) == 2):
                        if spin_cat == isKeyword('spin_constraint'):
                            # print('updated ANN propdict')
                            ANN_prop_dict.update({this_gene: this_prop})
                        elif isKeyword('spin_constraint') and get_metals()[metal].lower() == 'cr' and ox == 5:
                            print('Cr(V) does not exist in HS')
                            ANN_prop_dict.update({this_gene: 10000})
                        else:
                            # ANN_prop_dict.update({this_gene: 10000})
                            print('SKIPPED ANN PROP DICT ASSIGNMENT BC '+str(this_gene)+' '+str(spin))
                elif runtype in ['homo','gap']:
                    if (split_energy > 0 and int(spin) <= 3) or (split_energy < 0 and int(spin) > 3):
                        this_prop = float(ANN_dict[keys][runtype])
                        this_dist = float(ANN_dict[keys][runtype + '_dist'])
                        ANN_prop_dict.update({this_gene: this_prop})
                elif type(runtype) == list:
                    this_prop = []
                    this_dist = []
                    for run in runtype:
                        this_prop.append(float(ANN_dict[keys][run]))
                        this_dist.append(float(ANN_dict[keys][run + '_dist']))
                    if spin_cat == isKeyword('spin_constraint'): #Constraining this to a single spin state.
                        print('freq fitness multiple prop update')
                        ANN_prop_dict.update({this_gene: this_prop})
                    elif isKeyword('spin_constraint') and get_metals()[metal].lower() == 'cr' and ox == 5:
                        print('Cr(V) does not exist in HS')
                        ANN_prop_dict.update({this_gene: [10000, 10000]})
                    else:
                        # ANN_prop_dict.update({this_gene: 10000})
                        print('SKIPPED ANN PROP DICT ASSIGNMENT BC '+str(this_gene)+' IS SPIN '+str(spin))
                # elif not (this_gene in ANN_prop_dict.keys()):
                #     sardines
                #     this_prop = float(ANN_dict[keys][runtype])
                #     this_dist = float(ANN_dict[keys][runtype + '_dist'])
                #     ANN_prop_dict.update({this_gene: this_prop})
                ############################################################################################
                # Note: in the way this loop is currently set up, it loops over everything in ANN_dict.    #
                # Thus, if the ANN prediction predicts something for multiple spin states, you need to     #
                # Make sure that the property is assigned for the right spin state (the ground state)      #
                # See logic in the function _find_distances() in get_distances. In order for the property  #
                # to be assigned, certain criteria must be met as can be seen above. For spin splitting    #
                # Not an issue for spin splitting because both spins inherit the same splitting energy     #
                ############################################################################################
        print('NOW PRINTING ANN PROP DICT.........')
        print(ANN_prop_dict)

    generation = 0
    print('------------!!!!! ASSIGNED GENERATION '+str(generation)+' !!!!!-------------')
    while (generation == 0 or generation < lastgen):
        print('@@@@@@@@@@@@@@@@@@@@@@@@ GENERATION @@@@@@@@@@@@@@@@@@@@@@@@@@: '+str(generation))
        end_results = []
        current_gene_list = list()
        # print('_____________________________________________________________________Entered part 1.')
        # First find all unique genes and add to list from current_genes.csv.
        ## If gene is already present, increase its frequency by 1.
        base_path = isKeyword('rundir') + "statespace/"
        read_path = base_path + "gen_" + str(generation) + "/_current_genes.csv"
        fi = open(read_path, 'r')
        print("opened Gen: " + str(generation))
        for line in fi:
            # print(fi.readline())
            geneName = line.split(",")[1]
            if len(line.split(",")) == 5:
                dist_param = (",".join(line.split(",")[2:4])).strip('\n')
                fit_func = line.split(",")[4].strip('\n')
            else:
                dist_param = line.split(",")[2].strip('\n')
                fit_func = line.split(",")[3].strip('\n')
            geneName = geneName.strip('\n')
            # print('GET GENERAL GENENAME:',geneName)
            current_gene_list.append(geneName)
            index = _find_gene(geneName, end_results)
            if index >= 0:
                end_results[index].frequency += 1
            else:
                end_results.append(gene(geneName, 0, 1, generation, fit_func, dist_param))
        print('printing end results')
        print([val.name for val in end_results], len([val.name for val in end_results]))
        fi.close()
        print('got to this point 1')
        # Second, find fitness values of genes, add to list, and calculate mean fitness.
        # print('_____________________________________________________________________Entered part 2.')
        sumt = 0
        read_path = base_path + "gen_" + str(generation) + "/gene_fitness.csv"
        with open(read_path, 'r') as fi:
            list_of_lines = fi.readlines()
            for line in list_of_lines:
                # geneName, fitness = line.split(",",1)
                geneName, _ = line.split(",",1)
                prop = ANN_prop_dict[geneName]
                if 'hinge' in str(isKeyword('scoring_function')):
                    fitness = find_prop_hinge_fitness(prop, isKeyword('property_parameter'))
                else:
                    fitness = find_prop_fitness(prop, isKeyword('property_parameter'))

                # fitness = fitness.strip('\n')
                # print('FITNESS ASSIGNED', fitness)
                # print(end_results)
                index = _find_gene(geneName, end_results)
                if index >= 0:
                    temp = end_results[index]
                    sumt += temp.frequency * float(fitness)
                    temp.fitness = format(float(fitness), '.12f')
                    print(temp._long())
        print('got to this point 2')
        ## Write a summary file for the property, distance, and fitness.
        # print('_____________________________________________________________________Entered part 3.')
        sum_results = []
        read_path = base_path + "gen_" + str(generation) + "/gene_fitness.csv"
        fi = open(read_path, 'r')
        for line in fi:
            geneName = line.split(",")[0]
            geneName = geneName.strip('\n')
            index = _find_gene(geneName, sum_results)
            if index >= 0:
                sum_results[index].frequency += 1
            else:
                sum_results.append(gene(geneName, 0, 1, generation, None, None))
        fi.close()
        # Third, output the unique genes and their fitness values to .txt and .csv files.
        _write_all_txt(base_path, generation, end_results)
        _write_all_csv(base_path, generation, end_results)
        _gen_gene_fitness_csv(base_path, generation, end_results, sumt)
        _human_readable_csv(base_path, generation, end_results)
        print('here')
        _write_summary_csv(base_path, generation, sum_results)
        
        # Fourth, recover actual splitting energies only in ANN case
        if not isKeyword('DFT'):
            if type(runtype) == list:
                mean_prop = [0.0, 0.0]
            else:
                mean_prop = 0.0
            count = 0
            print('ANN keys are ')
            print(ANN_prop_dict.keys())
            for geneName in current_gene_list:
                if geneName in ANN_prop_dict.keys():
                    if ANN_prop_dict[geneName] == 10000:
                        continue
                    else:
                        if runtype == 'split':
                            mean_prop += abs(ANN_prop_dict[geneName])
                            count += 1
                        elif type(runtype) == list:
                            mean_prop[0] += ANN_prop_dict[geneName][0]
                            mean_prop[1] += ANN_prop_dict[geneName][1]
                        else:
                            mean_prop += ANN_prop_dict[geneName]
                            count += 1
                else:
                    print('Error: expected ' + geneName + ' to be in ANN results...')
                    pass
            if count > 0:
                if type(runtype) == list:
                    mean_prop[0] = mean_prop[0]/float(count)
                    mean_prop[1] = mean_prop[1]/float(count)
                else:
                    mean_prop = mean_prop / float(count)
            ## write
            _write_prop_csv(base_path, generation, mean_prop)
        generation += 1
        print('Moving to generation: ',generation)

def format_freqeuncies():
    lastgen, npool = _get_gen_npool(isKeyword('rundir'))
    print('lastgen', lastgen)
    _get_freq_fitness(lastgen, npool)
