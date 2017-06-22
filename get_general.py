import csv
from ga_tools import get_run_dir

class gene:
    def __init__(self,name,fitness,frequency,generation):
        self.name = name
        self.fitness = '{0:.5g}'.format(float(fitness))
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
        for i in xrange(len(end_results)):
            t = end_results[i]
            writer.writerow( (t.generation,t.name,t.fitness,t.frequency) )
    fo.close()

def _gen_gene_fitness_csv(base_path, generation, end_results,sum):
    mean_fitness = sum/float(10)
    csv_results_path = base_path + "_generation_meanFitness_diversity.csv"
    with open(csv_results_path,_write_mode(generation)) as fo:
        writer = csv.writer(fo)
        writer.writerow( (generation,mean_fitness,len(end_results)))
    fo.close()

##########################################################################################
#Find unique genes and their frequencies by name in current_genes and their fitness from gene_fitness. Output to a text file named results.txt
def _get_freq_fitness(lastgen, npool):
    for generation in xrange(lastgen+1):
        end_results = []

        #First find all unique genes and add to list from current_genes.csv.
        ## If gene is already present, increase its frequency by 1.
        base_path = get_run_dir()+"statespace/"
        read_path = base_path + "gen_"+str(generation)+"/current_genes.csv"
        fi = open(read_path,'r')
        print "opened Gen: " + str(generation)

        for h in xrange(npool):
            n, geneName = fi.readline().split(",")
            geneName = geneName.strip('\n')
            index = _find_gene(geneName, end_results)
            if index >= 0:
                end_results[index].frequency += 1
            else:
                end_results.append(gene(geneName, 0, 1,generation))

        fi.close()

        #Second, find fitness values of genes, add to list, and calculate mean fitness.
        sum = 0
        read_path = base_path + "gen_"+str(generation)+"/gene_fitness.csv"
        fi = open(read_path,'r')

        while True:
            line = fi.readline()
            if len(line) < 1:
                break
            geneName, fitness = line.split(",")
            fitness = fitness.strip('\n')
            index = _find_gene(geneName, end_results)
            if index >= 0:
                temp = end_results[index]
                sum += temp.frequency * float(fitness)
                temp.fitness = format(float(fitness), '.5f')
                print temp._long()

        fi.close()

        #Third, output the unique genes and their fitness values to .txt or .csv file.
        _write_all_txt(base_path, generation, end_results)
        _write_all_csv(base_path, generation, end_results)
        _gen_gene_fitness_csv(base_path, generation, end_results,sum)

def _get_gen_npool():
    with open(get_run_dir()+"statespace/current_status.csv",'r') as fi:
        reader = csv.reader(fi)
        for row in reader:
            if row[0] == "gen":
                lastgen = int(row[1].strip('\n'))
            elif row[0] == "npool":
                npool = int(row[1].strip('\n'))
    return lastgen,npool

lastgen,npool = _get_gen_npool()
_get_freq_fitness(lastgen, npool)
###!!!!Would it be good to write a separate file to run these functions?
