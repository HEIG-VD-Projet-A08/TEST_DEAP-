#!/usr/bin/python3
#    This file is part of DEAP.
#
#    DEAP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    DEAP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with DEAP. If not, see <http://www.gnu.org/licenses/>.


#    example which maximizes the sum of a list of integers
#    each of which can be 0 or 1

import random
from random import randint

from deap import base
from deap import creator
from deap import tools
import os

creator.create("FitnessMax", base.Fitness, weights=(0.5,1,10))
creator.create("Individual", list, fitness=creator.FitnessMax)

# CXPB  is the probability with which two individuals
#       are crossed
#
# MUTPB is the probability for mutating an individual
CXPB, MUTPB = 0.4, 0.3
proteinAlphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
target = "MLMPKKNRIAIHELLFKEGVMVAKKDVHMPKHPELADKNVPNLHVMKAMQSLKSRGCVKEQFAWRHFYWYLTNEGSQYLR"
minSizeWord = 10
maxSizeWord = 20
individualSize = 10
populationSize = 100
eValue=10
bestResult=0
g = 0
def generateText(intervalMin, intervalMax):
    size = randint(intervalMin, intervalMax)
    string = ""
    for i in range(0, size ):
        string += proteinAlphabet[randint(0, len(proteinAlphabet) - 1)]
    return string


def compareString(individual, target):
    counter = 0
    length =0
    for i in range(0, len(individual)):
        length += (len(individual[i])-minSizeWord)/(maxSizeWord-minSizeWord)
    return counter,


def switchChar(individual, indpb):
    """Shuffle the attributes of the input individual and return the mutant.
    The *individual* is expected to be a :term:`sequence`. The *indpb* argument is the
    probability of each attribute to be moved. Usually this mutation is applied on
    vector of indices.

    :param individual: Individual to be mutated.
    :param indpb: Independent probability for each attribute to be exchanged to
                  another position.
    :returns: A tuple of one individual.

    This function uses the :func:`~random.random` and :func:`~random.randint`
    functions from the python base :mod:`random` module.
    """
    size = len(individual)
    for i in range(size):
        if random.random() < indpb:
            swap_indx = random.randint(0, size - 2)
            if swap_indx >= i:
                swap_indx += 1
            tmp = individual[i]
            individual[i] = individual[swap_indx]
            individual[swap_indx]=tmp

    return individual


def mutateText(individual, indpb):
    for i in range(len(individual)):
        if random.random() < indpb:
            individual[i] = generateText(minSizeWord,maxSizeWord)
    return individual,

"""
    for i in range(len(individual)):
        list1 = list(individual[i])
        for j in range(len(individual[i])):
            if random.random() < indpb:
                list1[j] = proteinAlphabet[randint(0, len(proteinAlphabet) - 1)]
        individual[i] = ''.join(list1)
    individual=switchChar(individual,indpb)

    return individual,
"""


def crossoverChar(ind1, ind2,indpb=1):
    """Executes a two-point crossover on the input :term:`sequence`
    individuals. The two individuals are modified in place and both keep
    their original length.

    :param ind1: The first individual participating in the crossover.
    :param ind2: The second individual participating in the crossover.
    :returns: A tuple of two individuals.

    This function uses the :func:`~random.randint` function from the Python
    base :mod:`random` module.
    """
    size = min(len(ind1), len(ind2))
    cxpoint1 = random.randint(1, size)
    cxpoint2 = random.randint(1, size - 1)
    if cxpoint2 >= cxpoint1:
        cxpoint2 += 1
    else:  # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1
    for i in range(0,len(ind1)):
        if random.random() < indpb:

            tmp=ind2[i][:cxpoint1]
            tmp+=ind1[i][cxpoint1:cxpoint2]
            tmp += ind2[i][cxpoint2:]
            tmp2 = ind1[i][:cxpoint1]
            tmp2 += ind2[i][cxpoint1:cxpoint2]
            tmp2 += ind1[i][cxpoint2:]
            ind1[i]=tmp
            ind2[i]=tmp2

    return ind1, ind2
toolbox = base.Toolbox()

# Attribute generator
#                      define 'attr_bool' to be an attribute ('gene')
#                      which corresponds to integers sampled uniformly
#                      from the range [0,1] (i.e. 0 or 1 with equal
#                      probability)
toolbox.register("generateRandomTexte", generateText, minSizeWord, maxSizeWord)

# Structure initializers
#                         define 'individual' to be an individual
#                         consisting of 100 'attr_bool' elements ('genes')
toolbox.register("individual", tools.initRepeat, creator.Individual,
                 toolbox.generateRandomTexte, individualSize)

# define the population to be a list of individuals
toolbox.register("population", tools.initRepeat, list, toolbox.individual)


# the goal ('fitness') function to be maximized
def evalOneMax(individuals):
    global g, bestResult

    fastaFile = open("../ncbi-blast-2.10.0+/bin/toBlast.fasta","w")
    individualList =iter(individuals)
    differentElem=[]
    differentElem.append(next(individualList))
    fastaFile.write(">" + differentElem[-1] + "\n" + differentElem[-1] + "\n")
    fastaFile.close()
    for individual in individualList:
        tmpFile = open("../ncbi-blast-2.10.0+/bin/tmp.fasta","w")
        tmpFile.write(">" + individual + "\n" +individual + "\n")
        tmpFile.close()
        os.system("../ncbi-blast-2.10.0+/bin/blastp  -subject ../ncbi-blast-2.10.0+/bin/tmp.fasta -query ../ncbi-blast-2.10.0+/bin/toBlast.fasta -out resultSim -evalue 0.1 -word_size 2 -gapopen 10 -gapextend 1  -task blastp -outfmt '6 saccver  qacc ppos qseq sseq qcovs' -subject_besthit")
        result = open("resultSim", "r")
        m=0
        for line in result:
            elements = line.split('\t')
            m+=1
            if float(elements[2]) < 90:
                fastaFile = open("../ncbi-blast-2.10.0+/bin/toBlast.fasta", "a")
                fastaFile.write(">" + individual + "\n" + individual + "\n")
                fastaFile.close()
        result.close()
        if m ==0:
            fastaFile = open("../ncbi-blast-2.10.0+/bin/toBlast.fasta", "a")
            fastaFile.write(">" + individual + "\n" + individual + "\n")
            fastaFile.close()
    os.system("rm ../ncbi-blast-2.10.0+/bin/tmp.fasta")

    global eValue
    os.system("../ncbi-blast-2.10.0+/bin/blastp -db ../ncbi-blast-2.10.0+/bin/in.fasta -query ../ncbi-blast-2.10.0+/bin/toBlast.fasta -out result -evalue "+ str(eValue)+" -word_size 2 -gapopen 10 -gapextend 1 -num_threads 4 -task blastp -outfmt '6 saccver  qacc ppos qseq sseq qcovs' -subject_besthit")


    sum = 0
    nbrElem=0
    average=0
    nbrDifferentString=[]

    if os.path.getsize("result")!=0:
        fastaFile =open("result","r")

        for line in fastaFile:
            elements = line.split('\t')
            sum += float(elements[2])
            nbrElem +=1
            if elements[1] not in nbrDifferentString:
                nbrDifferentString.append(elements[1])
        fastaFile.close()
        average = sum/nbrElem
    if sum > bestResult:
        os.system("cp result bestRatio/result"+str(g))
        bestResult = sum
    return sum,average,len(nbrDifferentString)


# ----------
# Operator registration
# ----------
# register the goal / fitness function
toolbox.register("evaluate", evalOneMax)

# register the crossover operator
toolbox.register("mate", tools.cxTwoPoint)

# register a mutation operator with a probability to
# flip each attribute/gene of 0.05
toolbox.register("mutate", mutateText, indpb=0.05)

# operator for selecting individuals for breeding the next
# generation: each individual of the current generation
# is replaced by the 'fittest' (best) of three individuals
# drawn randomly from the current generation.
toolbox.register("select", tools.selTournament, tournsize=3)


# ----------

def main():
    #random.seed(64)

    # create an initial population of 300 individuals (where
    # each individual is a list of integers)

    print("Start of evolution")
    pop = toolbox.population(n=populationSize)
    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate,pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    print("  Evaluated %i individuals" % len(pop))

    # Extracting all the fitnesses of
    fits = [ind.fitness.values[0] for ind in pop]

    # Variable keeping track of the number of generations

    global g
    # Begin the evolution
    while  g < 100:
        global eValue
        global bestResult
        bestResult =0
        if eValue >= 0.001:
            eValue = eValue / 2

        #eValue = eValue /1.2
        # A new generation
        g = g + 1
        print("-- Generation %i --" % g)

        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))

        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):

            # cross two individuals with probability CXPB
            if random.random() < CXPB:
                toolbox.mate(child1, child2)

                # fitness values of the children
                # must be recalculated later
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:

            # mutate an individual with probability MUTPB
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
            del mutant.fitness.values

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        print("  Evaluated %i individuals" % len(invalid_ind))

        # The population is entirely replaced by the offspring
        pop[:] = offspring

        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        ind1 =pop[0]
        for ind in pop:
            if ind1.fitness.values[0] < ind.fitness.values[0]:
                ind1 = ind

        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x * x for x in fits)
        std = abs(sum2 / length - mean ** 2) ** 0.5

        print("  Min %s" % min(fits))
        print("  ratio %s" % ind1.fitness.values[1])
        print("  Max %s" % max(fits))
        print("  Avg %s" % mean)
        print("  Std %s" % std)
        print("evalue %0.10f" %eValue)
        if max(fits)< 500:
            eValue= eValue*2
    print("-- End of (successful) evolution --")

    best_ind = tools.selBest(pop, 1)[0]
    worst_ind = tools.selWorst(pop,1)[0]
    print("Best individual is %s, %s\n" % (best_ind, best_ind.fitness.values))
    print("worst individual is %s, %s" % (worst_ind, worst_ind.fitness.values))


if __name__ == "__main__":
    main()
