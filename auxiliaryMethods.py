import random
from math import sqrt,sin,pi
from classes import Individual

def weightVectors(N): #Obtains the weight vector's
    #First, we are gonna do the weights. We need N vectors of two elements 
    poblation = [] #Return an array of individuals with their weights 
    equidistance = 1/N #The distance between the first element of the vectors
    j = 0 #Loop's initialization

    while j < N: #A loop for the vectors
        if(len(poblation) < 1): #First iteration
            weights = [0.0]
        else:
            weights = [(len(poblation))*equidistance] 

        weights.append(1 - weights[0]) #The y element of the weigths vector 
        ind = Individual(weights)
        poblation.append(ind)
        j = j + 1

    return poblation

def neighbors(individuals, T): #Individual vector's with their weights, Number of neighbors
    i = 0 #Iterator
    while i < len(individuals): #Going through all the individuals
        cVector = individuals[i] #Current vector
        distance = {} #Distance's dictionary for the current vector
        k = 0 #We are gonna use it to know the position of the neighbor in the array
        for vector in individuals: 
            if vector.id == cVector.id: 
                distance[k] = 0.0
            else:
                sum = 0 
                j = 0 #Iterator
                while j < len(vector.weights): #That's the inside of the sqrt
                    sum = sum + abs(cVector.weights[j] - vector.weights[j])**2
                    j = j + 1
                squareRoot = sqrt(sum)
                distance[k] = squareRoot
            k = k + 1
        
        keys = sorted(distance, key=distance.__getitem__) #Returns the keys ordered by the values
        l = 0 #Loop's iterator
        for k in keys:
            if(l >= T):
                break
            else:
                cVector.add_neighbor(k) #We add the id instead of the individual cause it could change
            l = l + 1
        i = i + 1 #Next iteration

    return individuals #With the neighbors inside

def poblation(individuals, el, min, max): #Obtains a random poblation of N size with el elements
    i = 0 #Loop's initialization
    while i < len(individuals): #It generates a element between the range and adds it to the return vector
        j = 0 #Loop's initialization
        ind = []
        while j < el:
            ind.append(random.uniform(min,max)) 
            j = j + 1
        individuals[i].add_chromosome(ind)
        i = i + 1 #Next iteration

    return individuals

def functionZDT3(individuals): #Obtains the poblation's fitness and returns it in an array
    i = 0 #Loop's initialization
    try:
        while i < len(individuals): #We have an array
            ind = individuals[i] #Current poblation element
            c = ind.chromosome
            f1 = c[0]
            j = 2 #Loop's initialization
            sum = 0 #Summation
            while j < len(c):
                sum = sum + c[j]
                j = j + 1

            g = 1 + ((9*sum)/(len(c)-1))
            h = 1 - sqrt(f1/g) - (f1/g)*sin(10*pi*f1)
            f2 = g*h
            ind.add_f1(f1)
            ind.add_f2(f2)
            i = i + 1
    except TypeError as identifier: #We have an individual
        ind = individuals
        c = ind.chromosome
        f1 = c[0]
        j = 2 #Loop's initialization
        sum = 0 #Summation
        while j < len(c):
            sum = sum + c[j]
            j = j + 1
        g = 1 + ((9*sum)/(len(c)-1))
        h = 1 - sqrt(f1/g) - (f1/g)*sin(10*pi*f1)
        f2 = g*h
        ind.add_f1(f1)
        ind.add_f2(f2)

    return individuals

def updateZ(z,individuals):
    if(len(z) < 1): #Initialization
        f1Vector = []
        f2Vector = []
        for ind in individuals:
            f1Vector.append(ind.f1)
            f2Vector.append(ind.f2)
        f1Vector.sort() #We make an ascend sorting
        z.append(f1Vector[0]) #Took the best f1 return 
        f2Vector.sort()
        z.append(f2Vector[0])
    else: #Evaluating f1 and f2 of the new Individual
        if(individuals.f1 < z[0]):
            z[0] = individuals.f1
        if(individuals.f2 < z[1]):
            z[1] = individuals.f2

    return z

def differentialEvolution(ind, F, GR, z, individuals, minValue, maxValue):
    #F -> Mutation rate [0,2]
    #GR -> Recombination rate (0,1)

    poblationIndex = ind.neighbors #Index of the neighbors in the individuals vector 
    poblation = []
    for p in poblationIndex:
        poblation.append(individuals[p])
    NP = len(poblation)

    #Mutation
    #Target vectors
    xa = random.choice(poblation)
    xb = random.choice(poblation)
    while xb == xa: 
        xb = random.choice(poblation)
    xc = random.choice(poblation)
    while xc == xb or xc == xa:
        xc = random.choice(poblation)

    i = 0 #Loop's iterator
    ngp = []
    while i < len(xa.chromosome):
        op = xc.chromosome[i] + F*(xa.chromosome[i]- xb.chromosome[i])
        if(op < minValue or op > maxValue): #We don't want negative values 
            op = random.uniform(minValue, maxValue)
        ngp.append(op) 
        i = i + 1

    #Recombination
    j = 0 #Loop's iterator
    tgp = []
    k = random.randint(0,len(poblation)-1) #We choose a neighbor randomly
    xgp = poblation[k]
    while j < len(xa.chromosome):
        rand = random.uniform(0.01, 1)
        if(rand > GR):
            tgp.append(xgp.chromosome[j])
        else:
            tgp.append(ngp[j])
        j = j + 1

    #Selection  
    l = 0 #Loop's iterator
    while l < len(poblation): #It is the moment to update the neighbors
        #We have tgp (the chromosome)
        cNeighbor = poblation[l]
        cWeights = cNeighbor.weights
        newInd = Individual(cWeights)
        newInd.add_chromosome(tgp)
        newInd.neighbors = cNeighbor.neighbors
        newInd = functionZDT3(newInd) #Obtains the fitness of the new element
        z = updateZ(z, newInd) #Update the z vector

        absxf1 = cNeighbor.f1 - z[0]
        xg1 = cWeights[0]*abs(absxf1)
        absxf2 = cNeighbor.f2 - z[1]
        xg2 = cWeights[1]*abs(absxf2)
        gtex = max(xg1, xg2) #Getting the Tchebychef function of the neighbor

        absyf1 = newInd.f1 - z[0] 
        yg1 = cWeights[0]*abs(absyf1)
        absyf2 = newInd.f2 - z[1]
        yg2 = cWeights[1]*abs(absyf2)
        gtey = max(yg1, yg2) #Getting the Tchebychef function of the newInd

        if(gtey <= gtex): 
            index = poblationIndex[l]
            individuals[index] = newInd
        l = l + 1

    return [individuals, z] #We can use it in the future 
