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
            weights = [equidistance]
        else:
            weights = [(len(poblation)+1)*equidistance] 

        weights.append(1 - weights[0]) #The y element of the weigths vector 
        ind = Individual(weights)
        poblation.append(ind)
        j = j + 1

    return poblation

def neighbors(individuals, T): #Indivual vector's with their weights, Number of neighbors
    i = 0 #Iterator
    while i < len(individuals): #Going through all the individuals
        cVector = individuals[i] #Current vector
        distance = {} #Distance's dictionary for the current vector
        for vector in individuals: 
            if vector.id == cVector.id: 
                distance[vector.id] = 0.0
            else:
                sum = 0 
                j = 0 #Iterator
                while j < len(vector.weights): #That's the inside of the sqrt
                    sum = sum + abs(cVector.weights[j] - vector.weights[j])**2
                    j = j + 1
                squareRoot = sqrt(sum)
                distance[vector.id] = squareRoot
        
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

def poblation(individuals, el, N, min, max): #Obtains a random poblation of N size with el elements
    i = 0 #Loop's initialization
    while i < N: #It generates a element between the range and adds it to the return vector
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
    while i < len(individuals): 
        ind = individuals[i] #Current poblation element
        c = ind.chromosome
        f1 = c[0]
        j = 2 #Loop's initialization
        sum = 0 #Summation
        while j < len(c):
            sum = sum + c[j]
            j = j + 1
        g = 1 + (9/(len(c)-1))*sum
        h = 1 - sqrt(f1/g) - (f1/g)*sin(10*pi*f1)
        f2 = g*h
        ind.f1 = f1
        ind.f2 = f2
        i = i + 1

    return individuals

def zDT3(z,individuals):
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
    else: #By now, I suppose individuals is only a one individual 
        if(individuals.f1 < z[0]):
            z[0] = individuals.f1
        if(individuals.f2 < z[0]):
            z[1] = individuals.f2

    return z

ind = weightVectors(5)
ind = neighbors(ind,2)
ind = poblation(ind,10,5,0.0,1.0)
ind = functionZDT3(ind)
z = zDT3([],ind)
#ALL SEEMS TO WORK 

def differentialEvolution(individuals, ind, F, GR):
    #F -> Mutation rate [0,2]
    #GR -> Recombination rate (0,1)

    #n -> tamaño de cada individuo de la poblacion
    #m -> iterador dentro de cada individuo de la poblacion
    #NP -> numero de individuos que componen la poblacion
    #p -> iterador de individuos de la poblacion
    #g -> generacion

    ##HAY QUE OBTENER A TODOS LOS VECINOS PARA PODER HACER ESTO. A VER COMO LO PLANTEAMOS

    noisyRandomVectors = []
    trialVectors = []
    selectionVectors = []
    poblationIndex = ind.neighbors
    NP = len(poblationIndex)
    
    #Mutation
    i = 0 #Loop's iterator
    while i < NP:
        #Target vectors
        xaI = random.choice(poblationIndex)
        xbI = random.choice(poblationIndex)
        while xbI == xaI: 
            xbI = random.choice(poblationIndex)
        xcI = random.choice(poblationIndex)
        while xcI == xbI or xcI == xaI:
            xcI = random.choice(poblationIndex)
        
        xa = filter(lambda x: (x.id == xaI), individuals)
        xb = filter(lambda y: (y.id == xbI), individuals)
        xc = filter(lambda z: (z.id == xcI), individuals)

        j = 0 #Loop's iterator
        ngp = []
        while j < len(xa.chromosome):
            op = xc.chromosome[j] + F*(xa.chromosome[j]- xb.chromosome[j])
            ngp.append(op)
            j = j + 1
        noisyRandomVectors.append(ngp)
        i = i + 1

    #Recombination
    k = 0 #Loop's iterator
    while k < NP:
        j = 0 #Loop's iterator
        tgp = []
        while j < len(xa.chromosome):
            rand = random.uniform(0.01, 1)
            if(rand > GR):
                xgp = poblation[k]
                tgp.append(xgp[j])
            else:
                ngp = noisyRandomVectors[k]
                tgp.append(ngp[j])
            j = j + 1
        trialVectors.append(tgp)
        k = k + 1
    
    #Selection
    f = functionZDT3(trialVectors)
    l = 0 #Loop's iterator
    while l < len(f):
        tfit = f[l]
        xfit = functions[l]
        if(zt[0] < zx[0] or [1] < xfit[1]):
            selectionVectors.append(trialVectors[l])
        else:
            selectionVectors.append(poblation[l])
        l = l + 1

    return [selectionVectors,f] #The new poblation and their fitness 

w = weightVectors(3)
n = neighbors(w,2)
p = poblation(3,3)
pn = obtainPoblationNeighbors(n,p)
f = functionZDT3(pn)

#pneighbors,neighbors,functions,F,GR
differentialEvolution(pn,n,f,0.3,0.3)