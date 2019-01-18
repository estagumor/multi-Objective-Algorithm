import random
from math import sqrt,sin,pi
from classes import Individual

def weightVectors(N): #Obtains the weight vector's
    #First, we are gonna do the weights. We need N vectors of two elements 
    vectors = [] #The final return
    poblation = []
    equidistance = 1/N #The distance between the first element of the vectors
    j = 0 #Loop's initialization

    while j < N: #A loop for the vectors

        if(len(vectors) < 1): #First iteration
            weights = [equidistance]
        else:
            weights = [(len(vectors)+1)*equidistance] 

        weights.append(1 - weights[0]) #The y element of the weigths vector 
        ind = Individual(weights)
        vectors.append(weights)
        poblation.append(ind)
        j = j + 1

    return poblation

def neighbors(vectors, T): #Weight vector's, Number of neighbors 
    neighbors = [] #Array of neighbors for every vector. Return of the function
    ind = [] #index of the neighbors for every vector.
    i = 0 #Iterator
    while i < len(vectors): #Going through all the vectors
        cVector = vectors[i] #Current vector
        distance = {} #Distance's dictionary for the current vector
        for vector in vectors: 
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
        ind.append(keys)
        l = 0 #Loop's iterator
        n = []
        while l < T: #We only choose the T better vectors of the array 
            for k in keys:
                #TE HAS QUEDADO POR AQUI
                n.append()
                l = l + 1
        neighbors.append(n)
        i = i + 1 #Next iteration

    return [ind,neighbors]

w = weightVectors(5)
neighbors(w,2)

def poblation(el,N,weights,neighbors): #Obtains a random poblation of N size with el elements
    i = 0 #Loop's initialization
    poblationVector = []
    while i < N: #It generates a element between the range and adds it to the return vector
        j = 0 #Loop's initialization
        ind = []
        while j < el:
            ind.append(random.uniform(0,1)) 
            j = j + 1
        individual = Individual(ind)

        poblationVector.append(individual)
        i = i + 1 #Next iteration

    return poblationVector

def obtainPoblationNeighbors(neighbors,poblation):
    poblationNeighbors = []
    ind = neighbors[0]
    for i in ind:
        n = []
        for j in i:
            n.append(poblation[j])
        poblationNeighbors.append(n)

    return poblationNeighbors

def functionZDT3(poblation): #Obtains the poblation's fitness and returns it in an array
    i = 0 #Loop's initialization
    f1Vector = [] 
    f2Vector = []
    while i < len(poblation): 
        ind = poblation[i] #Current poblation element
        f1 = ind[0]
        j = 2 #Loop's initialization
        sum = 0 #Summation
        while j < len(ind):
            sum = sum + ind[j]
            j = j + 1
        g = 1 + (9/(len(ind)-1))*sum
        h = 1 - sqrt(f1/g) - (f1/g)*sin(10*pi*f1)
        f2 = g*h
        f1Vector.append(f1)
        f2Vector.append(f2)
        i = i + 1

    return [f1Vector,f2Vector]

def zDT3(z,functions):
    f1Vector = functions[0]
    f1Vector.sort()
    f2Vector = functions[1]
    f2Vector.sort()
    f1 = f1Vector[0]
    f2 = f2Vector[0]
    if(len(z) < 1): #initialization
        z.append(f1)
        z.append(f2)
    else:
        if(z[0] > f1):
            z[0] = f1
        
        if(z[1] > f2):
            z[1] = f2
    return z

def differentialEvolution(pneighbors,neighbors,functions,F,GR):
    #F -> Mutation rate [0,2]
    #GR -> Recombination rate (0,1)

    #n -> tamaño de cada individuo de la poblacion
    #m -> iterador dentro de cada individuo de la poblacion
    #NP -> numero de individuos que componen la poblacion
    #p -> iterador de individuos de la poblacion
    #g -> generacion
    noisyRandomVectors = []
    trialVectors = []
    selectionVectors = []
    poblation = pneighbors
    NP = len(poblation)
    
    #Mutation
    i = 0 #Loop's iterator
    while i < NP:
        #Target vectors
        xa = random.choice(poblation)
        xb = random.choice(poblation)
        while xb == xa: 
            xb = random.choice(poblation)
        xc = random.choice(poblation)
        while xc == xb or xc == xa:
            xc = random.choice(poblation)
        
        j = 0 #Loop's iterator
        ngp = []
        while j < len(xa):
            op = xc[j] + F*(xa[j]- xb[j])
            ngp.append(op)
            j = j + 1
        noisyRandomVectors.append(ngp)
        i = i + 1

    #Recombination
    k = 0 #Loop's iterator
    while k < NP:
        j = 0 #Loop's iterator
        tgp = []
        while j < len(xa):
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