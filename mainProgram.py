from auxiliaryMethods import weightVectors,neighbors,poblation,functionZDT3,zDT3,differentialEvolution
#First we define the parametrers that we can change:
##GENERAL PARAMETERS
N = 10 #Number of subproblems
G = 0 #Number of generations

##EVOLUTION OPERATORS
F = 0.5 #Should be between [0,2]
GR = 0.5 #Should be between (0,1)
##OTHERS
T = 5 #Number of neighbors. Should be between the 10%/30% of the N size
SearchSpace = [0.0,1.0]

##MAIN
#1.-Create an array of individuals with each individual's weight inside 
array = weightVectors(N)
#2.-Adds to the array the index of the neighbors for each individual
array = neighbors(array, T)
#3.-Adds to the array the chromosome, f1 and f2 fitness for each individual 
array = poblation(array, 30, N, SearchSpace[0], SearchSpace[1])
array = functionZDT3(array)
#4.-Initialize the z vector
z = zDT3([],array)

IT = 0 #Loop's iterator
while IT < G: #when IT is equal to G, the algorithm will stop
    for individual in array:
        #ind, F, GR, z
        differentialEvolution(individual, F, GR, z)

