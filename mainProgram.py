from auxiliary import weightVectors,neighbors,poblation,functionZDT3,zDT3
#First we define the parametrers that we can change:
##GENERAL PARAMETERS
N = 10 #Number of subproblems
G = 0 #Number of generations

##EVOLUTION OPERATORS

##OTHERS
T = 5 #Number of neighbors
SearchSpace = [0.0,0.0]

##MAIN
#1.-Create the weight vectors 
wei = weightVectors(N)
#2.-Identify the neighbors
nei = neighbors(wei,T)
#3.-Initialize the poblation and obtain the functions
p = poblation(30,N)
f = functionZDT3(p)
#4.-Initialize the z vector
z = zDT3([],f)
