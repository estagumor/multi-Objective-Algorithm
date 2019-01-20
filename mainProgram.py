from auxiliaryMethods import weightVectors, neighbors, poblation
from auxiliaryMethods import functionZDT3, updateZ, differentialEvolution
import matplotlib.pyplot as plt
#First we define the parametrers that we can change:
##GENERAL PARAMETERS
N = 40 #Number of subproblems
G = 250 #Number of generations

##EVOLUTION OPERATORS
F = 0.5 #Should be between [0,2]
GR = 0.65 #Should be between (0,1)
##OTHERS
T = 3 #Number of neighbors. Should be between the 10%/30% of the N size
SearchSpace = [0.0,1.0]

##MAIN
#1.-Create an array of individuals with each individual's weight inside 
array = weightVectors(N)
#2.-Adds to the array the index of the neighbors for each individual
array = neighbors(array, T)
#3.-Adds to the array the chromosome, f1 and f2 fitness for each individual 
array = poblation(array, 30, SearchSpace[0], SearchSpace[1])
array = functionZDT3(array)
#4.-Initialize the z vector
z = updateZ([],array)

IT = 0 #Loop's iterator
f1Vector = [] #We're gonna save the f1 of every individual and every generation
f2Vector = [] #We're gonna save the f2 of every individual and every generation
writeFile = open('outFilesZDT3/40x250.out', 'wb') #To save the results (metrics)
while IT < G: #when IT is equal to G, the algorithm will stop
    i = 0 #Loop's iterator
    while i < len(array):
        ret = differentialEvolution(array[i], F, GR, z, array, SearchSpace[0], SearchSpace[1])
        array = ret[0]
        z = ret[1]
        i = i + 1
        for r in ret[0]: #Saving the f1 and f2 functions of each element in the poblation
            f1Vector.append(r.f1)
            f2Vector.append(r.f2)

    for r in ret[0]:
        #Saving the data in a file
        escrib = str(str(r.f1) + "\t" + str(r.f2) + "\t" + str(0.0) + "\n").encode()
        writeFile.write(escrib)
    IT = IT + 1

writeFile.close()
print("Iteration: " + str(IT) + ". Actual z: " + str(ret[1]))

##PLOTTING
f1Ideal= []
f2Ideal= []
f = open("outFilesZDT3Ideal/40x250.out", "r", -1, "UTF-8") #NSGAII results 
for line in f.readlines():
    split = line.split("\t")
    f1Ideal.append(float(split[0]))
    split1 = split[1].split("\n")
    f2Ideal.append(float(split1[0]))
plt.plot(f1Vector, f2Vector, 'ro', f1Ideal, f2Ideal, 'bs', z[0], z[1], 'g^')
plt.ylabel('f1/f2 comparation')
plt.show()
