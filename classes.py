class Individual:
    next_id = 0
    """Individual abstraction"""
    def __init__(self, weights):
        self.weights = weights #An array with the weights of the individual
        self.chromosome = []
        self.neighbors = []
        self.f1 = 0.0
        self.f2 = 0.0
        self.id = Individual.next_id #We'll need it in the future 
        Individual.next_id += 1
        print('Individual: ' + str(self.id) + '. Weights: ' + str(weights))

    def add_neighbor(self, n): #An array with the index of the neighbors 
        self.neighbors.append(n)
    
    def add_chromosome(self, c): #An array with the chromosome of the individual
        self.chromosome = c

    def add_f1(self, f1): #A float that contains the fitness of the first objective function
        self.f1 = f1

    def add_f2(self, f2): #A float that contains the fitness of the second objective function
        self.f2 = f2