class Individual:
    next_id = 0
    """Individual abstraction"""
    def __init__(self, weights):
        self.weights = weights
        self.chromosome = []
        self.neighbors = []
        self.id = Individual.next_id
        Individual.next_id += 1
        print('Individual: ' + str(self.id) + '. Weights: ' + str(weights))

    def add_neighbor(self, n):
        self.neighbors.append(n)