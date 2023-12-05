import itertools
import sys
from queue import PriorityQueue

from StarAlignment import StarTreeAlignment

from utils import tuple_sum

INT_MAX = sys.maxsize

class CarrilloLipman:
    def __init__(self, sequences, delta, z):

        self. z = z 
        self.sequences = sequences
        self.delta = delta
        self.k = len(sequences.keys())
        self.id2key = {}
        self.key2id = {}

        self.lengths = []

        self.pQueue = PriorityQueue()


        for id, key in enumerate(sequences.keys()):

            self.id2key[id] =  key

            self.key2id[key] = id

            self.lengths.append(len(sequences[key]))


        self.alignment = {}

        
        self.edit_graph = {}

        self.init_edit_graph()

    
    def SP_cost(self, current_vertex, steps):

        gaps = 0

        total_cost = 0

        chars = []

        for idx ,step in enumerate(steps):

            if step == 0:
                gaps += 1
            elif step == 1:
                seq_pos = current_vertex[idx]+ 1

                new_char = self.sequences[self.id2key[idx]][seq_pos]

                if len(chars) > 0:
                    for char in chars:
                        total_cost += self.delta[char][new_char]

                chars.append(new_char)

        total_cost = gaps * delta['-'][new_char]

        return total_cost

    


    def find_neighbours(self, vertex):
        """
        Find all neighbors for a given vertex in the context of edit distance.
        
        :param vertex: A tuple representing the current position in each sequence.
        :return: A list of tuples representing neighboring vertices.
        """
        neighbours = []
        # Generate all combinations of sequences in which to move forward
        for r in range(1, self.k+1):
            for seq_indices in itertools.combinations(range(self.k), r):
                new_vertex = [0]*self.k
                valid_move = True

                for i in seq_indices:
                    if vertex[i] < self.lengths[i]:
                        new_vertex[i] += 1
                    else:
                        valid_move = False
                        break

                if valid_move:
                    neighbours.append(tuple(new_vertex))

        return neighbours

    def init_edit_graph(self):

        vertices = itertools.product(*(range(length + 1) for length in self.lengths))
        # print(list(vertices))
        for vertex in vertices:
        # Each vertex maps to a dictionary of its neighbors and their edge weights
        # Initially, neighbors can be empty or set to default values based on the problem
            if vertex == tuple([0]*self.k):

                self.edit_graph[vertex] = (0, False, None)
            else:
                self.edit_graph[vertex] = (INT_MAX, False, None)


    def shortest_path(self):

        self.pQueue.put((0, tuple([0]*self.k)))
        
        
        breakpoint()



if __name__ == "__main__":

    sequences = {"v1":"TGGGAGCGA",
             "v2":"TGCCAGGGA",
             "v3":"TGCCGGA",
             "v4":"AGCCGGGAA"}
    
    
    alphabet = ['A', 'C', 'G', 'T', '-']

    delta = {}
    for i in range(len(alphabet)):
        delta[alphabet[i]] = {k : v for (k,v)
                            in zip(alphabet, [0 if alphabet[i] == alphabet[j]  else 1
                                    for j in range(len(alphabet))]
                            )}

    aligner = CarrilloLipman(sequences, delta, z = 5)

    star = StarTreeAlignment(sequences, delta)

    a = star.align()
    print(aligner.lengths)

    # print(aligner.find_neighbours((9,9,6,6)))

    print(aligner.shortest_path())
    # print(SP_cost(a,delta))ß