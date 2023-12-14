import itertools
import sys
from queue import PriorityQueue

from StarAlignment import StarTreeAlignment


import time



from utils import tuple_sum, tuple_diff, global_align

from config import INT_MAX
from tqdm import tqdm

class CarrilloLipman:
    def __init__(self, sequences, delta):

        centerStar = StarTreeAlignment(sequences, delta)
        start_time = time.time()

        _, self.z = centerStar.align()
    
        end_time = time.time()
        print(f"Execution time Star Alignment: {round(end_time - start_time,2)} seconds")


        breakpoint()
        self.sequences = sequences
        self.delta = delta
        self.k = len(sequences.keys())

        self.id2key = {}
        self.key2id = {}
        self.lengths = []
        self.source = tuple([0]*self.k)
        self.pQueue = PriorityQueue()
        self.pQueue.put((0, self.source))

        self.num_visited = 0

        for id, key in enumerate(sequences.keys()):

            self.id2key[id] =  key

            self.key2id[key] = id

            self.lengths.append(len(sequences[key]))

        self.total_nodes  = 1
        for length in self.lengths:
            self.total_nodes *= (length + 1)


        self.alignment = {key : "" for key in self.sequences}

        
        self.edit_graph = {}


        # print("Starting initializing graph")
        
        # start_time = time.time()
        # self.init_edit_graph()
        # end_time = time.time()

        # print("Finished initializing graph")
        # print(f"Execution time: {end_time - start_time} seconds")

    def SP_step_cost(self, current_vertex, steps):

        gaps = 0

        total_cost = 0

        chars = []

        for idx ,step in enumerate(steps):

            if step == 0:
                gaps += 1
            elif step == 1:
                seq_pos = current_vertex[idx]

                new_char = self.sequences[self.id2key[idx]][seq_pos]

                if len(chars) > 0:
                    for char in chars:
                        total_cost += self.delta(char,new_char)

                chars.append(new_char)

        total_cost += gaps * (self.k - gaps) * delta('-',new_char)

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

    # def init_edit_graph(self):

    #     vertices = itertools.product(*(range(length + 1) for length in self.lengths))
    #     # print(list(vertices))
    #     for vertex in vertices:
    #     # Each vertex maps to a dictionary of its neighbors and their edge weights
    #     # Initially, neighbors can be empty or set to default values based on the problem
    #         if vertex == self.source:

    #             self.edit_graph[vertex] = (0, False, None) # (priority, visited, predecessor)
    #         else:
    #             self.edit_graph[vertex] = (INT_MAX, False, None)

    def optimal_cost(self, current_vertex):

        suffix_cost = 0
        
        for idx1, pos1 in enumerate(current_vertex):
            for idx2, pos2 in enumerate(current_vertex):

                if idx1 < idx2:
                    true_pos1 = pos1 
                    true_pos2 = pos2 

                    # print((self.id2key[idx1], self.sequences[self.id2key[idx1]][true_pos1:]))
                    # print((self.id2key[idx2], self.sequences[self.id2key[idx2]][true_pos2:]))
                    pairwise_optimal_cost, _ = global_align((self.id2key[idx1], self.sequences[self.id2key[idx1]][true_pos1:]),(self.id2key[idx2], self.sequences[self.id2key[idx2]][true_pos2:]), self.delta)
                    suffix_cost += pairwise_optimal_cost
        optimal_cost = self.edit_graph[current_vertex][0] + suffix_cost
        # breakpoint()
        return optimal_cost
        
    def shortest_path(self):

        self.num_visited = 1

        with tqdm() as iteration:

            while not self.pQueue.empty():
                
                current_priority,current_vertex  = self.pQueue.get()

                neighbours_step= self.find_neighbours(current_vertex)

                
                self.edit_graph[current_vertex]  (0, True, None)

                if any(index < 1 for index in current_vertex) or (self.optimal_cost(current_vertex) <= self.z):

                    for neighbour_step in neighbours_step:
        
                        new_priority = current_priority + self.SP_step_cost(current_vertex, neighbour_step)

                        neighbour_vertex = tuple_sum(current_vertex, neighbour_step)

                        if neighbour_vertex not in self.edit_graph.keys():
                           
                            self.edit_graph[neighbour_vertex] = (new_priority, False, current_vertex)                            
                            self.pQueue.put((new_priority,neighbour_vertex)) 
                            self.num_visited += 1

                        elif new_priority < self.edit_graph[neighbour_vertex][0]:


                            self.edit_graph[neighbour_vertex] = (new_priority, self.edit_graph[neighbour_vertex][1], current_vertex)

                            if not self.edit_graph[neighbour_vertex][1]:
                                self.pQueue.put((self.edit_graph[neighbour_vertex][0],neighbour_vertex))
                                self.num_visited += 1

                iteration.update()
                    
    def back_trace(self):

        current_vertex = tuple(self.lengths)

        while True:

            predeccesor_vertex = self.edit_graph[current_vertex][2]

            if predeccesor_vertex is None:
                break

            step = tuple_diff(current_vertex, predeccesor_vertex)

            for idx, s in enumerate(step):
                if s == 0:
                    self.alignment[self.id2key[idx]] = '-' + self.alignment[self.id2key[idx]]
                elif s == 1:
                    self.alignment[self.id2key[idx]] = self.sequences[self.id2key[idx]][predeccesor_vertex[idx]] + self.alignment[self.id2key[idx]]

            current_vertex = predeccesor_vertex

        return self.alignment 

    def align(self):
        self.shortest_path()

        return self.back_trace(), self.edit_graph[tuple(self.lengths)][0]


    def print_shortest_path(self):

        last_vertex = tuple(self.lengths)

        current_vertex = last_vertex

        print(current_vertex, self.edit_graph[current_vertex])
        while True:

            current_vertex = self.edit_graph[current_vertex][2]

            if current_vertex is None:
                break
            
            print(current_vertex,self.edit_graph[current_vertex])


if __name__ == "__main__":

    sequences = {"v1":"TGGGAGCGA",
             "v2":"TGCCAGGGA",
             "v3":"TGCCGGA",
              "v4":"AGCCGGGAA" }
    # ,
    #          "v4":"AGCCGGGAA"}
    
    
    alphabet = ['A', 'C', 'G', 'T', '-']

    def delta(i , j):
        if i != j:
            return 1
        else:
            return 0

    aligner = CarrilloLipman(sequences, delta)
  
    # print(aligner.find_neighbours((9,9,6,6)))

    print(aligner.align())




    pruned = 1 - aligner.num_visited/aligner.total_nodes

    print(pruned)


    print(aligner.z)

    # for current_vertex, value in aligner.print_shortest_path():
    #     print(current_vertex, value)
   