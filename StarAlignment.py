from Bio import pairwise2
from Bio.pairwise2 import format_alignment

class StarTreeAlignment:
    def __init__(self, sequences, delta):
    
        self.sequences = sequences
        self.delta = delta

        self.UP = (-1,0)
        self.LEFT = (0, -1)
        self.TOPLEFT = (-1, -1)
        self.ORIGIN = (0, 0)

        self.pariwise_alignments = self.generate_pairwise_alignments()

        self.MSA = None


    def generate_pairwise_alignments(self):

        pairwise_alignments = {}

        for i, seq_id1 in enumerate(self.sequences):
            pairwise_alignments[seq_id1] = {}

        for i, seq_id1 in enumerate(self.sequences):
            for j, seq_id2 in enumerate(self.sequences):
                if i < j:
                    pairwise_alignments[seq_id1][seq_id2] = self.global_align((seq_id1,self.sequences[seq_id1]), (seq_id2,self.sequences[seq_id2]))
                    pairwise_alignments[seq_id2][seq_id1] = pairwise_alignments[seq_id1][seq_id2]
        return pairwise_alignments


    def traceback_global(self, v_tuple, w_tuple, pointers):

        v_id, v = v_tuple

        w_id, w = w_tuple

        i,j = len(v), len(w)
        new_v = []
        new_w = []
        while True:
            di, dj = pointers[i][j]
            if (di,dj) == self.LEFT:
                new_v.append('-')
                new_w.append(w[j-1])
            elif (di,dj) == self.UP:
                new_v.append(v[i-1])
                new_w.append('-')
            elif (di,dj) == self.TOPLEFT:
                new_v.append(v[i-1])
                new_w.append(w[j-1])
            i, j = i + di, j + dj
            if (i <= 0 and j <= 0):
                break
        return {v_id: ''.join(new_v[::-1]),
            w_id: ''.join(new_w[::-1])}

    def global_align(self, v_tuple, w_tuple, delta = None):
        """
        Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as
        computed by traceback_global.

        :param: v
        :param: w
        :param: delta
        """
        if delta is None:
            delta = self.delta
        
        v_id, v = v_tuple

        w_id, w = w_tuple

        def reccurence(i,j):

            scores = []

        
            if i > 0:

                scores.append((M[i-1][j] + delta[v[i-1]]['-'], self.UP))

            if j > 0:
                scores.append((M[i][j-1] + delta['-'][w[j-1]], self.LEFT))

            if i > 0 and j > 0:
                scores.append((M[i-1][j-1] + delta[v[i-1]][w[j-1]], self.TOPLEFT))

            if i == 0 and j == 0:
                pass
            else:
                M[i][j], pointers[i][j] =  min(scores)



        M = [[0 for _ in range(len(w)+1)] for _ in range(len(v)+1)]
        pointers = [[self.ORIGIN for _ in range(len(w)+1)] for _ in range(len(v)+1)]
        score, alignment = None, None

        for i in range(len(v)+1):
            for j in range(len(w)+1):
                reccurence(i,j)

        score = M[len(v)][len(w)]
        alignment = self.traceback_global(v_tuple,w_tuple, pointers)

        return score, alignment





    def find_center_sequence(self):
        """
        Find center of a star that minimizes the distance 
        """
        min_distance = float('inf')
        center_sequence_id = None

        for seq_id in self.sequences:
            total_distance = sum(self.pariwise_alignments[seq_id][other_id][0] for other_id in self.pariwise_alignments[seq_id])
            
            if total_distance < min_distance:
                min_distance = total_distance
                center_sequence_id = seq_id
        print(min_distance)

        return center_sequence_id

    def update_MSA(curr_center_sequence, new_center_sequence):
        pass

    def align(self):
        """
        Perform the star alignment, aligning all sequences to the center sequence.
        """
        center_sequence_id = self.find_center_sequence()
        center_sequence = self.sequences[center_sequence_id]
        aligned_sequences = {center_sequence_id: center_sequence}

        for seq_tuple in self.sequences.items():
            if seq_tuple[0] != center_sequence_id:
                score, alignment = self.global_align((center_sequence_id, aligned_sequences[center_sequence_id]),seq_tuple)

                gapped_center_sequence = alignment[center_sequence_id]

                best_alignment = max(alignment, key=lambda x: x.score)
                aligned_sequences[seq_id] = best_alignment[1]  # Aligned sequence to the center


        self.MSA = aligned_sequences
        return aligned_sequences

if __name__ == "__main__":
    pass