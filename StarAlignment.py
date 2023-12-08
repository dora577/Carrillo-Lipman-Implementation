
from utils import global_align, INT_MAX, SP_cost

class StarTreeAlignment:
    def __init__(self, sequences, delta):
    
        self.sequences = sequences
        self.delta = delta
        self.pariwise_alignments = self.generate_pairwise_alignments()

        self.center_id = None
        self.alignment = None

    def generate_pairwise_alignments(self):

        pairwise_alignments = {}

        for i, seq_id1 in enumerate(self.sequences):
            pairwise_alignments[seq_id1] = {}

        for i, seq_id1 in enumerate(self.sequences):
            for j, seq_id2 in enumerate(self.sequences):
                if i < j:
                    pairwise_alignments[seq_id1][seq_id2] = global_align((seq_id1,self.sequences[seq_id1]), (seq_id2,self.sequences[seq_id2]),self.delta)
                    pairwise_alignments[seq_id2][seq_id1] = pairwise_alignments[seq_id1][seq_id2]
        return pairwise_alignments

    def find_center_sequence(self):
        """
        Find center of a star that minimizes the distance 
        """
        min_distance = INT_MAX   
        center_sequence_id = None

        for seq_id in self.sequences:
            total_distance = sum(self.pariwise_alignments[seq_id][other_id][0] for other_id in self.pariwise_alignments[seq_id])
            
            if total_distance < min_distance:
                min_distance = total_distance
                center_sequence_id = seq_id

        return center_sequence_id
  
    def find_new_gaps(self, original, new):
        # Initialize pointers for both sequences
        orig_pointer, new_pointer = 0, 0
        # List to hold the positions of new gaps
        new_gap_positions = []

        # Loop until the end of the new sequence
        while new_pointer < len(new):
            # If both pointers have a non-gap character, move both pointers
            if orig_pointer < len(original) and new_pointer < len(new) and original[orig_pointer] == new[new_pointer]:
                orig_pointer += 1
                new_pointer += 1
            # If the new has a gap, and orginal someything else move the original pointer
            elif new_pointer < len(new) and new[new_pointer] == '-':
                new_gap_positions.append(new_pointer)
                new_pointer += 1
            # If the original has a gap and new one something else OR they both have different non-gap letters then something is wrong
            else:
                raise ValueError

        return new_gap_positions

    def update_alignment(self, new_center_sequence):

        old_center_sequence = self.alignment[self.center_id]

        newGappedPositions =  self.find_new_gaps(old_center_sequence, new_center_sequence)

        for id in self.alignment.keys():

            for position in newGappedPositions:

                self.alignment[id] = self.alignment[id][:position] + "-" + self.alignment[id][position:]
        

    def align(self):

        """
        Perform the star alignment, aligning all sequences to the center sequence.
        """
        self.center_id = self.find_center_sequence()
        center_sequence = self.sequences[self.center_id]
        self.alignment = {self.center_id: center_sequence}
        for seq_tuple in self.sequences.items():
            if seq_tuple[0] != self.center_id:

                _, alignment = global_align((self.center_id, self.alignment[self.center_id]),seq_tuple, self.delta)

                new_center_sequence = alignment[self.center_id]

                self.update_alignment(new_center_sequence)

                self.alignment[seq_tuple[0]] = alignment[seq_tuple[0]]

        return self.alignment, SP_cost(self.alignment,self.delta)

if __name__ == "__main__":
    pass