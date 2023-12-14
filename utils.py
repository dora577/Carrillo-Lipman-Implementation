from config import *


from Bio import AlignIO


def read_MSA(file_path):

    reference_alignment_dict = {}
    reference_alignment = AlignIO.read(file_path, "stockholm")

    for record in reference_alignment:
        reference_alignment_dict[record.id] = str(record.seq) 

    return reference_alignment_dict


def to_sequences(alignment):

    original_sequences = {}

    for seq_id, sequence in alignment.items():
        original_sequences[seq_id] = sequence.replace('-', '')
    return original_sequences




def generate_pairs(alignment):
    sequence_ids = sorted(list(alignment.keys()))    
    pairs = set()
    for i, seq_id1 in enumerate(sequence_ids):
        for j, seq_id2  in enumerate(sequence_ids):
            if i < j:
                for idx1, l_1 in enumerate(alignment[seq_id1]):
                    for idx2, l_2 in enumerate(alignment[seq_id2]):
                        if l_1 == l_2:
                            pairs.add((seq_id1, idx1, seq_id2,idx2))
    return pairs

def get_precision(pairs_ref, pairs):

    TP = pairs.intersection(pairs_ref)

    return len(TP)/len(pairs)
    
def get_recall(pairs_ref, pairs):

    TP = pairs.intersection(pairs_ref)

    FN = pairs_ref.difference(TP)

    return len(TP)/(len(TP)+ len(FN))

def get_F_score(precision, recall, b = 1):

    return (1 + b**2) * (precision * recall)/(b**2 * precision + recall) 


def evaluate(alignment_ref, alignment):

    pairs_ref = generate_pairs(alignment_ref)

    pairs = generate_pairs(alignment)

    precision = get_precision(pairs_ref, pairs)

    recall = get_recall(pairs_ref, pairs)

    F_1_score = get_F_score(precision, recall)


    return precision, recall, F_1_score


def tuple_sum(tuple1, tuple2):

    if len(tuple1) != len(tuple2):
        raise ValueError("summing two tuples of different sizes")
    summation = []
    for i in range(len(tuple1)):
        summation.append(tuple1[i]+tuple2[i])

    return tuple(summation)

def tuple_diff(tuple1, tuple2):
    if len(tuple1) != len(tuple2):
        raise ValueError("summing two tuples of different sizes")
    difference = []
    for i in range(len(tuple1)):
        difference.append(tuple1[i]-tuple2[i])

    return tuple(difference)

def pairwise_cost(seq1, seq2, delta):
        """ Compute pairwise cost of aligning two sequences. """
        cost = 0

        if len(seq1) != len(seq2):
            raise ValueError("two sequences are not of the same size")
        
        for i in range(len(seq1)):
            cost += delta(seq1[i], seq2[i])
        return cost

def SP_cost(alignment, delta):
    total_cost = 0
    for i, key_i in enumerate(alignment.keys()):
        for j, key_j in enumerate(alignment.keys()):
            if i < j:
                total_cost+= pairwise_cost(alignment[key_i], alignment[key_j],delta)

    return total_cost

def traceback_global(v_tuple, w_tuple, pointers):

    v_id, v = v_tuple

    w_id, w = w_tuple

    i,j = len(v), len(w)
    new_v = []
    new_w = []
    while True:
        di, dj = pointers[i][j]
        if (di,dj) == LEFT:
            new_v.append('-')
            new_w.append(w[j-1])
        elif (di,dj) == UP:
            new_v.append(v[i-1])
            new_w.append('-')
        elif (di,dj) == TOPLEFT:
            new_v.append(v[i-1])
            new_w.append(w[j-1])
        i, j = i + di, j + dj
        if (i <= 0 and j <= 0):
            break
    return {v_id: ''.join(new_v[::-1]),
        w_id: ''.join(new_w[::-1])}

def global_align(v_tuple, w_tuple, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as
    computed by traceback_global.

    :param: v
    :param: w
    :param: delta
    """
    v_id, v = v_tuple

    w_id, w = w_tuple

    def reccurence(i,j):

        scores = []

    
        if i > 0:

            scores.append((M[i-1][j] + delta(v[i-1], '-'), UP))

        if j > 0:
            scores.append((M[i][j-1] + delta('-', w[j-1]), LEFT))

        if i > 0 and j > 0:
            scores.append((M[i-1][j-1] + delta(v[i-1], w[j-1]), TOPLEFT))

        if i == 0 and j == 0:
            pass
        else:
            M[i][j], pointers[i][j] =  min(scores)



    M = [[0 for _ in range(len(w)+1)] for _ in range(len(v)+1)]
    pointers = [[ORIGIN for _ in range(len(w)+1)] for _ in range(len(v)+1)]
    score, alignment = None, None

    for i in range(len(v)+1):
        for j in range(len(w)+1):
            reccurence(i,j)

    score = M[len(v)][len(w)]
    alignment = traceback_global(v_tuple,w_tuple, pointers)

    return score, alignment


