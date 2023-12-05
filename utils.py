from config import *


def tuple_sum(tuple1, tuple2):

    if len(tuple1) != len(tuple2):
        raise ValueError("summing two tuples of different sizes")
    summation = []
    for i in range(len(tuple1)):
        summation.append(tuple1[i]+tuple2[i])

    return tuple(summation)

def pairwise_cost(seq1, seq2, delta):
        """ Compute pairwise cost of aligning two sequences. """
        cost = 0

        if len(seq1) != len(seq2):
            raise ValueError("two sequences are not of the same size")
        
        for i in range(len(seq1)):
            cost += delta[seq1[i]][seq2[i]]
        return cost

def SP_cost(alignment, delta):
    total_cost = 0
    breakpoint()
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
    if delta is None:
        delta = delta
    
    v_id, v = v_tuple

    w_id, w = w_tuple

    def reccurence(i,j):

        scores = []

    
        if i > 0:

            scores.append((M[i-1][j] + delta[v[i-1]]['-'], UP))

        if j > 0:
            scores.append((M[i][j-1] + delta['-'][w[j-1]], LEFT))

        if i > 0 and j > 0:
            scores.append((M[i-1][j-1] + delta[v[i-1]][w[j-1]], TOPLEFT))

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


