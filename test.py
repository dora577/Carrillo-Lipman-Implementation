from StarAlignment import StarTreeAlignment

import itertools
from utils import tuple_sum
alphabet = ['A', 'C', 'G', 'T', '-']

delta = {}
for i in range(len(alphabet)):
    delta[alphabet[i]] = {k : v for (k,v)
                          in zip(alphabet, [0 if alphabet[i] == alphabet[j]  else 1
                                  for j in range(len(alphabet))]
                         )}


sequences = {"v1":"TGGGAGCGA",
             "v2":"TGCCAGGGA",
             "v3":"TGCCGGA",
             "v4":"AGCCGGGAA"}

print(tuple_sum((1,2,3,5),(2,3,4,5,6)))


# ST = StarTreeAlignment(sequences, delta)

# print(ST.find_center_sequence())

# print(ST.align())


# for t in sequences.items():
#     print(t)






