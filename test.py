from StarAlignment import StarTreeAlignment

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



ST = StarTreeAlignment(sequences, delta)

print(ST.find_center_sequence())


for t in sequences.items():
    print(t)