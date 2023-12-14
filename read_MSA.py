from Bio import AlignIO



# Path to your alignment file
file_path = 'data/PF00005.alignment.seed'

# Read the alignment
alignment = AlignIO.read(file_path, "stockholm")

# Search for an entry with an unknown amino acid ('X')
unknown_entry = None
for record in alignment:
    if 'X' in record.seq:
        unknown_entry = record
        break

if unknown_entry:
    print("Found an entry with unknown amino acid:", unknown_entry.id)
else:
    print("No entry with unknown amino acid found")


print("Alignment length %i" % alignment.get_alignment_length())

all_leters = set()
lengths = set()
for record in alignment:
    lengths.add(len(str(record.seq)))
    for letter in str(record.seq):
        all_leters.add(letter)

    # print(record.seq)# + " id: " + record.id)

print(lengths)
print(all_leters)

original_sequences, lengths = transformToInput(alignment)


print(lengths)