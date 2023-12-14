from CarrilloLipman import CarrilloLipman

from utils import to_sequences, evaluate, SP_cost, read_MSA

from datetime import date, datetime
import argparse
import csv
import os


def delta(i , j):
    if i != j:
        return 1
    else:
        return 0

# Create the parser
parser = argparse.ArgumentParser(description='Evaluate Carrillo-Lipman performance.')

# Add an argument for the parameter file location
parser.add_argument('input_file', type=str, help='The location of the PFAM seed alignment file')

# Parse the arguments
args = parser.parse_args()


reference_alignment = read_MSA(args.input_file)

input_sequences = to_sequences(reference_alignment)

reference_SP_score = SP_cost(reference_alignment, delta)

aligner = CarrilloLipman(input_sequences, delta)



alignment, SP_score = aligner.align()


precision, recall, F_1_score = evaluate(alignment_ref= reference_alignment, alignment= alignment)


results_path = "results/results.csv"

file_exists = os.path.isfile(results_path)

pruned = 1 - aligner.num_visited/aligner.total_nodes

average_length = sum(aligner.lengths)/len(aligner.lengths)


columns = ['Alignment Name','Date & Time', 'reference SP cost', 'output SP cost', 'z', 'Precision','Recall','F-1 score','Pruning ratio','Total Nodes', 'Average length']
current_datetime = date.today().strftime("%m/%d") + " "+ datetime.now().strftime("%H:%M:%S")


alignment_name = args.input_file.split('/')[-1]

# Open the file in append mode
with open(results_path, mode='a', newline='') as csvfile:
    # Create a writer object
    writer = csv.DictWriter(csvfile, fieldnames=columns)

    # If the file does not exist, write the header
    if not file_exists:
        writer.writeheader()

    row_data = [alignment_name, current_datetime, reference_SP_score, SP_score, aligner.z, precision, recall, F_1_score,pruned, aligner.total_nodes, average_length]
    # Write the row of data
    writer.writerow(dict(zip(columns, row_data)))






