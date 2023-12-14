from CarrilloLipman import CarrilloLipman

from utils import to_sequences, evaluate, SP_cost, read_MSA

from datetime import date, datetime
import argparse
import csv
import os





# Check if the file exists






def delta(i , j):
    if i != j:
        return 1
    else:
        return 0

# Create the parser
parser = argparse.ArgumentParser(description='Evaluate Carrillo-Lipman performance.')

# Add an argument for the parameter file location
parser.add_argument('seed_file', type=str, help='The location of the PFAM seed alignment file')

# Parse the arguments
args = parser.parse_args()


reference_alignment = read_MSA(args.seed_file)


input_sequences = to_sequences(reference_alignment)

reference_SP_score = SP_cost(reference_alignment, delta)

aligner = CarrilloLipman(input_sequences, delta)



alignment, SP_score = aligner.align()


precision, recall, F_1_score = evaluate(alignment_ref= reference_alignment, alignment= alignment)


results_path = "results/results.txt"

file_exists = os.path.isfile(results_path)

columns = ['Date & Time', 'reference SP cost', 'output SP cost', 'Precision','Recall','F-1 score']

current_datetime = date.today().strftime("%m/%d") + " "+ datetime.now().strftime("%H:%M:%S")



# Open the file in append mode
with open(results_path, mode='a', newline='') as csvfile:
    # Create a writer object
    writer = csv.DictWriter(csvfile, fieldnames=columns)

    # If the file does not exist, write the header
    if not file_exists:
        writer.writeheader()

    row_data = [current_datetime, reference_SP_score, SP_score, precision, recall, F_1_score]
    # Write the row of data
    writer.writerow(dict(zip(columns, row_data)))






