# Carrillo-Lipman Python Implementation & Center Star Alignment

## Project Overview

This repository contains the implementation of the Carrillo-Lipman algorithm for multiple sequence alignment (MSA). Our focus is on the minimum edit distance problem, aligning a collection of sequences within a k-dimensional space to optimize the sum-of-pairs (SP) cost. This approach is particularly effective in bioinformatics for aligning protein sequences.

### Initial Sequence Processing
In this stage, we prepared the dataset for alignment. This involved processing protein sequences from various species and genes (e.g., AHSG, APOA1, HBAZ). To manage time and computational demands, sequence lengths were reduced from 120-250 to 60-70 amino acids. The reference multiple alignments can be accessed [here](https://drive.google.com/file/d/1S_QzhT34_zmEmEoUo4SxTpVnLb1KLgrw/view?usp=share_link).


## How to Run Carrillo-Lipman

### Setup Instructions
1. **Environment Setup:**
   - The `environment.yaml` file lists all the necessary Python packages.
   - Create the required environment using the command: `conda env create -f environment.yml`.

2. **Running the Algorithm:**
   - Execute the Carrillo-Lipman algorithm with the command: `python <path to your reference alignment file>`.
   - The alignment file should be in "stockholm" format.

3. **Reproducing Experiments:**
   - To replicate our experiments, download the data provided and use the `run.evaluation.zsh` script.

### Output
The output is csv file with precision, recall, F-1 scores for the corresponding alignment. In addition there is information on the ratio of the pruned nodes while executing MSA.
