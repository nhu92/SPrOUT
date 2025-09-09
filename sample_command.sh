#!/bin/bash
#SBATCH -J RENAME
#SBATCH -p nocona
#SBATCH -o log/%x.out
#SBATCH -e log/%x.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

# Step 1: Sequence Assembly
python 01_exons_assembly.py -t 64 -r1 RENAME.R1.fastq -r2 RENAME.R2.fastq \
	-p RENAME -g gene.list.txt -m 50targetfiles.fasta \
	--output_hyb RENAME_hyb \
	--output_exon RENAME_exon

# Step 2: Exon Tree Creation
python 02_exon_trees.py -t 64 -p RENAME \
	-e RENAME_exon \
	-r order106_refs_50genes \
	--output_dir RENAME_phylo

# Step 3: Distance Matrix Calculation
python 03_distance_matrices.py -t 64 -p RENAME --threshold 1 \
	--input_dir RENAME_phylo \
	--output_dir RENAME_matrix

# Step 4: Prediction and Identification into Order
## We perform 3 different order level predictions at z-score=-0.1/0.5/0.9
python 04_prediction.py -i RENAME.cumulative_dist.csv -o RENAME.predictions.csv -tl o -z -0.1 -to RENAME.order_candidates_NPV95.txt

python 04_prediction.py -i RENAME.cumulative_dist.csv -o RENAME.predictions.csv -tl o -z 0.5 -to RENAME.order_candidates_NPV90.txt

python 04_prediction.py -i RENAME.cumulative_dist.csv -o RENAME.predictions.csv -tl o -z 0.9 -to RENAME.order_candidates_NPV85.txt

# ----------------
# If only need the order level, you can ignore the following commands

# Family Level Predictions
## We run 3 tests here for different NPV levels.
## For NPV=90% (z-score=0.5)
# Step 1: Select families from predicted orders
mkdir RENAME_ref_NPV90
while read line; do python pick_match_list.py family298_refs_50genes/${line} RENAME_ref_NPV90/${line} RENAME.order_candidates_NPV90.txt; done < gene.list.txt

# Step 2: Use the selected families as reference to reconstruct phylogeny
python 02_exon_trees.py -t 64 -p RENAME -r RENAME_ref_NPV90 \
        -e RENAME_exon \
        --output_dir RENAME_fam90_phylo

# Step 3: Distance Matrix Calculation
python 03_distance_matrices.py -t 64 -p RENAME --threshold 1 \
        --input_dir RENAME_fam90_phylo \
        --output_dir RENAME_fam90_matrix

# Step 4: Prediction and Identification into Order
python 04_prediction.py -i RENAME.cumulative_dist.csv -o RENAME.predictions_fam90.csv -tl f -z 0 -to RENAME.family90_candidates.txt

# ---
## For NPV=95% (z-score=-0.1)
# Step 1: Select families from predicted orders
mkdir RENAME_ref_NPV95
while read line; do python pick_match_list.py family298_refs_50genes/${line} RENAME_ref_NPV95/${line} RENAME.order_candidates_NPV95.txt; done < gene.list.txt

# Step 2: Use the selected families as reference to reconstruct phylogeny
python 02_exon_trees.py -t 64 -p RENAME -r RENAME_ref_NPV95 \
        -e RENAME_exon \
        --output_dir RENAME_fam95_phylo

# Step 3: Distance Matrix Calculation
python 03_distance_matrices.py -t 64 -p RENAME --threshold 1 \
        --input_dir RENAME_fam95_phylo \
        --output_dir RENAME_fam95_matrix

# Step 4: Prediction and Identification into Order
python 04_prediction.py -i RENAME.cumulative_dist.csv -o RENAME.predictions_fam95.csv -tl f -z 0 -to RENAME.family95_candidates.txt

# ---
## For NPV=85% (z-score=0.9)
# Step 1: Select families from predicted orders
mkdir RENAME_ref_NPV85
while read line; do python pick_match_list.py family298_refs_50genes/${line} RENAME_ref_NPV85/${line} RENAME.order_candidates_NPV85.txt; done < gene.list.txt

# Step 2: Use the selected families as reference to reconstruct phylogeny
python 02_exon_trees.py -t 64 -p RENAME -r RENAME_ref_NPV85 \
        -e RENAME_exon \
        --output_dir RENAME_fam85_phylo

# Step 3: Distance Matrix Calculation
python 03_distance_matrices.py -t 64 -p RENAME --threshold 1 \
        --input_dir RENAME_fam85_phylo \
        --output_dir RENAME_fam85_matrix

# Step 4: Prediction and Identification into Order
python 04_prediction.py -i RENAME.cumulative_dist.csv -o RENAME.predictions_fam85.csv -tl f -z 0 -to RENAME.family85_candidates.txt
