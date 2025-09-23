#!/bin/bash
#SBATCH -J 01x02x03
#SBATCH -p nocona
#SBATCH -o log/%x.out
#SBATCH -e log/%x.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=3G

# Please adjust the submitting argument according to your institute partition

# Step 1: Sequence Assembly
python 01_exons_assembly.py -t 64 -r1 SPrOut_test/01x02x03.R1.fastq -r2 SPrOut_test/01x02x03.R2.fastq \
	-p 01x02x03 -g sample_data/gene.list.txt -m sample_data/50targetfiles.fasta \
	--output_hyb 01x02x03_hyb \
	--output_exon 01x02x03_exon

# Step 2: Exon Tree Creation
python 02_exon_trees.py -t 64 -p 01x02x03 \
	-e 01x02x03_exon \
	-r sample_data/order106_refs_50genes \
	--output_phylo 01x02x03_phylo

# Step 3: Distance Matrix Calculation
python 03_distance_matrices.py -t 64 -p 01x02x03 --threshold 1 \
	--input_phylo 01x02x03_phylo \
	--output_tree 01x02x03_matrix

# Step 4: Prediction and Identification into Order
python 04_prediction.py -i 01x02x03.cumulative_dist.csv -o 01x02x03.predictions.csv -tl o -z 0.5 -to 01x02x03.order_candidates.txt

# ----------------
# If only need the order level, you can ignore the following commands

# Family Level Predictions
# Step 1: Select families from predicted orders
mkdir 01x02x03_ref_family
while read line; do python misc/pick_match_list.py family298_refs_50genes/${line} 01x02x03_ref_family/${line} 01x02x03.order_candidates.txt; done < gene.list.txt

# Step 2: Use the selected families as reference to reconstruct phylogeny
python 02_exon_trees.py -t 64 -p 01x02x03 -r 01x02x03_ref \
        -e 01x02x03_exon \
        --output_phylo 01x02x03_fam_phylo

# Step 3: Distance Matrix Calculation
python 03_distance_matrices.py -t 64 -p 01x02x03 --threshold 1 \
        --input_phylo 01x02x03_fam_phylo \
        --output_tree 01x02x03_fam_matrix

# Step 4: Prediction and Identification into Order
python 04_prediction.py -i 01x02x03.cumulative_dist.csv -o 01x02x03.predictions_fam.csv -tl f -z 0 -to 01x02x03.family_candidates.txt



