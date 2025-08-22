# Tinamou phylogenomics
This repository contains scripts and code for Musher et al. (in review)

Folder details:

ABBA-BABA-windows contains data needed to run and the results of the 100 kb window introgression analysis, including an r-script to generate Figure 4B.

"ASTRAL" contains files with all collapsed gene trees in a single file for each dataset used as input for astral

"phylonet" contains input and output files for phylonet along with breakpoint analysis for determining the optimal m-value

"RF.distances" contains a script for replicating several statistics and figures in the manuscript, "PIS_rf_dist_collapsed.R", output tables, and alignment/genetree files for each dataset. The R-script also contains code for looking at node heights of phylogenetic triplets.

"MSCquartets" contains data and a script that will replicate quartet analyses for assessing the effect of ILS

"beast" contains input and results form the Beast analysis 

"Filter_loci_for_beast" contains autsomal UCEs and gene trees for filtering prior to beast analysis

"Trees" contains all concatenated and astral phylogenies for each dataset. Those with prefix "concord" contain gene concordance values at all nodes. Those with prefix "crypt" show only relationships for clade A.

TESS contains a script for replicating the TESS analysis

Table S1 is the supplementary table 1 referenced in the article text containing a list of samples used in the study

Table S6 is a table of clade ages

Tinamou_assembly_pipeline.txt is a list of commands used for bioinformatic assembly of tinamou whole genomes
