# BacterialCheating

Code and data for the study:

**"Bacterial cheating is a widespread resource acquisition strategy in soils"**  

This repository contains scripts and processed datasets used to analyze microbial resource acquisition strategies across large collections of bacterial genomes and metagenome-assembled genomes (MAGs).

---

## Study Overview

To test whether microbial cheating strategies extend beyond controlled litter decomposition experiments, we applied a genomic framework developed from **litter MAGs in a Mediterranean biome (n = 431)** to three additional genome collections:

- Soil MAGs from global biomes 
- Soil isolate genomes
- Environmental and host-associated MAGs from GlobDB 

All four genome datasets were combined to compute a **ranked cheating index**, which ranks genomes based on their raw cheating index across the pooled dataset and scales ranks between 0 and 1. 

---


# Repository Organization

`Scripts`  
Contains R scripts used to conduct statistical analyses and generate figures.

`Data`  
Contains raw and processed datasets used in the analyses.

`Plots`  
Contains intermediate outputs, phylogenetic trees, statistical summaries, and final figures generated during the analyses.

---

# How to Reproduce the Analyses

All analyses were performed in **R (v4.4.2)**.

Scripts can be executed directly in **R or RStudio**. The scripts in the `Scripts` directory contain the full workflow used to generate the analyses and figures presented in the manuscript.

Most statistical analyses and figure generation can be reproduced on a standard desktop computer.

---

# Contact

**Yingyi Fu**  
University of Edinburgh  
yingyi.fu@ed.ac.uk  

**Ashish Malik**  
University of Edinburgh  
ashish.malik@ed.ac.uk

