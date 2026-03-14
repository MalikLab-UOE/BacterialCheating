# BacterialCheating

Code and data for the study:

**"Bacterial cheating is a widespread resource acquisition strategy in soils"**  
(Manuscript in preparation)

This repository contains scripts and processed datasets used to analyze microbial resource acquisition strategies across large collections of bacterial genomes and metagenome-assembled genomes (MAGs).

---

## Study Overview

To test whether microbial cheating strategies extend beyond controlled litter decomposition experiments, we applied a genomic framework developed from **litter MAGs in a Mediterranean biome (n = 431)** to three additional genome collections:

- Soil MAGs from global biomes (n = 28,932)
- Soil isolate genomes (n = 5,356)
- Environmental and host-associated MAGs from GlobDB (n = 249,311)

All four genome datasets were combined to compute a **ranked cheating index**, which ranks genomes based on their raw cheating index across the pooled dataset and scales ranks between 0 and 1. Genomes with ranked cheating index > 0.9 (top 10%) were classified as **cheaters**.

---

## Repository Structure

Data  
Raw and processed datasets used in the analyses.

Scripts  
R scripts used for statistical analyses and figure generation.

Plots  
Intermediate outputs, phylogenetic trees, statistical results, and final figures.

---

## Software

All analyses were performed in **R (v4.4.2)**.  
Scripts can be executed directly in **R or RStudio**.

---

## Contact

University of Edinburgh  

Yingyi Fu  
yingyi.fu@ed.ac.uk  

Ashish Malik  
ashish.malik@ed.ac.uk

