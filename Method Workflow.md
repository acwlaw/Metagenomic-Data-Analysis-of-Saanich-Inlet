# Project 2: Metagenomic Data Analysis 
## Reconstruction of Nitrogen Cycle at 100m Depth of the Saanich Inlet

Project Members: Alex Law, Bing Fan, Hajar El Bakkouri, Hannah Generoso, Tiffany Leung, Tommy (He) Huan

Link to the [full report.](https://docs.google.com/document/d/1pEF2VYwSx9Uh5TSB2BYKabZZf5Ln8tIQLf5F-qeYLQU/edit?usp=sharing)

## Method Workflow
### Step 1: MEGAHIT
MEGAHIT is used to assemble the metagenomes using succinct de Brujin graphs to output a single FASTA file contianin contigs to be used as input for binning.

```
megahit -1 /home/micb405/data/project_2/SI072_LV_100m_DNA_R1.fastq.gz -2 /home/micb405/data/project_2/SI072_L_100m_DNA_R2.fastq.gz \
--k-min 27 --k-max 147 --k-step 20 \
--min-contig-len 1000 -m 0.07 -t 2 \
--out-dir /home/micb405/Group10/Project2/Megahit/SI072_LV_100
```

### Step 2: MaxBin 2.0
MaxBin 2.0 
