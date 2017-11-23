# Project 2: Metagenomic Data Analysis 
## Reconstruction of Nitrogen Cycle at 100m Depth of the Saanich Inlet

Project Members: Alex Law, Bing Fan, Hajar El Bakkouri, Hannah Generoso, Tiffany Leung, Tommy (He) Huan

Link to the [full report.](https://docs.google.com/document/d/1pEF2VYwSx9Uh5TSB2BYKabZZf5Ln8tIQLf5F-qeYLQU/edit?usp=sharing)

## Method Workflow
### Step 1: MEGAHIT - Genome Assembly
MEGAHIT is used to assemble the metagenomes using succinct de Brujin graphs to output a single FASTA file containing contigs to be used as input for binning.

```
megahit -1 /home/micb405/data/project_2/SI072_LV_100m_DNA_R1.fastq.gz -2 /home/micb405/data/project_2/SI072_L_100m_DNA_R2.fastq.gz \
--k-min 27 --k-max 147 --k-step 20 \
--min-contig-len 1000 -m 0.07 -t 2 \
--out-dir /home/micb405/Group10/Project2/Megahit/SI072_LV_100
```

### Step 2: MaxBin 2.0 - Assigning MAGs
MaxBin 2.0 is used to cluster the metagenomic contigs into appropriate MAGs (metagenome assembled genomes), partitioning each bin to holds contigs consisting of the same species.

### Step 3: checkM - Identifying Best MAGs
checkM identifies high quality MAGs by estimating the completeness and contamination of a genome using marker genes that are specific to a genome's inferred lineage within a reference genome tree. The code was ran by Connor Morgan-Lang

``` 
checkm lineage_wf --tab_table -x .fasta --threads 4 --pplacer_threads 4 $BIN_DIR \
/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkm_output/ >/mnt/nfs/sharknado/Connor_MICB405_sandbox/ProcessedData/checkM/Reference/$sid\_checkM_stdout.tsv
```
The results were then filtered based on the following criteria:
* Completeness (13th column) > 10
* Contamination (12th column) < 5

```
awk -F"\t" '{ if ($12>10 && $13<5) print $0 }' *checkM_stdout.tsv >GT10Complete_LT5Contam_MAGs_checkM.tsv
```

### Step 4.1: Mash - Identify Closest Genomic Relative
Mash is used to identify the closest genomic relative in RefSeq by creating a compressed sketch representation of the sequence and determining the distance of query sequences against the provided input references. Results were outputted with a maximum P-Value of 0.01.

Mash was run against three different databases:
* refseq.genomes.k21s100.msh
* Saanich_365_taxassigned_SAGs_k21s1000.sig.msh
* Saanich_QCd_SAGs_k21s1000.sig.msh 

```
mash dist -v 0.01 /home/micb405/resources/project_2/refseq.genomes.k21s100.msh *.fasta > mash_ref.tab

mash dist -v 0.01 /home/micb405/resources/project_2/Saanich_365_taxassigned_SAGs_k21s1000.sig.msh *.fasta > mash_365.tab

mash dist -v 0.01 /home/micb405/resources/project_2/Saanich_QCd_SAGs_k21s1000.sig.msh *.fasta > mash_QCd.tab 
```
The following command was used to sort the results based on the top hits and the p-values:
```
sort -gk3 distances.tab | head
```

### Step 4.2: LAST - Identify Closest Genomic Relative
LAST is used annotate unassigned hits from Mash using an alignment-based method for quering marker genes that finds similar regions between the sequences. The tool was run against a SILVA 128 database.

### Step 5: Prokka - Annotate MAGs
Prokka is used to annotate bacterial genomes and is thus able to identify features of the gene.

### Step 6: RPKM/BWA - Abundance Estimation


