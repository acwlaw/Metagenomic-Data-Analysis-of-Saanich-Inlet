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

```
/home/micb405/Group10/Project2$ perl5.26.0 
/home/micb405/resources/project_2/MaxBin-2.2.4/run_MaxBin.pl -contig ../megahit_out/final.contigs.fa -out maxBin -reads 
/home/micb405/data/project_2/SI072_LV_100m_DNA_R1.fastq.gz -reads2 
/home/micb405/data/project_2/SI072_LV_100m_DNA_R2.fastq.gz
```

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

```
while read line; do bin=$( echo $line | awk '{ print $1 }'); sid=$( echo $bin | awk -F. '{ print $1 }'); if [ -f MaxBin/$bin.fasta ]; then best_hit=$(lastal -f TAB -P 4 /home/micb405/resources/project_2/db_SILVA_128_SSURef_tax_silva MaxBin/$bin.fasta | grep -v "^#" | head -1); echo $bin,$sid,$best_hit | sed 's/,\| /\t/g'; fi; done<GT10Complete_LT5Contam_MAGs_checkM.tsv >LAST_SILVA_alignments.BEST.tsv


while read line; do accession=$( echo $line | awk '{ print $4 }'); bin=$( echo $line | awk '{ print $1 }' ); if [ ! -z $accession ]; then last_hit=$( grep "$accession" /home/micb405/resources/project_2/SILVA_128_SSURef_taxa_headers.txt | awk '{ $1=""; print $0 }'); echo $bin,$last_hit; fi; done<LAST_SILVA_alignments.BEST.tsv >LAST_SILVA_classifications.BEST.csv
```

### Step 5: Prokka - Annotate MAGs
Prokka is used to annotate bacterial genomes and is thus able to identify features of the gene.
The generalized command for prokka. A prefix of the bin was assigned for clarity.
```
prokka maxBin.*.fasta --outdir /home/micb405/Group10/Project2/prokka_out/maxBin* --prefix maxBin*
```

### Step 6: BWA/RPKM - Abundance Estimation
BWA is used to align a reference sequence against the FASTQ files of the genome sequence at 100m depth (both forward and reverse reads). The reference sequence used is a collection of nitrogen cycle genes that were a result of a grep search against the Prokka output .tsv file from all of our high quality MAGs.
RPKM (Reads per Kilobase per Million) is then used to normalize the differences in sequence length of contigs and number of contigs (sequence depth) allowing us to visualize the abundance of nitrogen cycle genes.

Indexing the reference sequence:
```
bwa index /home/micb405/Group10/Project2/MEGAHIT/SI072_LV_100m/SI072_LV_100m.contigs.fa
```
Alignment:
```
nohup bwa mem -t 4 /home/micb405/Group10/Project2/MEGAHIT/SI072_LV_100m/SI072_LV_100m.contigs.fa \
/home/micb405/data/project_2/SI072_LV_100m_DNA_R1.fastq.gz /home/micb405/data/project_2/SI072_LV_100m_DNA_R2.fastq.gz \
1>/home/micb405/Group10/Project2/bwa/SI072_LV_100m_DNA.sam 2>/home/micb405/Group10/Project2/bwa/SI072_LV_100m_DNA.bwa.stderr &
```
RPKM Calculations:
```
/home/micb405/resources/project_2/rpkm -c /home/micb405/Group10/Project2/MEGAHIT/SI072_LV_100m/SI072_LV_100m.contigs.fa \
-a /home/micb405/Group10/Project2/bwa/SI072_LV_100m_DNA.sam -o /home/micb405/Group10/Project2/RPKM/SI072_LV_100m_DNA_RPKM.csv
```
Mean RPKM of each MAG:
```
/home/micb405/resources/project_2/find_mag_rpkm_average.py
```
Generating list of all the bins to generate the desired .csv file:
```
ls maxBin/highquality/*fasta >mag_list.txt

/home/micb405/resources/project_2/find_mag_rpkm_average.py -l mag_list.txt \
-r /home/micb405/Group10/Project2/RPKM/SI072_LV_100m_DNA_RPKM.csv \
-o /home/micb405/Group10/Project2/RPKM/SI072_LV_100m_DNA_RPKM.csv
```
### Step 7: Estimation of N-cycling Gene Abundance within the Microbial Community
This step is to create a normalized abundance for our nitrogen genes of interest and will be outputted in a .csv file with the Bin identifer, Gene, Prokka gene ID, and RPKM.

**Creating a reference set:**

Extracting annotated N-cycling genes outputted by prokka:
```
while read line; 
do grep $line Prokka/Bins/SI072_LV_*/*tsv >>bin_nitrogen_cycler_genes.txt; 
done<nitrogen_cyclers.txt

while read line
do ffn=$( echo $line | awk -F':' '{ print $1 }' | sed 's/.tsv/.ffn/g' )
prefix=$( echo $line | awk '{ print $1 }' | awk -F':' '{ print $2 }' )
grep "$prefix" $ffn; done<bin_nitrogen_cycler_genes.txt >bin_nitrogen_cycler_headers.txt
```

**Alignment using BWA**
```
bwa index bin_nitrogen_cycler_genes.ffn
```
```
nohup bwa mem -t 12 bin_nitrogen_cycler_genes.ffn \
/home/micb405/data/project_2/SI072_LV_165m_DNA_R1.fastq.gz \
/home/micb405/data/project_2/SI072_LV_165m_DNA_R2.fastq.gz \
1>bin_nitrogen_cycler_genes_165m.sam 2>bin_nitrogen_cycler_genes.bwa.stderr &
```
**Abundance Calculation**
```
/home/micb405/resources/project_2/rpkm -c bin_nitrogen_cycler_genes.ffn \
-a bin_nitrogen_cycler_genes_165m.sam \
-o bin_nitrogen_cycler_genes_165m_RPKM.csv \
--verbose
```
