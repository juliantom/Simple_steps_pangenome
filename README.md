# Simple Pipeline on How to Create a Pangenome with Anvi'o
The folowing file contains **simple steps** to download genomes and create a **pangenome** using anvi'o platform. 

### Why did I made this space?
I understand that sometimes we just want to have the code and get the result without worrying about the whys. Although this way of thinking is inadvisable from both a computational and a biological perspective, some of you might find it useful. If you want to understand the reasons for choosing each program and parameter, you can find an extensive description on anvi'o's documentation website [Anvi'o pangenome workflow](https://merenlab.org/2016/11/08/pangenomics-v2/) or dig even deeper into each program reference.

### Citation
* If you use this pipeline please acknowledge it by including the link to the repository "https://github.com/juliantom/Simple_steps_pangenome".
* Also, cite Anvi'o and third party programs accordingly.

### Pipepline
1. Check if Anvi'o and NCBI *datasets* are installed in your system. If not, you can get them from here: [Anvi'o install](https://anvio.org/install/) and [NCBI datasets install](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).
```bash
# ANVIO
# Activate anvio and check version. I am using anvio v8.
conda activate anvio-8
anvi-self-test -v

# Anvi'o .......................................: marie (v8)
# Python .......................................: 3.10.13
# Profile database .............................: 38
# Contigs database .............................: 21
# Pan database .................................: 16
# Genome data storage ..........................: 7
# Auxiliary data storage .......................: 2
# Structure database ...........................: 2
# Metabolic modules database ...................: 4
# tRNA-seq database ............................: 2

# DATASETS
# Easy install (if necessary)
conda install -c conda-forge ncbi-datasets-cli
# Check version. datasets version: 16.9.0.
datasets --version
```
2. Prepare the workspace where the magic will happen.
```bash
# Create a folder and move there.
mkdir  ~/Desktop/my_anvio_pangenomes && cd ~/Desktop/my_anvio_pangenomes

# Create a working folder and subfolders.
mkdir -p easy_pangenome/{00-download_genomes,01-genomes_raw,02-genomes_edited,03-contigs_db,04-pangenome,05-pan_summary,99-data}

# Change to the working folder.
cd easy_pangenome
```
3. Copy your genomes/MAGs/SAGs/bins in fasta format (e.g. my_genome_1.fa; my_mag_A.fasta) to the subfolder 01-genomes_raw and skip to Step 4. **Optionally** you can download genomes from NCBI.
```bash
#######################################
#                                     #
# I already have genomic fasta files. #
#                                     #
#######################################
# Copy the fasta files to 01-genomes_raw
# NOTE: The filenames should only contain alphanumeric and underscore (_) characters.
cp path/to/my_genomic_files/*.fa $PWD/easy_pangenome/01-genomes_raw

# Create a list with the filenames.
ls 01-genomes_raw | sed -e 's/\.fa//' > genome_ids.txt

#######################################
#                                     #
#     I want to download genomes.     #
#                                     #
#######################################
# I will use Buchnera as example and download all RefSeq genomes assembled to chromosome level (n=26; date March 21,2024)
datasets download genome taxon "Buchnera" --include genome --assembly-level chromosome --assembly-source 'RefSeq' --filename 00-download_genomes/genomic_file-Buchnera.zip

# Unzip file
unzip 00-download._genomes/genomic_file-Buchnera.zip -d 00-download_genomes

# Create list with the Assembly IDs.
ls -1 -d 00-download_genomes/ncbi_dataset/data/*/ | sed -e 's/.*data//' |sed -e 's/\///g' | sed -e 's/\./_/' > genome_ids.txt

# Rename files.
for g_id in `cat genome_ids.txt`
do
assembly_id=$( echo "$g_id" | rev | sed -e 's/_/\./' | rev ) # Variable with original Assembly id
old_genome_file_name=$( ls 00-download_genomes/ncbi_dataset/data/${assembly_id}/*_genomic.fna ) # Variable to store old name
new_genome_file_name=${g_id}.fa # Variable to store new name
# Copy the genomic sequence and change name
cp ${old_genome_file_name} 01-genomes_raw/$new_genome_file_name
done

```
4. Prepare genomes and contigs databases
```bash
# Reformat fasta files.
# EXPLANATION. This program will take your input fasta file and reformat it. 1) The headers will be simplified using as prefix the genome ID, 2) contigs smaller than 2000 nucelotides will be removed, 3) non-canonical bases other than {A,T,C,G}, will be subtituted for N's, and 4) a file reporting the changes to the headers will be produced.
for g_id in `cat genome_ids.txt`
do
anvi-script-reformat-fasta 01-genomes_raw/${g_id}.fa \
    -o 02-genomes_edited/${g_id}-edited.fa \
    --report-file 02-genomes_edited/${g_id}-edited.fa-report.txt \
    --simplify-names \
    --prefix ${g_id} \
    --min-len 2000 \
    --seq-type NT
done

# Create contigs DB for each genome.
for g_id in `cat genome_ids.txt`
do
exec_th=4 # Number of threads (modify accordingly)
anvi-gen-contigs-database -f 02-genomes_edited/${g_id}-edited.fa \
    -n ${g_id} \
    --num-threads $exec_th \
    --output-db-path 03-contigs_db/${g_id}-contigs.db
done
```
5. Genome annotation is **optional** but gratifying. First, check if databases are setup and install if necessary.
```bash
# MARKER GENE ANNOTATION
for g_id in `cat genome_ids.txt`
do
exec_th=4 # Number of threads (modify accordingly)
contigs_db=03-contigs_db/${g_id}-contigs.db
anvi-run-hmms -c $contigs_db --num-threads $exec_th --also-scan-trnas # Marker genes: Bac 71, Archeae, BUSCO, 16S rRNA genes, and scan for tRNAs
anvi-run-scg-taxonomy -c $contigs_db --num-threads $exec_th # Single-copy core gene taxonomy (GTDB)
anvi-run-trna-taxonomy -c $contigs_db --num-threads $exec_th # tRNA gene taxonomy (GTDB)
done

# FUNCTIONAL ANNOTATION
for g_id in `cat genome_ids.txt`
do
exec_th=4 # Number of threads (modify accordingly)
contigs_db=03-contigs_db/${g_id}-contigs.db
anvi-run-cazymes -c $contigs_db --num-threads $exec_th # CAZymes
anvi-run-interacdome -c 03-contigs_db/${g_id}-contigs.db --num-threads $exec_th --output-file-prefix 99-data/${g_id}-interacdome # interaction domains
anvi-run-kegg-kofams -c 03-contigs_db/${g_id}-contigs.db --num-threads $exec_th # KOfams KEGG
anvi-run-pfams -c 03-contigs_db/${g_id}-contigs.db --num-threads $exec_th # Pfams
anvi-run-ncbi-cogs -c 03-contigs_db/${g_id}-contigs.db --num-threads $exec_th # COG20 NCBI COGs
done
```
6. Create pangenome.
```bash
# Create an "external" ini file. This is a two-column TSV file. The first column is the name you want to give to each genome (here, I will use the genome ID, but you could give them a different name - alhpanumeric only). The second column is the full path to its contigs database.
echo -e "name\tcontigs_db_path" > 04-pangenome/external-ini.txt # Do not change headers.
for g_id in `cat genome_ids.txt`
do
echo -e "$g_id\t$PWD/03-contigs_db/$g_id-contigs.db" >> 04-pangenome/external-ini.txt
done

# Store gene calls, annotations of all genomes into a single file. The output must have '-GENOMES.db' as ending.
anvi-gen-genomes-storage -e 04-pangenome/external-ini.txt -o 04-pangenome/Buchnera_ref-GENOMES.db

# Run the pangenome. 
# This command will 1) compare every pair of genes (amino acid sequences) from your genomes using ncbi blast, 2) identify genes that are similar using a bitscore (minbit = 0.5), 3) cluster the genes by implementing mcl clustering algorithm (mcl=10, 10 used for not closlely related organisms such as genus or family level), 4) perform multiple alignment for each gene cluster. Enforce hierarchical clustering, if more than 20,000 gene clusters are produced. This happens when genomes are highly diverse.
anvi-pan-genome -g 04-pangenome/Buchnera_ref-GENOMES.db \
    --use-ncbi-blast \
    --minbit 0.5 \
    --mcl-inflation 10 \
    --align-with muscle \
    --project-name "Buchnera_ref" \
    --output-dir 04-pangenome/pan_db \
    --num-threads 30 \
    --enforce-hierarchical-clustering

# [Optional but exciting] Compute average nucelotide identity (ANI) and add it the pangenome (results will also be stored in a separete folder).
anvi-compute-genome-similarity -e 04-pangenome/external-ini.txt \
    -o 04-pangenome/genome_similarity \
    --pan-db 04-pangenome/pan_db/Buchnera_ref-PAN.db \
    --program pyANI \
    --method ANIb \
    --num-threads 10

```
7. Display pangenome, edit figure and find gene clusters of intereset
```bash
# Display pangenome, edit figure, store bins and save to a collection. Save edits as a state.
anvi-display-pan -g 04-pangenome/Buchnera_ref-GENOMES.db -p 04-pangenome/pan_db/Buchnera_ref-PAN.db

# I categorized gene into three groups (Core, Accessory, and Singletons) and store this collection as "cas3". Core= gene cluster present in all (100%) genomes; Accessory=gene cluster present in more than one and less than all genomes; Singletons= gene cluster present at most in one genome. 
# The display was stored as "simple_view"
# I generated an on-the-fly phylogeny by concatenating all single-copy core genes and running "generate phylogeny". The display was stored as "simple_view_with_phylo"
```
8. Produce summary tables for the pangenome. By including the collection, the table will contain gene-level information linking them to their gene-clusters, genomes, sequence, functions, and collection.
```bash
# Summarize pangenome. Use collection.
anvi-summarize -g 04-pangenome/Buchnera_ref-GENOMES.db \
    -p 04-pangenome/pan_db/Buchnera_ref-PAN.db \
    --collection-name cas3 \
    --output-dir 05-pan_summary/Buchnera_ref-pan_sum-cas3
```
