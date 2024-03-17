# Simple Pipeline on How to Create a Pangenome with Anvi'o
The folowing lines describe **simple steps** to create a **pangenome** using anvi'o platform. For an extensive description of the workflow and insights into each program please visit (and read) anvi'o documentation website [Anvi'o pangenome workflow](https://merenlab.org/2016/11/08/pangenomics-v2/). <br> If you use this pipeline please reference the repository and cite Anvi'o and third party programs accordingly.

### Pipepline
1. Install [Anvi'o](https://anvio.org/install/) and [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) if not in your system.
```bash
# Install ncbi-download genomes
conda install -c conda-forge ncbi-datasets-cli
```
2. Prepare the workspace where the magic will happen.
```bash
# Change to a folder of preference.
# In this example, I will create everythin in the "~/Downloads/" folder
cd ~/Downloads/pangenome_anvio_simple_scripts

# Create working folder and subfolders
mkdir -p pangenome_simple/{00-download_genomes,01-genomes_raw,02-genomes_edited,03-contigs_db,04-pangenome,05-pan_summary,99-data}

# Change to working directory 
cd pangenome_simple
```
3. Move your genomes/MAGs/SAGs/bins ( fasta format; e.g. my_genome_1.fa; my_mag_A.fasta) to the subfolder 01-genomes_raw and skip to Step 4. **Optionally** you can download genomes from ncbi ().
```bash
#######################################
#                                     #
# I already have genomic fasta files. #
#                                     #
#######################################
# Alternatively copy genomic sequences (files.fa) to a new folder (01-genomes_raw)
# NOTE: The fasta filenames should start with a letter and contain only alphanumeric and underscore (_) characters
cp path/to/my_genomic_files/*.fa $PWD/pangenome_simple/01-genomes_raw

# Create list of fasta filenames (these will be used as IDs for the genomes and can be changed in when constructiing the pangenome)
ls 01-genomes_raw | sed -e 's/\.fa//' > genome_ids.txt

###############################
#                             #
# I want to download genomes. #
#                             #
###############################
# Here I will download all RefSeq Buchnera genomes that are assembled to chromosome level (n=26; date March 14,2024)
datasets download genome taxon "Buchnera" --include genome --assembly-level chromosome --assembly-source 'RefSeq' --filename 00-download_genomes/genomic_file-Buchnera.zip

# Unzip file
unzip 00-download_genomes/genomic_file-Buchnera.zip -d 00-download_genomes

# Create list of Assembly IDs (these will be used as IDs for the genomes)
ls -1 -d 00-download_genomes/ncbi_dataset/data/*/ | sed -e 's/.*data//' |sed -e 's/\///g'  | sed -e 's/\./_/' > genome_ids.txt

# Rename files (should not contain any special character such as dots '.')
for g_id in `cat genome_ids.txt`
do
assembly_id=$( echo "$g_id" | rev | sed -e 's/_/\./' | rev ) # Reformat to original assembly id
old_genome_file_name=$( ls 00-download_genomes/ncbi_dataset/data/${assembly_id}/ ) # Store old name
new_genome_file_name=${g_id}.fa # Store new name
# Create a copy of the genomic sequence 
cp 00-download_genomes/ncbi_dataset/data/${assembly_id}/${old_genome_file_name} 01-genomes_raw/$new_genome_file_name
done

```
4. Prepare genomes and contigs databases
```bash
# Reformat fasta files. This program will: simplyfy header names, remove contigs smaller than 2000 nucelotides,  substitute non-canonical nuceoltides for N's, and create a report for the changes. In this example, the new header will be the same as the assembly ID + number.
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
anvi-gen-contigs-database -f 02-genomes_edited/${g_id}-edited.fa \
    -n ${g_id} \
    --num-threads 10 \
    --output-db-path 03-contigs_db/${g_id}-contigs.db
done
```
5. Annotate genomes.
```bash
# MARKER GENE ANNOTATION
for g_id in `cat genome_ids.txt`
do
contigs_db=03-contigs_db/${g_id}-contigs.db
anvi-run-hmms -c $contigs_db --num-threads 10 --also-scan-trnas # Marker genes: Bac 71, Archeae, BUSCO, 16S rRNA genes, and scan for tRNAs
anvi-run-scg-taxonomy -c $contigs_db --num-threads 10 # Single-copy core gene taxonomy (GTDB)
anvi-run-trna-taxonomy -c $contigs_db --num-threads 10 # tRNA gene taxonomy (GTDB)
done

# FUNCTIONAL ANNOTATION (check databases are installed and programs available and adjust if necessary)
for g_id in `cat genome_ids.txt`
do
contigs_db=03-contigs_db/${g_id}-contigs.db
threads=30
anvi-run-cazymes -c $contigs_db --num-threads $threads # cazymes
anvi-run-interacdome -c 03-contigs_db/${g_id}-contigs.db --num-threads $threads --output-file-prefix 02-genomes_edited/${g_id}-interacdome # interaction domains
anvi-run-kegg-kofams -c 03-contigs_db/${g_id}-contigs.db --num-threads $threads # KOfams KEGG
anvi-run-pfams -c 03-contigs_db/${g_id}-contigs.db --num-threads $threads # Pfams
anvi-run-ncbi-cogs -c 03-contigs_db/${g_id}-contigs.db --num-threads $threads # COG20 NCBI COGs
done
```
6. Create pangenome.
```bash
# First you will need to create an "external" ini file. This is a two-column TSV file. The first column is the ID of the genome you want to give it in your pangenome, and the second column is the full path to its contigs database. You can change the genome names accordingly. I prefere to keep the same name as in the contigs database.
echo -e "name\tcontigs_db_path" > 04-pangenome/external-ini.txt # Do not change headers.
for g_id in `cat genome_ids.txt`
do
echo -e "$g_id\t$PWD/03-contigs_db/$g_id-contigs.db" >> 04-pangenome/external-ini.txt
done

# Create a storage database for the genomes. This will store the information (gene calls, annotation, name) of all your genomes in one single database. The output must have '-GENOMES.db' as ending.
anvi-gen-genomes-storage -e 04-pangenome/external-ini.txt -o 04-pangenome/Buchnera_ref-GENOMES.db

# Run the pangenome. This command will 1) compare every pair of genes (amino acid sequences) from your genomes using ncbi blast, 2) identify genes that are similar using a bitscore (minbit = 0.5), 3) cluster the genes by implementing mcl clustering algorithm (mcl=10, 10 used for noot closlely related organisms such as genus or family level), 4) perform multiple alignment for each gene cluster. If the genomes are highly diverse, a large amount of gene clusters will be produced; anvi'o will not perform hierarchical clustering if more than 20,000 gene clusters are produced unless it is enforced to do so (But it might take a long time).
anvi-pan-genome -g 04-pangenome/Buchnera_ref-GENOMES.db \
    --use-ncbi-blast \
    --minbit 0.5 \
    --mcl-inflation 10 \
    --align-with muscle \
    --project-name "Buchnera_ref" \
    --output-dir 04-pangenome/pan_db \
    --num-threads 30 \
    --enforce-hierarchical-clustering

# [OPTIONAL] Compute genome similarity, average nucelotide identity (ANI), using the program pyANI and ANIb as method).
# This will compute the percentage of identity among all genomic sequences in your pangenome. Results will be stored in a dedicated folder and appended to the pangenome database for easy display. There are multiple parameters that can be furhter adjusted.
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
    --output-dir 05-pan_summary/Buchnera_ref-cas3-SUMMARY
```
