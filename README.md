

# Example commands and workflows for data analysis

Datasets can be found at 10.5281/zenodo.14861210 and under BioProject PRJNA1165788

R figure generation described in https://github.com/carleton-envbiotech/Consortium_meta_omics/blob/main/Nitrifying_consortium_data_analyses_codebook_2025_02_15.Rmd with input files in https://github.com/carleton-envbiotech/Consortium_meta_omics/tree/main/Raw_data_for_figures

# Metagenome Assembly

## PacBio HiFi Long Reads

Iterative assemblies were performed using metaMDBG v. 0.3 and hifiasm_meta v.03-r063.2:
```
#!/bin/bash -l
#SBATCH --mem=464000
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --comment=""
#SBATCH --partition=standard
#SBATCH --job-name=metaMDBG
#SBATCH --account=
#SBATCH --array=1-3   # Adjust the array range based on the number of input FASTQ files

conda activate metaMDBG

# Define an array of input FASTQ files
input_files=(
    "PB127C_MBSox_4ml_fresh_rep3.hifi_reads.fastq.gz"
    "PB127BC_combined_Qps_4ml_frozen_rep1.hifi_reads.fastq.gz"
    "PB127BC_combined_MBSox_2ml_frozen_test_rep2.hifi_reads.fastq.gz"
    # Add more input FASTQ files as needed
)

# Get the input FASTQ file corresponding to the current array index
input_file=${input_files[$SLURM_ARRAY_TASK_ID - 1]}
# Extract the base name of the input FASTQ file (without extension)
input_file_basename=$(basename "$input_file" .fastq.gz)

# Run metaMDBG for the current input FASTQ file
metaMDBG asm "./metaMDBG_${input_file_basename}" "$input_file" -t $SLURM_CPUS_PER_TASK
```
Hifiasm_meta can be performed as above with the following command:

```
conda activate hifiasmmeta
hifiasm_meta -t $SLURM_CPUS_PER_TASK -o NitrifyingCombined ../NitrifyingCombined.fastq.gz  2>combinedNitrifying.log
```

### Assembly

Assemblies from each iterative run from the Pacbio HiFi read sets and combined readset were ran through the HiFi-MAG-Pipeline v. 2.0.2 for binning and QC
https://github.com/PacificBiosciences/pb-metagenomics-tools/tree/v2.0.2

```
#!/bin/bash -l
#SBATCH --mem=456000
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --comment="tmpfs_size=30G"
#SBATCH --partition=standard
#SBATCH --job-name=hifi-mag
#SBATCH --account=

conda activate snakemake

snakemake --snakefile Snakefile-hifimags.smk --configfile configs/Sample-Config.yaml --rerun-incomplete --use-conda -c 64 \
# --workflow-profile /HiFi-MAG-Pipeline/profile/ --rerun-incomplete --use-conda
```
### dRep v 3.5.0
```
dRep dereplicate dRepPB -g *.fa -p 48 -sa 0.98
```
## Illumina
### Nextflow  nf-core-mag Assembly Parameters
```
nextflow run/nf-core-mag_2.5.4/2_5_4 -profile singularity --outdir /ShortReads/Nitrifying/MAGoutGroupcorrected --input /ShortReads/Nitrifying/samplesheet.csv --gtdb_db /release214/ --kraken2_db /Databases/MAG/k2_pluspf_20230314.tar.gz --checkm_db /Databases/MAG/checkm -resume --gtdbtk_max_contamination 10 --gtdbtk_min_completeness 70 --krona_db /taxonomy --skip_spades --skip_concoct --busco_db /Databases/MAG/BUSCO/busco-data.ezlab.org/v5/data --gtdbtk_pplacer_scratch TRUE -c nextflowst.config --metaeuk_db /Databases/Nextflow '--megahit_options=--min-contig-len 1000' --postbinning_input refined_bins_only --refine_bins_dastool --binqc_tool checkm --bin_domain_classification --coassemble_group --skip_quast 
```
### Metagenome QC

Workflow adapted from pb-metagenomics-tools/HiFi-MAG-Pipeline at master · PacificBiosciences/pb-metagenomics-tools to filter short-read MAGs  output from the binning step of the nf-core-mag workflow with CheckM2 and subsequently classify using GTDB-TK R220 and on MAGs generated via short-read assembly for QC before replication.

See https://github.com/carleton-envbiotech/Consortium_meta_omics/tree/main/Genomics/ModifiedPBMetaWorkflow-MAG_QC_Filtering for the modified Snakemake workflow for Step 5 and on processing and to update the code to be compatible with GTDB-TK R220. 

# Final Metagenome Dereplication and Profiling

## Combined dataset 

Filtered short read MAGs and dereplicated PacBio HiFi Long Reads consensus MAGs from Step 1
### dRep v 3.5.0P
```
dRep dereplicate dRepNFGroupsPB -g *.fa -p 48 -sa 0.98
```
## Sylph Profiling

Example commands: 

```
sylph sketch -1 ../*_R1_*.fastq.gz -2 ../*_R2_*.fastq.gz -d Nitrifying -t 12 &&

sylph sketch /dRepNFGroupsPB/dereplicated_genomes/*.fa -o database

sylph profile /dRepNFGroupsPB/dereplicated_genomes/*.syldb Nitrifying061/*sylsp > Nitrifying061/all-to-all-profile_NF_PBExpMAGs_derep.tsv
```

Final output was NF_PBGroupMAGs_derep_Sylph.tsv after using sylph utility scripts to profile and merge. No taxonomy database was assigned as the reads were profiled against the MAGs alone.
```
python sylph_to_taxprof.py - -s all-to-all-profile_NF_PBExpMAGs_derep.tsv -o NF_PBGroupMAGs/

python merge_sylph_taxprof.py NF_PBGroupMAGs/*.sylphmpa --column sequence_abundance -o  NF_PBGroupMAGs_derep_Sylph.tsv
```
# SqueezeMeta Analysis

## Reads processed via fastp v. 0.23.4: array script
```
#!/bin/bash -l
#SBATCH --mem=20000
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-28%10  # Adjust this to your number of read pairs and %10 to control how many run simultaneously
#SBATCH --job-name=fastp
#SBATCH --account=
#SBATCH --partition=standard
#SBATCH --comment=""
#SBATCH --cluster=


source /conda/etc/profile.d/conda.sh
conda activate base
# Define the directory containing reads and results
READS_DIR="/Nitrifying/DNA"
RESULTS_DIR="/Nitrifying/DNA/cleaned"

R1_BASE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" r1_reads.txt)
R1="${R1_BASE}_R1_001.fastq.gz"
R2="${R1_BASE}_R2_001.fastq.gz"
#example read name in r1_reads.txt "NitrifyingPelletDNA_Day0_DNARNAkit_rep1_S1"
# Output files
OUT_R1="${RESULTS_DIR}/${R1_BASE}_cleaned_R1.fastq.gz"
OUT_R2="${RESULTS_DIR}/${R1_BASE}_cleaned_R2.fastq.gz"
JSON="${RESULTS_DIR}/${R1_BASE}_fastp.json"
HTML="${RESULTS_DIR}/${R1_BASE}_fastp.html"

# Run fastp
fastp \
    -i "${READS_DIR}/${R1}" \
    -o "${OUT_R1}" \
    -I "${READS_DIR}/${R2}" \
    -O "${OUT_R2}" \
    --json "${JSON}" \
    --html "${HTML}" \
    -w ${SLURM_CPUS_PER_TASK}

echo "fastp processing complete for ${R1_BASE}"

```

## Against fastp processed DNA and RNA reads (SqueezeMeta v. 1.6.3):
```
SqueezeMeta.pl -p Nitrifying/Squeeze/SqueezeCoassemblyNF -taxbinmode s+c -m coassembly -f /Nitrifying/Squeeze -b 50 -s /Squeeze/samples.tsv -t $SLURM_CPUS_PER_TASK-c 500 -binners maxbin,metabat2,concoct
```
##For 98 MAGs provided as external bins (SqueezeMeta v. 1.7.0b8):
```
SqueezeMeta.pl -p /Nitrifying/Squeeze/SqueezeExtBins -m coassembly -extbins /dRepGroupsPB-Final/MAGs \
 -f /Nitrifying/Squeeze -b 25 \
 -s /Nitrifying/Squeeze/samples.tsv -t $SLURM_CPUS_PER_TASK 
```
## SQMtools analysis
```
library(SQMtools)

NFBin <- loadSQM("SqueezeExtBins/")
NFCo <- loadSQM("SqueezeCoassemblyNF/")

rnadna <- c("Nitrifying_Day0_rep1_DNA", "Nitrifying_Day0_rep2_DNA", "Nitrifying_Day0_rep3_DNA", "Nitrifying_Week4_Aerobic_rep1_DNA", "Nitrifying_Week4_Aerobic_rep2_DNA", "Nitrifying_Week4_Aerobic_rep3_DNA", "Nitrifying_Week4_Ironreduction_rep1_DNA", "Nitrifying_Week4_Ironreduction_rep2_DNA", "Nitrifying_Week4_Ironreduction_rep3_DNA", "Nitrifying_Week4_Nitratereduction_rep1_DNA", "Nitrifying_Week4_Nitratereduction_rep2_DNA", "Nitrifying_Week4_Sulphatereduction_rep1_DNA", "Nitrifying_Week4_Sulphatereduction_rep2_DNA", "Nitrifying_Week4_Sulphatereduction_rep3_DNA", "Nitrifying_Week8_Aerobic_rep1_DNA", "Nitrifying_Week8_Aerobic_rep2_DNA", "Nitrifying_Week8_Aerobic_rep3_DNA", "Nitrifying_Week8_Ironreduction_rep1_DNA", "Nitrifying_Week8_Ironreduction_rep2_DNA", "Nitrifying_Week8_Ironreduction_rep3_DNA", "Nitrifying_Week8_Nitratereduction_rep1_DNA", "Nitrifying_Week8_Nitratereduction_rep2_DNA", "Nitrifying_Week8_Nitratereduction_rep3_DNA", "Nitrifying_Week8_Sulphatereduction_rep1_DNA", "Nitrifying_Week8_Sulphatereduction_rep2_DNA", "Nitrifying_Week8_Sulphatereduction_rep3_DNA", "Nitrifying_Day0_rep1_RNA", "Nitrifying_Day0_rep2_RNA", "Nitrifying_Day0_rep3_RNA", "Nitrifying_Week4_Aerobic_rep1_RNA", "Nitrifying_Week4_Aerobic_rep2_RNA", "Nitrifying_Week4_Aerobic_rep3_RNA", "Nitrifying_Week4_Ironreduction_rep1_RNA", "Nitrifying_Week4_Ironreduction_rep2_RNA", "Nitrifying_Week4_Ironreduction_rep3_RNA", "Nitrifying_Week4_Nitratereduction_rep1_RNA", "Nitrifying_Week4_Nitratereduction_rep2_RNA", "Nitrifying_Week4_Sulphatereduction_rep1_RNA", "Nitrifying_Week4_Sulphatereduction_rep2_RNA", "Nitrifying_Week8_Aerobic_rep1_RNA", "Nitrifying_Week8_Aerobic_rep2_RNA", "Nitrifying_Week8_Aerobic_rep3_RNA", "Nitrifying_Week8_Ironreduction_rep3_RNA", "Nitrifying_Week8_Nitratereduction_rep1_RNA", "Nitrifying_Week8_Nitratereduction_rep2_RNA", "Nitrifying_Week8_Sulphatereduction_rep1_RNA", "Nitrifying_Week8_Sulphatereduction_rep2_RNA", "Nitrifying_Week8_Sulphatereduction_rep3_RNA")
```
### Amino acid subset
```
combinedAA <- "K00022|K00128|K00135|K00137|K00149|K00252|K00290|K00293|K00306|K00382|K00468|K00471|K00473|K00474|K00626|K00658|K00824|K00825|K01034|K01035|K01041|K01506|K01582|K01692|K01782|K01825|K01843|K01844|K04462|K06101|K07250|K07511|K07514|K07515|K09186|K09187|K09188|K09189|K09251|K11419|K11420|K11421|K11422|K11423|K11424|K11425|K11426|K11427|K11428|K11429|K11430|K11431|K11432|K11433|K11703|K13609|K13645|K13646|K13647|K14085|K14157|K14268|K14959|K15588|K15736|K15737|K15791|K17451|K18011|K18012|K18013|K18014|K18201|K18202|K18494|K18804|K18826|K18854|K19743|K20795|K20796|K21672|K22410|K22748|K23385|K23700|K24034|K24405|K24406|K24407|K25316|K25317|K25986|K26061|K26062|K26063|K26064|K26065|K00020|K00140|K00157|K00166|K00167|K00186|K00187|K00188|K00189|K00248|K00249|K00253|K00263|K00271|K00632|K00822|K00826|K00827|K01027|K01028|K01029|K01640|K01641|K01847|K01848|K01849|K01907|K01964|K01965|K01966|K01968|K01969|K03334|K05605|K05606|K05607|K07508|K07509|K07513|K08683|K09478|K09699|K11263|K11381|K11410|K11538|K13524|K13766|K15036|K15037|K18472|K18660|K18661|K19312|K22568|K23146|K27094|K27095|K00139|K00259|K00260|K00261|K00262|K00264|K00265|K00266|K00272|K00278|K00294|K00608|K00609|K00610|K00764|K00811|K00812|K00813|K00814|K00820|K00823|K00830|K01424|K01425|K01437|K01580|K01744|K01755|K01756|K01779|K01914|K01915|K01939|K01940|K01948|K01953|K01954|K01955|K01956|K05597|K08324|K08590|K09758|K11358|K11540|K11541|K13051|K13566|K13821|K14260|K14272|K14454|K14455|K14592|K14681|K15371|K16871|K17761|K18309|K18310|K18311|K19244|K22457|K23265"

AANFco <- subsetFun(NFCo, fun=combinedAA)

AAnfcoALL_plot <- plotFunctions(AANFco, samples =rnadna, N = 200)

This table was filtered for RNA and top 25 hits, then unique identifiers were used for Fig S7.

write.table(AAnfcoALL_plot$data, file = "AAplotALL_Nfco.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
```
### Nitrogen metabolism functions

```
Nitrogen_kos <- "K10944|K10945|K10946|K10535|K00371|K00370|K00368|K15864|K15876|K03385|K04561|K00376|K02567|K02568"

nitrokoplot_NFco <- plotFunctions(NFCo, samples = rnadna, fun = nitrogen_kos)

write.table(nitrokoplot_NFco$data, file = "nitrokoplot_Nfco.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
```
### Mapping RNA reads to external bins 
```
rna <- c( "Nitrifying_Day0_rep1_RNA", "Nitrifying_Day0_rep2_RNA", "Nitrifying_Day0_rep3_RNA", "Nitrifying_Week4_Aerobic_rep1_RNA", "Nitrifying_Week4_Aerobic_rep2_RNA", "Nitrifying_Week4_Aerobic_rep3_RNA", "Nitrifying_Week4_Ironreduction_rep1_RNA", "Nitrifying_Week4_Ironreduction_rep2_RNA", "Nitrifying_Week4_Ironreduction_rep3_RNA", "Nitrifying_Week4_Nitratereduction_rep1_RNA", "Nitrifying_Week4_Nitratereduction_rep2_RNA", "Nitrifying_Week4_Sulphatereduction_rep1_RNA", "Nitrifying_Week4_Sulphatereduction_rep2_RNA", "Nitrifying_Week8_Aerobic_rep1_RNA", "Nitrifying_Week8_Aerobic_rep2_RNA", "Nitrifying_Week8_Aerobic_rep3_RNA", "Nitrifying_Week8_Ironreduction_rep3_RNA", "Nitrifying_Week8_Nitratereduction_rep1_RNA", "Nitrifying_Week8_Nitratereduction_rep2_RNA", "Nitrifying_Week8_Sulphatereduction_rep1_RNA", "Nitrifying_Week8_Sulphatereduction_rep2_RNA", "Nitrifying_Week8_Sulphatereduction_rep3_RNA")

nfbintax <- plotBins(NFBin, samples =rna)

write.table(nfbintax$data, file = "nfbins_RNA.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
```
## MAG Phylogenetics and Annotation

### Abricate

```
#!/bin/bash
#SBATCH --job-name=abricate
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --comment=""
#SBATCH --partition=standard
#SBATCH --account=
#SBATCH --cluster 

# Activate Conda environment
echo "Activating Conda environment..."
source /condaforge/etc/profile.d/conda.sh
conda activate /condaforge/abricate

# Run

abricate -db card --threads 32 --minid 60 MAGs/*.fa > cardmags.out
abricate -db vfdb --threads 32 --minid 60 MAGs/*.fa > vfdbmags.out
```
### Phlyogenetics

Acquire reference genomes:
```
esearch -query '"Legionellales"[ORGN] AND "latest refseq"[filter] AND "reference genome"[filter]' -db assembly | esummary | xtract -pattern DocumentSummary -def "NA" -element AssemblyAccession > Legionellales-refseq-reference-accs.txt
```
RefSeq accession for Legionellales accessed 2025/07/14 – removed GCF_002776555.1
```
esearch -query '"Achromobacter"[ORGN] AND "latest refseq"[filter] AND "reference genome"[filter]' -db assembly | esummary | xtract -pattern DocumentSummary -def "NA" -element AssemblyAccession > Achromobacter-refseq-reference-accs.txt
```
Refseq accession for Achromobacter accessed 2025/07/11
 ```   
esearch -query '"Mycobacteriaceae"[ORGN] AND "latest refseq"[filter] AND "reference genome"[filter]' -db assembly | esummary | xtract -pattern DocumentSummary -def "NA" -element AssemblyAccession > Mycobacteriaceae-refseq-reference-accs.txt
 ```
 RefSeq accession for Mycobacteriaceae accessed 2025/07/11

```
GToTree -a Legionellales-refseq-reference-ref-accs-GCF_002776555.1-removed.txt -g genbank_files.txt -f fasta_files.txt -H Gammaproteobacteria -t -L Species,Strain -j 12 -o Legionellales_GCF_002776555.1removed 
```
Rooted on Pseudomonas aeruginosa PAO1 GCF_000006765.1 - retrieved 2024/04/27
```
GToTree -a Achromobacter-refseq-reference-accs.txt -g genbank_files.txt -f fasta_files.txt -H Gammaproteobacteria -t -L Species,Strain -j 48 -o Achromobacter -F
```
Rooted on Pseudomonas aeruginosa PAO1 GCF_000006765.1 - retrieved 2024/04/27
```
GToTree -a Mycobacteriaceae-refseq-reference-accs.txt -g genbank_files.txt -f fasta_files.txt -H Actinobacteria -t -L Species,Strain -j 48 -o Myco_phylo
```
Rooted on Corynebacterium diphtheriae NCTC 13129 GCF_000195815.1 retrieved  2025/07/11

### Bakta Annotation
```
!/bin/bash
#SBATCH --job-name=bakta
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-98%20  # Adjust the 1-100 to match the number of files and %10 for the number of simultaneous tasks
#SBATCH --output=bakta-array-%A_%a.out
#SBATCH --error=bakta-array-%A_%a.err
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --comment=""
#SBATCH --partition=standard
#SBATCH --account=
#SBATCH --cluster 

# Activate Conda environment

# Set BAKTA_DB environment variable
export BAKTA_DB=/Databases/Bakta/db

# Define path to the input files and output directory
INPUT_DIR="/FinalMAGS"
OUTPUT_DIR="/FinalMAGS/baktaarray"

# Use SLURM_ARRAY_TASK_ID to get the specific file for this array job
# Assuming input files are named like 'sample1.fa', 'sample2.fa', etc.
FILE=$(ls $INPUT_DIR/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")
FILENAME=$(basename $FILE .fa)
set -e

# Verbose output
set -x

# Activate Conda environment
echo "Activating Conda environment..."
source /conda/etc/profile.d/conda.sh
conda activate /conda/bakta

# Run Bakta
bakta --db $BAKTA_DB --output $OUTPUT_DIR/${FILENAME}_bakta --prefix $FILENAME --threads $SLURM_CPUS_PER_TASK --verbose $FILE
```
### Clinker Plots
```
clinker *.gb -p plot.html -o alignments.csv -dl "," -dc 4 --force -r region_[435000_555000]:37000-97000 region_[44500_146900]:29000-89000 region_[30000_130000]:20000-60000
```
# Probe Design

Approach was adapted from methods described in "Rational probe design for efficient rRNA depletion and improved metatranscriptomic analysis of human microbiomes" doi: 10.1186/s12866-023-03037-y

See https://github.com/carleton-envbiotech/Consortium_meta_omics/tree/main/Genomics/SupplementalProbeFiles for input/outputs and 10.5281/zenodo.14861210 for .fastq and MAGs used for design

## Barrnap v. 0.9  rRNA detection 

For bacterial MAGs:
```
barrnap --threads 32 --kingdom bac --outseq complete.10.fasta --reject 0.5 complete.10.fa
```
rRNA sequences were combined, duplicates were removed using 
```
#awk '/^>/{if(header[$0]++)next; printf "%s\n", $0; getline; print}' rRNA_nitrifyingMAGS.fasta >rRNA_nitrifyingMAGS.fasta
```
Sequences were clustered at 99% using  CD-Hit v. 4.8.1
```
 cd_hit_clustering rRNA_nitrifyingMAGS.fasta rRNA_nitrifyingMAGS._99pct.fasta 0.99
```
## sortmerna v. 4.3 - ribosome outputs were concatenated and used for downstream analysis
```
sortmerna -ref /sortmerna/smr_v4.3_default_db.fasta --reads NitrifyingRNA_Day0_RiboZeroPlusKit_rep1_S1_R1_001.fastq.gz --reads NitrifyingRNA_Day0_RiboZeroPlusKit_rep1_S1_R2_001.fastq.gz --fastx --threads 14 --aligned NitrifyingRNA_Day0_RiboZeroPlusKit_rep1_S1_ribosomes --workdir /sortmerna/NitrifyingRNA_Day0_RiboZeroPlusKit_rep1_S1 --paired_in --other NitrifyingRNA_Day0_RiboZeroPlusKit_rep1_S1_R1 --out2 NitrifyingRNA_Day0_RiboZeroPlusKit_rep1_S1_R2 -v 2

sortmerna -ref /home/smithd/sortmerna/smr_v4.3_default_db.fasta --reads NitrifyingRNA_Day0_RiboZeroPlusKit_rep2_S2_R1_001.fastq.gz --reads NitrifyingRNA_Day0_RiboZeroPlusKit_rep2_S2_R2_001.fastq.gz --fastx --threads 14 --aligned NitrifyingRNA_Day0_RiboZeroPlusKit_rep2_S2_ribosomes --workdir /sortmerna/NitrifyingRNA_Day0_RiboZeroPlusKit_rep2_S2 --paired_in --other NitrifyingRNA_Day0_RiboZeroPlusKit_rep2_S2_R1 --out2 NitrifyingRNA_Day0_RiboZeroPlusKit_rep2_S2_R2 -v 2
```

Both ribosome files were concatenated into NitrifyingRNA_Day0_RiboZeroPlusKit_combined_R1.fastq.gz and NitrifyingRNA_Day0_RiboZeroPlusKit_combined_R2.fastq.gz

## Ribosomal reads were then mapped against the rRNAs identified in MAGs

Initially, we ran the bowtie alignment against the MAG rRNAs to generate the unmapped reads, those were then ran through rnaspades and the transcripts file was ran through barrnap for arc, euk, and bac, concatenated, and then deduplicated and concatenated with the MAG derived rRNA that were also deduplicated as identical fragments were annotated in many binned contigs
```
bowtie2-build --threads $threads  rRNA_nitrifyingMAGS._99pct.fasta 0.99 rRNA_index

bowtie2 -p $threads -x rRNA_index -1 NitrifyingRNA_Day0_RiboZeroPlusKit_combined_R1.fastq.gz -2 NitrifyingRNA_Day0_RiboZeroPlusKit_combined_R2.fastq.gz -S mapped_reads.sam --un-conc unaligned_reads.fastq 2> original_bowtie2_stats.txt
```
## Assembling unaligned reads using  rnaSPAdes v 3.13.0 
```
rnaspades.py -1	unaligned_reads.1.fastq -2 unaligned_reads.2.fastq --threads 32 	-o unaligned
```
Assembled rRNAs from reads were detected using:
```
barrnap --threads 32 --kingdom arc -outseq Archaea.fasta --reject 0.5 transcripts.fasta
barrnap --threads 32 --kingdom euk --outseq Bacteria.fasta --reject 0.5 transcripts.fasta
barrnap --threads 32 --kingdom bac --outseq Eukaryote.fasta  --reject 0.5 transcripts.fasta
```
rRNA sequences were combined, duplicates were removed using 
```
#awk '/^>/{if(header[$0]++)next; printf "%s\n", $0; getline; print}' Spades_Assembled_rrna.fasta > Spades_Assembled_rrna.fasta
```
Sequences were clustered at 99% using  CD-Hit v. 4.8.1
```
 cd_hit_clustering Spades_Assembled_rrna.fasta Spades_Assembled_rrna_99pct.fasta 0.99
```

## Iterating and clustering targets 
Findtargets.sh
```
#!/bin/bash

# Define input files
fastq1="/probedesign/NitrifyingRNA_Day0_RiboZeroPlusKit_combined_R1.fastq.gz"
fastq2="/probedesign/NitrifyingRNA_Day0_RiboZeroPlusKit_combined_R2.fastq.gz"
fasta_file="Assembled_nitrifyingMAGs_99pct.fasta"
# Define parameters
threads=16
min_coverage=500
output_file1="top100_abundant_regions.fasta"
output_file2="top50_abundant_regions.fasta"
output_file3="top25_abundant_regions.fasta"
# Index the rRNA sequences
bowtie2-build --threads $threads $fasta_file rRNA_index

bowtie2 -p $threads -x rRNA_index -1 $fastq1 -2 $fastq2 -S mapped_reads.sam --un-conc unaligned_reads.fastq 2> original_bowtie2_stats.txt

# Convert SAM to BAM and sort
samtools view -@ $threads -bS mapped_reads.sam | samtools sort -@ $threads -o sorted_reads.bam

# Index sorted BAM file
samtools index -@ $threads sorted_reads.bam

# Merge regions in close proximity (1–3 nt) using bedtools
bedtools bamtobed -i sorted_reads.bam | bedtools merge -d 3 -c 4 -o count,collapse -i - > merged_regions.bed

# Filter merged regions based on coverage
awk -v min_cov=$min_coverage '$4 >= min_cov' merged_regions.bed > high_coverage_regions.bed

# Extract the top 100 and top 50 most abundant regions
sort -k 4,4nr high_coverage_regions.bed | head -n 100 > top100_regions.bed
sort -k 4,4nr high_coverage_regions.bed | head -n 50 > top50_regions.bed
sort -k 4,4nr high_coverage_regions.bed | head -n 25 > top25_regions.bed
# Extract sequences of the top 100, top 50, and top 25 regions from the original fasta file
bedtools getfasta -fi $fasta_file -bed top100_regions.bed -fo $output_file1
bedtools getfasta -fi $fasta_file -bed top50_regions.bed -fo $output_file2
bedtools getfasta -fi $fasta_file -bed top25_regions.bed -fo $output_file3
# Clean up intermediate files
rm rRNA_index* mapped_reads.sam sorted_reads.bam sorted_reads.bam.bai merged_regions.bed high_coverage_regions.bed top*_regions.bed

echo "Current working directory: $(pwd)"

# CD-HIT clustering with different parameters
cd_hit_clustering() {
    local input_file="$1"
    local output_prefix="$2"
    local clustering_rate="$3"

    cd-hit -i $input_file -o ${output_prefix}_${clustering_rate}.fasta -c $clustering_rate -n 5 -M 4000 -T $threads
}

# Clustering with different parameters and capture Bowtie2 stats
for clustering_rate in 0.80 0.85 0.90 0.95; do
    for top_regions_file in $output_file1 $output_file2 $output_file3; do
        output_prefix=$(basename $top_regions_file .fasta)
        echo "Clustering file: clustered_sequences_${output_prefix}_${clustering_rate}"
        cd_hit_clustering $top_regions_file clustered_sequences_${output_prefix} $clustering_rate

        echo "Building Bowtie2 index: clustered_index_${output_prefix}_${clustering_rate}"
        bowtie2-build --threads $threads clustered_sequences_${output_prefix}_${clustering_rate}.fasta clustered_index_${output_prefix}_${clustering_rate}

        echo "Mapping reads with Bowtie2: remapped_reads_${output_prefix}_${clustering_rate}.sam"
        bowtie2 -p $threads -x clustered_index_${output_prefix}_${clustering_rate} -1 $fastq1 -2 $fastq2 -S remapped_reads_${output_prefix}_${clustering_rate}.sam --un-conc unaligned_reads_after_clustering_${output_prefix}_${clustering_rate}#.fastq 2> bowtie2_stats_${output_prefix}_${clustering_rate}.txt
        rm remapped_reads_${output_prefix}_${clustering_rate}.sam
    done
done


# Combine all output to AlignmentStats.txt
{
  cat original_bowtie2_stats.txt
  for clustering_rate in 0.80 0.85 0.90 0.95; do
    for top_regions_file in $output_file1 $output_file2 $output_file3; do
        output_prefix=$(basename $top_regions_file .fasta)
        alignment_rate=$(awk '/^([0-9]+([.][0-9]+)?)% overall alignment rate$/{print $1}' bowtie2_stats_${output_prefix}_${clustering_rate}.txt)
        echo "Clustering file: clustered_sequences_${output_prefix}_${clustering_rate}.fasta"
        echo "Overall alignment rate: $alignment_rate"
    done
  done
} > AlignmentStats.txt

# Clean up intermediate files
# rm original_bowtie2_stats.txt

# Combine all header counts into a single file with filenames
combined_headers_count="combined_headers_count.txt"

# Function to count headers in a FASTA file and append to the combined file
count_fasta_headers() {
    local fasta_file="$1"
    local output_file="$2"

    echo "Counting headers in $fasta_file..."
    echo  "$(basename $fasta_file)\t$(grep -c "^>" $fasta_file)" >> $output_file
}

# Count headers in the original FASTA file and append to the combined file
count_fasta_headers $fasta_file $combined_headers_count

# Count headers in the clustered sequences files and append to the combined file
for clustering_rate in 0.80 0.85 0.90 0.95; do
    for top_regions_file in $output_file1 $output_file2 $output_file3; do
        output_prefix=$(basename $top_regions_file .fasta)
        count_fasta_headers clustered_sequences_${output_prefix}_${clustering_rate}.fasta $combined_headers_count
    done
done


pigz -p $threads *.fastq
echo "Script completed successfully!"
```
