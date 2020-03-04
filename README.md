# amelia

<!-- MarkdownTOC autolink="true" levels="1,2,3,4" -->

- [Set software and directory paths](#set-software-and-directory-paths)
  - [Software](#software)
  - [Directories](#directories)
  - [Create directories](#create-directories)
  - [Set up reference files](#set-up-reference-files)
    - [Download D. melanogaster, wMel, and SINV reference files](#download-d-melanogaster-wmel-and-sinv-reference-files)
    - [Create FASTA for the spike sequence](#create-fasta-for-the-spike-sequence)
    - [Create D. melanogaster and wMel combined reference files](#create-d-melanogaster-and-wmel-combined-reference-files)
    - [Create D. melanogaster, wMel, SINV, and spike combined reference files](#create-d-melanogaster-wmel-sinv-and-spike-combined-reference-files)
- [Count the number of reads in each sample that mapped to D. melanogaster, wMel, SINV, and the spike sequence](#count-the-number-of-reads-in-each-sample-that-mapped-to-d-melanogaster-wmel-sinv-and-the-spike-sequence)
  - [Count the number of total reads in each sample](#count-the-number-of-total-reads-in-each-sample)
  - [Align long reads to combined D. melanogaster, wMel, SINV, and spike combined reference](#align-long-reads-to-combined-d-melanogaster-wmel-sinv-and-spike-combined-reference)
  - [Find number of reads that map to D. melanogaster, wMel, SINV, and spike](#find-number-of-reads-that-map-to-d-melanogaster-wmel-sinv-and-spike)
    - [D. melanogaster](#d-melanogaster)
    - [wMel](#wmel)
    - [SINV](#sinv)
    - [spike](#spike)
- [Identify differential isoform expression between treated D. melanogaster versus non-treated D. melanogaster](#identify-differential-isoform-expression-between-treated-d-melanogaster-versus-non-treated-d-melanogaster)
  - [Align long reads to combined D. melanogaster and wMel reference](#align-long-reads-to-combined-d-melanogaster-and-wmel-reference)
  - [Sort and index BAM files](#sort-and-index-bam-files)
  - [Create annotations for the non-treated and tet-treated samples from long-reads](#create-annotations-for-the-non-treated-and-tet-treated-samples-from-long-reads)
    - [Run Stringtie as the annotator](#run-stringtie-as-the-annotator)
      - [Use Stringtie to create a new annotation from the long read data](#use-stringtie-to-create-a-new-annotation-from-the-long-read-data)
    - [Combine Stringtie annotations made from the separate long read samples](#combine-stringtie-annotations-made-from-the-separate-long-read-samples)
  - [Run TALON as the annotator](#run-talon-as-the-annotator)
    - [Initialize TALON database with existing annotation](#initialize-talon-database-with-existing-annotation)
    - [Create SAM files](#create-sam-files)
    - [Create TALON config files](#create-talon-config-files)
    - [Use TALON to create a new annotation from the long read data](#use-talon-to-create-a-new-annotation-from-the-long-read-data)
    - [Filter TALON transcripts](#filter-talon-transcripts)
    - [Assess TALON summary file](#assess-talon-summary-file)
    - [Generate TALON abundance file](#generate-talon-abundance-file)
    - [Create TALON GTF file](#create-talon-gtf-file)
- [Identify differential modification in SINV between tet-treated D. melanogaster versus non-treated D. melanogaster](#identify-differential-modification-in-sinv-between-tet-treated-d-melanogaster-versus-non-treated-d-melanogaster)
  - [Resquiggle the MinION FAST5 files](#resquiggle-the-minion-fast5-files)
  - [Detect base modifications in the tet-treated and non-treated samples using Tombo](#detect-base-modifications-in-the-tet-treated-and-non-treated-samples-using-tombo)
    - [Run Tombo in de novo mode for detecting modifications](#run-tombo-in-de-novo-mode-for-detecting-modifications)
    - [Detect modifications using Tombo de novo mode](#detect-modifications-using-tombo-de-novo-mode)
    - [Output dampened fraction values for each position from Tombo de novo mode](#output-dampened-fraction-values-for-each-position-from-tombo-de-novo-mode)
  - [Run Tombo in 5mC alternative model mode](#run-tombo-in-5mc-alternative-model-mode)
    - [Detect modifications using Tombo 5mC alternative model mode](#detect-modifications-using-tombo-5mc-alternative-model-mode)
    - [Output dampened fraction values for each position from Tombo 5mC alternative model mode](#output-dampened-fraction-values-for-each-position-from-tombo-5mc-alternative-model-mode)
  - [Run Tombo in model_sample_compare mode](#run-tombo-in-model_sample_compare-mode)
    - [Detect modifications using Tombo model_sample_compare mode](#detect-modifications-using-tombo-model_sample_compare-mode)
    - [Output dampened fraction values for each position from Tombo model_sample_compare mode](#output-dampened-fraction-values-for-each-position-from-tombo-model_sample_compare-mode)
  - [Print out the C positions in the SINV genome](#print-out-the-c-positions-in-the-sinv-genome)
  - [Plot Tombo results from each run mode to identify modified positions](#plot-tombo-results-from-each-run-mode-to-identify-modified-positions)
    - [Set R inputs](#set-r-inputs)
  - [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo)
  - [Create a plot for each dampened fraction file](#create-a-plot-for-each-dampened-fraction-file)
  - [Create a plot for only the C positions in the dampened fraction file for model_sample_compare](#create-a-plot-for-only-the-c-positions-in-the-dampened-fraction-file-for-model_sample_compare)
  - [Create a plot for each dampened fraction file](#create-a-plot-for-each-dampened-fraction-file-1)
  - [aaaa](#aaaa)

<!-- /MarkdownTOC -->


# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software

```{bash, eval = F}
PYTHON_LIB_PATH=/usr/local/packages/python-3.5/lib

GFFCOMPARE_BIN_DIR=/usr/local/packages/gffcompare-0.10.5/bin
GFFREAD_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/gffread_v0.10.4
MINIMAP2_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/minimap2-2.17_x64-linux
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
STRINGTIE_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/stringtie_v2.1.1
TALON_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin/
TOMBO_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin

SCRIPTS_BIN_DIR=/home/mattchung/scripts
```

## Directories

```{bash, eval = F}
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/amelia
```

## Create directories

```{bash, eval = F}
mkdir -p "$WORKING_DIR"/references
mkdir -p "$WORKING_DIR"/bam
mkdir -p "$WORKING_DIR"/stringtie

mkdir -p "$WORKING_DIR"/plots

mkdir -p "$WORKING_DIR"/talon
```

```{bash, eval = F}
## in vitro (20190215)
/usr/local/projects/RDMIN/SEQUENCE/20190215/SINV_IVT/20190215_1549_MN21969_FAJ05343_c9ab6447/fastq_pass/FAJ05343_207c601fce7a926411ae726282c35aed37ce5e1f_0.fastq

## non-treated (20181024 + 20181026)
/local/aberdeen2rw/julie/Matt_dir/amelia/nontreated_native.fastq

## tet-treated (20190307)
/usr/local/projects/RDMIN/SEQUENCE/20190307/JW18_Tet_1/20190307_1511_MN21969_FAK36034_cc9722d6/fastq_pass/FAK36034_1c11934bf425d689c4359d9929334e23470cdcdd_0.fastq
```

## Set up reference files

### Download D. melanogaster, wMel, and SINV reference files

##### Commands:
```{bash, eval = F}
## D. melanogaster
wget -O "$WORKING_DIR"/references/dmelanogaster.fna.gz  ftp://ftp.flybase.net/releases/FB2020_01/dmel_r6.32/fasta/dmel-all-chromosome-r6.32.fasta.gz
wget -O "$WORKING_DIR"/references/dmelanogaster.gff.gz  ftp://ftp.flybase.net/releases/FB2020_01/dmel_r6.32/gff/dmel-all-filtered-r6.32.gff.gz
wget -O "$WORKING_DIR"/references/dmelanogaster.gtf.gz  ftp://ftp.flybase.net/releases/FB2020_01/dmel_r6.32/gtf/dmel-all-r6.32.gtf.gz
gunzip "$WORKING_DIR"/references/dmelanogaster.fna.gz
gunzip "$WORKING_DIR"/references/dmelanogaster.gff.gz
gunzip "$WORKING_DIR"/references/dmelanogaster.gtf.gz

## wMel
wget -O "$WORKING_DIR"/references/wMel.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/025/GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/wMel.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/025/GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/wMel.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/025/GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.gtf.gz
gunzip "$WORKING_DIR"/references/wMel.fna.gz
gunzip "$WORKING_DIR"/references/wMel.gff.gz
gunzip "$WORKING_DIR"/references/wMel.gtf.gz

## SINV
wget -O "$WORKING_DIR"/references/GCA_000860545.1_ViralProj15316_genomic.fna.gz  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/860/545/GCA_000860545.1_ViralProj15316/GCA_000860545.1_ViralProj15316_genomic.fna.gz
gunzip "$WORKING_DIR"/references/GCA_000860545.1_ViralProj15316_genomic.fna.gz 
```

### Create FASTA for the spike sequence

##### Commands:
```{bash, eval = F}
vim "$WORKING_DIR"/references/yeast_enolase.fasta
```

```{bash, eval = F}
>YHR174W ENO2 SGDID:S000001217, Chr VIII from 451327-452640, Genome Release 64-2-1, Verified ORF, "Enolase II, a phosphopyruvate hydratase; catalyzes conversion of 2-phosphoglycerate to phosphoenolpyruvate during glycolysis and the reverse reaction during gluconeogenesis; expression induced in response to glucose; ENO2 has a paralog, ENO1, that arose from the whole genome duplication"
ATGGCTGTCTCTAAAGTTTACGCTAGATCCGTCTACGACTCCCGTGGTAACCCAACCGTCGAAGTCGAATTAACCACCGAAAAGGGTGTTTTCAGATCCATTGTTCCATCTGGTGCCTCCACCGGTGTCCACGAAGCTTTGGAAATGAGAGATGAAGACAAATCCAAGTGGATGGGTAAGGGTGTTATGAACGCTGTCAACAACGTCAACAACGTCATTGCTGCTGCTTTCGTCAAGGCCAACCTAGATGTTAAGGACCAAAAGGCCGTCGATGACTTCTTGTTGTCTTTGGATGGTACCGCCAACAAGTCCAAGTTGGGTGCTAACGCTATCTTGGGTGTCTCCATGGCCGCTGCTAGAGCCGCTGCTGCTGAAAAGAACGTCCCATTGTACCAACATTTGGCTGACTTGTCTAAGTCCAAGACCTCTCCATACGTTTTGCCAGTTCCATTCTTGAACGTTTTGAACGGTGGTTCCCACGCTGGTGGTGCTTTGGCTTTGCAAGAATTCATGATTGCTCCAACTGGTGCTAAGACCTTCGCTGAAGCCATGAGAATTGGTTCCGAAGTTTACCACAACTTGAAGTCTTTGACCAAGAAGAGATACGGTGCTTCTGCCGGTAACGTCGGTGACGAAGGTGGTGTTGCTCCAAACATTCAAACCGCTGAAGAAGCTTTGGACTTGATTGTTGACGCTATCAAGGCTGCTGGTCACGACGGTAAGGTCAAGATCGGTTTGGACTGTGCTTCCTCTGAATTCTTCAAGGACGGTAAGTACGACTTGGACTTCAAGAACCCAGAATCTGACAAATCCAAGTGGTTGACTGGTGTCGAATTAGCTGACATGTACCACTCCTTGATGAAGAGATACCCAATTGTCTCCATCGAAGATCCATTTGCTGAAGATGACTGGGAAGCTTGGTCTCACTTCTTCAAGACCGCTGGTATCCAAATTGTTGCTGATGACTTGACTGTCACCAACCCAGCTAGAATTGCTACCGCCATCGAAAAGAAGGCTGCTGACGCTTTGTTGTTGAAGGTTAACCAAATCGGTACCTTGTCTGAATCCATCAAGGCTGCTCAAGACTCTTTCGCTGCCAACTGGGGTGTTATGGTTTCCCACAGATCTGGTGAAACTGAAGACACTTTCATTGCTGACTTGGTTGTCGGTTTGAGAACTGGTCAAATCAAGACTGGTGCTCCAGCTAGATCCGAAAGATTGGCTAAGTTGAACCAATTGTTGAGAATCGAAGAAGAATTGGGTGACAAGGCTGTCTACGCCGGTGAAAACTTCCACCACGGTGACAAGTTGTAA
```

### Create D. melanogaster and wMel combined reference files

##### Commands:
```{bash, eval = F}
cat "$WORKING_DIR"/references/dmelanogaster.fna "$WORKING_DIR"/references/wMel.fna > "$WORKING_DIR"/references/combined_dmelanogaster_wMel.fna
cat "$WORKING_DIR"/references/dmelanogaster.gff "$WORKING_DIR"/references/wMel.gff | grep -v "#" > "$WORKING_DIR"/references/combined_dmelanogaster_wMel.gff
cat "$WORKING_DIR"/references/dmelanogaster.gff "$WORKING_DIR"/references/wMel.gtf | grep -v "#" > "$WORKING_DIR"/references/combined_dmelanogaster_wMel.gtf


awk '$3 == "gene" || $3 == "exon" {print $0}' "$WORKING_DIR"/references/combined_dmelanogaster_wMel.gff > "$WORKING_DIR"/references/temp.gff
mv "$WORKING_DIR"/references/temp.gff "$WORKING_DIR"/references/combined_dmelanogaster_wMel.gff
```

### Create D. melanogaster, wMel, SINV, and spike combined reference files


```{bash, eval = F}
cat "$WORKING_DIR"/references/dmelanogaster.fna "$WORKING_DIR"/references/wMel.fna "$WORKING_DIR"/references/GCA_000860545.1_ViralProj15316_genomic.fna "$WORKING_DIR"/references/yeast_enolase.fasta > "$WORKING_DIR"/references/combined.fna
```

# Count the number of reads in each sample that mapped to D. melanogaster, wMel, SINV, and the spike sequence

## Count the number of total reads in each sample

##### Input Sets:
```{bash, eval = F}
## in vitro
FASTQ=/usr/local/projects/RDMIN/SEQUENCE/20190215/SINV_IVT/20190215_1549_MN21969_FAJ05343_c9ab6447/fastq_pass/FAJ05343_207c601fce7a926411ae726282c35aed37ce5e1f_0.fastq

## non-treated
FASTQ=/local/aberdeen2rw/julie/Matt_dir/amelia/nontreated_native.fastq

## tet-treated
FASTQ=/usr/local/projects/RDMIN/SEQUENCE/20190307/JW18_Tet_1/20190307_1511_MN21969_FAK36034_cc9722d6/fastq_pass/FAK36034_1c11934bf425d689c4359d9929334e23470cdcdd_0.fastq
```

##### Commands:
```{bash, eval = F}
echo $(($(wc -l $FASTQ | awk '{print $1}')/4))
```

in vitro: 12811 
non-treated: 448747   
tet-treated: 304310  

## Align long reads to combined D. melanogaster, wMel, SINV, and spike combined reference

##### Input Sets:
```{bash, eval = F}
REF_GENE_FNA=/local/projects-t3/EBMAL/mchung_dir/amelia/references/combined.fna
THREADS=4

## in vitro
OUTPUT_PREFIX="$WORKING_DIR"/bam/invitro_all
FASTQ=/usr/local/projects/RDMIN/SEQUENCE/20190215/SINV_IVT/20190215_1549_MN21969_FAJ05343_c9ab6447/fastq_pass/FAJ05343_207c601fce7a926411ae726282c35aed37ce5e1f_0.fastq

## non-treated
OUTPUT_PREFIX="$WORKING_DIR"/bam/nontreated_all
FASTQ=/local/aberdeen2rw/julie/Matt_dir/amelia/nontreated_native.fastq

## tet-treated
OUTPUT_PREFIX="$WORKING_DIR"/bam/tettreated_all
FASTQ=/usr/local/projects/RDMIN/SEQUENCE/20190307/JW18_Tet_1/20190307_1511_MN21969_FAK36034_cc9722d6/fastq_pass/FAK36034_1c11934bf425d689c4359d9929334e23470cdcdd_0.fastq
```

##### Commands:
```{bash, eval = F}
echo -e ""$MINIMAP2_BIN_DIR"/minimap2 -ax splice -uf -k14 -t "$THREADS" "$REF_GENE_FNA" "$FASTQ" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$OUTPUT_PREFIX".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N minimap2 -wd "$WORKING_DIR"
```

## Find number of reads that map to D. melanogaster, wMel, SINV, and spike

##### Input Sets:
```{bash, eval = F}
## in vitro
BAM="$WORKING_DIR"/bam/invitro_all.bam

## non-treated
BAM="$WORKING_DIR"/bam/nontreated_all.bam

## tet-treated
BAM="$WORKING_DIR"/bam/tettreated_all.bam
```

### D. melanogaster
##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools view -F4 "$BAM" | grep -v NC_002978.6 | grep -v J02363.1 | grep -v YHR174W | cut -f1 | sort -n | uniq | wc -l
```

in vitro: 0  
non-treated: 419057  
tet-treated: 264036  

### wMel
##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools view -F4 "$BAM" | grep NC_002978.6 | cut -f1 | sort -n | uniq | wc -l
```

in vitro: 0  
non-treated: 219  
tet-treated: 135  

### SINV
##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools view -F4 "$BAM" | grep J02363.1 | cut -f1 | sort -n | uniq | wc -l
```

in vitro: 521  
non-treated: 540  
tet-treated: 1  

### spike
##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools view -F4 "$BAM" | grep YHR174W | cut -f1 | sort -n | uniq | wc -l
```

in vitro: 11512  
non-treated: 10023  
tet-treated: 4008   


##### Input Sets:
```{bash, eval = F}
REF_GENE_FNA=/local/projects-t3/EBMAL/mchung_dir/amelia/references/combined.fna
THREADS=4

## non-treated
OUTPUT_PREFIX="$WORKING_DIR"/bam/invitro_all
FASTQ=/usr/local/projects/RDMIN/SEQUENCE/20190215/SINV_IVT/20190215_1549_MN21969_FAJ05343_c9ab6447/fastq_pass/FAJ05343_207c601fce7a926411ae726282c35aed37ce5e1f_0.fastq

## non-treated
OUTPUT_PREFIX="$WORKING_DIR"/bam/nontreated_all
FASTQ=/local/aberdeen2rw/julie/Matt_dir/amelia/nontreated_native.fastq

## tet-treated
OUTPUT_PREFIX="$WORKING_DIR"/bam/tettreated_all
FASTQ=/usr/local/projects/RDMIN/SEQUENCE/20190307/JW18_Tet_1/20190307_1511_MN21969_FAK36034_cc9722d6/fastq_pass/FAK36034_1c11934bf425d689c4359d9929334e23470cdcdd_0.fastq
```

# Identify differential isoform expression between treated D. melanogaster versus non-treated D. melanogaster

## Align long reads to combined D. melanogaster and wMel reference

##### Input Sets:

```{bash, eval = F}
REF_FNA="$WORKING_DIR"/references/combined_dmelanogaster_wMel.fna
THREADS=16

## non-treated
OUTPUT_PREFIX="$WORKING_DIR"/bam/nontreated
FASTQ=/local/aberdeen2rw/julie/Matt_dir/amelia/nontreated_native.fastq

## tet-treated
OUTPUT_PREFIX="$WORKING_DIR"/bam/tettreated
FASTQ=/usr/local/projects/RDMIN/SEQUENCE/20190307/JW18_Tet_1/20190307_1511_MN21969_FAK36034_cc9722d6/fastq_pass/FAK36034_1c11934bf425d689c4359d9929334e23470cdcdd_0.fastq
```

```{bash, eval = F}
echo -e ""$MINIMAP2_BIN_DIR"/minimap2 -ax splice --MD -uf -k14 -t "$THREADS" --secondary=yes "$REF_FNA" "$FASTQ" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$OUTPUT_PREFIX".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N minimap2 -wd "$WORKING_DIR"
```

## Sort and index BAM files

##### Input Sets:
```{bash, eval = F}
THREADS=4

## non-treated
BAM="$WORKING_DIR"/bam/nontreated.bam

## tet-treated
BAM="$WORKING_DIR"/bam/tettreated.bam
```

##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools sort -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/g")" "$BAM"
"$SAMTOOLS_BIN_DIR"/samtools index -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/g")"
```

## Create annotations for the non-treated and tet-treated samples from long-reads

### Run Stringtie as the annotator

#### Use Stringtie to create a new annotation from the long read data

##### Input Sets:
```{bash, eval = F}
THREADS=4

## non-treated
FEATURE_PREFIX=non-treated
BAM="$WORKING_DIR"/bam/nontreated.sortedbyposition.bam

## tet-treated
FEATURE_PREFIX=tet-treated
BAM="$WORKING_DIR"/bam/tettreated.sortedbyposition.bam
```

##### Commands:
```{bash, eval = F}
"$STRINGTIE_BIN_DIR"/stringtie "$BAM" -l "$FEATURE_PREFIX" -o "$WORKING_DIR"/stringtie/"$(basename "$BAM" | sed "s/[.]bam$/.gtf/g")" -p "$THREADS"
```

### Combine Stringtie annotations made from the separate long read samples

##### Input Sets:
```{bash, eval = F}
OUTPUT_PREFIX="$WORKING_DIR"/stringtie/gffcmp
REF_GFF="$WORKING_DIR"/references/combined_dmelanogaster_wMel.gff
GFF1="$WORKING_DIR"/stringtie/nontreated.sortedbyposition.gtf
GFF2="$WORKING_DIR"/stringtie/tettreated.sortedbyposition.gtf
```

##### Commands:
```{bash, eval = F}
"$GFFCOMPARE_BIN_DIR"/gffcompare -C "$GFF1" "$GFF2" -R "$REF_GFF" -o "$OUTPUT_PREFIX"
```

## Run TALON as the annotator

### Initialize TALON database with existing annotation

##### Inputs
```{bash, eval = F}
ANNOTATION_NAME=dmelanogaster
GENOME_BUILD=r6.32
GTF="$WORKING_DIR"/references/dmelanogaster.gtf
OUTPUT_DB="$WORKING_DIR"/talon/dmelanogaster
```

##### Commands
```{bash, eval = F}
"$TALON_BIN_DIR"/talon_initialize_database --f="$GTF" --g="$GENOME_BUILD" --a="$ANNOTATION_NAME" --o="$OUTPUT_DB"
```

### Create SAM files 

##### Input Sets
```{bash, eval = F}
THREADS=4

## non-treated
BAM=/local/projects-t3/EBMAL/mchung_dir/amelia/bam/nontreated.bam

## tet-treated
BAM=/local/projects-t3/EBMAL/mchung_dir/amelia/bam/nontreated.bam
```

##### Commands
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools view -@ 4 -ho "$(echo "$BAM" | sed "s/[.]sam$//g")" "$BAM"
```

### Create TALON config files

TALON config file is a csv consisting of the sample ID, sample type, sequencing instrument, and SAM path

##### Inputs
```{bash, eval = F}
vim "$WORKING_DIR"/talon/config.txt
```

##### Commands
```{bash, eval = F}
nontreated,nontreated,Oxford_Nanopore_MinION,/local/projects-t3/EBMAL/mchung_dir/amelia/bam/nontreated.sam
tettreated,tettreated,Oxford_Nanopore_MinION,/local/projects-t3/EBMAL/mchung_dir/amelia/bam/tettreated.sam
```

### Use TALON to create a new annotation from the long read data

##### Inputs
```{bash, eval = F}
THREADS=4
CONFIG="$WORKING_DIR"/talon/config.txt
GENOME_BUILD=r6.32
DB="$WORKING_DIR"/talon/dmelanogaster.db
OUTPUT_PREFIX="$WORKING_DIR"/talon/talon
```

##### Commands
```{bash, eval = F}
qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N talon -wd "$(dirname "$OUTPUT_PREFIX")" -b y "$TALON_BIN_DIR"/talon --f="$CONFIG" --db="$DB" --build="$GENOME_BUILD" --t="$THREADS" --o "$OUTPUT_PREFIX"
```

### Filter TALON transcripts 


```{bash, eval = F}
DB="$WORKING_DIR"/talon/dmelanogaster.db
ANNOTATION_NAME=dmelanogaster
GENOME_BUILD=r6.32
PAIRINGS=/local/projects-t3/EBMAL/mchung_dir/amelia/talon/pairings.csv
OUTPUT_PREFIX="$WORKING_DIR"/talon/talon
```

```{bash, eval = F}
"$TALON_BIN_DIR"/talon_filtered_transcripts --db="$DB" -a "$ANNOTATION_NAME" --b="$GENOME_BUILD" --p="$PAIRINGS" --o="$OUTPUT_PREFIX".csv
```

### Assess TALON summary file

### Generate TALON abundance file

### Create TALON GTF file

# Identify differential modification in SINV between tet-treated D. melanogaster versus non-treated D. melanogaster

## Resquiggle the MinION FAST5 files

##### Input Sets:
```{bash, eval = F}
THREADS=16
REF_TRANSCRIPT_FNA="$WORKING_DIR"/references/GCA_000860545.1_ViralProj15316_genomic.fna

## in vitro
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/invitro_fast5
FASTQ_FILE=/usr/local/projects/RDMIN/SEQUENCE/20190215/SINV_IVT/20190215_1549_MN21969_FAJ05343_c9ab6447/fastq_pass/FAJ05343_207c601fce7a926411ae726282c35aed37ce5e1f_0.fastq

## non-treated
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/native_fast5
FASTQ_FILE=/local/aberdeen2rw/julie/Matt_dir/amelia/nontreated_native.fastq

## tet-treated
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/tet_native_fast5
FASTQ_FILE=/usr/local/projects/RDMIN/SEQUENCE/20190307/JW18_Tet_1/20190307_1511_MN21969_FAK36034_cc9722d6/fastq_pass/FAK36034_1c11934bf425d689c4359d9929334e23470cdcdd_0.fastq
```

##### Commands:
```{bash, eval = F}
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_PATH":"$LD_LIBRARY_PATH"\n"$TOMBO_BIN_DIR"/tombo resquiggle --overwrite --rna "$FAST5_DIR" "$REF_TRANSCRIPT_FNA" --processes "$THREADS" --num-most-common-errors 5" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N tombo_resquiggle -wd "$FAST5_DIR"
```

in vitro:
```{bash, eval = F}
[15:45:22] Final unsuccessful reads summary (96.3% reads unsuccessfully processed; 12342 total reads):
    96.1% (  12312 reads) : Alignment not produced
     0.2% (     28 reads) : Poor raw to expected signal matching (revert with `tombo filter clear_filters`)
     0.0% (      1 reads) : Read event to sequence alignment extends beyond bandwidth
     0.0% (      1 reads) : Base calls not found in FAST5 (see `tombo preprocess`)
```

non-treated:
```{bash, eval = F}
[21:33:51] Final unsuccessful reads summary (99.9% reads unsuccessfully processed; 486771 total reads):
    99.9% ( 486681 reads) : Alignment not produced
     0.0% (     84 reads) : Poor raw to expected signal matching (revert with `tombo filter clear_filters`)
     0.0% (      6 reads) : Read event to sequence alignment extends beyond bandwidth
```

tet-treated
```{bash, eval = F}
[18:45:55] Final unsuccessful reads summary (100.0% reads unsuccessfully processed; 304309 total reads):
   100.0% ( 304309 reads) : Alignment not produced
```

## Detect base modifications in the tet-treated and non-treated samples using Tombo 

Tombo was only run for the in vitro and non-treated data sets because the tet-treated data set had only 1 read mapping to SINV.

### Run Tombo in de novo mode for detecting modifications

### Detect modifications using Tombo de novo mode

##### Input Sets:
```{bash, eval = F}
THREADS=16

## in vitro
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/invitro_fast5
FASTQ_FILE="$WORKING_DIR"/tombo/invitro

## non-treated
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/native_fast5
STATS_PREFIX="$WORKING_DIR"/tombo/nontreated
```

##### Commands:
```{bash, eval = F}
echo -e ""$TOMBO_BIN_DIR"/tombo detect_modifications de_novo \
   --fast5-basedirs "$FAST5_DIR" \
   --per-read-statistics-basename "$STATS_PREFIX".denovo_base_detection \
   --statistics-file-basename "$STATS_PREFIX".denovo_base_detection \
   --processes "$THREADS"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N tombo_detect_modifications_de_novo -wd "$FAST5_DIR"
```

### Output dampened fraction values for each position from Tombo de novo mode

##### Input Sets:
```{bash, eval = F}
THREADS=16

## in vitro
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/invitro_fast5
FASTQ_FILE="$WORKING_DIR"/tombo/invitro

## non-treated
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/native_fast5
STATS_PREFIX="$WORKING_DIR"/tombo/nontreated
```

##### Commands:
```{bash, eval = F}
echo -e ""$TOMBO_BIN_DIR"/tombo text_output browser_files --fast5-basedirs "$FAST5_DIR" --statistics-filename "$STATS_PREFIX".denovo_base_detection.tombo.stats --file-types dampened_fraction statistic --browser-file-basename "$STATS_PREFIX".denovo_base_detection" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=50G -N tombo_text_output_browser_files -wd "$(dirname "$STATS_PREFIX")"
```

## Run Tombo in 5mC alternative model mode

### Detect modifications using Tombo 5mC alternative model mode

##### Input Sets:
```{bash, eval = F}
THREADS=16

## in vitro
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/invitro_fast5
FASTQ_FILE="$WORKING_DIR"/tombo/invitro

## non-treated
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/native_fast5
STATS_PREFIX="$WORKING_DIR"/tombo/nontreated
```

##### Commands:
```{bash, eval = F}
echo -e ""$TOMBO_BIN_DIR"/tombo detect_modifications alternative_model --fast5-basedirs "$FAST5_DIR" \
  --statistics-file-basename "$STATS_PREFIX".5mC_base_detection \
  --per-read-statistics-basename "$STATS_PREFIX".5mC_base_detection \
  --alternate-bases 5mC \
  --processes "$THREADS"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N tombo_detect_modifications -wd "$FAST5_DIR"
```

### Output dampened fraction values for each position from Tombo 5mC alternative model mode

##### Input Sets:
```{bash, eval = F}
THREADS=16

## in vitro
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/invitro_fast5
FASTQ_FILE="$WORKING_DIR"/tombo/invitro

## non-treated
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/native_fast5
STATS_PREFIX="$WORKING_DIR"/tombo/nontreated
```

##### Commands
```{bash, eval = F}
echo -e ""$TOMBO_BIN_DIR"/tombo text_output browser_files --fast5-basedirs "$FAST5_DIR" --statistics-filename "$STATS_PREFIX".5mC_base_detection.5mC.tombo.stats --file-types dampened_fraction statistic --browser-file-basename "$STATS_PREFIX".5mC_base_detection" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=50G -N tombo_text_output_browser_files -wd "$(dirname "$STATS_PREFIX")"
```

## Run Tombo in model_sample_compare mode

### Detect modifications using Tombo model_sample_compare mode

##### Input Sets:
```{bash, eval = F}
THREADS=16

## non-treated
FAST5_DIR1=/local/aberdeen2rw/julie/Matt_dir/amelia/native_fast5
FAST5_DIR2=/local/aberdeen2rw/julie/Matt_dir/amelia/invitro_fast5

STATS_PREFIX="$WORKING_DIR"/tombo/model_sample_compare
```

##### Commands:
```{bash, eval = F}
echo -e ""$TOMBO_BIN_DIR"/tombo detect_modifications model_sample_compare --fast5-basedirs "$FAST5_DIR1" --control-fast5-basedirs "$FAST5_DIR2" \
  --statistics-file-basename "$STATS_PREFIX" \
  --processes "$THREADS"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N tombo_detect_modifications -wd "$(dirname "$STATS_PREFIX")"
```

### Output dampened fraction values for each position from Tombo model_sample_compare mode

##### Input Sets:
```{bash, eval = F}
STATS_PREFIX="$WORKING_DIR"/tombo/model_sample_compare
```
##### Commands:
```{bash, eval = F}
echo -e ""$TOMBO_BIN_DIR"/tombo text_output browser_files --statistics-filename "$STATS_PREFIX".tombo.stats --file-types dampened_fraction statistic --browser-file-basename "$STATS_PREFIX"" | qsub -P jdhotopp-lab -l mem_free=50G -N tombo_text_output_browser_files -wd "$(dirname "$STATS_PREFIX")"
```

## Print out the C positions in the SINV genome

##### Inputs:
```{bash, eval = F}
REF_FNA="$WORKING_DIR"/references/combined.fna
BAM="$WORKING_DIR"/bam/nontreatd_all.bam
```

##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools mpileup -aa -d 1000000 -f "$REF_FNA" --reference "$REF_FNA" "$BAM" | grep J02363.1 | cut -f2,3 | sort -n | awk '$2 == "C" {print $0}' > "$WORKING_DIR"/ref_c_positions.list
```

## Plot Tombo results from each run mode to identify modified positions

### Set R inputs
```{R}
WORKING_DIR="Z:/EBMAL/mchung_dir/amelia/"

DAMPENED_FRACTION.PATH.LIST <- c(paste0(WORKING_DIR,"/tombo/nontreated.denovo_base_detection.dampened_fraction_modified_reads.plus.wig"),
                                 paste0(WORKING_DIR,"/tombo/nontreated.5mC_base_detection.dampened_fraction_modified_reads.plus.wig"),
                                 paste0(WORKING_DIR,"/tombo/invitro.denovo_base_detection.dampened_fraction_modified_reads.plus.wig"),
                                 paste0(WORKING_DIR,"/tombo/invitro.5mC_base_detection.dampened_fraction_modified_reads.plus.wig"),
                                 paste0(WORKING_DIR,"/tombo/model_sample_compare.dampened_fraction_modified_reads.plus.wig"))
REF_C_POSITIONS.PATH <- "Z:/EBMAL/mchung_dir/amelia/ref_c_positions.list"
```

## Load R packages and view sessionInfo

```{R}
library(cowplot)
library(ggplot2)
sessionInfo()
```

```{R, eval = F}
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] cowplot_0.9.4 ggplot2_3.2.1

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1       withr_2.1.2      assertthat_0.2.1 dplyr_0.8.3      crayon_1.3.4     R6_2.4.0         grid_3.5.1       gtable_0.3.0    
 [9] magrittr_1.5     scales_1.0.0     pillar_1.3.1     rlang_0.4.0      lazyeval_0.2.2   rstudioapi_0.10  tools_3.5.1      glue_1.3.1      
[17] purrr_0.3.2      munsell_0.5.0    xfun_0.6         yaml_2.2.0       compiler_3.5.1   pkgconfig_2.0.2  colorspace_1.4-1 tidyselect_0.2.5
[25] knitr_1.22       tibble_2.1.1
```

## Create a plot for each dampened fraction file

```{R}
ref_c_position <- read.delim(REF_C_POSITIONS.PATH)

plot.list <- list()
for(i in 1:length(DAMPENED_FRACTION.PATH.LIST)){
  dampened_fraction <- read.delim(DAMPENED_FRACTION.PATH.LIST[i], sep = " ")
  
  xcoords <- seq(1,11703)
  ycoords <- as.numeric(as.character(dampened_fraction[match(xcoords,dampened_fraction[,1]),2]))
  ycoords[is.na(ycoords)] <- 0
  
  plot.list[[i]] <- ggplot(mapping=aes_string(x=xcoords,y=ycoords))+
    geom_line()+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_cartesian(ylim=c(0.75,1))+
    labs(x="SINV position", y="dampened fraction")+
    theme_bw()
}
```

## Create a plot for only the C positions in the dampened fraction file for model_sample_compare

```{R}
ref_c_position <- read.delim(REF_C_POSITIONS.PATH)
ycoords[!(seq(1,11703) %in% ref_c_position[,1])] <- 0
plot.list[[6]] <- ggplot(mapping=aes_string(x=xcoords,y=ycoords))+
  geom_line()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0.75,1))+
  labs(x="SINV position", y="dampened fraction")+
  theme_bw()
```

## Create a plot for each dampened fraction file
```{R, fig.height = 8, fig.width = 10}
pdf(paste0(WORKING_DIR,"/plots/tombo_plot_raw.pdf"),
    height=8,
    width=10)
plot_grid(plotlist = plot.list, nrow = 3, ncol = 2, labels = c("A","B","C","D","E","F"))
dev.off()

png(paste0(WORKING_DIR,"/plots/tombo_plot_raw.png"),
    height=8,
    width=10,
    units = "in",res=300)
plot_grid(plotlist = plot.list, nrow = 3, ncol = 2, labels = c("A","B","C","D","E","F"))
dev.off()

plot_grid(plotlist = plot.list, nrow = 3, ncol = 2, labels = c("A","B","C","D","E","F"))
```

![Image description](/images/tombo_plot_raw.png)

A: nontreated - de novo base detection  
B: nontreated - 5mC alternative model base detection  
C: in vitro - de novo base detection  
D: in vitro - 5mC alternative model base detection  
E: nontreated v. in vitro - model sample compare  
F: nontreated v. in vitro - model sample compare C positions  

## aaaa





```{bash, eval = F}
"$TALON_BIN_DIR"/talon_summarize --db="$DB" --v --o test
```

```{bash, eval = F}
---------------nontreated---------------
Number of reads: 202829
Known genes: 6705
Novel genes: 1086
----Antisense genes: 556
----Other novel genes: 530
Known transcripts: 5355
Novel transcripts: 21111
----ISM transcripts: 8828
--------Prefix ISM transcripts: 913
--------Suffix ISM transcripts: 5967
----NIC transcripts: 2434
----NNC transcripts: 5167
----antisense transcripts: 524
----genomic transcripts: 3674
---------------tettreated---------------
Number of reads: 190386
Known genes: 6651
Novel genes: 1108
----Antisense genes: 570
----Other novel genes: 538
Known transcripts: 4731
Novel transcripts: 19174
----ISM transcripts: 8059
--------Prefix ISM transcripts: 701
--------Suffix ISM transcripts: 5141
----NIC transcripts: 2043
----NNC transcripts: 4313
----antisense transcripts: 526
----genomic transcripts: 3754
```






```{bash, eval = F}
"$TALON_BIN_DIR"/talon_abundance --db="$DB" -a "$ANNOTATION_NAME" --b="$GENOME_BUILD" --whitelist "$OUTPUT_PREFIX".csv --o="$OUTPUT_PREFIX"
```

```{bash, eval = F}
"$TALON_BIN_DIR"/talon_create_GTF --db="$DB" -a "$ANNOTATION_NAME" --b="$GENOME_BUILD" --whitelist "$OUTPUT_PREFIX".csv --o="$OUTPUT_PREFIX"
```

"$SAMTOOLS_BIN_DIR"/samtools view -F16 -o "$(echo "$BAM" | sed "s/[.]bam$/.f.bam/g")" "$BAM"
"$SAMTOOLS_BIN_DIR"/samtools view -f16 -o "$(echo "$BAM" | sed "s/[.]bam$/.r.bam/g")" "$BAM"


awk -F "\t" '$3 == "gene" {print $1"\t"$4"\t"$5"\t"$9}' "$GTF" > "$WORKING_DIR"/splice_index_analysis/talon_gene.bed
"$SAMTOOLS_BIN_DIR"/samtools depth -aa -d 1000000 -b "$WORKING_DIR"/splice_index_analysis/talon_gene.bed -f "$WORKING_DIR"/splice_index_analysis/bam.list > "$WORKING_DIR"/splice_index_analysis/gene.depth



WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/amelia
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin

GTF="$WORKING_DIR"/talon/talon_talon.gtf

rm "$WORKING_DIR"/splice_index_analysis/avg_exon2.depth
awk -F "\t" '$3 == "exon" {print $9}' "$GTF" | sed "s/.* exon_id .//g" | sed "s/.;.*//g" | sort -n | uniq | while read exon
do
  chr="$(grep " exon_id .$exon.;" "$GTF" | awk -F "\t" '$3 == "exon" {print $1}' | uniq)"
  start="$(grep " exon_id .$exon.;" "$GTF" | awk -F "\t" '$3 == "exon" {print $4}' | uniq)"
  stop="$(grep " exon_id .$exon.;" "$GTF" | awk -F "\t" '$3 == "exon" {print $5}' | uniq)"

  "$SAMTOOLS_BIN_DIR"/samtools depth -aa -d 1000000 -r "$chr":"$start"-"$stop" -f "$WORKING_DIR"/splice_index_analysis/bam.list > "$WORKING_DIR"/splice_index_analysis/gene.subset2.depth

  echo -e ""$exon"\t"$(awk '{ sum += $3 } END { if (NR > 0) print sum / NR }' "$WORKING_DIR"/splice_index_analysis/gene.subset2.depth)"""\t"$(awk '{ sum += $4 } END { if (NR > 0) print sum / NR }' "$WORKING_DIR"/splice_index_analysis/gene.subset2.depth)"""\t"$(awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' "$WORKING_DIR"/splice_index_analysis/gene.subset2.depth)"""\t"$(awk '{ sum += $6 } END { if (NR > 0) print sum / NR }' "$WORKING_DIR"/splice_index_analysis/gene.subset2.depth)"""" >> "$WORKING_DIR"/splice_index_analysis/avg_exon2.depth
done
qsub -P jdhotopp-lab -l mem_free=5G -N splice_index_analysis -wd "$WORKING_DIR"/splice_index_analysis -b y "$WORKING_DIR"/splice_index_analysis/calc_exon_depth.sh

qsub -P jdhotopp-lab -l mem_free=5G -N splice_index_analysis -wd "$WORKING_DIR"/splice_index_analysis -b y "$WORKING_DIR"/splice_index_analysis/calc_exon_depth2.sh

