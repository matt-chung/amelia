# mansonella

<!-- MarkdownTOC autolink="true" levels="1,2,3,4" -->

- [Set software and directory paths](#set-software-and-directory-paths)
  - [Software](#software)
  - [Directories](#directories)
  - [Create directories](#create-directories)
- [Identify differential isoform expression between treated D. melanogaster versus non-treated D. melanogaster](#identify-differential-isoform-expression-between-treated-d-melanogaster-versus-non-treated-d-melanogaster)
  - [Set up reference files](#set-up-reference-files)
    - [Download D. melanogaster and wMel reference files](#download-d-melanogaster-and-wmel-reference-files)
    - [Create D. melanogaster and wMel combined reference files](#create-d-melanogaster-and-wmel-combined-reference-files)
  - [Align long reads to combined D. melanogaster and wMel reference](#align-long-reads-to-combined-d-melanogaster-and-wmel-reference)
    - [Align long reads from non-treated and tet-treated](#align-long-reads-from-non-treated-and-tet-treated)
    - [Sort and index BAM files](#sort-and-index-bam-files)
  - [Create annotations for the non-treated and tet-treated samples from long-reads](#create-annotations-for-the-non-treated-and-tet-treated-samples-from-long-reads)
  - [Combine annotations made from the separate long read samples](#combine-annotations-made-from-the-separate-long-read-samples)
  - [Quantitate long-read FASTQs using combined GFF](#quantitate-long-read-fastqs-using-combined-gff)
    - [Create transcript FASTA from combined GFF](#create-transcript-fasta-from-combined-gff)
    - [Align long reads to transcript FASTA](#align-long-reads-to-transcript-fasta)
    - [Quantitate transcripts using Salmon](#quantitate-transcripts-using-salmon)

<!-- /MarkdownTOC -->


# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software

```{bash, eval = F}
GFFCOMPARE_BIN_DIR=/usr/local/packages/gffcompare-0.10.5/bin
GFFREAD_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/gffread_v0.10.4
MINIMAP2_BIN_DIR=/usr/local/packages/minimap2-2.9/bin
SALMON_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/salmon_v1.1.0/bin
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
STRINGTIE_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/stringtie_v2.1.1

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

mkdir -p "$WORKING_DIR"/salmon/nontreated
mkdir -p "$WORKING_DIR"/salmon/tettreated

```

```{bash, eval = F}
## in vitro (20190215)
/usr/local/projects/RDMIN/SEQUENCE/20190215/SINV_IVT/20190215_1549_MN21969_FAJ05343_c9ab6447/fastq_pass/FAJ05343_207c601fce7a926411ae726282c35aed37ce5e1f_0.fastq

## non-treated (20191024 + 20191026)
/local/aberdeen2rw/julie/Matt_dir/amelia/nontreated_native.fastq

## tet-treated (20190307)
/usr/local/projects/RDMIN/SEQUENCE/20190307/JW18_Tet_1/20190307_1511_MN21969_FAK36034_cc9722d6/fastq_pass/FAK36034_1c11934bf425d689c4359d9929334e23470cdcdd_0.fastq
```

# Identify differential isoform expression between treated D. melanogaster versus non-treated D. melanogaster

## Set up reference files

### Download D. melanogaster and wMel reference files

##### Commands:
```{bash, eval = F}
## D. melanogaster
wget -O "$WORKING_DIR"/references/dmelanogaster.fna.gz  ftp://ftp.flybase.net/releases/FB2020_01/dmel_r6.32/fasta/dmel-all-chromosome-r6.32.fasta.gz
wget -O "$WORKING_DIR"/references/dmelanogaster.gff.gz  ftp://ftp.flybase.net/releases/FB2020_01/dmel_r6.32/gff/dmel-all-filtered-r6.32.gff.gz
gunzip "$WORKING_DIR"/references/dmelanogaster.fna.gz
gunzip "$WORKING_DIR"/references/dmelanogaster.gff.gz

## wMel
wget -O "$WORKING_DIR"/references/wMel.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/025/GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/wMel.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/025/GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.gff.gz
gunzip "$WORKING_DIR"/references/wMel.fna.gz
gunzip "$WORKING_DIR"/references/wMel.gff.gz

```

### Create D. melanogaster and wMel combined reference files

##### Commands:
```{bash, eval = F}
cat "$WORKING_DIR"/references/dmelanogaster.fna "$WORKING_DIR"/references/wMel.fna > "$WORKING_DIR"/references/combined_dmelanogaster_wMel.fna
cat "$WORKING_DIR"/references/dmelanogaster.gff "$WORKING_DIR"/references/wMel.gff | grep -v "#" > "$WORKING_DIR"/references/combined_dmelanogaster_wMel.gff

awk '$3 == "gene" || $3 == "exon" {print $0}' "$WORKING_DIR"/references/combined_dmelanogaster_wMel.gff > "$WORKING_DIR"/references/temp.gff
mv "$WORKING_DIR"/references/temp.gff "$WORKING_DIR"/references/combined_dmelanogaster_wMel.gff
```

## Align long reads to combined D. melanogaster and wMel reference

### Align long reads from non-treated and tet-treated
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
echo -e ""$MINIMAP2_BIN_DIR"/minimap2 -ax splice -uf -k14 -t "$THREADS" "$REF_FNA" "$FASTQ" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$OUTPUT_PREFIX".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N minimap2 -wd "$WORKING_DIR"
```

### Sort and index BAM files

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

## Combine annotations made from the separate long read samples

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

## Quantitate long-read FASTQs using combined GFF

### Create transcript FASTA from combined GFF

##### Input Sets:
```{bash, eval = F}
REF_FNA="$WORKING_DIR"/references/combined_dmelanogaster_wMel.fna
GTF="$WORKING_DIR"/stringtie/gffcmp.combined.gtf
```

##### Commands:
```{bash, eval = F}
"$GFFREAD_BIN_DIR"/gffread -w "$(echo "$REF_FNA" | sed "s/[.]fna$/.gene.fna/g")" -g "$REF_FNA" "$WORKING_DIR"/stringtie/gffcmp.combined.gtf
```

### Align long reads to transcript FASTA

##### Input Sets:
```{bash, eval = F}
REF_GENE_FNA="$WORKING_DIR"/references/combined_dmelanogaster_wMel.gene.fna
THREADS=4

## non-treated
OUTPUT_PREFIX="$WORKING_DIR"/bam/nontreated
FASTQ=/local/aberdeen2rw/julie/Matt_dir/amelia/nontreated_native.fastq

## tet-treated
OUTPUT_PREFIX="$WORKING_DIR"/bam/tettreated
FASTQ=/usr/local/projects/RDMIN/SEQUENCE/20190307/JW18_Tet_1/20190307_1511_MN21969_FAK36034_cc9722d6/fastq_pass/FAK36034_1c11934bf425d689c4359d9929334e23470cdcdd_0.fastq
```

##### Commands:
```{bash, eval = F}
echo -e ""$MINIMAP2_BIN_DIR"/minimap2 -ax splice -uf -k14 -t "$THREADS" "$REF_GENE_FNA" "$FASTQ" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$OUTPUT_PREFIX".gene.bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N minimap2 -wd "$WORKING_DIR"
```

### Quantitate transcripts using Salmon

##### Input Sets:
```{bash, eval = F}
REF_GENE_FNA="$WORKING_DIR"/references/combined_dmelanogaster_wMel.gene.fna
THREADS=4

## non-treated
OUTPUT_DIR="$WORKING_DIR"/salmon/nontreated
BAM="$WORKING_DIR"/bam/nontreated.gene.bam

## tet-treated
OUTPUT_DIR="$WORKING_DIR"/salmon/tettreated
BAM="$WORKING_DIR"/bam/tettreated.gene.bam
```

##### Commands:
```{bash, eval = F}
echo -e ""$SALMON_BIN_DIR"/salmon quant -t "$REF_GENE_FNA" --libType A -a "$BAM" -p "$THREADS" -o "$OUTPUT_DIR"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$OUTPUT_DIR"
```
