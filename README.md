# amelia

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
- [Identify differential modification in SINV between tet-treated D. melanogaster versus non-treated D. melanogaster](#identify-differential-modification-in-sinv-between-tet-treated-d-melanogaster-versus-non-treated-d-melanogaster)
  - [Set up reference files](#set-up-reference-files-1)
  - [Resquiggle the MinION FAST5 files](#resquiggle-the-minion-fast5-files)
  - [Detect base modifications in the tet-treated and non-treated samples using Tombo de novo mode](#detect-base-modifications-in-the-tet-treated-and-non-treated-samples-using-tombo-de-novo-mode)
  - [Output dampened fraction values for each position from Tombo de novo mode](#output-dampened-fraction-values-for-each-position-from-tombo-de-novo-mode)
  - [Detect base modifications in the tet-treated and non-treated samples using Tombo 5mC alternative model mode](#detect-base-modifications-in-the-tet-treated-and-non-treated-samples-using-tombo-5mc-alternative-model-mode)
  - [Output dampened fraction values for each position from Tombo 5mC alternative model mode](#output-dampened-fraction-values-for-each-position-from-tombo-5mc-alternative-model-mode)

<!-- /MarkdownTOC -->


# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software

```{bash, eval = F}
PYTHON_LIB_PATH=/usr/local/packages/python-3.5/lib

GFFCOMPARE_BIN_DIR=/usr/local/packages/gffcompare-0.10.5/bin
GFFREAD_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/gffread_v0.10.4
MINIMAP2_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/minimap2-2.17_x64-linux
SALMON_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/salmon_v1.1.0/bin
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

mkdir -p "$WORKING_DIR"/salmon/nontreated
mkdir -p "$WORKING_DIR"/salmon/tettreated

```

```{bash, eval = F}
## in vitro (20190215)
/usr/local/projects/RDMIN/SEQUENCE/20190215/SINV_IVT/20190215_1549_MN21969_FAJ05343_c9ab6447/fastq_pass/FAJ05343_207c601fce7a926411ae726282c35aed37ce5e1f_0.fastq

## non-treated (20181024 + 20181026)
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
wget -O "$WORKING_DIR"/references/dmelanogaster.gtf.gz  ftp://ftp.flybase.net/releases/FB2020_01/dmel_r6.32/gtf/dmel-all-r6.32.gtf.gz
gunzip "$WORKING_DIR"/references/dmelanogaster.fna.gz
gunzip "$WORKING_DIR"/references/dmelanogaster.gff.gz
gunzip "$WORKING_DIR"/references/dmelanogaster.gtf.gz

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
echo -e ""$MINIMAP2_BIN_DIR"/minimap2 -ax splice --MD -uf -k14 -t "$THREADS" --secondary=yes "$REF_FNA" "$FASTQ" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$OUTPUT_PREFIX".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N minimap2 -wd "$WORKING_DIR"
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

########3

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

# Identify differential modification in SINV between tet-treated D. melanogaster versus non-treated D. melanogaster

## Set up reference files

## Resquiggle the MinION FAST5 files

##### Input Sets:
```{bash, eval = F}
THREADS=16
REF_TRANSCRIPT_FNA=/local/projects-t3/EBMAL/mchung_dir/amelia/references/GCA_000860545.1_ViralProj15316_genomic.fna

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

## Detect base modifications in the tet-treated and non-treated samples using Tombo de novo mode

##### Input Sets:
```{bash, eval = F}
THREADS=16

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

## Output dampened fraction values for each position from Tombo de novo mode

##### Input Sets:
```{bash, eval = F}
THREADS=16

## non-treated
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/native_fast5
STATS_PREFIX="$WORKING_DIR"/tombo/nontreated
```

##### Commands:
```{bash, eval = F}
echo -e ""$TOMBO_BIN_DIR"/tombo text_output browser_files --fast5-basedirs "$FAST5_DIR" --statistics-filename "$STATS_PREFIX".denovo_base_detection.tombo.stats --file-types dampened_fraction --browser-file-basename "$STATS_PREFIX".denovo_base_detection" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=50G -N tombo_text_output_browser_files -wd "$(dirname "$STATS_PREFIX")"
```

## Detect base modifications in the tet-treated and non-treated samples using Tombo 5mC alternative model mode

##### Input Sets:
```{bash, eval = F}
THREADS=16

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

## Output dampened fraction values for each position from Tombo 5mC alternative model mode

##### Input Sets:
```{bash, eval = F}
THREADS=16

## non-treated
FAST5_DIR=/local/aberdeen2rw/julie/Matt_dir/amelia/native_fast5
STATS_PREFIX="$WORKING_DIR"/tombo/nontreated
```

##### Commands
```{bash, eval = F}
echo -e ""$TOMBO_BIN_DIR"/tombo text_output browser_files --fast5-basedirs "$FAST5_DIR" --statistics-filename "$STATS_PREFIX".5mC_base_detection.5mC.tombo.stats --file-types dampened_fraction --browser-file-basename "$STATS_PREFIX".5mC_base_detection" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=50G -N tombo_text_output_browser_files -wd "$(dirname "$STATS_PREFIX")"
```


```{bash, eval = F}
ANNOTATION_NAME=dmelanogaster_r6.32
GENOME_BUILD=r6.32
GTF="$WORKING_DIR"/references/dmelanogaster.gtf
OUTPUT_DB="$WORKING_DIR"/talon/dmelanogaster_r6.32
```

```{bash, eval = F}
"$TALON_BIN_DIR"/talon_reformat_gtf -gtf "$GTF"

"$TALON_BIN_DIR"/talon_initialize_database --f="$GTF" --g="$GENOME_BUILD" --a="$ANNOTATION_NAME" --o="$OUTPUT_DB"


```

```{bash, eval = F}
CONFIG="$WORKING_DIR"/talon/config.txt
GENOME_BUILD=r6.32
THREADS=4
DB="$WORKING_DIR"/talon/dmelanogaster_r6.32.db
OUTPUT_PREFIX="$WORKING_DIR"/talon/talon

ANNOTATION_NAME="$WORKING_DIR"/references/dmelanogaster_r6.32

GTF="$WORKING_DIR"/references/dmelanogaster.gtf
```

```{bash, eval = F}
qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N talon -wd "$(dirname "$OUTPUT_PREFIX")" -b y "$TALON_BIN_DIR"/talon --f="$CONFIG" --db="$DB" --build="$GENOME_BUILD" --t="$THREADS" --o "$OUTPUT_PREFIX"



```





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
DB="$WORKING_DIR"/talon/dmelanogaster_r6.32.db
ANNOTATION_NAME=dmelanogaster_r6.32
GENOME_BUILD=r6.32
OUTPUT_PREFIX="$WORKING_DIR"/talon/talon

PAIRINGS=/local/projects-t3/EBMAL/mchung_dir/amelia/talon/pairings.csv
```

```{bash, eval = F}
"$TALON_BIN_DIR"/talon_filtered_transcripts --db="$DB" -a "$ANNOTATION_NAME" --b="$GENOME_BUILD" --p="$PAIRINGS" --o="$OUTPUT_PREFIX".csv


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

