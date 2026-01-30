## Female and Male m6A peak analysis and visualization
## Including transcriptomic distribution, genomic annotation, and motif logo plotting

rm(list = ls())
options(scipen = 999)
options(stringsAsFactors = FALSE)

## ----------------------------
## Load required packages
## ----------------------------
library(Guitar)
library(GenomicFeatures)
library(ChIPseeker)
library(clusterProfiler)
library(motifStack)
library(data.table)

## ----------------------------
## Build TxDb from GFF3 annotation
## ----------------------------
txdb <- makeTxDbFromGFF(
  file = "schistosoma_japonicum.annotations.gff3",
  format = "gff3",
  dataSource="my",
  organism="Schistosoma japonicum",
  taxonomyId=NA,
  chrominfo=NULL,
  miRBaseBuild=NA,
  metadata=NULL,
)

## ----------------------------
## Load exomePeak2 peak results
## ----------------------------
female_peak_raw <- read.csv("female_peaks.csv")
male_peak_raw   <- read.csv("male_peaks.csv")

## Keep valid chromosomes only
female_peak_raw <- female_peak_raw[female_peak_raw$chr %in% seqlevels(txdb), ]
male_peak_raw   <- male_peak_raw[male_peak_raw$chr %in% seqlevels(txdb), ]

## Filter high-confidence m6A peaks
female_peak_sig <- subset(
  female_peak_raw,
  RPM.IP >= 1 & fdr <= 1e-10 & log2FC >= 1
)

male_peak_sig <- subset(
  male_peak_raw,
  RPM.IP >= 1 & fdr <= 1e-10 & log2FC >= 1
)

female_peak_sig <- na.omit(female_peak_sig)
male_peak_sig   <- na.omit(male_peak_sig)

## ----------------------------
## Export BED3 files for visualization
## ----------------------------
female_bed <- female_peak_sig[, c("chr", "chromStart", "chromEnd")]
male_bed   <- male_peak_sig[, c("chr", "chromStart", "chromEnd")]

write.table(
  female_bed, "female_sig_peaks.bed",
  row.names = FALSE, col.names = FALSE,
  quote = FALSE, sep = "\t"
)

write.table(
  male_bed, "male_sig_peaks.bed",
  row.names = FALSE, col.names = FALSE,
  quote = FALSE, sep = "\t"
)

## ----------------------------
## Transcriptomic distribution (Guitar)
## ----------------------------
GuitarPlot(
  txTxdb = txdb,
  stBedFiles = list("female_sig_peaks.bed", "male_sig_peaks.bed"),
  headOrtail = TRUE,
  enableCI = FALSE,
  mapFilterTranscript = TRUE,
  pltTxType = "mrna",
  stGroupName = c("Female_m6A", "Male_m6A")
) +
  ggtitle("Distribution of m6A peaks on mRNA") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

## ----------------------------
## Genomic annotation (ChIPseeker)
## ----------------------------
female_gr <- readPeakFile("female_sig_peaks.bed")
male_gr   <- readPeakFile("male_sig_peaks.bed")

options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)

female_anno <- annotatePeak(
  female_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron", "Intergenic")
)

male_anno <- annotatePeak(
  male_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron", "Intergenic")
)

plotAnnoPie(female_anno, main = "Female m6A peak annotation")
plotAnnoPie(male_anno,   main = "Male m6A peak annotation")

## ----------------------------
## Export BED4 files for motif analysis
## ----------------------------
female_motif_bed <- female_peak_sig[, c("chr", "chromStart", "chromEnd", "strand")]
female_motif_bed$chromStart <- female_motif_bed$chromStart - 5
female_motif_bed$chromEnd   <- female_motif_bed$chromEnd + 5
female_motif_bed <- na.omit(female_motif_bed)

male_motif_bed <- male_peak_sig[, c("chr", "chromStart", "chromEnd", "strand")]
male_motif_bed$chromStart <- male_motif_bed$chromStart - 5
male_motif_bed$chromEnd   <- male_motif_bed$chromEnd + 5
male_motif_bed <- na.omit(male_motif_bed)

write.table(
  female_motif_bed, "female_sig_peaks_motif.bed",
  row.names = FALSE, col.names = FALSE,
  quote = FALSE, sep = "\t"
)

write.table(
  male_motif_bed, "male_sig_peaks_motif.bed",
  row.names = FALSE, col.names = FALSE,
  quote = FALSE, sep = "\t"
)

## ============================================================
## Motif logo plotting from FASTA sequences
## ============================================================

plot_motif <- function(sequence_file, label) {
  
  sequences <- fread(sequence_file, header = FALSE)[[1]]
  seq_length <- nchar(sequences[1])
  
  ## Initialize count matrix
  count_matrix <- matrix(
    0,
    nrow = 4,
    ncol = seq_length,
    dimnames = list(c("A", "C", "G", "T"), NULL)
  )
  
  ## Count base frequencies
  for (seq in sequences) {
    for (i in seq_len(seq_length)) {
      base <- substr(seq, i, i)
      if (base %in% rownames(count_matrix)) {
        count_matrix[base, i] <- count_matrix[base, i] + 1
      }
    }
  }
  
  ## Remove flanking positions
  trimmed_matrix <- count_matrix[, 3:(ncol(count_matrix) - 2)]
  
  ## Base colors
  base_colors <- c(
    A = "#018001",
    C = "#0000FF",
    G = "#FFA400",
    T = "#FF0000"
  )
  
  motif <- new(
    "pcm",
    mat = trimmed_matrix,
    name = label,
    color = base_colors
  )
  
  plot(motif, ic.scale = FALSE, ylab = "Probability", font = "serif")
  return(trimmed_matrix)
}

## ----------------------------
## Plot motif logos
## ----------------------------
plot_motif("female_sig_peaks.fasta", "Female m6A motif")
plot_motif("male_sig_peaks.fasta",   "Male m6A motif")
