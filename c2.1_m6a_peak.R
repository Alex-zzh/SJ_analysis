## m6A peak calling using exomePeak2
## Example workflow for Schistosoma japonicum
## This script is designed for public sharing

rm(list = ls())

## ----------------------------
## Build TxDb from GFF3
## ----------------------------
library(GenomicFeatures)

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
## Load BSgenome
## ----------------------------
library(BSgenome)
library(BSgenome.Schistosoma.japonicum.WormBase.FuDan)

## ----------------------------
## BAM files (female samples as example)
## Replace filenames with your own BAM files
## ----------------------------
FEMALE_IP_BAM <- c(
  "IP_1.bam",
  "IP_2.bam",
  "IP_3.bam"
)

FEMALE_INPUT_BAM <- c(
  "Input_1.bam",
  "Input_2.bam",
  "Input_3.bam"
)

## ----------------------------
## Peak calling with exomePeak2
## ----------------------------
library(exomePeak2)

female_peak <- exomePeak2(
  bam_ip = FEMALE_IP_BAM,
  bam_input = FEMALE_INPUT_BAM,
  txdb = txdb,
  genome = BSgenome.Schistosoma.japonicum.WormBase.FuDan,
  strandness = "unstrand",
  fragment_length = 100,
  p_cutoff = 1e-10,
  parallel = 4,
  mode = "exon",
  bin_size = 25,
  step_size = 25,
  experiment_name = "m6A_peak",
  motif_based = TRUE,
  motif_sequence = "NNNNNANNNNN"
)
