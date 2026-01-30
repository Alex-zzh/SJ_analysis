## Ortholog identity heatmap based on BioMart results
## Example for m6A-related genes in Schistosoma japonicum

rm(list = ls())
options(stringsAsFactors = FALSE)

## ----------------------------
## Load required packages
## ----------------------------
library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(tibble)

## ----------------------------
## Load BioMart ortholog table
## ----------------------------
bio_res <- fread(
  "schistosoma_japonicum.PRJNA520774.WBPS18.orthologs.tsv"
)

## ----------------------------
## Select human orthologs (writers)
## ----------------------------
human_data <- bio_res[grep("Human", bio_res[[2]], ignore.case = TRUE), ]

## Human m6A writers: METTL3, METTL14, WTAP, VIRMA
target_ids <- c(
  "ENSG00000165819",  # METTL3
  "ENSG00000145388",  # METTL14
  "ENSG00000146457",  # WTAP
  "ENSG00000164944"   # VIRMA
)

writers <- human_data[human_data[[3]] %in% target_ids]
writers <- writers[, .SD[which.max(query_identity)], by = ortholog_gene_id]

## Selected S. japonicum genes
selected_sj <- unlist(writers$gene_id)

## ----------------------------
## Species selection
## ----------------------------
target_species <- c(
  "Human", "Mouse", "Drosophila melanogaster",
  "Schistosoma haematobium", "Schistosoma mattheei",
  "Schistosoma bovis", "Schistosoma mansoni",
  "Fasciolopsis buski", "Clonorchis sinensis",
  "Paragonimus westermani", "Schmidtea mediterranea"
)

species_pattern <- paste0("\\b", target_species, "\\b", collapse = "|")

ts_data <- bio_res %>%
  mutate(
    species_match = grepl(species_pattern, ortholog_species_name, ignore.case = TRUE)
  ) %>%
  filter(gene_id %in% selected_sj, species_match) %>%
  select(-species_match)

## Keep best hit per species
ts_data <- ts_data %>%
  mutate(
    gene_id = as.character(gene_id),
    ortholog_species_name = as.character(ortholog_species_name),
    query_identity = as.numeric(query_identity)
  ) %>%
  group_by(gene_id, ortholog_species_name) %>%
  arrange(desc(query_identity), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

## ----------------------------
## Add S. japonicum self-alignment (100%)
## ----------------------------
sj_rows <- data.frame(
  gene_id = c(
    "EWB00_000332",
    "EWB00_007032",
    "EWB00_007521",
    "EWB00_004181"
  ),
  ortholog_species_name = "Schistosoma japonicum",
  query_identity = 100
)

for (col in setdiff(names(ts_data), names(sj_rows))) {
  sj_rows[[col]] <- NA
}

sj_rows <- sj_rows[, names(ts_data)]
ts_data <- rbind(ts_data, sj_rows)

## ----------------------------
## Normalize species names
## ----------------------------
ts_data$ortholog_species_name <- sapply(
  strsplit(ts_data$ortholog_species_name, " (", fixed = TRUE),
  `[`, 1
)

species_map <- c(
  "Human" = "Homo sapiens",
  "Mouse" = "Mus musculus"
)

ts_data <- ts_data %>%
  mutate(
    ortholog_species_name = recode(
      ortholog_species_name,
      !!!species_map
    )
  )

## ----------------------------
## Add manually curated readers (Excel)
## ----------------------------
hnrnpc <- read_excel("hnrnpc.xlsx")
ythdf2 <- read_excel("ythdf2.xlsx")

ts_data <- rbind(ts_data, hnrnpc, ythdf2)

## ----------------------------
## Heatmap ordering
## ----------------------------
row_order <- c(
  "EWB00_007032",
  "EWB00_007521",
  "EWB00_000332",
  "EWB00_004181",
  "EWB00_004879",
  "EWB00_002818"
)

col_order <- c(
  "Schistosoma japonicum",
  "Schistosoma mansoni",
  "Schistosoma bovis",
  "Schistosoma mattheei",
  "Schistosoma haematobium",
  "Paragonimus westermani",
  "Clonorchis sinensis",
  "Fasciolopsis buski",
  "Schmidtea mediterranea",
  "Drosophila melanogaster",
  "Mus musculus",
  "Homo sapiens"
)

## ----------------------------
## Convert to heatmap matrix
## ----------------------------
heatmap_data <- ts_data %>%
  select(gene_id, ortholog_species_name, query_identity) %>%
  filter(gene_id %in% row_order) %>%
  mutate(
    gene_id = factor(gene_id, levels = row_order),
    ortholog_species_name = factor(ortholog_species_name, levels = col_order)
  ) %>%
  group_by(gene_id, ortholog_species_name) %>%
  summarise(query_identity = max(query_identity, na.rm = TRUE), .groups = "drop") %>%
  complete(gene_id, ortholog_species_name) %>%
  pivot_wider(
    names_from = ortholog_species_name,
    values_from = query_identity
  ) %>%
  arrange(match(gene_id, row_order)) %>%
  select(gene_id, all_of(col_order)) %>%
  mutate(across(-gene_id, ~ ifelse(is.na(.), 0, .)))

heatmap_matrix <- heatmap_data %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

## biomart identity for this pair is inaccurate;
## the value was replaced by a manually curated BLASTP result
heatmap_matrix[3, 4] <- 29.39

## ----------------------------
## Vertical heatmap visualization
## ----------------------------
ggplot(
  data = as.data.frame(heatmap_matrix) %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(
      cols = -gene_id,
      names_to = "species",
      values_to = "identity"
    ) %>%
    mutate(
      gene_id = factor(gene_id, levels = row_order),
      species = factor(species, levels = rev(col_order))
    ),
  aes(x = gene_id, y = species, fill = identity)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("white", "#FFE2D8", "#E15844", "#650311"),
    values = rescale(c(0, 20, 25, 100)),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Query Identity (%)",
    guide = guide_colorbar(
      barheight = unit(5, "cm"),
      barwidth = unit(0.6, "cm"),
      frame.colour = "black"
    )
  ) +
  geom_text(
    data = function(d) d[d$identity > 30, ],
    aes(label = round(identity, 1)),
    color = "white",
    size = 3.5,
    fontface = "bold"
  ) +
  geom_text(
    data = function(d) d[d$identity <= 30, ],
    aes(label = round(identity, 1)),
    color = "black",
    size = 3.5,
    fontface = "bold"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14, color = "black"),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 1)
