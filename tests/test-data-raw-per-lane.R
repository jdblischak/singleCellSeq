library("testit")
library("data.table")
suppressPackageStartupMessages(library("dplyr"))

# Test existence of raw data files ---------------------------------------------

assert(
  "Raw data files exist",
  file.exists("../data/reads-raw-bulk-per-lane.txt"),
  file.exists("../data/reads-raw-single-per-lane.txt"),
  file.exists("../data/molecules-raw-single-per-lane.txt")
)

# Import raw data files --------------------------------------------------------

reads_raw <- fread("../data/reads-raw-single-per-lane.txt")
setDF(reads_raw)

molecules_raw <- fread("../data/molecules-raw-single-per-lane.txt")
setDF(molecules_raw)

reads_bulk_raw <- fread("../data/reads-raw-bulk-per-lane.txt")
setDF(reads_bulk_raw)

# Test dimensions --------------------------------------------------------------

# exected number of lanes: 3 individuals * 3 reps * 96 wells * 3 lanes = 2,592
assert(
  "Reads and molecules for single cells contain 2,592 lanes (rows)",
  nrow(reads_raw) == 2592,
  nrow(molecules_raw) == 2592)

# expected number of lanes: 3 individuals * 3 reps * 2 indexes * 4 lanes = 72
assert(
  "Bulk reads data has 72 lanes (rows)",
  nrow(reads_bulk_raw) == 72)

assert(
  "Per lane files have the same number of columns",
  ncol(reads_bulk_raw) == ncol(reads_raw),
  ncol(reads_raw) == ncol(molecules_raw)
)

# Test column names ------------------------------------------------------------

assert(
  "Per lane files contain the same column names",
  colnames(reads_bulk_raw) == colnames(reads_raw),
  colnames(reads_raw) == colnames(molecules_raw)
)

assert(
  "Columns include Ensembl gene names",
  any(grep("ENSG", colnames(reads_raw)))
)

assert(
  "Columns include ERCC gene names",
  any(grep("ERCC", colnames(reads_raw)))
)

exons <- fread("../data/exons.saf")
setDF(exons)
gene_names_uniq <- unique(exons$GeneID)

assert(
  "All genes in exons.saf present",
  gene_names_uniq %in% colnames(reads_raw)
)

# Test annotation---------------------------------------------------------------

assert(
  "Single cell samples are ordered the same in reads and molecules files",
  reads_raw$individual == molecules_raw$individual,
  reads_raw$replicate == molecules_raw$replicate,
  reads_raw$well == molecules_raw$well
)

# Transpose data to gene-by-sample ---------------------------------------------

lane_id_single <- paste(reads_raw$individual, reads_raw$replicate,
                        reads_raw$well, reads_raw$index, reads_raw$lane,
                        reads_raw$flow_cell, sep = ".")

lane_id_bulk <- paste(reads_bulk_raw$individual, reads_bulk_raw$replicate,
                        reads_bulk_raw$well, reads_bulk_raw$index, reads_bulk_raw$lane,
                        reads_bulk_raw$flow_cell, sep = ".")

reads <- reads_raw %>%
  select(starts_with("ENSG"), starts_with("ERCC")) %>%
  t
colnames(reads) <- lane_id_single

molecules <- molecules_raw %>%
  select(starts_with("ENSG"), starts_with("ERCC")) %>%
  t
colnames(molecules) <- lane_id_single

reads_bulk <- reads_bulk_raw %>%
  select(starts_with("ENSG"), starts_with("ERCC")) %>%
  t
colnames(reads_bulk) <- lane_id_bulk

# Test counts ------------------------------------------------------------------

assert(
  "Every sample has at least one count",
  colSums(reads_bulk) > 0,
  colSums(reads) > 0,
  colSums(molecules) > 0
)

reads_bulk_rowsums <- rowSums(reads_bulk)
reads_rowsums <- rowSums(reads)
molecules_rowsums <- rowSums(molecules)

assert(
  "The files includes some genes with count data",
  any(reads_bulk_rowsums > 0),
  any(reads_rowsums > 0),
  any(molecules_rowsums > 0)
)

assert(
  "The files includes some genes without count data",
  any(reads_bulk_rowsums == 0),
  any(reads_rowsums == 0),
  any(molecules_rowsums == 0)
)

assert(
  "The counts should be higher in reads data compared to molecules data",
  reads_rowsums >= molecules_rowsums,
  reads >= molecules
)

assert(
  "Genes with zero reads have zero molecules",
  molecules[reads == 0] == 0
)

assert(
  "Genes with at least one molecule have at least one read",
  reads[molecules > 0] > 0
)


# fails
assert(
  "Genes with at least one read have at least one molecule",
  molecules[reads > 0] > 0
)

# fails
assert(
  "Genes with zero molecules have zero reads",
  reads[molecules == 0] == 0
)

# fails
assert(
  "Genes with zero counts should match in reads and molecules data",
  (reads == 0) == (molecules == 0)
)
