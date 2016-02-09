library("testit")
library("data.table")
suppressPackageStartupMessages(library("dplyr"))

# Bulk reads -------------------------------------------------------------------

# The bulk reads per sample should be the sum of the bulk reads per lane.

reads_bulk_raw <- fread("../data/reads-raw-bulk-per-sample.txt")
setDF(reads_bulk_raw)

reads_bulk_raw_lane <- fread("../data/reads-raw-bulk-per-lane.txt")
setDF(reads_bulk_raw_lane)

reads_bulk_raw_lane_summed <- reads_bulk_raw_lane %>%
  select(individual, replicate, well, starts_with("ENSG"), starts_with("ERCC")) %>%
  group_by(individual, replicate, well) %>%
  summarise_each(funs(sum)) %>%
  arrange(individual, replicate, well) %>%
  ungroup

assert(
  "The bulk reads per sample are the sum of the bulk reads per lane",
  dim(reads_bulk_raw) == dim(reads_bulk_raw_lane_summed),
  # The following checks both the metadata and count columns
  reads_bulk_raw == reads_bulk_raw_lane_summed
)

# Single cell reads -------------------------------------------------------------------

# The single cell reads per sample should be the sum of the single cell reads
# per lane.

reads_raw <- fread("../data/reads-raw-single-per-sample.txt")
setDF(reads_raw)

reads_raw_lane <- fread("../data/reads-raw-single-per-lane.txt")
setDF(reads_raw_lane)

reads_raw_lane_summed <- reads_raw_lane %>%
  select(individual, replicate, well, starts_with("ENSG"), starts_with("ERCC")) %>%
  group_by(individual, replicate, well) %>%
  summarise_each(funs(sum)) %>%
  arrange(individual, replicate, well) %>%
  ungroup

assert(
  "The single cell reads per sample are the sum of the single cell reads per lane",
  dim(reads_raw) == dim(reads_raw_lane_summed),
  # The following checks both the metadata and count columns
  reads_raw == reads_raw_lane_summed
)

# Single cell molecules -------------------------------------------------------------------

# The single cell molecules per sample should be less than or equal to the sum
# of the single cell molecules per lane. This is because the per sample counts
# are obtained by first combining all the reads before removing duplicates using
# the UMI information. In contrast, for the per lane counts the same molecule
# can be independently counted in each lane that it is sequenced.

molecules_raw <- fread("../data/molecules-raw-single-per-sample.txt")
setDF(molecules_raw)

molecules_raw_lane <- fread("../data/molecules-raw-single-per-lane.txt")
setDF(molecules_raw_lane)

molecules_raw_lane_summed <- molecules_raw_lane %>%
  select(individual, replicate, well, starts_with("ENSG"), starts_with("ERCC")) %>%
  group_by(individual, replicate, well) %>%
  summarise_each(funs(sum)) %>%
  arrange(individual, replicate, well) %>%
  ungroup

molecules_metadata_cols <- colnames(molecules_raw) %in%
  c("individual", "replicate", "well")
molecules_ensg_cols <- grep("ENSG", colnames(molecules_raw))
molecules_ercc_cols <- grep("ERCC", colnames(molecules_raw))

# fails
assert(
  "The single cell molecules per sample less than or equal to the sum of the single cell molecules per lane",
  dim(molecules_raw) == dim(molecules_raw_lane_summed),
  molecules_raw[, molecules_metadata_cols] == molecules_raw_lane_summed[, molecules_metadata_cols],
  molecules_raw[, molecules_ensg_cols] <= molecules_raw_lane_summed[, molecules_ensg_cols],
  molecules_raw[, molecules_ercc_cols] <= molecules_raw_lane_summed[, molecules_ercc_cols]
)
