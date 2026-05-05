# scAbsolute

```bash
Schneider, Michael P., et al. 2024.
scAbsolute: Measuring Single-Cell Ploidy and Replication Status.
Genome Biology 25 (1): 62.
https://doi.org/10.1186/s13059-024-03204-y.

```

Measuring single-cell ploidy and cell cycle status in single-cell DNA sequencing data.
Please see the [manuscript](https://doi.org/10.1186/s13059-024-03204-y) for details.

This repository contains the source code for the *scAbsolute* software and scripts to replicate the figures in the manuscript.

In order to run the code, to see examples of how to use the package and an easy-to-use workflow, we recommend using the lab's [single-cell sequencing pipeline](https://github.com/markowetzlab/scDNAseq-workflow).

## Per-cell failure detection

When running scAbsolute on many cells (e.g. via the scDNAseq-workflow Snakemake pipeline), some cells may fail the scaling step. Failed cells are now retained in the output object with a `failure_reason` field in their metadata, rather than crashing the process and producing an empty file.

### Inspecting failures after merging

```r
merged <- readRDS("results/500/MySample_500.rds")

# Summary of failure reasons (NA = cell passed)
table(pData(merged)$failure_reason, useNA = "ifany")

# Full details for failed cells only
library(dplyr)
pData(merged) %>%
  filter(!is.na(failure_reason)) %>%
  select(name, failure_reason, used.reads, rpc)
```

### Failure reasons

| `failure_reason` | Cause |
|---|---|
| `too_few_reads` | Average reads per bin < 0.5 — insufficient coverage |
| `zero_valid_bins_in_ploidy_region` | No bins survive blacklisting/filtering in the ploidy region (chr1–22) |
| `all_solutions_filtered` | DPGMM produced solutions but all were rejected by ploidy constraints (`minPloidy`/`maxPloidy`) |
| `process_crash` | R process crashed before output was written (recorded in `failed_cells.csv` only) |

### failed_cells.csv

The `merge.R` script in the workflow writes `failed_cells.csv` alongside the merged RDS, listing cells that either crashed (`process_crash`) or failed with a recorded reason. This file has two columns: `name` and `failure_reason`.
