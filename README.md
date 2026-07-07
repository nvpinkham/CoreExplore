# CoreExplore

> Explore the core microbiome of a sample-by-feature count table by prevalence.

`CoreExplore` provides `bin_by_prevalence()`, which groups the features (OTUs /
ASVs) of a microbiome count table by how many samples they occur in, then draws
a stacked bar plot of each taxonomic group's mean relative abundance within each
prevalence bin. Features can be aggregated to a taxonomic level (e.g. family)
before binning, making it easy to see which taxa form the stable "core" of a
community versus the rare tail.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("nvpinkham/CoreExplore")
```

## Quick start

```r
library(CoreExplore)

# Example data ships with the package
otu <- read.csv(system.file("extdata", "otu.csv", package = "CoreExplore"),
                row.names = 1)
tax <- read.csv(system.file("extdata", "tax.csv", package = "CoreExplore"))
map <- read.csv(system.file("extdata", "map.csv", package = "CoreExplore"))

# Bin every sample at the family level
bin_by_prevalence(mat = otu, tax_level = tax$family, seed = 1)
```

Split by host to compare cores between groups:

```r
# Human host samples
human <- map$host.species == "human"
bin_by_prevalence(mat = otu[human, ], tax_level = tax$family, seed = 1)
title(unique(map$host.species[human]))

# Mouse host samples
mouse <- map$host.species == "mouse"
bin_by_prevalence(mat = otu[mouse, ], tax_level = tax$family, seed = 1)
title(unique(map$host.species[mouse]))
```

## How it works

Prevalence is the proportion of samples (rows of `mat`) in which a feature is
present (count > 0). Each feature is assigned to one of ten prevalence bins —
`0-10%`, `10-20%`, …, `90-100%` — and the mean relative abundance of each
taxonomic group is plotted per bin. The rightmost bins hold the taxa found in
nearly every sample: the community core.

## Arguments

| Argument | Description | Default |
|---|---|---|
| `mat` | Count table, samples in rows, features in columns (names required) | — |
| `tax_level` | Taxonomic label per column of `mat`; length must equal `ncol(mat)` | — |
| `bin_tax_first` | Aggregate to `tax_level` before binning | `TRUE` |
| `top.X` | Max number of groups to plot, by total abundance | `20` |
| `only_when_present` | Average abundance only over samples where a group is present | `TRUE` |
| `seed` | Optional integer for reproducible plot colours | `NULL` |

Returns (invisibly) the taxon-to-colour key used in the plot.

## Example data

The bundled dataset (`inst/extdata/`) is a demo gut-microbiome study:

- `otu.csv` — sample-by-OTU count table (119 samples × 243 OTUs)
- `tax.csv` — per-OTU taxonomy (`family`, `phylum`)
- `map.csv` — per-sample metadata (`host.species`, `participant`, alpha
  diversity, etc.)

## License

MIT © Nick Pinkham
