## Demo for the CoreExplore package
## Run after installing:  library(CoreExplore); source(system.file("examples", "demo.R", package = "CoreExplore"))

library(CoreExplore)

otu <- read.csv(system.file("extdata", "otu.csv", package = "CoreExplore"),
                row.names = 1)
tax <- read.csv(system.file("extdata", "tax.csv", package = "CoreExplore"))
map <- read.csv(system.file("extdata", "map.csv", package = "CoreExplore"))

dim(otu)
dim(map)
dim(tax)

# All samples
bin_by_prevalence(mat = otu, tax_level = tax$family, seed = 1)

# Human host samples
p <- map$host.species == "human"
bin_by_prevalence(mat = otu[p, ], tax_level = tax$family, seed = 1)
title(unique(map$host.species[p]))

# Mouse host samples
p <- map$host.species == "mouse"
bin_by_prevalence(mat = otu[p, ], tax_level = tax$family, seed = 1)
title(unique(map$host.species[p]))
