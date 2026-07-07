#' Bin taxa by prevalence and plot their abundance
#'
#' Bins the features (columns) of a sample-by-feature count matrix according to
#' how prevalent each feature is across samples, then draws a stacked bar plot of
#' the mean relative abundance of each taxonomic group within each prevalence
#' bin. Features can optionally be aggregated to a supplied taxonomic level
#' before binning.
#'
#' Prevalence is defined as the proportion of samples (rows of \code{mat}) in
#' which a feature is present (count greater than zero). Features are assigned to
#' one of ten prevalence bins: \code{"0-10\%"}, \code{"10-20\%"}, ...,
#' \code{"90-100\%"}.
#'
#' @param mat A numeric matrix or data frame of counts with samples in rows and
#'   features (e.g. OTUs / ASVs) in columns. Row names are sample IDs and column
#'   names are feature IDs.
#' @param tax_level A vector giving the taxonomic label for each column of
#'   \code{mat}. Its length must equal \code{ncol(mat)}. \code{NA} values are
#'   relabelled \code{"Unknown"}.
#' @param bin_tax_first Logical. If \code{TRUE} (default) features are aggregated
#'   (summed) to \code{tax_level} before prevalence is computed, so binning is
#'   done at the taxonomic level. If \code{FALSE} binning is done per feature.
#' @param top.X Integer. Maximum number of taxonomic groups to display, ordered
#'   by total abundance across bins. Defaults to \code{20}.
#' @param only_when_present Logical. If \code{TRUE} (default) abundances are
#'   averaged only over samples in which the group is present (zeros treated as
#'   missing); the y-axis then reads "\% of population when present". If
#'   \code{FALSE} zeros are included in the mean.
#' @param seed Optional integer. If supplied, \code{set.seed(seed)} is called
#'   before colours are sampled so the plot is reproducible. Defaults to
#'   \code{NULL} (colours vary between calls, matching the original behaviour).
#'
#' @return Invisibly, a two-column character matrix (\code{col.key}) mapping each
#'   displayed taxonomic group to the colour used for it in the plot. Called for
#'   its side effect of drawing a bar plot.
#'
#' @examples
#' # Load the bundled example data (shipped in inst/extdata)
#' otu <- read.csv(system.file("extdata", "otu.csv", package = "CoreExplore"),
#'                 row.names = 1)
#' tax <- read.csv(system.file("extdata", "tax.csv", package = "CoreExplore"))
#' map <- read.csv(system.file("extdata", "map.csv", package = "CoreExplore"))
#'
#' # All samples, binned at family level, reproducible colours
#' bin_by_prevalence(mat = otu, tax_level = tax$family, seed = 1)
#'
#' # Restrict to human host samples
#' human <- map$host.species == "human"
#' bin_by_prevalence(mat = otu[human, ], tax_level = tax$family, seed = 1)
#'
#' @importFrom stats aggregate
#' @importFrom graphics barplot legend
#' @export
bin_by_prevalence <- function(mat,
                              tax_level,
                              bin_tax_first = TRUE,
                              top.X = 20,
                              only_when_present = TRUE,
                              seed = NULL) {

  if (ncol(mat) != length(tax_level)) {
    stop("mismatched taxonomy: ncol(mat) must equal length(tax_level)")
  } else {
    message("taxonomy matches")
  }

  tax_level[is.na(tax_level)] <- "Unknown"

  if (bin_tax_first) {

    mat.t.agg <- aggregate(t(mat), list(tax_level), sum)
    rownames(mat.t.agg) <- mat.t.agg$Group.1

    tax_level <- mat.t.agg$Group.1
    mat.t.agg$Group.1 <- NULL
    mat <- t(mat.t.agg)
  }

  mat.bi <- mat > 0

  bins <- seq(0, nrow(mat), length.out = 11)
  bins <- bins[-11]
  names(bins) <- paste0(0:9 * 10, "-", 1:10 * 10, "%")

  binning <- NULL

  for (i in seq_len(ncol(mat.bi))) {
    binning[i] <- names(rev(bins[sum(mat.bi[, i]) >= bins])[1])
  }

  print(table(binning))

  key <- as.data.frame(cbind(tax_level, binning))

  levels <- sort(unique(key$binning))

  tax_level.unique <- unique(tax_level)

  res <- as.data.frame(matrix(ncol = length(levels),
                              nrow = length(tax_level.unique)))

  colnames(res) <- levels
  rownames(res) <- tax_level.unique

  for (i in seq_along(levels)) {

    key.i <- key[key$binning == levels[i], ]
    mat.i <- mat[, key$binning == levels[i]]

    mat.i.agg <- aggregate(t(mat.i), list(key.i$tax_level), sum)
    rownames(mat.i.agg) <- mat.i.agg$Group.1
    mat.i.agg$Group.1 <- NULL

    mat.i.agg <- t(mat.i.agg)

    if (only_when_present) {
      mat.i.agg[mat.i.agg == 0] <- NA
      lab <- "% of population when present"
    } else {
      lab <- "% of population"
    }

    res.i <- colMeans(mat.i.agg, na.rm = TRUE)
    m <- match(names(res.i), rownames(res))

    res[m, i] <- res.i
  }

  color <- c("pink1", "green", "mediumpurple1", "slateblue1",
             "gold", "orchid", "turquoise2", "skyblue", "steelblue",
             "tan2", "navyblue", "orange", "orangered", "coral2",
             "palevioletred", "violetred", "darkred", "springgreen2",
             "yellowgreen", "palegreen4", "wheat2", "tan", "magenta",
             "tan3", "brown", "yellow", "snow2", "blue")

  o <- order(rowSums(res, na.rm = TRUE), decreasing = TRUE)
  res <- res[o, ]

  res[is.na(res)] <- 0

  top.X <- min(nrow(res), top.X)

  res <- res[1:top.X, ]

  if (!is.null(seed)) {
    set.seed(seed)
  }

  cols <- sample(color, nrow(res))

  col.key <- cbind(rownames(res), cols)

  p <- col.key[, 1] == "Unknown"
  col.key[p, 2] <- "grey50"

  res <- res / mean(rowSums(mat)) * 100

  barplot(as.matrix(res), width = 1, beside = FALSE,
          space = 0.1, col = col.key[, 2],
          xlim = c(0, ncol(res) * 2.1),
          las = 2,
          ylab = lab,
          font = 4, font.lab = 4,
          cex.axis = 1, cex.names = 1)

  legend("right",
         bty = "n",
         text.font = 4,
         fill = rev(col.key[, 2]),
         legend = rev(col.key[, 1]))

  invisible(col.key)
}
