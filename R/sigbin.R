###
## Signal-based binning.
##
## Author: Jiri Hon <jiri.hon@gmail.com>
## Date: 2016/11/03
## Package: sigbin
##

#' Perform binning based on Hjorth parameters.
#'
#' @param hp Data frame of Hjorth parameters.
#' @param k Number of clusters.
#'
#' @return Input data frame extended by binning annotation.
#' @examples
#' file <- system.file("extdata", "reads.fasta", package="sigbin")
#' hp <- fasta_hjorth_params(file)
#' binning <- sigbin(hp)
#'
sigbin <- function(hp, k = NULL) {
  dst_m <- dist(hp, method = "euclidean")
  tree <- hclust(dst_m, method = "average")
  if (is.null(k)) {
    bin <- cutreeDynamicTree(tree, deepSplit = FALSE)
  } else {
    bin <- cutree(tree, k)
  }
  return(data.frame(
    activity = hp$activity,
    mobility = hp$mobility,
    complexity = hp$complexity,
    bin = factor(bin)
  ))
}

#' Plot signal based binning in 3D.
#'
#' @param sb Data frame with binning annotation.
#' @examples
#' file <- system.file("extdata", "reads.fasta", package="sigbin")
#' hp <- fasta_hjorth_params(file)
#' binning <- sigbin(hp)
#' plot_sigbin_3d(binning)
#'
plot_sigbin_3d <- function(sb) {
  colors <- tim.colors(length(levels(sb$bin)))
  bin_colors <- sb$bin
  levels(bin_colors) <- colors

  open3d()
  plot3d(sb$activity, sb$mobility, sb$complexity, col = bin_colors,
    xlab = 'Activity', ylab = 'Mobility', zlab = 'Complexity',
    main = 'Hjorth phase descriptors')
}

# ref_csv <- read.csv("../feature_select/features.csv")
# hp <- data.frame(activity = ref_csv$HP1_aktivita, mobility = ref_csv$HP2_mobilita, complexity = ref_csv$HP3_komplexita)
# binning <- sigbin(hp)
#
# colors <- tim.colors(length(levels(binning$bin)))
# bin_colors <- binning$bin
# levels(bin_colors) <- colors
#
# open3d()
# plot3d(binning$activity, binning$mobility, binning$complexity, col = bin_colors,
#   xlab = 'Activity', ylab = 'Mobility', zlab = 'Complexity',
#   main = 'Hjorth descriptors from phase')
#
# for (i in 1:nrow(binning)) {
#   print(i)
#   points3d(binning$activity[i], binning$mobility[i], binning$complexity[i], col = colors[as.integer(binning$bin[i])+1])
#   # if (i == 1) {
#   #   plot3d(x = binning$activity[i], y = binning$mobility[i], z = binning$complexity[i],
#   #          col = colors[as.integer(binning$bin[i])+1], size=2.5,
#   #          xlab = 'Activity', ylab = 'Mobility', zlab = 'Complexity', axes = TRUE, box = TRUE, top = TRUE,
#   #          main = 'Hjorth descriptors from phase')
#   # } else {
#   #   plot3d(x = binning$activity[i], y = binning$mobility[i], z = binning$complexity[i],
#   #          col = colors[as.integer(binning$bin[i])+1], size=2.5)
#   # }
# }
#
# i = 1
# print_point <- function(binning, i) {
#   points3d(binning$activity[i], binning$mobility[i], binning$complexity[i], col = colors[as.integer(binning$bin[i])+1])
# }
