## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(clustur)
set.seed(19760620)

## -----------------------------------------------------------------------------
full_count_table <- read_count(example_path("amazon.full.count_table"))
sparse_count_table <- read_count(example_path("amazon.sparse.count_table"))

## -----------------------------------------------------------------------------
column_distance <- read_dist(example_path("amazon_column.dist"), full_count_table, cutoff = 0.03)

## -----------------------------------------------------------------------------
phylip_distance <- read_dist(example_path("amazon_phylip.dist"), full_count_table, cutoff = 0.03)

## -----------------------------------------------------------------------------
get_distance_df(column_distance)
get_distance_df(phylip_distance)

## -----------------------------------------------------------------------------
cutoff <- 0.03
cluster_data <- cluster(column_distance, cutoff)

## -----------------------------------------------------------------------------

cluster_data <- cluster(column_distance, cutoff, method = "furthest")
cluster_data <- cluster(column_distance, cutoff, method = "nearest")
cluster_data <- cluster(column_distance, cutoff, method = "average")
cluster_data <- cluster(column_distance, cutoff, method = "weighted")

## ----R.options=list(max.print=15)---------------------------------------------
clusters <- cluster(column_distance, cutoff, method = "opticlust")
get_cutoff(clusters)
get_bins(clusters)
get_abundance(clusters)
get_metrics(clusters)

