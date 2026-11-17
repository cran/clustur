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

## -----------------------------------------------------------------------------
  amazon_count <- read_count(example_path("amazon.full.count_table"))
  amazon_dist <- read_dist(example_path("amazon_phylip.dist"), amazon_count, 0.03)

  esophagus_fasta <- read_fasta(example_path("esophagus.fasta"))
  esophagus_count <- read_count(example_path("esophagus.count_table"))
  amazon_fasta <- read_fasta(example_path("amazon.align"))
  esophagus_list <- read_mothur_list(example_path("esophagus.list"))


  denovo <- optifit(amazon_dist, 0.03)
  denovo_open <- optifit(amazon_dist, 0.03, close = FALSE)
  accnos <- optifit(amazon_dist, 0.03, accnos = amazon_count$Representative_Sequence[1:40])
  accnos_open <- optifit(amazon_dist, 0.03, accnos = amazon_count$Representative_Sequence[1:40],
                         closed=FALSE)

  ref_cluster <- optifit(amazon_dist, 0.03, "feature", "bin", accnos = NULL,
                         ref_fasta = esophagus_fasta, ref_count = esophagus_count,
                         ref_list = esophagus_list, fit_fasta = amazon_fasta,
                         number_of_threads = parallel::detectCores(), print_ref = TRUE, closed = FALSE)


## -----------------------------------------------------------------------------
  cutoff <- 0.03
  count_table <- read_count(example_path("amazon.full.count_table"))
  fasta <- read_fasta(example_path("amazon.align"))
  taxonomy <- read_mothur_taxonomy(example_path("amazon.taxonomy"))
  df <- cluster_split(fasta = fasta, taxonomy = taxonomy, count_table = count_table, cluster_method = "opti", cutoff = cutoff, taxonomy_cutoff = 3, number_of_threads = 1)

## ----R.options=list(max.print=15)---------------------------------------------
clusters <- cluster(column_distance, cutoff, method = "opticlust")
get_cutoff(clusters)
get_bins(clusters)
get_abundance(clusters)
get_metrics(clusters)

