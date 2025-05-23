% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.R
\name{cluster}
\alias{cluster}
\title{Cluster entities together}
\usage{
cluster(
  distance_object,
  cutoff,
  method = "opticlust",
  feature_column_name_to = "feature",
  bin_column_name_to = "bin",
  random_seed = 123
)
}
\arguments{
\item{distance_object}{The distance object that
was created using the `read_dist()` function.}

\item{cutoff}{The cutoff you want to cluster towards.}

\item{method}{The method of clustering to be performed: opticlust (default),
furthest, nearest, average, or weighted.}

\item{feature_column_name_to}{Set the name of the column in the cluster
dataframe that contains the sequence names.}

\item{bin_column_name_to}{Set the name of the column in the
cluster dataframe that contains the name of the group of sequence names.}

\item{random_seed}{the random seed to use, (default = 123).}
}
\value{
A list of `data.frames` that contain abundance, and clustering
results. If you used `method = opticlust`, it will also return clustering
performance metrics.
}
\description{
Clusters entities represented in a distance matrix and count table using
one of several algorithms and outputs information about the composition and
abundance of each cluster
}
\examples{

 cutoff <- 0.03
 count_table <- read_count(example_path("amazon.full.count_table"))
 distance_data <- read_dist(example_path("amazon_column.dist"),
                            count_table, cutoff)

 cluster_results <- cluster(distance_data,
                            cutoff, method = "opticlust",
                            feature_column_name_to = "sequence",
                            bin_column_name_to = "omu")
 cluster_results <- cluster(distance_data,
                            cutoff, method = "furthest")
 cluster_results <- cluster(distance_data,
                            cutoff, method = "nearest")
 cluster_results <- cluster(distance_data,
                            cutoff, method = "average")
 cluster_results <- cluster(distance_data,
                            cutoff, method = "weighted")


}
