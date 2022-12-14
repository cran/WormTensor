% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/worm_visualize.R
\name{worm_visualize,WormTensor-method}
\alias{worm_visualize,WormTensor-method}
\title{Plots evaluation result
A visualization result is generated from a WormTensor object.}
\usage{
\S4method{worm_visualize}{WormTensor}(
  object,
  out.dir = tempdir(),
  algorithm = c("tSNE", "UMAP"),
  seed = 1234,
  tsne.dims = 2,
  tsne.perplexity = 15,
  tsne.verbose = FALSE,
  tsne.max_iter = 1000,
  umap.n_neighbors = 15,
  umap.n_components = 2,
  silhouette.summary = FALSE
)
}
\arguments{
\item{object}{WormTensor object with a result of worm_evaluate}

\item{out.dir}{Output directory (default: tempdir())}

\item{algorithm}{Dimensional reduction methods}

\item{seed}{Arguments passed to set.seed (default: 1234)}

\item{tsne.dims}{Output dimensionality (default: 2)}

\item{tsne.perplexity}{Perplexity parameter (default: 15)}

\item{tsne.verbose}{logical; Whether progress updates should be printed
(default: TRUE)}

\item{tsne.max_iter}{Number of iterations (default: 1000)}

\item{umap.n_neighbors}{The size of the local neighborhood (default: 15)}

\item{umap.n_components}{The dimension of the space to embed into
(default: 2)}

\item{silhouette.summary}{logical; If true a summary of cluster silhouettes
are printed.}
}
\value{
Silhouette plots. ARI with a merge result and each animal(with MCMI).
Dimensional reduction plots colored by cluster, no. of identified cells,
consistency(with labels), Class_label(with labels).
}
\description{
Plots evaluation result
A visualization result is generated from a WormTensor object.
}
\examples{
\donttest{
    # Temporary directory to save figures
    out.dir <- tempdir()

    # Labels
    worm_download("mSBD", qc = "PASS")$Ds |>
        as_worm_tensor() |>
            worm_membership(k = 6) |>
                worm_clustering() -> object
    Ds_mSBD <- worm_download("mSBD", qc = "PASS")
    labels <- list(
        label1 = replace(
            Ds_mSBD$labels$Class,
            which(is.na(Ds_mSBD$labels$Class)),
            "NA"
        ),
        label2 = sample(4, length(object@clustering), replace = TRUE),
        label3 = sample(5, length(object@clustering), replace = TRUE)
    )

    # Pipe Operation (without Labels)
    worm_download("mSBD", qc = "PASS")$Ds |>
        as_worm_tensor() |>
            worm_membership(k = 6) |>
                worm_clustering() |>
                    worm_evaluate() |>
                        worm_visualize(out.dir) -> object_no_labels

    # Pipe Operation (with Labels)
    worm_download("mSBD", qc = "PASS")$Ds |>
        as_worm_tensor() |>
            worm_membership(k = 6) |>
                worm_clustering() |>
                    worm_evaluate(labels) |>
                        worm_visualize(out.dir) -> object_labels
}
}
\references{
The .dist_nn function is quoted from dist_nn
(not exported function) in package uwot(\url{https://github.com/jlmelville/uwot/tree/f467185c8cbcd158feb60dde608c9da153ed10d7}).
}
