% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/worm_download.R
\name{worm_download}
\alias{worm_download}
\title{Downloads distance matrices
28 animals' data including 24 normal and 4 noisy are retrieved from figshare.}
\usage{
worm_download(distance = c("mSBD", "Euclid"), qc = c("PASS", "WARN", "FAIL"))
}
\arguments{
\item{distance}{"mSBD" or "Euclid" can be specified}

\item{qc}{"PASS" or "WARN" or "FAIL" can be specified. "PASS" downloads
24 data except 4 noisy data. "WARN" downloads 27 data except 1 noisy data.
"FAIL" downloads all 28 data.}
}
\value{
A List of containing distance matrices. The list also includes
metadata for each animals.
}
\description{
Downloads distance matrices
28 animals' data including 24 normal and 4 noisy are retrieved from figshare.
}
\examples{
\donttest{
    Ds_Euclid <- worm_download("Euclid", qc = "WARN")
    Ds_mSBD <- worm_download("mSBD", qc = "PASS")
}
}
