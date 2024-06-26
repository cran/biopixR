% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edgeDetection.R
\name{edgeDetection}
\alias{edgeDetection}
\title{Canny edge detector}
\usage{
edgeDetection(img, t1, t2, alpha = 1, sigma = 2)
}
\arguments{
\item{img}{image (import by \code{\link[biopixR]{importImage}})}

\item{t1}{threshold for weak edges (if missing, both thresholds are
determined automatically)}

\item{t2}{threshold for strong edges}

\item{alpha}{threshold adjustment factor (default 1)}

\item{sigma}{smoothing (default 2)}
}
\value{
Object of class 'cimg', displaying detected edges.
}
\description{
Adapted code from the 'imager' \code{\link[imager]{cannyEdges}} function
without the usage of 'dplyr' and 'purrr'. If the threshold parameters are
missing, they are determined automatically using a k-means heuristic. Use
the alpha parameter to adjust the automatic thresholds up or down. The
thresholds are returned as attributes. The edge detection is based on a
smoothed image gradient with a degree of smoothing set by the sigma
parameter.
}
\examples{
edgeDetection(beads, alpha = 0.5, sigma = 0.5) |> plot()
}
\references{
https://CRAN.R-project.org/package=imager
}
