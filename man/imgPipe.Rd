% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imgPipe.R
\name{imgPipe}
\alias{imgPipe}
\title{Image analysis pipeline}
\usage{
imgPipe(
  img1 = img,
  color1 = "color1",
  img2 = NULL,
  color2 = "color2",
  img3 = NULL,
  color3 = "color3",
  method = "edge",
  alpha = 1,
  sigma = 2,
  sizeFilter = FALSE,
  upperlimit = "auto",
  lowerlimit = "auto",
  proximityFilter = FALSE,
  radius = "auto"
)
}
\arguments{
\item{img1}{image (import by \code{\link[biopixR]{importImage}})}

\item{color1}{name of color in img1}

\item{img2}{image (import by \code{\link[biopixR]{importImage}})}

\item{color2}{name of color in img2}

\item{img3}{image (import by \code{\link[biopixR]{importImage}})}

\item{color3}{name of color in img3}

\item{method}{choose method for object detection ('edge' / 'threshold')
(from \code{\link[biopixR]{objectDetection}})}

\item{alpha}{threshold adjustment factor (numeric / 'static' / 'interactive' / 'gaussian')
(from \code{\link[biopixR]{objectDetection}})}

\item{sigma}{smoothing (numeric / 'static' / 'interactive' / 'gaussian')
(from \code{\link[biopixR]{objectDetection}})}

\item{sizeFilter}{applying \code{\link[biopixR]{sizeFilter}} function (default - FALSE)}

\item{upperlimit}{highest accepted object size (numeric / 'auto')
(only needed if sizeFilter = TRUE)}

\item{lowerlimit}{smallest accepted object size (numeric / 'auto')
(only needed if sizeFilter = TRUE)}

\item{proximityFilter}{applying \code{\link[biopixR]{proximityFilter}} function (default - FALSE)}

\item{radius}{distance from one object in which no other centers
are allowed (in pixels) (only needed if proximityFilter = TRUE)}
}
\value{
list of 2 to 3 objects:
\itemize{
\item Summary of all the objects in the image.
\item Detailed information about every single object.
\item (optional) Result for every individual color.
}
}
\description{
This function serves as a pipeline that integrates tools for complete
start-to-finish image analysis. It enables the handling of images from
different channels, for example the analysis of dual-color micro particles.
This approach simplifies the workflow, providing a straightforward method to
analyze complex image data.
}
\examples{
result <- imgPipe(
  beads,
  alpha = 1,
  sigma = 2,
  sizeFilter = TRUE,
  upperlimit = 150,
  lowerlimit = 50
  )

# Highlight remaining microparticles
plot(beads)
with(
  result$detailed,
  points(
    result$detailed$x,
    result$detailed$y,
    col = "darkgreen",
    pch = 19
    )
  )
}
\seealso{
\code{\link[=objectDetection]{objectDetection()}}, \code{\link[=sizeFilter]{sizeFilter()}}, \code{\link[=proximityFilter]{proximityFilter()}}, \code{\link[=resultAnalytics]{resultAnalytics()}}
}
