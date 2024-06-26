% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objectDetection.R
\name{objectDetection}
\alias{objectDetection}
\title{Object detection}
\usage{
objectDetection(img, method = "edge", alpha = 1, sigma = 2, vis = TRUE)
}
\arguments{
\item{img}{image (import by \code{\link[biopixR]{importImage}})}

\item{method}{choose method for object detection ('edge' / 'threshold')}

\item{alpha}{threshold adjustment factor (numeric / 'static' / 'interactive' / 'gaussian') (only needed for 'edge')}

\item{sigma}{smoothing (numeric / 'static' / 'interactive' / 'gaussian') (only needed for 'edge')}

\item{vis}{creates image were object edges/coordinates (purple) and detected centers (green) are highlighted (TRUE | FALSE)}
}
\value{
list of 3 objects:
\itemize{
\item \code{data.frame} of labeled regions with the central coordinates (including size information).
\item All coordinates that are in labeled regions.
\item Image where object edges/coordinates (purple) and detected centers (green) are colored.
}
}
\description{
This function identifies objects in an image using either edge detection or
thresholding methods. It gathers the coordinates and centers of the
identified objects, highlighting the edges or overall coordinates for easy
recognition.
}
\details{
The \code{\link[biopixR]{objectDetection}} function provides several methods
for calculating the alpha and sigma parameters, which are critical for edge
detection:
\enumerate{
\item \strong{Input of a Numeric Value:}
\itemize{
\item Users can directly input numeric values for alpha and sigma, allowing for precise control over the edge detection parameters.
}
\item \strong{Static Scanning:}
\itemize{
\item When both alpha and sigma are set to "static", the function systematically tests all possible combinations of these parameters within the range (alpha: 0.1 - 1.5, sigma: 0 - 2). This exhaustive search helps identify the optimal parameter values for the given image. (Note: takes a lot of time)
}
\item \strong{Interactive Selection:}
\itemize{
\item Setting the alpha and sigma values to "interactive" initiates a Tcl/Tk graphical user interface (GUI). This interface allows users to adjust the parameters interactively, based on visual feedback. To achieve optimal results, the user must input the necessary adjustments to align the parameters with the specific requirements of the image. The user can also switch between the methods through the interface.
}
\item \strong{Multi-Objective Optimization:}
\itemize{
\item For advanced parameter optimization, the function \code{\link[GPareto]{easyGParetoptim}} will be utilized for multi-objective optimization using Gaussian process models. This method leverages the 'GPareto' package to perform the optimization. It involves building Gaussian Process models for each objective and running the optimization to find the best parameter values.
}
}
}
\examples{
res_objectDetection <- objectDetection(beads,
                                       method = 'edge',
                                       alpha = 1,
                                       sigma = 0)
res_objectDetection$marked_objects |> plot()

res_objectDetection <- objectDetection(beads,
                                       method = 'threshold')
res_objectDetection$marked_objects |> plot()
}
