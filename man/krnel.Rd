% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/krnel.R
\name{krnel}
\alias{krnel}
\title{main function to identify features in an image}
\usage{
krnel(
  img,
  crop = NULL,
  resizw = NULL,
  watershed = F,
  huethres,
  vthres = NULL,
  minsize,
  maxsize,
  save.outline = F,
  img.name = NULL,
  blackbg = F,
  whitebg = F,
  ws.avg = F,
  bw = F,
  color.erode = F,
  colerode.rad.ratio = 0.75,
  open = FALSE,
  open.brush.size = 5
)
}
\arguments{
\item{img}{the path to the image file to analyze or an Image object as defined by the EBImage package}

\item{crop}{a numeric vector in the form c(xmin, xmax, ymin, ymax) to define cropping area}

\item{resizw}{the width of resized image to use for feature detection. decreasing the size of the image will improve performance}

\item{huethres}{a numeric vector in the form c(Huemin, Huemax) to define Hue threshold. Can be easily determined using the get_huethr_values function}

\item{vthres}{an optional v value threshold under which pixels are considered as part of features (override huesthres values)}

\item{minsize}{the minimum size of features to be considered, in pixel, as determined on the resized image}

\item{maxsize}{the maximum size of features to be considered, in pixel, as determined on the resized image}

\item{save.outline}{boolean, should an outline image be saved in the same directory?}

\item{img.name}{name of the image. This is useful in the case img is an Image object and save.outline is TRUE. This will be used for the name of the outline file.}

\item{blackbg}{behave differently in case the background is black}

\item{whitebg}{behave differently in case the background is white}

\item{ws.avg}{watershed average : this allows to perform a watershed after a feature detection with no watershed, and to compute for each feature bbox width/height average and sum on the different sub-features obtained through watershed}

\item{bw}{boolean, is it a black and white (greyscale) image? In this case huethres is not required.}

\item{color.erode}{boolean, apply morphological erosion to each feature before getting feature's color. This is way to estimate color from the center part of each feature, to avoid border effects.}

\item{colerode.rad.ratio}{if color.erode is true, this is the size of the round brush used for morphological erosion expressed as a ratio (>0, <=1) of the average radius of all feaures detected on the image}

\item{open}{apply an open operation (erode+dilate) to suppress possible artefact attached to the detected feature. If watershed=T open is applied before watershed.}

\item{open.brush.size}{size of the brush to use for the open operation. Brush shape is hardcoded as 'disc'.}
}
\value{
the function will return krnel object, ie. a list with three components
\itemize{
\item features : a data.frame with as many rows as the number of detected features and description variables as returned by \link[EBImage:computeFeatures]{EBImage::computeFeatures} It also includes : (i) bounding box width and height aka Feret min and max diameter (bbox.width, bbox.height), (ii) pole of inaccessibility coordinates and longest distance to 'coastline' (poi.x,poi.y,poi.dist). In complex shapes, poi.dist might be a good measure of object width.
\item contours : a list with as many elements as the number of detected features. Each element is a matrix with the coordinates of each feature.
\item bbox : a list with as many elements as the number of detected features. Each element is a list with 4 components: $pts that contains the coordinates of each corner of the bounding box, $width, $height, and $angle
\item params : a list with the analysis parameters. can be used for further plotting
\item ws.contours : if ws.avg=T. a list with as many elements as the number of detected main features. Each element is a list of matrices with contours coordinates of each sub-features.
\item ws.bbox : if ws.avg=T. The same with bounding boxes
\item ws.pois : if ws.avg=T. The same with point of inaccessibility
}
}
\description{
main function to identify features in an image
}
