#' plot features of a krnel object
#'
#' @param x a krnel object as returned by the krnel function
#' @param label boolean, plot label (ie number)
#' @param label.cex label size
#' @param label.col label color
#' @param plot.bbox whether to plot bounding boxes or not
#' @param bbox.col color of bounding box
#'
#' @return nothing
#' @export
#'
#' @examples
plot.krnel<-function(x, label=TRUE, label.cex=0.6, label.col="red", plot.bbox=T, bbox.col="blue"){
  plot(0,0, type="n",xlim=c(0,x$params$w),ylim=c(0,x$params$h), xlab="", ylab="", asp = 1)

  for (g in 1:length(x$contours)){
    polygon(x$contours[[g]][,1],x$params$h-x$contours[[g]][,2],col = x$features$cols[g], border = NA)
  }
  if (label){
    text(x = x$features$m.cx, y = x$params$h-x$features$m.cy, label=names(x$contours), cex=label.cex, col = label.col)
  }
  if (plot.bbox){
    if (!is.null(x$bbox)){
      noret <- lapply(x$bbox, function(a) polygon(a$pts[,1], x$params$h-a$pts[,2], border = bbox.col))
    }
  }
}
