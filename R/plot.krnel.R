#' plot features of a krnel object
#'
#' @param x a krnel object as returned by the krnel function
#' @param label boolean, plot label (ie number)
#' @param label.cex label size
#' @param label.col label color
#'
#' @return nothing
#' @export
#'
#' @examples
plot.krnel<-function(x, label=TRUE, label.cex=0.6, label.col="red"){
  plot(0,0, type="n",xlim=c(0,x$params$w),ylim=c(0,x$params$h), xlab="", ylab="", asp = 1)

  for (g in 1:length(x$contours)){
    polygon(x$contours[[g]][,1],x$params$h-x$contours[[g]][,2],col = x$features$cols[g], border = NA)
  }
  if (label){
    text(x = x$features$m.cx, y = x$params$h-x$features$m.cy, label=1:nrow(x$features), cex=label.cex, col = label.col)
  }
}
