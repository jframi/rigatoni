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
plot.krnel<-function(x, label=TRUE, label.cex=0.6, label.col="red", plot.bbox=T, bbox.col="blue", title=T, plot.radii=T, plot.ellipse=T, plot.ws.bbox=F, ws.bbox.col="pink"){
  if (title){
    plot(0,0, type="n",xlim=c(0,x$params$w),ylim=c(0,x$params$h), xlab="", ylab="", asp = 1, main=x$params$image.name)
  }else{
    plot(0,0, type="n",xlim=c(0,x$params$w),ylim=c(0,x$params$h), xlab="", ylab="", asp = 1)
  }

  for (g in 1:length(x$contours)){
    polygon(x$contours[[g]][,1],x$params$h-x$contours[[g]][,2],col = x$features$cols[g], border = NA)
  }
  if (label){
    text(x = x$features$m.cx, y = x$params$h-x$features$m.cy, label=x$features$id, cex=label.cex, col = label.col)
  }
  if (plot.bbox){
    if (!is.null(x$bbox)){
      noret <- lapply(x$bbox, function(a) polygon(a$pts[,1], x$params$h-a$pts[,2], border = bbox.col))
    }
  }
  if (plot.radii){
    segments(x0 = x$features$m.cx,
             x1 = x$features$m.cx + x$features$s.radius.max*cos(-x$features$m.theta),
             y0 = x$params$h-x$features$m.cy,
             y1 = x$params$h-x$features$m.cy + x$features$s.radius.max*sin(-x$features$m.theta),col = "red")

    segments(x0 = x$features$m.cx,
             x1 = x$features$m.cx + x$features$s.radius.min*cos(-x$features$m.theta+pi/2),
             y0 = x$params$h-x$features$m.cy,
             y1 = x$params$h-x$features$m.cy + x$features$s.radius.min*sin(-x$features$m.theta+pi/2),col = "green")
  }
  if (plot.ellipse){
    plotrix::draw.ellipse(x = x$features$m.cx,
                          y = x$params$h-x$features$m.cy,
                          a = x$features$m.majoraxis/2,
                          b = sqrt(x$features$m.majoraxis^2*(1 - x$features$m.eccentricity^2))/2,
                          angle = -x$features$m.theta, deg = F, border = "blue")

  }
  if (plot.ws.bbox){
      if (!is.null(x$ws.bbox)){
        noret <- lapply(x$ws.bbox, function(b) lapply(b, function(a) polygon(a$pts[,1], x$params$h-a$pts[,2], border = ws.bbox.col)))
      }

  }
}
