#' Plot the histogram of hue values of image to determine the hue range of the background
#'
#' @param img
#' @param return.values
#'
#' @return
#' @export
#'
#' @examples
get_huethr_values<-function(img, return.values=TRUE){
  img<-resize(img,w = 500)
  rgb<-matrix(c(c(img[,,1]),c(img[,,2]),c(img[,,3])),nrow = 3, byrow = T)
  # convert to hsv
  x<-rgb2hsv(rgb)
  # hist h component
  hist(x[1,],xlim=c(0,1),xlab="Hue value", main="")
  if (return.values){
    bounds<-locator(2)
    return(round(255*bounds$x,0))
  }
}


get_athr_values<-function(img){
  img<-resize(img,w = 500)
  rgb<-t(matrix(c(c(img[,,1]),c(img[,,2]),c(img[,,3])),nrow = 3, byrow = T))
  # convert to Lab
  x <- farver::convert_colour(rgb, from = "rgb", to = "lab")
  # get a component
  xa<-matrix(x[,2],nrow = 500)
  xht<-xa-min(xa)
  # hist h component
  hist(xht,xlab="a value")
  return(locator(1)[[1]][1])
}
