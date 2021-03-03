#' main function to identify features in an image
#'
#' @param img the path to the image file to analyze or an Image object as defined by the EBImage package
#' @param crop a numeric vector in the form c(xmin, xmax, ymin, ymax) to define cropping area
#' @param resizw the width of resized image to use for feature detection. decreasing the size of the image will improve performance
#' @param huethres a numeric vector in the form c(Huemin, Huemax) to define Hue threshold. Can be easily determined using the get_huethr_values function
#' @param minsize the minimum size of features to be considered, in pixel, as determined on the resized image
#' @param maxsize the maximum size of features to be considered, in pixel, as determined on the resized image
#' @param save.outline boolean, should an outline image be saved in the same directory
#'
#' @return the function will return krnel object, ie. a list with three components
#' * features : a data.frame with as many rows as the number of detected features and description variables as returned by EBImage::computeFeatures
#' * contours : a list with as many elements as the number of detected features Each element is a matrix with the coordinates of each feature.
#' * params : a list with the analysis parameters. can be used for further plotting
#' @md
#' @export
#'
#' @examples
krnel<- function(img, crop, resizw, huethres, minsize, maxsize, save.outline=F){

  mf<-match.call()
  if(class(eval(mf$img))=="Image"){
    img.name<-as.character(mf$img)
  }else{
    img.name<-eval(mf$img)
  }

  if (!class(img)=="Image"){
    img = readImage(img)
  }
  #crop
  img <- img[crop[1]:crop[2], crop[3]:crop[4],1:3]
  #resize
  img<-resize(img,w = resizw)
  # Make rgb matrix
  #rgb<-rbind(c(img[,,1]),c(img[,,2]),c(img[,,3]))
  rgb<-matrix(c(c(img[,,1]),c(img[,,2]),c(img[,,3])),nrow = 3, byrow = T)
  # convert to hsv
  x<-rgb2hsv(rgb)
  # get h component
  xh<-matrix(x[1,],nrow = resizw)
  xht<-xh
  # hue threshold
  minhue<-huethres[1]
  maxhue<-huethres[2]
  xht[xht<=(minhue/255)]<-0
  xht[xht>=(maxhue/255)]<-0
  xht[xht>(minhue/255) & xht<(maxhue/255)]<-1
  # make it an image
  xhti<-as.Image(xht)

  # get mask
  nmask = fillHull(1-xhti)
  nmask = bwlabel(nmask)

  # filter on minsize
  nmask[!nmask%in%(which(table(c(imageData(nmask)))>minsize)-1)]<-0
  # filter on maxsize
  nmask[!nmask%in%(names(which(table(c(imageData(nmask)))<maxsize)))]<-0
  # get features
  nmask.shp<-computeFeatures.shape(nmask)
  nmask.mom<-computeFeatures.moment(nmask)
  nmask.bas<-computeFeatures.basic(nmask,img)
  #nmask.basr<-computeFeatures.basic(nmask,channel(img,"r"))
  #nmask.basg<-computeFeatures.basic(nmask,channel(img,"g"))
  #nmask.basb<-computeFeatures.basic(nmask,channel(img,"b"))
  nmask.basr<-computeFeatures.basic(nmask,img[,,1])
  nmask.basg<-computeFeatures.basic(nmask,img[,,2])
  nmask.basb<-computeFeatures.basic(nmask,img[,,3])
  nmask.cont<-ocontour(nmask)
  cols<-rgb(nmask.basr[,1],nmask.basg[,1],nmask.basb[,1])

  ret<-list(features= data.table(nmask.shp,
                                 nmask.mom,
                                 nmask.bas,
                                 cols),
            contours=nmask.cont,
            params=list(crops=crop,
                        w=dim(img)[1],
                        h=dim(img)[2],
                        huethres=huethres,
                        minsize=minsize,
                        image.name=img.name))
  if (save.outline){
    writeImage(paintObjects(nmask,img,thick = resizw>1500),
               files=paste0(img.name,".outline.jpg"),
               type="jpeg",
               quality = 80)
  }
  class(ret) <- c(class(ret), "krnel")
  return(ret)
}

#' lab space version of krnel function (experimental)
#'
#' @param img
#' @param crop
#' @param resizw
#' @param athres
#' @param minsize
#' @param maxsize
#' @param save.outline
#'
#' @return
#' @export
#' @examples
krnellab<- function(img, crop=NULL, resizw=NULL, athres="auto", minsize,maxsize, save.outline=F){

  mf<-match.call()
  if(class(eval(mf$img))=="Image"){
    img.name<-as.character(mf$img)
  }else{
    img.name<-eval(mf$img)
  }

  if (!class(img)=="Image"){
    img = readImage(img)
  }
  #crop
  if (!is.null(crop)) img <- img[crop[1]:crop[2], crop[3]:crop[4],1:3]
  #resize
  if (is.null(resizw)) resizw <- dim(img)[1]
  img<-resize(img,w = resizw)

  # Make rgb matrix
  #rgb<-rbind(c(img[,,1]),c(img[,,2]),c(img[,,3]))
  rgb<-t(matrix(c(c(img[,,1]),c(img[,,2]),c(img[,,3])),nrow = 3, byrow = T))
  # convert to lab
  x <- farver::convert_colour(rgb, from = "rgb", to = "lab")
  # get a component
  xa<-matrix(x[,2],nrow = resizw)
  xht<-xa-min(xa)

  # a threshold
  if (athres=="auto") athres <- otsu(xht, range = c(min(xht),max(xht)))
  xht[xht<=athres]<-0
  xht[xht>athres]<-1

  # make it an image
  xhti<-as.Image(xht)

  # get mask
  nmask = fillHull(xhti)
  nmask = bwlabel(nmask)

  # filter on minsize
  #nmask[!nmask%in%(which(table(c(imageData(nmask)))>minsize)-1)]<-0
  nmask[!nmask%in%data.table(c(imageData(nmask)))[,.N,V1][N>minsize,V1]]<-0

  # filter on maxsize
  #nmask[!nmask%in%(names(which(table(c(imageData(nmask)))<maxsize)))]<-0
  nmask[!nmask%in%data.table(c(imageData(nmask)))[,.N,V1][N<maxsize,V1]]<-0

    # get contours
  nmask.cont<-ocontour(nmask)
  # get features
  nmasks <- EBImage:::splitObjects(nmask)
  nmask.shp <- cbind(s.area = sapply(nmasks, length),
                     do.call(rbind, lapply(nmask.cont, function(z) {
                       cz = apply(z, 2, mean)
                       radius = sqrt(rowSums((z - rep(cz, each = nrow(z)))^2))
                       radius.mean = mean(radius)
                       radius.sd = sqrt(mean((radius - radius.mean)^2))
                       c(s.perimeter = nrow(z), s.radius.mean = radius.mean,
                         s.radius.sd = radius.sd, s.radius.min = min(radius),
                         s.radius.max = max(radius))
                       }))
  )

  nmask.mom<-computeFeatures.moment(nmask)
  nmask.bas<-computeFeatures.basic(nmask,img)
  #nmask.basr<-computeFeatures.basic(nmask,channel(img,"r"))
  #nmask.basg<-computeFeatures.basic(nmask,channel(img,"g"))
  #nmask.basb<-computeFeatures.basic(nmask,channel(img,"b"))
  nmask.basr<-computeFeatures.basic(nmask,img[,,1])
  nmask.basg<-computeFeatures.basic(nmask,img[,,2])
  nmask.basb<-computeFeatures.basic(nmask,img[,,3])
  cols<-rgb(nmask.basr[,1],nmask.basg[,1],nmask.basb[,1])

  ret<-list(features= data.table(nmask.shp,
                                 nmask.mom,
                                 nmask.bas,
                                 cols),
            contours=nmask.cont,
            params=list(crops=crop,
                        w=dim(img)[1],
                        h=dim(img)[2],
                        athres=athres,
                        minsize=minsize,
                        image.name=img.name))
  if (save.outline){
    writeImage(paintObjects(nmask,img,thick = resizw>1500),
               files=paste0(img.name,".outline.jpg"),
               type="jpeg",
               quality = 80)
  }
  class(ret) <- c(class(ret), "krnel")
  return(ret)
}
