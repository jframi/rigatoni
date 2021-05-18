#' main function to identify features in an image
#'
#' @param img the path to the image file to analyze or an Image object as defined by the EBImage package
#' @param crop a numeric vector in the form c(xmin, xmax, ymin, ymax) to define cropping area
#' @param resizw the width of resized image to use for feature detection. decreasing the size of the image will improve performance
#' @param huethres a numeric vector in the form c(Huemin, Huemax) to define Hue threshold. Can be easily determined using the get_huethr_values function
#' @param minsize the minimum size of features to be considered, in pixel, as determined on the resized image
#' @param maxsize the maximum size of features to be considered, in pixel, as determined on the resized image
#' @param save.outline boolean, should an outline image be saved in the same directory
#' @param img.name name of the image. This is useful in the case img is an Image object and save.outline is TRUE. This will be used for the name of the outline file.
#' @param blackbg behave differently in case the background is black
#' @param ws.avg watershed average : this allows to perform a watershed after a feature detection with no watershed, and to compute for each feature bbox width/height average and sum on the different sub-features obtained through watershed
#'
#' @return the function will return krnel object, ie. a list with three components
#' * features : a data.frame with as many rows as the number of detected features and description variables as returned by EBImage::computeFeatures. It also includes bounding box width and height aka Feret min and max diameter.
#' * contours : a list with as many elements as the number of detected features. Each element is a matrix with the coordinates of each feature.
#' * bbox : a list with as many elements as the number of detected features. Each element is a list with 4 components: $pts that contains the coordinates of each corner of the bounding box, $width, $height, and $angle
#' * params : a list with the analysis parameters. can be used for further plotting
#' * ws.contours : if ws.avg=T. a list as many elements as the number of detected main features. Each elementis a list of matrix with the coordinates of each sub-features.
#' * ws.bbox : if ws.avg=T. The same with bounding boxes
#' * ws.pois : if ws.avg=T. The same with point of inaccessibility
#' @md
#' @importFrom polylabelr poi
#' @export
#' @examples
krnel<- function(img, crop=NULL, resizw=NULL, watershed=F, huethres, minsize, maxsize, save.outline=F, img.name=NULL, blackbg=F, ws.avg=F){

  #mf<-match.call()
  #if(class(eval(mf$img))=="Image"){
  #  img.name<-as.character(mf$img)
  #}else{
  #  img.name<-eval(mf$img)
  #}

  if (!class(img)=="Image"){
    img.name <- img
    img = suppressWarnings(readImage(img))
  } else {
    if (save.outline==T & is.null(img.name)) stop("If img is an EBImage image object and save.outline=T you should provide an image.name for the output file")
  }
  #crop
  if (!is.null(crop)){
    img <- img[crop[1]:crop[2], crop[3]:crop[4],1:3]
  }
  #resize
  orig.w <- dim(img)[1]
  orig.h <- dim(img)[2]
  if (!is.null(resizw)){
    img<-resize(img,w = resizw)
  } else {
    resizw <- dim(img)[1]
  }

  if (blackbg == TRUE){
    img@colormode <- Grayscale
    xhti <- (1-(img>otsu(img)))
    if (numberOfFrames(img, type="total")>1) {
      xhti <- xhti[,,1]
    }
  }else{
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
  }

  # get mask
  nmask = fillHull(1-xhti)
  nmask = bwlabel(nmask)

  # filter on minsize
  #nmask[!nmask%in%(which(table(c(imageData(nmask)))>minsize)-1)]<-0
  nmask <- rmObjects(nmask,which(computeFeatures.shape(nmask)[,"s.area"]<minsize))
  ## filter on maxsize
  #nmask[!nmask%in%(names(which(table(c(imageData(nmask)))<maxsize)))]<-0

  if (watershed){
    dmap = distmap(nmask)
    nmask = watershed(dmap)
    # filter on minsize
    #nmask[!nmask%in%(which(table(c(imageData(nmask)))>minsize)-1)]<-0
    nmask <- rmObjects(nmask,which(computeFeatures.shape(nmask)[,"s.area"]<minsize))
  }
  # filter on maxsize
  #nmask[!nmask%in%(names(which(table(c(imageData(nmask)))<maxsize)))]<-0
  nmask <- rmObjects(nmask,which(computeFeatures.shape(nmask)[,"s.area"]>maxsize))

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

  # Compute bounding boxes (bbox)
  nmask.bbox <- lapply(nmask.cont, getBBox)
  # bbox width should always be smaller than height
  nmask.bbox.wh <- do.call(rbind,lapply(nmask.bbox, function(a) data.table(bbox.width=ifelse(a$width<a$height, a$width, a$height), bbox.height=ifelse(a$width<a$height, a$height, a$width))))

  # Compute the pole of inaccessibility to get the Maximum inscribed Circle
  pois<-setnames(data.table(do.call(rbind,lapply(lapply(nmask.cont, polylabelr::poi), unlist))), c("poi.x","poi.y","poi.dist"))

  # build the krnel object to return
  ret<-list(features= data.table(id=as.numeric(names(nmask.cont)),
                                 nmask.shp,
                                 nmask.bbox.wh,
                                 pois,
                                 nmask.mom,
                                 nmask.bas,
                                 cols),
            contours=nmask.cont,
            bbox = nmask.bbox,
            params=list(crops=crop,
                        w=dim(img)[1],
                        h=dim(img)[2],
                        orig.w=orig.w,
                        orig.h=orig.h,
                        huethres=huethres,
                        minsize=minsize,
                        maxsize=maxsize,
                        watershed=watershed,
                        save.outline=save.outline,
                        image.name=img.name,
                        blackbg=blackbg,
                        ws.avg=ws.avg))

  # If ws.avg is true, redo an analysis with watershed, link sub-features to features detected withous watershed and compute stats
  if (ws.avg){
    # watershed
    dmap = distmap(nmask)
    nmask.ws = watershed(dmap)
    # filter on minsize after watershed
    nmask.ws <- rmObjects(nmask.ws,which(computeFeatures.shape(nmask.ws)[,"s.area"]<minsize))
    # do measurements on watershed sub-features
    nmask.ws.mom<-computeFeatures.moment(nmask.ws)
    nmask.ws.cont<-ocontour(nmask.ws)
    # Compute bounding boxes (bbox)
    nmask.ws.bbox <- lapply(nmask.ws.cont, getBBox)
    nmask.ws.bbox.wh <- do.call(rbind,lapply(nmask.ws.bbox, function(a) data.table(bbox.width=ifelse(a$width<a$height, a$width, a$height), bbox.height=ifelse(a$width<a$height, a$height, a$width))))
    # Compute the pole of inaccessibility to get the Maximum inscribed Circle
    pois.ws<-setnames(data.table(do.call(rbind,lapply(lapply(nmask.ws.cont, polylabelr::poi), unlist))), c("poi.x","poi.y","poi.dist"))
    # Construct a data.table with sub-features measurements
    ret.ws <- data.table(nmask.ws.mom,nmask.ws.bbox.wh,pois.ws)
    # Find the parent feature using bounding boxes of features detected before watershed and the is_p_in_rectangle applied to centroids of sub-features
    ret.ws[, parent:=apply(ret.ws[,.(m.cx,m.cy)], 1, function(a) which(unlist(lapply(ret$bbox, function(b) is_p_in_rectangle(pt=a, r=b$pts)))))]
    # Compute stats group by parent feature
    ret.ws2 <- ret.ws[, .(bbox.width.ws.avg=mean(bbox.width),
                          bbox.height.ws.avg=mean(bbox.height),
                          bbox.width.ws.sum=sum(bbox.width),
                          bbox.height.ws.sum=sum(bbox.height),
                          poi.dist.max=max(poi.dist),
                          poi.dist.min=min(poi.dist),
                          poi.dist.avg=mean(poi.dist)), parent]
    # Join the new parameters with the ones obtained previously
    #ret$features <- ret$features[ret.ws2, on=c(id="parent")]
    ret$features <- ret.ws2[ret$features, on=c(parent="id")]
    setnames(ret$features, old="parent", new="id")
    ret$features[is.na(bbox.width.ws.avg), `:=`(bbox.width.ws.avg=bbox.width,
                                                bbox.height.ws.avg=bbox.height,
                                                bbox.width.ws.sum=bbox.width,
                                                bbox.height.ws.sum=bbox.height,
                                                poi.dist.max=poi.dist,
                                                poi.dist.min=poi.dist,
                                                poi.dist.avg=poi.dist)]
    # Output the contours of the ws sub-features in final krnel object
    cont.ids <- ret.ws[,.(id=1:.N,parent)][,.(ids=list(id)),parent][order(parent)]
    ret$ws.contours <- lapply(1:nrow(cont.ids), function(a) nmask.ws.cont[unlist(cont.ids[a, ids])])
    # Output the bounding boxes of the ws sub-features in final krnel object
    ret$ws.bbox <- lapply(1:nrow(cont.ids), function(a) nmask.ws.bbox[unlist(cont.ids[a, ids])])
    # Output the points of inaccessibility of the ws sub-features in final krnel object
    ret$ws.pois <- ret.ws[,.(parent, poi.x,poi.y,poi.dist)]
  }
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
#' @noRd
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
