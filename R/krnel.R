#' main function to identify features in an image
#'
#' @param img the path to the image file to analyze or an Image object as defined by the EBImage package
#' @param crop a numeric vector in the form c(xmin, xmax, ymin, ymax) to define cropping area
#' @param resizw the width of resized image to use for feature detection. decreasing the size of the image will improve performance
#' @param huethres a numeric vector in the form c(Huemin, Huemax) to define Hue threshold. Can be easily determined using the get_huethr_values function
#' @param vthres an optional v value threshold under which pixels are considered as part of features (override huesthres values)
#' @param minsize the minimum size of features to be considered, in pixel, as determined on the resized image
#' @param maxsize the maximum size of features to be considered, in pixel, as determined on the resized image
#' @param save.outline boolean, should an outline image be saved in the same directory?
#' @param img.name name of the image. This is useful in the case img is an Image object and save.outline is TRUE. This will be used for the name of the outline file.
#' @param blackbg behave differently in case the background is black
#' @param whitebg behave differently in case the background is white
#' @param bw boolean, is it a black and white (greyscale) image? In this case huethres is not required.
#' @param ws.avg watershed average : this allows to perform a watershed after a feature detection with no watershed, and to compute for each feature bbox width/height average and sum on the different sub-features obtained through watershed
#' @param color.erode boolean, apply morphological erosion to each feature before getting feature's color. This is way to estimate color from the center part of each feature, to avoid border effects.
#' @param colerode.rad.ratio if color.erode is true, this is the size of the round brush used for morphological erosion expressed as a ratio (>0, <=1) of the average radius of all feaures detected on the image
#' @param open apply an open operation (erode+dilate) to suppress possible artefact attached to the detected feature. If watershed=T open is applied before watershed.
#' @param open.brush.size size of the brush to use for the open operation. Brush shape is hardcoded as 'disc'.
#'
#' @return the function will return krnel object, ie. a list with three components
#' * features : a data.frame with as many rows as the number of detected features and description variables as returned by [EBImage::computeFeatures] It also includes : (i) bounding box width and height aka Feret min and max diameter (bbox.width, bbox.height), (ii) pole of inaccessibility coordinates and longest distance to 'coastline' (poi.x,poi.y,poi.dist). In complex shapes, poi.dist might be a good measure of object width.
#' * contours : a list with as many elements as the number of detected features. Each element is a matrix with the coordinates of each feature.
#' * bbox : a list with as many elements as the number of detected features. Each element is a list with 4 components: $pts that contains the coordinates of each corner of the bounding box, $width, $height, and $angle
#' * params : a list with the analysis parameters. can be used for further plotting
#' * ws.contours : if ws.avg=T. a list with as many elements as the number of detected main features. Each element is a list of matrices with contours coordinates of each sub-features.
#' * ws.bbox : if ws.avg=T. The same with bounding boxes
#' * ws.pois : if ws.avg=T. The same with point of inaccessibility
#' @md
#' @importFrom polylabelr poi
#' @importFrom polyclip pointinpolygon
#' @export
#' @examples
krnel<- function(img, crop=NULL, resizw=NULL, watershed=F, huethres, vthres=NULL, minsize, maxsize, save.outline=F, img.name=NULL, blackbg=F, whitebg=F, ws.avg=F, bw=F, color.erode=F, colerode.rad.ratio=0.75, open = FALSE, open.brush.size=5){

  #mf<-match.call()
  #if(class(eval(mf$img))=="Image"){
  #  img.name<-as.character(mf$img)
  #}else{
  #  img.name<-eval(mf$img)
  #}
  if ((bw | blackbg | whitebg) & missing(huethres)){
    huethres <- c(0,0)
  }
  if (color.erode==TRUE & colerode.rad.ratio==0){
    color.erode <- FALSE
  }

  if (!class(img)=="Image"){
    img.name <- img
    img = suppressWarnings(readImage(img))
  } else {
    if (save.outline==T & is.null(img.name)) stop("If img is an EBImage image object and save.outline=T you should provide an image.name for the output file")
  }
  #crop
  if (!is.null(crop)){
    if (numberOfFrames(img)==1){
      img <- img[crop[1]:crop[2], crop[3]:crop[4]]
    } else{
      img <- img[crop[1]:crop[2], crop[3]:crop[4],1:dim(img)[3]]
    }
  }
  #resize
  orig.w <- dim(img)[1]
  orig.h <- dim(img)[2]
  if (!is.null(resizw)){
    img<-resize(img,w = resizw)
  } else {
    resizw <- dim(img)[1]
  }
  if (bw){
      img@colormode <- Grayscale
      if (blackbg == TRUE){
        xhti <- (1-(img>otsu(img)[1]))
      } else {
        xhti <- img>otsu(img)[1]
      }
      if (numberOfFrames(img, type="total")>1) {
        xhti <- xhti[,,1]
      }
    } else {
    if (blackbg | whitebg){
      img2 <- img
      img2@colormode <- Grayscale
      if (whitebg) {
        xhti <- ((img2>otsu(img2)[1]))
      } else {
        xhti <- (1 - (img2>otsu(img2)[1]))
      }
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
      if (!is.null(vthres)){
        xv<-matrix(x[3,],nrow = resizw)
        xht[xv<vthres]<-0
      }
      # make it an image
      xhti<-as.Image(xht)
    }
  }

  # get mask
  nmask = fillHull(1-xhti)
  nmask = bwlabel(nmask)
  nmask.shp <- computeFeatures.shape(nmask)

  # filter on minsize
  #nmask[!nmask%in%(which(table(c(imageData(nmask)))>minsize)-1)]<-0
  nmask <- rmObjects(nmask,which(nmask.shp[,"s.area"]<minsize))
  #nmask <- rmObjects(nmask,data.table(o=c(nmask@.Data))[,.N,o][N<minsize]$o)
  if (open){
    nmask <- dilate(erode(nmask, makeBrush(open.brush.size, shape='disc')), makeBrush(open.brush.size, shape='disc'))
    nmask = bwlabel(nmask)
  }

  if (watershed){
    dmap = distmap(nmask)
    nmask = watershed(dmap)
    # filter on minsize
    #nmask[!nmask%in%(which(table(c(imageData(nmask)))>minsize)-1)]<-0
    nmask <- rmObjects(nmask,which(computeFeatures.shape(nmask)[,"s.area"]<minsize))
    #nmask <- rmObjects(nmask,data.table(o=c(nmask@.Data))[,.N,o][N<minsize]$o)

  }
  # filter on maxsize
  #nmask[!nmask%in%(names(which(table(c(imageData(nmask)))<maxsize)))]<-0
  nmask <- rmObjects(nmask,which(computeFeatures.shape(nmask)[,"s.area"]>maxsize))
  #nmask <- rmObjects(nmask,data.table(o=c(nmask@.Data))[,.N,o][N>maxsize]$o)

  nmask.cont<-ocontour(nmask)

  # get features
  nmask.shp<-computeFeatures.shape(nmask)
  nmask.mom<-computeFeatures.moment(nmask)
  nmask.bas<-computeFeatures.basic(nmask,img)

  if (numberOfFrames(img)>=3){
    if (!color.erode){
      nmask.basr<-computeFeatures.basic(nmask,img[,,1])
      nmask.basg<-computeFeatures.basic(nmask,img[,,2])
      nmask.basb<-computeFeatures.basic(nmask,img[,,3])
      cols<-rgb(nmask.basr[,1],nmask.basg[,1],nmask.basb[,1])

    } else {
      # Erode the nmask image
      # compute the brushsize
      brushsize <- 2*round(min(nmask.shp[,which(colnames(nmask.shp)=="s.radius.mean")])*colerode.rad.ratio/2,0)-1
      # apply morphological erosion to the nmask image
      nmask.er <- bwlabel(erode(nmask, makeBrush(brushsize, shape = "disc")))
      if (watershed){
        dmap = distmap(nmask.er)
        nmask.er = watershed(dmap)
      }

          # Compute features on the eroded image
          nmask.er.shp <- computeFeatures.shape(nmask.er)
          nmask.er.mom <- computeFeatures.moment(nmask.er)
          nmask.basr<-computeFeatures.basic(nmask.er,img[,,1])
          nmask.basg<-computeFeatures.basic(nmask.er,img[,,2])
          nmask.basb<-computeFeatures.basic(nmask.er,img[,,3])

          #xi <- as.Image(x)
          #nmask.bash<-computeFeatures.basic(nmask.er,xi[,,1])
          #nmask.bass<-computeFeatures.basic(nmask.er,xi[,,2])
          #nmask.basv<-computeFeatures.basic(nmask.er,xi[,,3])

          # We need to match the features detected on the eroded image
          # to the ones detected on the non-eroded image
          # we use the centers of the features detected on the eroded image
          # and use the pointinpolygon function and the contours or features detected on non-eroded image

          nmask.er.mom.dt <- data.table(nmask.er.mom)
          nmask.er.shp.dt <- data.table(nmask.er.shp)
          nmask.er.shp.dt[,id:=1:.N]

          # pip contains the result of pointinpolygon function
          # for all contours on the original image
          pip <- lapply(1:length(nmask.cont),
                 function(a) polyclip::pointinpolygon(P=list(x=nmask.er.mom.dt$m.cx,
                                                             y=nmask.er.mom.dt$m.cy),
                                                      A=list(x=nmask.cont[[a]][,1],
                                                             y=nmask.cont[[a]][,2])))
          # pip.dt is a datatable that as all eroded features attributed to each orginal contour, NA if none.
          pip.dt <- do.call(rbind,Map(function(p,n) data.table(eroded_id=ifelse(any(p==1),which(p==1),NA),id=n), pip, 1:length(pip)))
          # pip.dt.j brings the size of eroded features together with pip.dt
          pip.dtj <- setnames(nmask.er.shp.dt[pip.dt, on=.(id=eroded_id)], c("id","i.id"),c("eroded_id","id"))
          # pip.dtjf is the previous one filtered to keep the largest eroded feature for each original contour
          # nrow(pip.dtjf) should always == length(nmask.cont)
          pip.dtjf <- pip.dtj[,.(.N, eroded_id=eroded_id[which.max(s.area)]),id]

      # compute hex colors of all eroded features
      cols<-rgb(nmask.basr[,1],
                nmask.basg[,1],
                nmask.basb[,1])
      # and keep only the ones that have a corresponding id in the orginal image,
      # and sort them according to the order of features in the original image
      cols <- cols[pip.dtjf[order(id)]$eroded_id]

      #colshsv <- hsv(nmask.bash[,1],
      #               nmask.bass[,1],
      #               nmask.basv[,1])
      #colshsv <- colshsv[pip.dtjf[order(id)]$eroded_id]

      if (any(is.na(cols))) warning("Some colors are. Erosion was too string Consider decreasing colerode.rad.ratio")
    }
  } else {
    cols <- rgb(nmask.bas[,1],nmask.bas[,1],nmask.bas[,1])
  }

  # Compute bounding boxes (bbox)
  nmask.bbox <- lapply(nmask.cont, Rigatoni:::getBBox)
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
                                 #colshsv),
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
                        bw=bw,
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
    # ret.ws[, parent:=apply(ret.ws[,.(m.cx,m.cy)], 1, function(a) which(unlist(lapply(ret$bbox, function(b) is_p_in_rectangle(pt=a, r=b$pts)))))]
    whoswhere <- do.call(rbind,
                          Map(function(n,a) data.table(parent=n,id2=a),
                              as.numeric(names(ret$contours)),
                              lapply(ret$contours,
                                     function(a) which(pointinpolygon(c(ret.ws[,.(x=m.cx,y=m.cy)]),
                                                                      list(x=a[,1],y=a[,2]))==1)
                                     )
                              )
                          )
    ret.ws[,id2:=1:.N]
    ret.ws <- whoswhere[ret.ws, on=.(id2)]
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
    if (color.erode){
      outline.img <- paintObjects(nmask.er,paintObjects(nmask,img,thick = resizw>1500), col=c("white",NA), thick = resizw>1500)
    } else {
      outline.img <- paintObjects(nmask,img,thick = resizw>1500)
    }
    writeImage(outline.img,
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
