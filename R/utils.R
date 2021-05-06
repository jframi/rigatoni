#' Rotating calipers algorithm
#'
#' @description
#' Calculates the most rectangle bounding box using the
#' rotating calipers algorithm.
#' Derived from the [flightplanning::getMinBBox] function from the flightplanning package <https://github.com/caiohamamura/flightplanning-R>
#'
#' @param xy A matrix of xy values from which to calculate the bounding box.
#' @noRd
#' @importFrom grDevices chull
getBBox <- function(xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)

  ## rotating calipers algorithm using the convex hull
  H    <- grDevices::chull(xy)      ## hull indices, vertices ordered clockwise
  n    <- length(H)      ## number of hull vertices
  hull <- xy[H, ]        ## hull vertices

  ## unit basis vectors for all subspaces spanned by the hull edges
  hDir  <- diff(rbind(hull, hull[1, ])) ## hull vertices are circular
  hLens <- sqrt(rowSums(hDir^2))        ## length of basis vectors
  huDir <- diag(1/hLens) %*% hDir       ## scaled to unit length

  ## unit basis vectors for the orthogonal subspaces
  ## rotation by 90 deg -> y' = x, x' = -y
  ouDir <- cbind(-huDir[ , 2], huDir[ , 1])

  ## project hull vertices on the subspaces spanned by the hull edges, and on
  ## the subspaces spanned by their orthogonal complements - in subspace coords
  projMat <- rbind(huDir, ouDir) %*% t(hull)

  ## range of projections and corresponding width/height of bounding rectangle
  rangeH  <- matrix(numeric(n*2), ncol=2)  ## hull edge
  rangeO  <- matrix(numeric(n*2), ncol=2)  ## orthogonal subspace
  widths  <- numeric(n)
  heights <- numeric(n)

  for(i in seq(along=numeric(n))) {
    rangeH[i, ] <- range(projMat[  i, ])

    ## the orthogonal subspace is in the 2nd half of the matrix
    rangeO[i, ] <- range(projMat[n+i, ])
    widths[i]   <- abs(diff(rangeH[i, ]))
    heights[i]  <- abs(diff(rangeO[i, ]))
  }

  ## extreme projections for min-area rect in subspace coordinates
  ### hull edge leading to minimum-area
  #eMin  <- which.min(widths*heights)
  ## hull edge leading to maximum difference between width and height
  eMin  <- which.max(abs(widths-heights))
  hProj <- rbind(   rangeH[eMin, ], 0)
  oProj <- rbind(0, rangeO[eMin, ])

  ## move projections to rectangle corners
  hPts <- sweep(hProj, 1, oProj[ , 1], "+")
  oPts <- sweep(hProj, 1, oProj[ , 2], "+")

  ## corners in standard coordinates, rows = x,y, columns = corners
  ## in combined (4x2)-matrix: reverse point order to be usable in polygon()
  ## basis formed by hull edge and orthogonal subspace
  basis <- cbind(huDir[eMin, ], ouDir[eMin, ])
  hCorn <- basis %*% hPts
  oCorn <- basis %*% oPts
  pts   <- t(cbind(hCorn, oCorn[ , c(2, 1)]))

  ## angle of longer edge pointing up
  dPts <- diff(pts)
  e    <- dPts[which.max(rowSums(dPts^2)), ] ## one of the longer edges
  eUp  <- e * sign(e[2])       ## rotate upwards 180 deg if necessary
  deg  <- atan2(eUp[2], eUp[1])*180 / pi     ## angle in degrees

  return(list(pts=pts, width=heights[eMin], height=widths[eMin], angle=deg))
}


#' Convert a krnel object into a coo list from Momocs2 package
#'
#' @param x a krnel object as returned by krnel function
#'
#' @return a coo list
#'
#' @export
#' @import Momocs2
#' @import tibble
#' @examples
kernl_to_coo_list <- function(x){
  coo_list(lapply(x$contours, function(a) new_coo_single(tibble(x=a[,1],y=a[,2]))))
}



#' Check if a point is in a rectangle
#'
#' @param pt
#' @param r
#'
#' @return
#' @noRd
#'
#' @examples
is_p_in_rectangle <- function(pt, r){
  xA <- round(r[1,1],0)
  xB <- round(r[2,1],0)
  xC <- round(r[3,1],0)
  xD <- round(r[4,1],0)
  yA <- round(r[1,2],0)
  yB <- round(r[2,2],0)
  yC <- round(r[3,2],0)
  yD <- round(r[4,2],0)
  xP <- round(pt[1],0)
  yP <- round(pt[2],0)

  ABCD = 0.5 * abs((yA - yC)*(xD - xB) + (yB - yD)*(xA - xC))

  ABP = 0.5 * abs(xA*(yB - yP) + xB*(yP - yA) + xP*(yA - yB))
  BCP = 0.5 * abs(xB*(yC - yP) + xC*(yP - yB) + xP*(yB - yC))
  CDP = 0.5 * abs(xC*(yD - yP) + xD*(yP - yC) + xP*(yC - yD))
  DAP = 0.5 * abs(xD*(yA - yP) + xA*(yP - yD) + xP*(yD - yA))
  #return(round(ABCD,3) == round(ABP + BCP + CDP + DAP, 3))
  #return(abs(ABCD - (ABP + BCP + CDP + DAP)) <2)
  return(ABCD == (ABP + BCP + CDP + DAP))
}
