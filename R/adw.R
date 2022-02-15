#' @title adw function, Angular distance weighting interpolation
#'
#' Gridding the irregularly-spaced data onto a regular latitude-longitude grid
#' by weighting each station according to its distance and angle from the centre
#' of a search radius
#'
#' @param dd is the dataframe
#' @param xmin is the minimum of the longitude
#' @param xmax is the minimum of the longitude
#' @param ymin is the minimum of the latitude
#' @param ymax is the maximum of the latitude
#' @param gridSize is the grid resolution
#' @param cdd is the correlation decay distance, unit: meter
#' @param m is used to adjust the weighting function further
#' @return a dataframe
#' @export
#' @importFrom sf st_as_sf st_buffer st_coordinates st_distance st_geometry
#' @importFrom magrittr %>%
#' @importFrom geosphere bearing
#' @importFrom stats na.omit

adw <- function(dd, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
                gridSize = 1, cdd = 1000000, m = 4) {
  # library(sf)
  # library(geosphere)
  # library(magrittr)
  requireNamespace("sf")
  requireNamespace("geosphere")
  requireNamespace("magrittr")
  if(is.null(xmin)) {
    xmin <- min(dd$lon)
    xmax <- max(dd$lon)
    ymin <- min(dd$lat)
    ymax <- max(dd$lat)
  }
  # create grids
  lon <- seq(xmin, xmax, gridSize)    # grid center
  lat <- seq(ymin, ymax, gridSize)    # grid center
  nlon <- length(lon)
  nlat <- length(lat)
  dg <- data.frame(lon = rep(lon, each = nlat),
                   lat = rep(lat, times = nlon)) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  dg[, c("lon", "lat")] <- st_coordinates(dg)
  dg[, "value"] <- NA
  ngrds <- nrow(dg)
  ds <- na.omit(dd) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)
  for (j in 1:ngrds) {
    circle <- st_buffer(dg[j, ], dist = cdd) %>% st_geometry()
    dx <- ds[circle,]
    npts <- nrow(dx)  # station points numbers in the searh radius
    if (npts >= 3) {
      distance <- st_distance(dg[j,], dx) %>% as.numeric()
      r <- exp(1)^(-distance/cdd)
      f <- r^m
      theta <- bearing(c(dg$lon[j], dg$lat[j]), st_coordinates(dx))
      theta <- theta * pi / 180
      alpha <- rep(NA, npts)
      for (k in 1:npts) {
        diffTheta <- theta[-k] - theta[k]
        alpha[k] <- sum(f[-k] * (1 - cos(diffTheta))) / sum(f[-k])
      }
      w <- f * (1 + alpha)
      dg[j, "value"] <- sum(dx$value * w) / sum(w)
    }
  }
  dg <- as.data.frame(dg)
  dg <- dg[, c("lon", "lat", "value")]
  return(dg)
}

