#' @title Angular Distance Weighting for land
#' @description
#' The irregularly-spaced data are interpolated onto regular latitude-longitude grids by weighting each station according to its distance and angle from the center of a search radius.
#' @param dd a input dataframe which contains column names of lon, lat, value
#' @param xmin the minimum longitude of the rectangular mesh
#' @param xmax the maximum longitude of the rectangular mesh
#' @param ymin the minimum latitude of the rectangular mesh
#' @param ymax the maximum latitude of the rectangular mesh
#' @param gridSize the grid resolution
#' @param cdd the correlation decay distance, unit: meter
#' @param m is used to adjust the weighting function further
#' @return a regular latitude-longitude grid dataframe
#' @examples
#' set.seed(123)
#' dd <- data.frame(lon = runif(100, min = 110, max = 117),
#'                  lat = runif(100, min = 31, max = 37),
#'                  value = runif(100, min = -10, max = 10))
#' head(dd)
#' dg <- adw(dd, gridSize = 1, cdd = 1e5)
#' # dg is the dataframe of grid (mesh)
#' head(dg)
#' @importFrom sf st_as_sf st_buffer st_coordinates st_distance st_geometry sf_use_s2
#' @importFrom magrittr %>%
#' @importFrom geosphere bearing
#' @importFrom stats na.omit
#' @importFrom rnaturalearth ne_countries
#' @export

# mask ocean
adw_land <- function(dd, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
                gridSize = 1, cdd = 1000000, m = 4) {
  requireNamespace("sf")
  requireNamespace("geosphere")
  requireNamespace("magrittr")
  requireNamespace("rnaturalearth")
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
  # sea mask
  sf_use_s2(FALSE) # Spherical geometry (s2) switched off
  wld_land <- ne_countries(returnclass = "sf")
  dg <- dg[wld_land,]  # land sea mask, removing the gridboxes in the sea
  sf_use_s2(TRUE)  # Spherical geometry (s2) switched on
  ngrds <- nrow(dg)
  ds <- na.omit(dd) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)
  for (j in 1:ngrds) {
    circle <- st_buffer(dg[j, ], dist = cdd) %>% st_geometry()
    dx <- ds[circle,]
    npts <- nrow(dx)  # station points numbers in the searh radius
    if (npts > 10) {
      dx[, "distance"] <- st_distance(dg[j,], dx) %>% as.numeric() # Units: [m]
      dx <- dx[order(dx$distance), ]
      dx <- dx[1:10, ]
      r <- exp(1)^(-dx$distance/cdd)
      f <- r^m
      theta <- bearing(c(dg$lon[j], dg$lat[j]), st_coordinates(dx))
      theta <- theta * pi / 180
      alpha <- rep(NA, 10)
      
      for (k in 1:10) {
        diffTheta <- theta[-k] - theta[k]
        alpha[k] <- sum(f[-k] * (1 - cos(diffTheta))) / sum(f[-k])
      }
      w <- f * (1 + alpha)
      dg[j, "value"] <- sum(dx$value * w) / sum(w)
    }
    else if (npts >= 3 & npts <= 10) {
      dx[, "distance"] <- st_distance(dg[j,], dx) %>% as.numeric() # Units: [m]
      r <- exp(1)^(-dx$distance/cdd)
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
