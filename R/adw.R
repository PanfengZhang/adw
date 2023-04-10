
# degree is converted to radian
deg2rad <- function(degree) {
  return(degree * pi / 180.0)
}

# great circle distance (http://www.movable-type.co.uk/scripts/latlong.html)
# This uses the ‘haversine’ formula to calculate the great-circle distance between two points.
# According to Yoder (1995), the radius of the earth is 6371.01 km. s2_earth_radius_meters()
# Yoder, C.F. 1995. "Astrometric and Geodetic Properties of Earth and the Solar System" in Global Earth Physics, A Handbook of Physical Constants, AGU Reference Shelf 1, American Geophysical Union, Table 2. doi:10.1029/RF001p0001
distance_gcd <- function(lon1, lat1, lon2, lat2) {
  radius <- 6371.01 # earth radius, unit: km
  # pi <- 3.1415926
  rlon1 <- deg2rad(lon1)
  rlat1 <- deg2rad(lat1)
  rlon2 <- deg2rad(lon2)
  rlat2 <- deg2rad(lat2)
  
  dLon = abs(rlon1 - rlon2)
  dLat = abs(rlat1 - rlat2)
  
  a = sin(dLat/2.0)^2.0 + cos(rlat1) * cos(rlat2) * sin(dLon/2.0)^2.0
  cc = 2.0 * atan2(sqrt(a), sqrt(1-a))
  # cc = 2.0 * asin(sqrt(a))
  return(radius * cc)
}

# bearing angle (http://www.movable-type.co.uk/scripts/latlong.html)
# input: degree; output: rad
bearingR <- function(lon1, lat1, lon2, lat2) {
  rlon1 = deg2rad(lon1)
  rlat1 = deg2rad(lat1)
  rlon2 = deg2rad(lon2)
  rlat2 = deg2rad(lat2)
  y = sin(rlon2-rlon1) * cos(rlat2)
  x = cos(rlat1)*sin(rlat2) - sin(rlat1) * cos(rlat2) * cos(rlon2-rlon1)
  return(atan2(y,x))
}

# calculating weight
weightCal <- function(stnDistances, corDecayDistance, m,
                      grdLon, grdLat, stnLon, stnLat) {
  r <- exp(-stnDistances/corDecayDistance)
  f <- r^m
  theta <- bearingR(grdLon, grdLat, stnLon, stnLat)
  npts <- length(stnDistances)
  alpha <- rep(NA, npts)
  for (k in 1:npts) {
    diffTheta <- theta[-k] - theta[k]
    alpha[k] <- sum(f[-k] * (1 - cos(diffTheta))) / sum(f[-k])
  }
  w <- f * (1 + alpha)
  return(w)
}

#' @title Angular Distance Weighting Interpolation, adw_vector.
#' @description
#' The irregularly-spaced data are interpolated onto regular latitude-longitude 
#' grids by weighting each station according to its distance and angle from the 
#' center of a search radius.
#' @param ds a input dataframe which contains the column names of lon, lat, value.
#' @param extent a extent numeric vector of length 4 in the order c(xmin, xmax, ymin, ymax),
#' or a polygon object with class 'sf' (package 'sf'),  or a polygon object with
#' class 'SpatVector' (package 'terra'). Assume that the coordinate reference 
#' system is WGS1984 (EPSG: 4326).
#' @param gridsize the grid size, i.e. the grid resolution. units: degree.
#' @param cdd correlation decay distance, i.e. the maximum search radius. 
#' unit: kilometer. default value: 1000km.
#' @param m is used to adjust the weighting function further, higher values of m
#' increase the rate at which the weight decays with distance. default value 4.
#' @param nmin the minimum number of observation points required to interpolate 
#' a grid within the search radius (i.e. cdd); if the number of stations within 
#' the search ridius (cdd) is less than nmin, a missing value will be generated  
#' to fill this grid. default value 3.
#' @param nmax The number of nearest points within the search radius to use for 
#' interpolation. default value 10.
#' @return a regular latitude-longitude dataframe grid (interpoled values).
#' @references Caesar, J., L. Alexander, and R. Vose, 2006: Large-scale changes in observed daily maximum and minimum temperatures: Creation and analysis of a new gridded data set. Journal of Geophysical Research, 111, https://doi.org/10.1029/2005JD006280.
#' @examples
#' set.seed(2)
#' dd <- data.frame(lon = runif(100, min = 110, max = 117),
#'                  lat = runif(100, min = 31, max = 37),
#'                  value = runif(100, min = -10, max = 10))
#' head(dd)
#' # example
#' grd <- adw(dd, extent = c(110, 117, 31, 37), gridsize = 0.5, cdd = 500)
#' head(grd)
#' @export
adw_vector <- function(ds, extent, gridsize = 5, cdd = 1e3, m = 4, nmin = 3, 
                       nmax = 10) {
  xmin = extent[1]
  xmax = extent[2]
  ymin = extent[3]
  ymax = extent[4]
  
  dg <- expand.grid(lon = seq(xmin+gridsize/2, xmax, gridsize), 
                    lat = seq(ymin+gridsize/2, ymax, gridsize))
  dg[, "value"] <- NA
  ngrds <- nrow(dg)
  for (j in 1:ngrds) {
    # ds[, "distance"] <- geosphere::distHaversine(dg[j, c("lon", "lat")], ds[, c("lon", "lat")])
    ds[, "distance"] <- distance_gcd(dg$lon[j], dg$lat[j], ds$lon, ds$lat)
    dx <- ds[ds$distance < cdd, ]
    npts <- nrow(dx)  # station points numbers within the searh radius (cdd).
    if (npts > nmax) {
      dx <- dx[order(dx$distance), ]
      dx <- dx[1:nmax, ]
      w <- weightCal(dx$distance, cdd, m, dg$lon[j], dg$lat[j], dx$lon, dx$lat)
      dg[j, "value"] <- sum(dx$value * w) / sum(w)
      
    } else if (npts >= nmin & npts <= nmax) {
      w <- weightCal(dx$distance, cdd, m, dg$lon[j], dg$lat[j], dx$lon, dx$lat)
      dg[j, "value"] <- sum(dx$value * w) / sum(w)
    } else {
      next
    }
  }
  return(dg)
}

#' @title Angular Distance Weighting Interpolation, adw_sf.
#' @description
#' The irregularly-spaced data are interpolated onto regular latitude-longitude 
#' grids by weighting each station according to its distance and angle from the 
#' center of a search radius.
#' @param ds a input dataframe which contains the column names of lon, lat, value.
#' @param extent a polygon object with class 'sf' (package 'sf'). Assume that 
#' the coordinate reference system is WGS1984 (EPSG: 4326).
#' @param gridsize the grid size, i.e. the grid resolution. units: degree.
#' @param cdd correlation decay distance, i.e. the maximum search radius. 
#' unit: kilometer. default value: 1000km.
#' @param m is used to adjust the weighting function further, higher values of m
#' increase the rate at which the weight decays with distance. default value 4.
#' @param nmin the minimum number of observation points required to interpolate 
#' a grid within the search radius (i.e. cdd); if the number of stations within 
#' the search ridius (cdd) is less than nmin, a missing value will be generated  
#' to fill this grid. default value 3.
#' @param nmax The number of nearest points within the search radius to use for 
#' interpolation. default value 10.
#' @return a regular latitude-longitude dataframe grid (interpoled values).
#' @references Caesar, J., L. Alexander, and R. Vose, 2006: Large-scale changes in observed daily maximum and minimum temperatures: Creation and analysis of a new gridded data set. Journal of Geophysical Research, 111, https://doi.org/10.1029/2005JD006280.
#' @examples
#' set.seed(2)
#' dd <- data.frame(lon = runif(100, min = 110, max = 117),
#'                  lat = runif(100, min = 31, max = 37),
#'                  value = runif(100, min = -10, max = 10))
#' head(dd)
#' urlmap <- "https://geo.datav.aliyun.com/areas_v3/bound/410000.json"
#' hmap <- sf::read_sf(urlmap, as_tibble = FALSE) |> sf::st_make_valid() # return a 'sf' object.
#' grd <- adw_sf(dd, extent = hmap, gridsize = 0.5, cdd = 500)
#' head(grd)
#' @importFrom methods is
#' @importFrom sf st_bbox st_as_sf
#' @export
adw_sf <- function(ds, extent, gridsize = 5, cdd = 1e3, m = 4, nmin = 3, nmax = 10) {
  #require(sf)
  bbox <- st_bbox(extent)
  xmin = bbox['xmin']
  xmax = bbox['xmax']
  ymin = bbox['ymin']
  ymax = bbox['ymax']

  dg <- expand.grid(lon = seq(xmin+gridsize/2, xmax, gridsize), 
                    lat = seq(ymin+gridsize/2, ymax, gridsize)) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  dg <- dg[extent, ]
  dg[, "value"] <- NA
  dg <- as.data.frame(dg)
  dg <- dg[, c("lon", "lat", "value")]
  ngrds <- nrow(dg)
  for (j in 1:ngrds) {
    # ds[, "distance"] <- geosphere::distHaversine(dg[j, c("lon", "lat")], ds[, c("lon", "lat")])
    ds[, "distance"] <- distance_gcd(dg$lon[j], dg$lat[j], ds$lon, ds$lat)
    dx <- ds[ds$distance < cdd, ]
    npts <- nrow(dx)  # station points numbers within the searh radius (cdd).
    if (npts > nmax) {
      dx <- dx[order(dx$distance), ]
      dx <- dx[1:nmax, ]
      w <- weightCal(dx$distance, cdd, m, dg$lon[j], dg$lat[j], dx$lon, dx$lat)
      dg[j, "value"] <- sum(dx$value * w) / sum(w)
      
    } else if (npts >= nmin & npts <= nmax) {
      w <- weightCal(dx$distance, cdd, m, dg$lon[j], dg$lat[j], dx$lon, dx$lat)
      dg[j, "value"] <- sum(dx$value * w) / sum(w)
    } else {
      next
    }
  }
  return(dg)
}

#' @title Angular Distance Weighting Interpolation, adw_terra.
#' @description
#' The irregularly-spaced data are interpolated onto regular latitude-longitude 
#' grids by weighting each station according to its distance and angle from the 
#' center of a search radius.
#' @param ds a input dataframe which contains the column names of lon, lat, value.
#' @param extent a polygon object with class 'SpatVector' (package 'terra'). 
#' Assume that the coordinate reference system is WGS1984 (EPSG: 4326).
#' @param gridsize the grid size, i.e. the grid resolution. units: degree.
#' @param cdd correlation decay distance, i.e. the maximum search radius. 
#' unit: kilometer. default value: 1000km.
#' @param m is used to adjust the weighting function further, higher values of m
#' increase the rate at which the weight decays with distance. default value 4.
#' @param nmin the minimum number of observation points required to interpolate 
#' a grid within the search radius (i.e. cdd); if the number of stations within 
#' the search ridius (cdd) is less than nmin, a missing value will be generated  
#' to fill this grid. default value 3.
#' @param nmax The number of nearest points within the search radius to use for 
#' interpolation. default value 10.
#' @return a regular latitude-longitude dataframe grid (interpoled values).
#' @references Caesar, J., L. Alexander, and R. Vose, 2006: Large-scale changes in observed daily maximum and minimum temperatures: Creation and analysis of a new gridded data set. Journal of Geophysical Research, 111, https://doi.org/10.1029/2005JD006280.
#' @examples
#' set.seed(2)
#' dd <- data.frame(lon = runif(100, min = 110, max = 117),
#'                  lat = runif(100, min = 31, max = 37),
#'                  value = runif(100, min = -10, max = 10))
#' head(dd)
#' # example
#' urlmap <- "https://geo.datav.aliyun.com/areas_v3/bound/410000.json"
#' hmap <- terra::vect(urlmap) # return a 'SpatVector' object.
#' grd <- adw(dd, extent = hmap, gridsize = 0.5, cdd = 500)
#' head(grd)
#' @importFrom terra ext mask
#' @importFrom methods is
#' @export
adw_terra <- function(ds, extent, gridsize = 5, cdd = 1e3, m = 4, nmin = 3, nmax = 10) {
  #require(terra)
  bbox <- terra::ext(extent)
  xmin = bbox[1]
  xmax = bbox[2]
  ymin = bbox[3]
  ymax = bbox[4]
  dg <- expand.grid(lon = seq(xmin+gridsize/2, xmax, gridsize), 
                    lat = seq(ymin+gridsize/2, ymax, gridsize)) |>
    terra::vect(crs = "+proj=longlat +datum=WGS84", keepgeom = TRUE)
  dg <- terra::mask(dg, extent)
  dg[, "value"] <- NA
  dg <- as.data.frame(dg)
  ngrds <- nrow(dg)
  for (j in 1:ngrds) {
    # ds[, "distance"] <- geosphere::distHaversine(dg[j, c("lon", "lat")], ds[, c("lon", "lat")])
    ds[, "distance"] <- distance_gcd(dg$lon[j], dg$lat[j], ds$lon, ds$lat)
    dx <- ds[ds$distance < cdd, ]
    npts <- nrow(dx)  # station points numbers within the searh radius (cdd).
    if (npts > nmax) {
      dx <- dx[order(dx$distance), ]
      dx <- dx[1:nmax, ]
      w <- weightCal(dx$distance, cdd, m, dg$lon[j], dg$lat[j], dx$lon, dx$lat)
      dg[j, "value"] <- sum(dx$value * w) / sum(w)
      
    } else if (npts >= nmin & npts <= nmax) {
      w <- weightCal(dx$distance, cdd, m, dg$lon[j], dg$lat[j], dx$lon, dx$lat)
      dg[j, "value"] <- sum(dx$value * w) / sum(w)
    } else {
      next
    }
  }
  return(dg)
}


#' @title Angular Distance Weighting Interpolation.
#' @description
#' The irregularly-spaced data are interpolated onto regular latitude-longitude 
#' grids by weighting each station according to its distance and angle from the 
#' center of a search radius.
#' @param ds a input dataframe which contains the column names of lon, lat, value.
#' @param extent a extent numeric vector of length 4 in the order c(xmin, xmax, ymin, ymax),
#' or a polygon object with class 'sf' (package 'sf'),  or a polygon object with
#' class 'SpatVector' (package 'terra'). Assume that the coordinate reference 
#' system is WGS1984 (EPSG: 4326).
#' @param gridsize the grid size, i.e. the grid resolution. units: degree.
#' @param cdd correlation decay distance, i.e. the maximum search radius. 
#' unit: kilometer. default value: 1000km.
#' @param m is used to adjust the weighting function further, higher values of m
#' increase the rate at which the weight decays with distance. default value 4.
#' @param nmin the minimum number of observation points required to interpolate 
#' a grid within the search radius (i.e. cdd); if the number of stations within 
#' the search ridius (cdd) is less than nmin, a missing value will be generated  
#' to fill this grid. default value 3.
#' @param nmax The number of nearest points within the search radius to use for 
#' interpolation. default value 10.
#' @return a regular latitude-longitude dataframe grid (interpoled values).
#' @references Caesar, J., L. Alexander, and R. Vose, 2006: Large-scale changes in observed daily maximum and minimum temperatures: Creation and analysis of a new gridded data set. Journal of Geophysical Research, 111, https://doi.org/10.1029/2005JD006280.
#' @examples
#' set.seed(2)
#' dd <- data.frame(lon = runif(100, min = 110, max = 117),
#'                  lat = runif(100, min = 31, max = 37),
#'                  value = runif(100, min = -10, max = 10))
#' head(dd)
#' 
#' # example 1
#' grd <- adw(dd, extent = c(110, 117, 31, 37), gridsize = 0.5, cdd = 500)
#' head(grd)
#' 
#' # example 2
#' urlmap <- "https://geo.datav.aliyun.com/areas_v3/bound/410000.json"
#' hmap <- sf::read_sf(urlmap, as_tibble = FALSE) |> sf::st_make_valid() # return a 'sf' object.
#' grd <- adw_sf(dd, extent = hmap, gridsize = 0.5, cdd = 500)
#' head(grd)
#' 
#' # example 3
#' urlmap <- "https://geo.datav.aliyun.com/areas_v3/bound/410000.json"
#' hmap <- terra::vect(urlmap) # return a 'SpatVector' object.
#' grd <- adw(dd, extent = hmap, gridsize = 0.5, cdd = 500)
#' head(grd)
#' @importFrom methods is
#' @export
adw <- function(ds, extent, gridsize = 5, cdd = 1e3, m = 4, nmin = 3, nmax = 10) {
  if (methods::is(extent, "vector")) {
    grd <- adw_vector(ds, extent, gridsize, cdd, m, nmin, nmax)
  } else if (methods::is(extent, "sf")) {
    grd <- adw_sf(ds, extent, gridsize, cdd, m, nmin, nmax)
  } else if (methods::is(extent, "SpatVector")) {
    grd <- adw_terra(ds, extent, gridsize, cdd, m, nmin, nmax)
  }
  return(grd)
}