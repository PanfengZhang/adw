## Angular Distance Weighting Interpolation

The irregularly-spaced data are interpolated onto regular latitude-longitude grids by weighting each station according to its distance and angle from the center of a search radius.

## Reference

Caesar, J., L. Alexander, and R. Vose, 2006: Large-scale changes in observed daily maximum and minimum temperatures: Creation and analysis of a new gridded data set. Journal of Geophysical Research, 111, <https://doi.org/10.1029/2005JD006280>.

## Installation

Install the latest CRAN release via command:

```r
install.packages("adw")
```

## load data and plot scatter


```r
library(sf)
library(ggplot2)
library(adw)
head(tavg)
#>      lon   lat value
#> 1 113.82 36.07  26.5
#> 2 114.40 36.05  27.3
#> 3 112.92 35.12  23.5
#> 4 114.18 35.62  29.2
#> 5 112.63 35.08  23.6
#> 6 113.08 35.15  24.5
mapurl <- "https://geo.datav.aliyun.com/areas_v3/bound/410000.json"
hmap <- read_sf(mapurl, as_tibble = FALSE)
ggplot() +
  geom_point(data = tavg, aes(x = lon, y = lat, colour = value), 
             pch = 17, size = 2.5) +
  geom_sf(data = st_cast(hmap, 'MULTILINESTRING')) +
  scale_colour_fermenter(palette = "YlOrRd",
                         direction = 1,
                         breaks = seq(from = 25, to = 32, by = 1),
                         limits = c(0, 100),
                         name = expression("\u00B0C")) +
  ggtitle("The irregularly-spaced data") +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.key.height = unit(1.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 11))
```

![](README_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

## Interpolation
### Usage 1

The parameter *extent* in the **adw** function is a *sf* class (sf package), and the coordinate reference system of the object is WGS1984 (EPSG: 4326).


```r
library(adw)
mapurl <- "https://geo.datav.aliyun.com/areas_v3/bound/410000.json"
hmap_sf <- read_sf(mapurl) |> st_make_valid()
dg <- adw(tavg, extent = hmap_sf, gridsize = 0.1, cdd = 400)
head(dg)
#>          lon      lat    value
#> 50  115.3105 31.43345 32.49991
#> 107 114.7105 31.53345 32.12537
#> 108 114.8105 31.53345 32.20815
#> 109 114.9105 31.53345 32.34560
#> 110 115.0105 31.53345 32.44178
#> 113 115.3105 31.53345 32.53088
ggplot() +
  geom_tile(data = dg, aes(x = lon, y = lat, fill = value)) +
  geom_sf(data = st_cast(hmap_sf, 'MULTILINESTRING')) +
  scale_fill_fermenter(palette = "YlOrRd",
                       direction = 1,
                       breaks = seq(from = 25, to = 32, by = 1),
                       limits = c(0, 100),
                       name = expression("\u00B0C"),
                       na.value = "white") +
  ggtitle("Angular distance weighting interpolation") +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.key.height = unit(1.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 11))
```

![](README_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

### Usage 2

The parameter *extent* in the **adw** function is a *SpatVector* class (terra packag), and the coordinate reference system of the object is WGS1984 (EPSG: 4326).


```r
library(adw)
library(terra)
#> terra 1.7.18
mapurl <- "https://geo.datav.aliyun.com/areas_v3/bound/410000.json"
hmap_terra <- terra::vect(mapurl)
dg <- adw(tavg, extent = hmap_terra, gridsize = 0.1, cdd = 400)
head(dg)
#>        lon      lat    value
#> 1 115.3105 31.43345 32.49991
#> 2 114.7105 31.53345 32.12537
#> 3 114.8105 31.53345 32.20815
#> 4 114.9105 31.53345 32.34560
#> 5 115.0105 31.53345 32.44178
#> 6 115.3105 31.53345 32.53088
ggplot() +
  geom_tile(data = dg, aes(x = lon, y = lat, fill = value)) +
  geom_sf(data = st_cast(hmap_sf, 'MULTILINESTRING')) +
  scale_fill_fermenter(palette = "YlOrRd",
                       direction = 1,
                       breaks = seq(from = 25, to = 32, by = 1),
                       limits = c(0, 100),
                       name = expression("\u00B0C"),
                       na.value = "white") +
  ggtitle("Angular distance weighting interpolation") +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.key.height = unit(1.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 11))
```

![](README_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

### Usage 3

The parameter *extent* in the **adw** function is a extent *vector* of length 4 in the order [xmin, xmax, ymin, ymax]


```r
library(adw)
interpExtent <- c(110.36, 116.65, 31.38, 36.37) # [xmin, xmax, ymin, ymax]
dg <- adw(tavg, extent = interpExtent, gridsize = 0.1, cdd = 400)
head(dg)
#>      lon   lat    value
#> 1 110.41 31.43 30.91942
#> 2 110.51 31.43 30.93326
#> 3 110.61 31.43 30.94749
#> 4 110.71 31.43 31.02596
#> 5 110.81 31.43 31.04038
#> 6 110.91 31.43 31.05498
ggplot() +
  geom_tile(data = dg, aes(x = lon, y = lat, fill = value)) +
  geom_sf(data = st_cast(hmap_sf, 'MULTILINESTRING')) +
  scale_fill_fermenter(palette = "YlOrRd",
                       direction = 1,
                       breaks = seq(from = 25, to = 32, by = 1),
                       limits = c(0, 100),
                       name = expression("\u00B0C"),
                       na.value = "white") +
  ggtitle("Angular distance weighting interpolation") +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.key.height = unit(1.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 11))
```

![](README_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

