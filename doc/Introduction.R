## ---- include = FALSE, message=FALSE------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("adw")

## ---- message=FALSE-----------------------------------------------------------
library(sf)
library(ggplot2)
library(adw)
head(tavg)
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

## -----------------------------------------------------------------------------
library(adw)
mapurl <- "https://geo.datav.aliyun.com/areas_v3/bound/410000.json"
hmap_sf <- read_sf(mapurl) |> st_make_valid()
dg <- adw(tavg, extent = hmap_sf, gridsize = 0.1, cdd = 400)
head(dg)
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

## -----------------------------------------------------------------------------
library(adw)
library(terra)
mapurl <- "https://geo.datav.aliyun.com/areas_v3/bound/410000.json"
hmap_terra <- terra::vect(mapurl)
dg <- adw(tavg, extent = hmap_terra, gridsize = 0.1, cdd = 400)
head(dg)
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

## -----------------------------------------------------------------------------
library(adw)
interpExtent <- c(110.36, 116.65, 31.38, 36.37) # [xmin, xmax, ymin, ymax]
dg <- adw(tavg, extent = interpExtent, gridsize = 0.1, cdd = 400)
head(dg)
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

