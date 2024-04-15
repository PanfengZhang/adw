## ----include = FALSE, message=FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("adw")

## ----message=FALSE------------------------------------------------------------
library(sf)
library(ggplot2)
library(adw)
library(cnmap)

set.seed(1)
tavg <- data.frame(lon = runif(100, min = 110, max = 117),
                   lat = runif(100, min = 31, max = 37),
                   value = runif(100, min = 20, max = 35))

hmap <- getMap(name = "河南省", returnClass = "sf")
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
hmap_sf <- getMap(name = "河南省", returnClass = "sf") |> st_make_valid()
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
hmap_sv <- getMap(name = "河南省", returnClass = "sv")
dg <- adw(tavg, extent = hmap_sv, gridsize = 0.1, cdd = 400)
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

## -----------------------------------------------------------------------------
library(adw)
interpExtent <- c(110.36, 116.65, 31.38, 36.37) # [xmin, xmax, ymin, ymax]
dg <- points2grid(tavg, extent = interpExtent, gridsize = 0.5)
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
dg <- na.omit(dg)
awa(dg$value, dg$lat)

