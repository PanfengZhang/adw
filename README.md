# Angular distance weighting (adw)

The irregularly-spaced data are interpolated onto a regular
latitude-longitude grid by weighting each station according to its
distance and angle from the center of a search radius

# Reference

Caesar, J., L. Alexander, and R. Vose, 2006: Large-scale changes in
observed daily maximum and minimum temperatures: Creation and analysis
of a new gridded data set. J. Geophys. Res., 111, D05101.

## Installation

The **development** version can be installed from GitHub
(<https://github.com/PanfengZhang/adw>) using:

    # install.packages("remotes")
    remotes::install_github("PanfengZhang/adw")

## Usage

    library(adw)
    library(ggplot2)
    library(sf)

    ## Linking to GEOS 3.9.1, GDAL 3.2.1, PROJ 7.2.1; sf_use_s2() is TRUE

    set.seed(123)
    dd <- data.frame(lon = runif(50, min = 110, max = 117),
                     lat = runif(50, min = 31, max = 37),
                     value = runif(50, min = -10, max = 10))
    head(dd)

    ##        lon      lat      value
    ## 1 112.0130 31.27499  1.9997792
    ## 2 115.5181 33.65320 -3.3435292
    ## 3 112.8628 35.79355 -0.2277393
    ## 4 116.1811 31.73140  9.0894765
    ## 5 116.5833 34.36569 -0.3419521
    ## 6 110.3189 32.23919  7.8070044

    dg <- adw(dd, gridSize = 0.5, cdd = 1e5, m = 4) %>% na.omit()
    head(dg)

    ##         lon      lat     value
    ## 2  110.1723 31.50375 0.6793303
    ## 3  110.1723 32.00375 3.9473573
    ## 4  110.1723 32.50375 6.3224020
    ## 13 110.6723 31.00375 0.7495071
    ## 14 110.6723 31.50375 0.4942404
    ## 15 110.6723 32.00375 0.1070331

    urlmap <- "https://geo.datav.aliyun.com/areas_v3/bound/410000_full.json"
    cmap <- read_sf(urlmap) %>% st_cast('MULTILINESTRING')
    library(ggplot2)
    ggplot() +
      geom_tile(data = dg, aes(x = lon, y = lat, fill = value)) +
      geom_sf(data = cmap) +
      coord_sf(expand = FALSE) +
      guides(fill = guide_coloursteps(title.position = "right")) +
      ggtitle("adw interpolation") +
      scale_fill_fermenter(palette = "RdYlBu", na.value = "white", 
                           breaks = seq(-8, 8, 2),
                           limits = c(-11, 11),
                           name = expression("\u00B0C")) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.key.width = unit(1.5, "cm"),
            legend.key.height = unit(0.3, "cm"),
            axis.title = element_blank(),
            legend.title = element_text(face = "plain", size = 9))

![](README_files/figure-markdown_strict/unnamed-chunk-1-1.png)
