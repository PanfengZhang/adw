# Angular distance weighting (adw)

The irregularly-spaced data are gridded onto a regular
latitude-longitude grid by weighting each station according to its
distance and angle from the center of a search radius

# Reference

Caesar, J., L. Alexander, and R. Vose, 2006: Large-scale changes in
observed daily maximum and minimum temperatures: Creation and analysis
of a new gridded data set. J. Geophys. Res., 111, D05101,
<https://doi.org/10.1029/2005JD006280>.

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
    dd <- data.frame(lon = runif(100, min = 110, max = 117),
                     lat = runif(100, min = 31, max = 37),
                     value = runif(100, min = -10, max = 10))
    head(dd)

    ##        lon      lat      value
    ## 1 112.0130 34.59993 -5.2254795
    ## 2 115.5181 32.99694  9.2471787
    ## 3 112.8628 33.93168  2.0273145
    ## 4 116.1811 36.72684  0.3005945
    ## 5 116.5833 33.89741 -1.9485332
    ## 6 110.3189 36.34210  7.6049308

    dg <- adw(dd, gridSize = 0.5, cdd = 1e5, m = 4) %>% na.omit()
    head(dg)

    ##         lon     lat      value
    ## 4  110.0044 32.5628  0.5288647
    ## 5  110.0044 33.0628  4.5636123
    ## 7  110.0044 34.0628 -4.6583732
    ## 8  110.0044 34.5628 -5.1216463
    ## 9  110.0044 35.0628 -6.1269176
    ## 10 110.0044 35.5628 -6.1508663

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
