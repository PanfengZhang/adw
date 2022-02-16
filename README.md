# Angular distance weighting (adw)

Gridding the irregularly-spaced data onto a regular latitude-longitude
grid by weighting each station according to its distance and angle from
the center of a search radius

# reference

Caesar, J., L. Alexander, and R. Vose, 2006: Large-scale changes in
observed daily maximum and minimum temperatures: Creation and analysis
of a new gridded data set. J. Geophys. Res., 111, D05101,
<https://doi.org/10.1029/2005JD006280>.

## Installation

The **development** version can be installed from GitHub using:

    # install.packages("remotes")
    remotes::install_github("PanfengZhang/adw")

## Usage

    library(adw)
    library(ggplot2)
    library(sf)

    ## Linking to GEOS 3.9.1, GDAL 3.2.1, PROJ 7.2.1

    ds <- read.csv("C:/documents/test.csv")
    ds$value <- runif(nrow(ds), min = -10, max = 10)
    head(ds)

    ##       lon    lat       value
    ## 1 113.061 32.928 -3.08832563
    ## 2 114.310 33.653  2.99203032
    ## 3 111.267 33.506  8.34545434
    ## 4 111.196 34.112  2.21211022
    ## 5 114.337 33.435 -0.32206607
    ## 6 115.850 34.236 -0.05104718

    dg <- adw(ds, xmin = 110, xmax = 117, ymin = 31, ymax = 37, 
        gridSize = 0.5, cdd = 100000, m = 4)
    # dg is the grid (mesh) dataframe
    head(dg)

    ##   lon  lat     value
    ## 1 110 31.0        NA
    ## 2 110 31.5        NA
    ## 3 110 32.0        NA
    ## 4 110 32.5        NA
    ## 5 110 33.0        NA
    ## 6 110 33.5 -1.014796

    # plot
    urlmap <- "https://geo.datav.aliyun.com/areas_v3/bound/410000_full.json"
    cmap <- read_sf(urlmap) %>% st_cast('MULTILINESTRING')
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
