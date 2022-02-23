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

The **development** version can be installed from GitHub
(<https://github.com/PanfengZhang/adw>) using:

    # install.packages("remotes")
    remotes::install_github("PanfengZhang/adw")

## Usage

    library(adw)
    library(ggplot2)
    library(sf)

    ## Linking to GEOS 3.9.1, GDAL 3.2.1, PROJ 7.2.1; sf_use_s2() is TRUE

    ds <- read.csv("C:/documents/test/adw_test.csv")
    ds$value <- runif(nrow(ds), min = -10, max = 10)
    head(ds)

    ##       lon    lat      value
    ## 1 113.061 32.928 -2.4436383
    ## 2 114.310 33.653 -0.8338343
    ## 3 111.267 33.506  8.4525327
    ## 4 111.196 34.112  8.4404065
    ## 5 114.337 33.435  1.5890665
    ## 6 115.850 34.236 -5.6577436

    dg <- adw(ds, gridSize = 0.5, cdd = 100000, m = 4)
    # dg is the grid (mesh) dataframe
    head(dg)

    ##       lon    lat      value
    ## 1 110.397 31.453         NA
    ## 2 110.397 31.953         NA
    ## 3 110.397 32.453         NA
    ## 4 110.397 32.953  1.0175670
    ## 5 110.397 33.453 -0.6547725
    ## 6 110.397 33.953  2.7219923

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
