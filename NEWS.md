# adw 0.3.1
-   The help files for funchtions 'adw_vector', 'adw_sf' and 'adw_terra' are changed, and some unnecessary content has been removed. The unit of parameter 'extent' was converted from meter to kilometer since the version 0.3.1.

# adw 0.3.0
-   The adw interpolation function is rewritten, and the parameter 'extent' can be a class of 'sf', 'SpatVector', or 'vector'. The calculation speed will be several times faster than before.
-   deleter parameter of 'maskON'.

# adw 0.2.1

-   fix the BugReports website.
-   add README file

# adw 0.2.0

-   all of the functions of sf package were replaced by the functions of terra package. The calculation speed will be several times faster than before.
-   add parameters of extent, nmin, nmax, and maskON
-   parameter gridSize were changed to gridsize
-   delete parameters of xmin, xmax, ymin, ymax
-   delete adw_land function. The new version of the function 'adw' adds a parameter 'maskON' to implement the function of the mask.