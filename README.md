GNSS software receiver written in python
===============================================
Playground project. Just to understand how GNSS works.

Project Contents:
-----------------
* proto/     -- prototype of software receiver (mostly `python` + `numpy`)
    + coord/ -- tools for coordinate transition
    + helper/ -- misc helper functions
    + tests/ -- unit test
    + visualization/ -- visualization tools
* test_data/ -- supplementary data files for tests

What works:
----------
For now one can obtain rover's position by using single- or dual-frequency
 observations, with accuracy of 50-100 meters. For static observations with GPS.
 Something might work in kinematic mode and for GLONASS, partly implemented but not tested yet.
 (But you must understand how to run all this machinery.
 Sorry, guys, I have no time for UI at the moment.)

+ There are parsers for the following:
    * RINEX (both with navigation data and observations)
    * SP3 -- precise orbits and clocks
    * IONEX ionospheric maps
    * Vienna mapping function sources: hydrostatic and wet coefficients (see [here](http://unb-vmf1.gge.unb.ca/Products.html))
+ It is possible to visualize satellites from RINEX Nav file (3D model with rotation).
+ It's possible to visualize points on the map (depends on `mpl_toolkits.basemap`).
+ Vienna mapping function (aka `VMF1_HT`) implemented.
+ Ionospheric delays from IONEX maps implemented.
+ Klobuchar (broadcast) ionospheric model implemented.
+ [MILES algorithm](http://www.cs.mcgill.ca/~chang/software/MILES_Theory_Alg.pdf) reimplemented in python

Dependencies:
------------
* python 2.7
* numpy
* [optional] matplotlib.pyplot -- sats and rover visualization
* [optional] mpl_toolkits -- satellite visualization
* [optional] mpl_toolkits.basemap-data -- rover on the map

Maybe some day:
--------------
* GPS and GLONASS with 10 meter accuracy
* moving rover
* phase ambiguity
