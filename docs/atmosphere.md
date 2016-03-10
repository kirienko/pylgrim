What already works

tropospheric delays:
-------------
* Saastamoinen model: `saast()` - based on RTKlib code. If I were you I wouldn't use it.
* Vienna mapping function (aka VMF1_HT): `vmf()` - rewritten from Fortran source code
* parsing of source tropospheric files also implemented:
one can use `find_VMF_coeffs(ah_file, aw_file, coords)` for this purpose

ionospheric delays
------------
* Klobuchar model: `klobuchar()` - based on RTKlib
* ionofree pseudorange linear combination

IONEX maps parsing is to be implemented.