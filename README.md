# pyCenTESS: in-/out-of transit centroid analysis for candidate exoplanets observed by TESS

Requirements:
- numpy
- astropy
- astroquery
- matplotlib
- tess-point

How to "install":
cd pyCenTESS/pycentess/stackima/
f2py -c --f90flags=-O2 -m stackima stackima.f90
cd ../../

How to use:

python pyCenTESS.py RA DEC Period T0 T14 min_sector max_sector

where:
RA, DEC are Right Ascension and Declination of the target star in degree
Period is the orbital period of the candidate exoplanet in days
T0 is the central time of the first transit in BTJD
T14 is the duration of the transit in hours
min_sector/max_sector are the minimum and maximum sector over which the software will find TESS observations of the target

Outputs:
vetting_tic-TICID.pdf : finding chart of the candidate with the mean in-/out-of-transit centroids reported for each sector.



