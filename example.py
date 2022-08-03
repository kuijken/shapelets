import shapeletpsf as sp
from sys import argv
import os
import numpy as np
import matplotlib.pyplot as plt

# Read in a PSF map and construct the PSF at a given pixel position.

# arguments are the psf map file and the x,y position for the PSF (in pixels).
# If x,y are integers the pixel is centered on the middle plxel of the cutout
# if x,y are not integers than the PSF is shifted by a fraction of a pixel, with
# the centre of the PSF at (x-int(x), y-int(y))
# by default the cutouts run from -12 to +12 pixels (25 x 25)

mapfile=argv[1]
xmap=float(argv[2])
ymap=float(argv[3])

map=sp.shapeletmap(mapfile=mapfile)   # Read in the PSF map
psfsh=map.atxy(xmap,ymap)             # interpolate the shapelet coefficients to the XY position
z=psfsh.toimage(display=True,xc=xmap%1,yc=ymap%1)         # make a postage stamp of the PSF and display it

print('Integral of PSF, scale radius:', z.sum(),psfsh.beta)
print('Shapelet coefs: [order 00 10 01 20 11 02 30 21 12 03 ...]')
print(psfsh.coefs)

