# shapelets
Python and fortran77 routines for modelling astronomical sources and Point Spread Functions as shapelets

In KiDS, PSFs are modelled as 'psf maps', which are shapelets with coefficients that vary as polynomials across the image.
These codes allow these maps to be read and turned into postage stamps of the PSF at any (fractional) position in the image.


# python3 (shapeletpsf.py)

Import shapeletpsf and read the help; it explains the two classes shapeletseries and shapeletmap, how to read a PSF map file, and how to use it.
A simple way to use shapeletpsf is provided in the example.py file, which reads in a map (an example is provided) and shows the PSF at a given pixel position. E.g.,:

**python3 example.py KIDS_0.0_-31.2_i2.psf.map 6000.8 5000.1** *gives output*

Loaded (nsh=10, n=4) shapelet map from KIDS_0.0_-31.2_i2.psf.map

Integral of PSF, scale radius: 0.997741330420072 2.17293406

Shapelet coefs: [order 00 10 01 20 11 02 30 21 12 03 ...]

[ 1.28344960e-01  2.21663057e-12  5.16109804e-12 -1.11720752e-02
  5.48167864e-04 -1.94338421e-02  1.77388759e-03  6.64638510e-04
  1.21116663e-03  1.32725960e-04  9.73213814e-03  2.54315604e-04
  8.80680095e-03 -1.44193630e-04  1.29987809e-02 -4.84052869e-04
 -1.62782875e-05 -5.51669858e-04 -2.36994195e-04 -4.14755639e-04
 -1.24234739e-04 -4.47892729e-05  4.86989024e-05 -2.31523612e-04
 -1.01267575e-04 -5.93706773e-04 -8.80008448e-05 -2.24179700e-03
  6.72891792e-04  1.69272111e-04  4.96527443e-04  6.88674259e-05
  4.65378178e-04  2.04453936e-04  4.23654098e-04  1.85584768e-04
  2.24442809e-03 -4.39750109e-06  1.94487892e-03 -1.03624986e-04
  2.06652075e-03 -3.30578626e-05  2.20895029e-03 -6.95906245e-05
  3.10555825e-03 -1.28869187e-04 -2.48128813e-05 -1.89544861e-04
 -7.03491347e-05 -1.52831671e-04 -4.57944046e-06 -1.80286357e-04
 -6.81788222e-05 -2.19711691e-04 -3.43651096e-05  7.21402447e-04
 -4.25156175e-06  4.74756475e-04  1.73903474e-05  4.35291878e-04
  5.14887535e-05  4.08161682e-04  5.93954295e-05  2.87058123e-04
 -2.96846608e-05  1.40137115e-05]

![example_600 8_5000 1](https://user-images.githubusercontent.com/6078683/182646093-832cfb69-35ff-4252-a561-5f139fd219cb.png)


**python3 example.py KIDS_0.0_-31.2_i2.psf.map 17000 15000.5**  *gives output*

Loaded (nsh=10, n=4) shapelet map from KIDS_0.0_-31.2_i2.psf.map

Integral of PSF, scale radius: 0.9994863401595164 2.17293406

Shapelet coefs: [order 00 10 01 20 11 02 30 21 12 03 ...]

[ 1.30287211e-01  1.69371847e-12 -8.56441620e-12 -1.44900016e-02
 -3.43722894e-03 -1.97632552e-02 -9.54054877e-04  1.01236595e-03
 -3.72260524e-04  2.35350267e-04  1.14176052e-02  4.88776368e-04
  1.03302134e-02  1.03566574e-03  1.32986561e-02 -3.07793259e-04
  1.25505602e-04  1.39728927e-06 -2.73773921e-04  2.99577898e-04
 -3.32267252e-05 -8.17240347e-04 -6.47750737e-04 -8.91520220e-04
 -4.99929565e-04 -1.43278619e-03 -5.72661828e-04 -2.54672702e-03
  6.81167867e-05  1.28098521e-04  7.18087971e-05  7.46116188e-05
  9.51700111e-05  1.52357434e-04 -1.94788147e-05  8.96528826e-05
  2.46068553e-03  2.77258221e-04  2.25420563e-03  2.79363814e-04
  2.40855819e-03  2.99293317e-04  2.63126690e-03  1.15828132e-04
  3.23146735e-03  8.71841247e-05  4.15101545e-05 -2.11110363e-05
 -6.74120699e-06  6.33921923e-05  1.60631780e-04  7.04070570e-05
 -7.80463844e-06  9.37305541e-05 -2.31307441e-05  4.12570369e-04
 -9.05276902e-05  2.94388903e-04 -1.27179892e-04  2.70780758e-04
 -9.31432771e-05  2.35115097e-04 -1.33334693e-05  4.76841547e-06
 -8.05574154e-05 -3.65866075e-05]
                            
![example_17000_15000 5](https://user-images.githubusercontent.com/6078683/182646126-fd4af06b-0d0b-46aa-84b6-bf08c6cf124b.png)


# fortran77 shapeletpsf.F

The original subroutines that do the same as the python module above.

A wrapper programme is also provided: compile with 

gfortran -Dmm=16 -Dmsh=155 -Dmmpoly=6 -Dmmpoly2=28
-ffixed-line-length-0 psfmap2pixpsf.F shapeletpsf.F -o psfmap2pixpsf

and run as

(echo KIDS_0.0_-31.2_i2.psf.map; echo 12000 8000) | psfmap2pixpsf

(this simple wrapper simple writes out all elements of the postage
stamp array stamp(i,j) )
