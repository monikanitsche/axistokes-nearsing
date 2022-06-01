# axistokes-nearsing
Axisymmetric Stokes flow - near singular integrals
=================================================================
This repository contains code and supporting files to evaluate the single and double layer potential in 
axisymmetric Stokes flow without swirl, as used in the paper
"Corrected trapezoidal rule for near-singular integrals in axi-symmetric Stokes flow"   

The directory mathematica/ contains pdf files of the mathematica notebooks used to obtain the expansions of the integrands for each of the kernels.

The directory simpletest/ contains the code (Fortran 77) to evaluate the integrals for a range of target points on a grid, for a simple case of a spherical sheet with prescribed densities, using the corrections presented in the above mentioned paper. The code reproduces figure 6 and the figures in the supplementary materials to that paper.

--------------------------------------------------------------
To run the code: set the flags in file flags.h
THIRD for third order method
FOURTH for fourth order method
(other flags can be ignored, QUAD was used to test round off error issues, which are now resolved)

% make driver.out
% driver.out < datinfile

where datinfile is one of datinl100, datinl200, datinl400, datinl800, datinl1600
This creates dat.m files in the corresponding subdirectories.
These are read into matlab and save as dat.mat files
which are then plotted by the plotting routines in the given subdirectories.

--------------------------------------------------------------
Notes:
1) by searching all .F files for "ifdef THIRD" or "ifdef FOURTH" the reader will find the changes necessary to go form O(h^3) method to O(h^4) method

2) The variables in the code and the mathematica files, due to historical reasons, dont quite match those in the paper.
For example

Paper  Code

u     ->  v

v     ->  u

r     ->  r or x

z      -> z or y

Q111   -> Qv22

Q112   -> Qv21

Q222  ->  Qu11

M11   ->  M22

M12   ->  M21

f1    ->  f2

f2    ->  f1


3) the given code only resolves the flow in the left half of the sphere, as shown in the figure in paper
(the correction near the point that meets the axis at alf=pi is not included in this code)

4) the code uses the fact that we know the position of the sphere r,z, exactly, as well as the density functions f1,f2 (single layer) and u,v (double layer). for example, all projections and values at the projections are exact, to eliminate other sources of errors.

5) the code contains remnants of the full code used to compute the evolution of the double emulsion through the nozzle, which are somewhat redundant here (for example, compxyvol.F)

6) The file LAYOUT describes the dependency between all subroutines (but not quite up to date, needs to be updated)

--------------------------------------------------------------

Please contact the author for a copy of the paper, or any questions.
