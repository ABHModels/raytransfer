# RayTransfer

A general relativistic ray-tracing code for calculating the transfer functions in Non-Kerr metrics. This code can be used for generating the FITS file to be used with RELXILL_NK and NKBB models.    

Cite the code: [![DOI](https://zenodo.org/badge/596872782.svg)](https://zenodo.org/doi/10.5281/zenodo.12703491)    
    
When using RayTransfer please cite the following:    

* _Public Release of RELXILL_NK: A Relativistic Reflection Model for Testing Einstein's Gravity_; A.B. Abdikamalov, D. Ayzenberg, C. Bambi, T. Dauser, J.A. Garcia, & S. Nampalliwar; [Astrophys.J. 878 no.2, 91 (2019)](https://doi.org/10.3847/1538-4357/ab1f89).    
* _Testing the Kerr black hole hypothesis using X-ray reflection spectroscopy and a thin disk model with finite thickness_; A.B. Abdikamalov, D. Ayzenberg, C. Bambi, T. Dauser, J.A. Garcia, S. Nampalliwar, A. Tripathi, & M. Zhou; [Astrophys.J. 889 no.1, 80 (2020)](https://doi.org/10.3847/1538-4357/aba625).

## Usage

RUNNING WITH BUILT-IN JOHANNSEN METRIC    

**Only one non-zero deformation parameter possible**

1. Generate isco.dat file using isco.cpp. First must modify isco.cpp to use the deformation parameter and bounds you want.    

2. If not using alpha13 deformation parameter, modify main.cpp to use the deformation parameter you want.    

3. Generate jobs files using job.cpp. First must modify job.cpp to use the isco.dat file, deformation parameter, and accretion rate you want.    

4. Compile main.cpp using:    
        g++ main.cpp -O3 -o photons4trf    

5. Make sure a folder named "photons" exists in the same directory.    

6. Modify Transfer.py to use the isco.dat file, deformation parameter, and accretion rate you want.    

7. Newly created FITS file will be called Trf.fits    
 
8. Support contact: <relxill_nk@fudan.edu.cn>

## Changing Metric

**Only one non-zero deformation parameter possible.**    

1. Using the Maple worksheet "CodeOptimization.mw" or the Mathematica notebook "CodeOptimization.nb", enter covariant metric components of new metric and run to generate optimized code segments. Replace optimized code segments in metric.cpp, diffeqs.cpp, and redshift.cpp as labeled in worksheet/notebook. Must also modify deformation parameter declarations and assignments in main.cpp and def.h. *Note:* the Maple worksheet generates a much more optimized code so it is preferred over the Mathematica notebook.    

2. Follow RUNNING WITH BUILT-IN JOHANNSEN METRIC above.    

## File description

### Core ray-tracing code
* main.cpp - Contains parameter assignments, primary loops of the code, search loops for redshift values, and outputs to file.
* def.h - Contains include statements, function declarations, definitions of constants, and global variable declarations.
* metric.cpp - Contains functions to get components of metric and metric derivatives.
* findisco.cpp - Finds ISCO radius.
* gauleg.cpp - Sets disk radial values sought after and eventually used in the FITS file.
* rayprecise.cpp - Contains function with search loop for radial values. Also contains function to convert from polar to cartesian coordinates on the screen.
* raytracing.cpp - Ray-tracing code.
* diffeqs.cpp - Differential equations solved for by ray-tracing code.
* redshift.cpp - Calculates redshift value at the end of the ray-trace.

### Additional files
* CodeOptimization.mw/CodeOptimization.nb - Maple/Mathematica scripts used to optimized code related to the metric.
* OptimizeExpressionToC.m - Mathematica script used by CodeOptimization.nb
* isco.cpp - Generates isco.dat file used in other files.
* Transfer.py - Python script to create FITS file after ray-tracing code has finished.
