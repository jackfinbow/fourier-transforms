# Discrete Fourier Transforms
Code written in C that performs discrete Fourier transforms (DFTs) and inverse discrete Fourier transforms (IDFTs) of complex functions.

## Introduction
This was a second-year project as part of my MPhys Physics degree. It involved writing code in C to be able to handle complex functions, and perform DFTs and IDFTs to these functions. Two of these functions were defined in the code and sampled at discrete time intervals and a third was obtained from a provided text file via raw data values. Python was used for plotting graphs of the data and comparing the results.

## Mathematical Background
The functions `h1(t)` and `h2(t)` are given by

<img src="/Figures/functions.png" alt="functions" width="250">

They were first sampled at discrete time values and a DFT applied to both to produde `H1(omega)` and `H2(omega)`. The function `h1(t)` features a sinudoidal wave envelope containing another sinusoidal shape; applying the IDFT but skipping the low frequency term removes the sinusoidal pattern of the wave envelope. When applying the IDFT to `h2(t)`, the largest amplitude term was skipped, translated the function in the negative y direction.

<img src="/Figures/h1.png" alt="h1" width="410"> <img src="/Figures/h2.png" alt="h2" width="440">

The data obtained for `h3(t)` displays a high level of noise in both the real and imaginary parts. The IDFT was applied only to the four largest amplitude terms which occurred at low frequency, hencing skipping the high frequency terms. This had the effect of reducing noise in the raw data to produce a smoother curve.

<img src="/Figures/h3_re.png" alt="h3 real" width="410"> <img src="/Figures/h3_im.png" alt="h3 imag" width="410">

## File Structure
The C code is contained within FourierTransforms.c and the code used to plot the figures is contained within Graphs.py. The raw data files and data produced from the IDFTs can be found in the DataFiles folder. The figures plotted using Python are located in the Figures folder.
