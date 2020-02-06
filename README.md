# Collection of scripts to analyse data from openQCD-Fastsum

This repository contains scripts used to produce and analyse data from the FASTSUM
collaboration. It is intended to be a collection of independent folders. In each
folder, there should be a fully functional collection of scripts.

1. **concatenateOpenQCD:** Scripts to concatenate openQCD correlation function
    data into a single file. Each file will contain N_C x N_S statistics of the
    particle to be analysed.

2. **extractEffMass:** Collection of scripts to extract the effective mass solving
    the cosh function and plotting the data into a pdf using gnp2tex.

1. **gnp2tex:** Script to transform gnuplot files to standalone pdf files with
    LaTeX style.

2. **cppFitCodes:** This folder contains two subfolders:

    1. **cppFitCode:** Main code to fit functions to data using Nelder-Mead
       algorithm. It is implemented in C++ and uses MPI as a parallel implementation.

    2. **scanRunServer:** Set of scripts to automatise the fit of data from
       different channels and number of points in the time direction using a
       variational method in the window used to estimate the systematic error.

    3. **PCACMass:** Code branched from **cppFitCode** to calculate the PCAC
       mass from correlation functions data using bootstrap.

3. **pythonFitCodes:** This folder contains two subfolders:

    1. **fitCode:** Python implementation of Nelder-Mead algorithm.

    3. **getData:** Python implementation of the bootstrap estimation of
       correlation functions.
