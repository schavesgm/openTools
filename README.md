# Collection of scripts to analyse data from openQCD-Fastsum
---
This repository contains scripts used to produce and analyse data from the FASTSUM
collaboration. It is intended to be a collection of independent folders. In each
folder, there should be a fully functional collection of scripts.

Some folders are submodules, in order to fetch those repositories, you do need to
use,

```bash
git submodule init
git submodule update
```

---
## concatenate_hadspec -- Submodule
Scripts to concatenate Fastum hadspec correlation function data into a single file. 
Each file will contain ``N_C x N_S`` statistics of the particle to be analysed. The
scripts support mesons and baryons.

## run_scan -- Submodule
Collection of scripts to run jobs on the server using a better approach than the 
one inside ``cppFitCode/scanRunServer``. The code is way more mantainable and 
extendable. It is

## extractEffMass
Collection of scripts to extract the effective mass solving the cosh function and 
plotting the data into a pdf using ``gnp2tex``.

## gnp2tex -- Submodule
Program to transform gnuplot files to standalone pdf files with LaTeX style.

## cppFitCodes
This folder contains two subfolders:

1. **cppFitCode:** Main code to fit functions to data using Nelder-Mead
   algorithm. It is implemented in C++ and uses MPI as a parallel implementation.

2. **scanRunServer:** Set of scripts to automatise the fit of data from
   different channels and number of points in the time direction using a
   variational method in the window used to estimate the systematic error.

3. **PCACMass:** Code branched from **cppFitCode** to calculate the PCAC
   mass from correlation functions data using bootstrap.

## pythonFitCodes
This folder contains two subfolders:

1. **fitCode:** Python implementation of Nelder-Mead algorithm.

3. **getData:** Python implementation of the bootstrap estimation of
   correlation functions.
