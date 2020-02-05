# Set of scripts to fit correlation functions data in the server.

This README contains all the information about the chain of scripts
used in the fitting process of correlation functions data in the 
server.

* The main script is scwLAUNCHER.sh. It iterates over all the folders
contained in the previous directory with names containing __x32__. The 
script copies the needed files inside each folder and copies the data
inside a given folder whose name must be declared. It runs the script
__moveToFolder.sh__, which copies the data contained in the concMeson 
folders inside the __nameFolder__ directory. Moreover, it copies all the
needed files to fit data inside __nameFolder__. Then it calls the script
__automaticFit.sh__, which performs the fitting process and calls.

* The script __moveToFolder.sh__ moves the data contained inside 
concMesons with the given channel to the __nameFolder__ directory 
to prepare it for the fitting process.

* The script __launchFit.sh__ iterates over all the folders containing
the hadronic correlation functions and launches __automaticFit.sh__.

* The script __automaticFit.sh__ changes the input inside the hadron
folder to adjust the fit. It copies the files inside cppFitCode.tar.gz,
which contains __convergeFit.sh__. This script generates the estimation
of the correlation function using bootstrap and calculates the fitted 
parameters convergently. It calls SLURM mpirun to generate the data.

* In this newer version, the script scans all the possible time
windows from \tau = a_\tau to N_\tau - 2 * a_\tau.

* The script __cleanMass.dat__ produces a clean version of the data when
we call the organize script. It allows us to generate a file with all the
masses extracted at one.
