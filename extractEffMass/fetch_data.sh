#!/bin/bash

NT=( 128 64 56 48 40 36 32 28 24 20 16 )
chan="g0"

for nt in ${NT[@]}
do

    mkdir -p ${nt}x32.cc
    for type in ll ss
    do
        # Copy the file into the folder
        # locat="/home/s.989336/Fastsum-Lat/${nt}x32/concMeson.cc/${type}/${chan}_Folder/*"
        # mkdir -p ${nt}x32.cc/${type}
        # scp scw:${locat} ${nt}x32.cc/${type}

        # Then execute the script
        cp extract_cosh_mass.py ${nt}x32.cc/${type}/    # Copy the analysis file
        cd ${nt}x32.cc/${type}/                         # Move into directory
        python ./extract_cosh_mass.py $( ls ) effMass_${type}_$( ls )
        cd ../../

    done
done
