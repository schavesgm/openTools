#!/bin/bash

# Name of the folder inside x32 -- where results are located
channelFit=$1
type=$2

case "$channelFit" in
    'g5' )
        chan="pseudoscalar${typeCalc^^}"
        ;;
    'vec' )
        chan="vectorConc${typeCalc^^}"
        ;;
    'ax' )
        chan="axialConc${typeCalc^^}"
        ;;
    'g0' )
        chan="scalar${typeCalc^^}"
        ;;
    'g6' )
        chan="axialZero${typeCalc^^}"
        ;;
    'teto' )
        chan="tempTensor${typeCalc^^}"
        ;;
    'tesp' )
        chan="spatTensor${typeCalc^^}"
        ;;
    'axvec' )
        chan="ratioAxVec${typeCalc^^}"
        ;;
    'g0g5' )
        chan="ratioG0G5${typeCalc^^}"
        ;;
    'vecg5' )
        chan="ratioVecG5${typeCalc^^}"
        ;;
esac

# Name of folder to look for
nameChannel=${chan}${type^^}

# Number of bootstrap iterations 
numBootCalc="20000"

# Number of parameters fitted in calculation
dimParam="2"

# Bin size used
binSize="36"

# Name to differ between runs
diffFolds=${numBootCalc}_${dimParam}_${binSize}

## NAME OF THE FOLDER 'PSEUDOSCALAR'
nameFolder=$( echo "${nameChannel}" | sed 's/Data//g' )
nameFolder=$( echo "${nameFolder}" | tr a-z A-Z )

lookFile=${numBootCalc}_${dimParam}_${binSize}

# All the possible hadron combinations
HADRONS=( uu us uc ss sc cc )

cd .. && FOLDS=$( ls -d *x32 ) && cd -

# Create a folder for the channel
mkdir -p ../${chan^^}

# Create a folder inside the channel with the type
mkdir -p ../${chan^^}/${type}/${diffFolds}

for NT in ${FOLDS}
do

    timePoints=$( echo ${NT} | sed "s/x.*//" )
           
    # Generate a folder for each NT
    mkdir -p ../${chan^^}/${type}/${diffFolds}/${timePoints}
    
    # Now for each hadron -- Generate its own folder and copy data
    for had in ${HADRONS[@]}
    do
        # Create a folder to hold the hadron
        mkdir -p ../${chan^^}/${type}/${diffFolds}/${timePoints}/${had}
        
        # Copy the data to that folder
        cp ../${NT}/${nameChannel}/*${had}*/*resultsFile_${numBootCalc}_${dimParam}_${binSize}* \
           ../${chan^^}/${type}/${diffFolds}/${timePoints}/${had}
        
        # Location now to return later
        locateNow=$( pwd )
        
        # Copy the cleaning file to folder
        cp ./cleanDataMass.sh ../${chan^^}/${type}/${diffFolds}/${timePoints}/${had}

        # # Move into that folder and execute file
        cd ../${chan^^}/${type}/${diffFolds}/${timePoints}/${had}
        bash ./cleanDataMass.sh ${chan} 
        rm ./cleanDataMass.sh && cd ${locateNow}
        
    done

done
