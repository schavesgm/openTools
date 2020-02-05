#!/bin/bash

#SBATCH --job-name=128x32_teto_ss
#SBATCH --output=128x32_teto_ss
#SBATCH --error=error_%J.out
#SBATCH --ntasks=40
#SBATCH --tasks-per-node=40
#SBATCH --array=1

## This definition Folder* will calculate all the mesons, change to Folder*cc*
## to calculate just one of them

## Look inside each folder
calculateMeson="Folder*"

## Name of the channel contained in the name -- used to get a proper name
channelFit=teto

for i in $( ls -d ${calculateMeson} )
do
    sed -i "/foldName=/c\foldName='${i}'" ./automaticFit.sh
    myHadron=$( echo "${i}" | sed "s/^.*${channelFit}/${channelFit}/" )
    sed -i "/typeOfHadron=/c\typeOfHadron='${myHadron}'" ./automaticFit.sh
    echo "Estoy en el hadron ${myHadron}"
    bash ./automaticFit.sh
    echo "-----"
done

# rm Folder*/Gen2l* 
