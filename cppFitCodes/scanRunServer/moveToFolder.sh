#!/bin/bash

type=ss
## This is the name of the channel -- inside concMeson.xx
channelFit=teto

if [ ${channelFit} == 'g5' ]; then           # Pseudoscalar
    lookFolder="g5_Folder"
elif [ ${channelFit} == 'vec' ]; then        # Conc vector
    lookFolder="vecConc_Folder"
elif [ ${channelFit} == 'ax' ]; then         # Conc axial
    lookFolder="axConc_Folder"
elif [ ${channelFit} == 'g0' ]; then         # Scalar
    lookFolder="g0_Folder"
elif [ ${channelFit} == 'g6' ]; then         # Temporal axial
    lookFolder="g6_Folder"
elif [ ${channelFit} == 'axvec' ]; then      # Axial / Vector
    lookFolder="axvecRatio_Folder"
elif [ ${channelFit} == 'teto' ]; then      # Axial / Vector
    lookFolder="tetoConc_Folder"
elif [ ${channelFit} == 'tesp' ]; then      # Axial / Vector
    lookFolder="tespConc_Folder"
elif [ ${channelFit} == 'g0g5' ]; then       # Scalar / Pseudoscalar
    lookFolder="g0g5Ratio_Folder"
elif [ ${channelFit} == 'vecg5' ]; then      # Vector / Pseudoscalar
    lookFolder="vecg5Ratio_Folder"
fi

## This is the name of the folder to move -- pseudoscalar
nameFolder=tempTensorSS

for dir in $( ls -d conc* )
do

    # Move into the hadron with definite type
    cd ${dir}/${type}

    # Get the name of the file to generate its Folder_*
    myHadron=$( cd ${lookFolder} && ls | sed "s/^/Folder_/" )

    # Move back twice
    cd ../..

    # Generate a folder to hold the data
    mkdir -p ./${nameFolder}/${myHadron}

    # Cooy the files inside the folder relative to the channel
    cp ./${dir}/${type}/${lookFolder}/* ./${nameFolder}/${myHadron}

done
