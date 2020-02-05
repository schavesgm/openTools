#!/bin/bash -- login

# Script para lanzar un conjunto de jobs

# Choose from g0, g5, vec, ax, axvec, g0g5, vecg5, teto, tesp
channelFit=$1

# Choose from ll or ss
typeCalc=$2

# LOGIC TO MAKE THINGS AUTOMATIC -----------------------------------------
case "$channelFit" in
    'g5' )
        nameFolderSave="pseudoscalar${typeCalc^^}"
        ansatz='cosh' 
        ;;
    'vec' )
        nameFolderSave="vectorConc${typeCalc^^}"
        ansatz='cosh' 
        ;;
    'ax' )
        nameFolderSave="axialConc${typeCalc^^}"
        ansatz='cosh' 
        ;;
    'g0' )
        nameFolderSave="scalar${typeCalc^^}"
        ansatz='cosh'
        ;;
    'g6' )
        nameFolderSave="axialZero${typeCalc^^}"
        ansatz='cosh'
        ;;
    'teto' )
        nameFolderSave="tempTensor${typeCalc^^}"
        ansatz='cosh'
        ;;
    'tesp' )
        nameFolderSave="spatTensor${typeCalc^^}"
        ansatz='cosh'
        ;;
    'axvec' )
        nameFolderSave="ratioAxVec${typeCalc^^}"
        ansatz='exp'
        ;;
    'g0g5' )
        nameFolderSave="ratioG0G5${typeCalc^^}"
        ansatz='exp'
        ;;
    'vecg5' )
        nameFolderSave="ratioVecG5${typeCalc^^}"
        ansatz='exp'
        ;;
esac
        
## I would say 40 tasks is the best option
numProc=40
binSize=36
numBootBest=20000
numBootCalc=20000

# BE CAREFUL IN HERE  -- LL - resBool = 0, SS - resBool = 1ss
if [ ${typeCalc} == 'll' ]; then
    resBool=0       # local-local
else
    resBool=1       # smeared-smeared
    # Smearing in ratios does not need rescaling
    if [ ${channelFit} == 'axvec' ]; then
        resBool=0
    elif [ ${channelFit} == 'g0g5' ]; then
        resBool=0
    elif [ ${channelFit} == 'vecg5' ]; then
        resBool=0
    fi
fi
resFactor="1e9"

# Change inside the python script the window size
# getWindow="35-65"
dimParam=2

if [ ${channelFit} == 'g5' ]; then
    initGuess=( 1 0.5 )                     # Pseudoscalar
elif [ ${channelFit} == 'g0' ]; then
    initGuess=( 1 0.5 )
elif [ ${channelFit} == 'axvec' ]; then
    initGuess=( 1 0.2 )
elif [ ${channelFit} == 'g0g5' ]; then
    initGuess=( 1 0.2 )
elif [ ${channelFit} == 'vecg5' ]; then
    initGuess=( 1 0.2 )
else
    initGuess=( -0.2 0.5 )    # Axial - Vector
fi


maxChiSq=100000

thisFolder=$( pwd )

## Change things inside the other scripts
# Inside moveToFolder.sh
sed -i "/nameFolder=/c\nameFolder=${nameFolderSave}" moveToFolder.sh
sed -i "/channelFit=/c\channelFit=${channelFit}" moveToFolder.sh
sed -i "/type=/c\type=${typeCalc}" moveToFolder.sh

# Inside automaticFit.sh
sed -i "/numProc=/c\numProc=${numProc}" automaticFit.sh
sed -i "/binSize=/c\binSize=${binSize}" automaticFit.sh
sed -i "/nBootBest=/c\nBootBest=${numBootBest}" automaticFit.sh
sed -i "/nBootCalc=/c\nBootCalc=${numBootCalc}" automaticFit.sh 
sed -i "/resBool=/c\resBool=${resBool}" automaticFit.sh 
sed -i "/resFactor=/c\resFactor=${resFactor}" automaticFit.sh 
sed -i "/maxChiSq=/c\maxChiSq=${maxChiSq}" automaticFit.sh
sed -i "/initGuess=/c\initGuess=( ${initGuess[0]} ${initGuess[1]} ${initGuess[2]} )" automaticFit.sh 
sed -i "/dimParam=/c\dimParam=${dimParam}" automaticFit.sh
sed -i "/ansatz=/c\ansatz=${ansatz}" automaticFit.sh

# Inside launchFit.sh
sed -i "/channelFit=/c\channelFit=${channelFit}" launchFit.sh
sed -i "/#SBATCH --ntasks=/c\#SBATCH --ntasks=${numProc}" launchFit.sh
sed -i "/#SBATCH --tasks-per-node/c\#SBATCH --tasks-per-node=${numProc}" launchFit.sh

cd .. ## Move backwards so we can get the data directories
for dir in $( ls -d *x32 | sort -n )
do  

    echo $dir
    ## Move back to runServer
    cd $thisFolder

    ## Change job name inside launcFit.sh
    sed -i "/#SBATCH --job-name/c\#SBATCH --job-name=${dir}_${channelFit}_${typeCalc}" \
    launchFit.sh
    sed -i "/#SBATCH --output/c\#SBATCH --output=${dir}_${channelFit}_${typeCalc}" \
    launchFit.sh
    
    
    ## Change the time extent in each file
    timePoints=$( echo ${dir} | sed "s/x.*//" )
    sed -i "/timePoints=/c\timePoints=${timePoints}" automaticFit.sh 

    # ## Copy and run moveToFolder to the directory
    cp moveToFolder.sh ../${dir}/ && cd ../${dir}
    mkdir -p ${nameFolderSave}
    bash moveToFolder.sh && rm moveToFolder.sh
    cp ../scanRunServer/automaticFit.sh ./${nameFolderSave}
    cp ../scanRunServer/cppFitCode.tar.gz ./${nameFolderSave}
    cp ../scanRunServer/launchFit.sh ./${nameFolderSave}

    # Run the slurm job using launchFit.sh
    cd ./${nameFolderSave} && sbatch launchFit.sh

done

