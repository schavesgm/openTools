#!/bin/bash 

## Folder to fit -- Plan to do it automatically
foldName='Folder_Gen2l_128x32.meson.g5.cc'
## The type of hadron could be change to any other
typeOfHadron='g5.cc'

timePoints=128
colSize=3
maxChiSq=100000

## convergeFit.sh parameters
numProc=40
improvBreak=0.0001
maxIter=6

binSize=36
nBootBest=20000
nBootCalc=20000

# BE CAREFUL HERE IN THE CASE OF --- LOCAL-LOCAL
resBool=1
resFactor=1e9

dimParam=2
initGuess=( -0.2 0.5  )

tInit=1
tFini=$( echo "${timePoints}/2-1" | bc )

ansatz=cosh

## ---------------------------------------------------------------------
# CHANGE INPUT FILES ACCORDINGLY
resultFile="${typeOfHadron}_${foldName}_resultsFile_${nBootCalc}_${dimParam}_${binSize}.dat"
resultFile=$( echo "$resultFile" | sed 's/_Folder//g' )
resultFile=$( echo "$resultFile" | sed 's/ /-/g' )

echo "# --" > ./$foldName/$resultFile

## Directory to input file
dirInput="./${foldName}/input.in"

## Copy the fitting code inside the folder
tar -xvf cppFitCode.tar.gz

# Change the time extent of the function to fit
sed -i "/\\#define TIME_EXTENT/c\\\#define TIME_EXTENT ${timePoints}" \
       ./cppfitCode.tarball/fitMain.cpp

if [ ${ansatz} == 'cosh' ]; then
    sed -i "/\\#define ANSATZ/c\\\#define ANSATZ 0" \
           ./cppfitCode.tarball/fitMain.cpp
elif [ ${ansatz} == 'exp' ]; then
    sed -i "/\\#define ANSATZ/c\\\#define ANSATZ 1" \
           ./cppfitCode.tarball/fitMain.cpp
fi

## LOOP TO GENERATE DIFFERENT TIME WINDOWS
# Run convergeFit.sh as a function of the initial time t1
for i in $( seq  ${tInit} 1 ${tFini} )
do
    
    inTime=$i
    
    # The final time differs from cosh to exp   
    if [ ${ansatz} == 'cosh' ]; then
        fnTime=$( echo "${timePoints} - $i" | bc ) 
    elif [ ${ansatz} == 'exp' ]; then
        fnTime=$( echo "${timePoints}/2" | bc ) 
    fi

    echo "I am using TW = $inTime - $fnTime" 

    # Copy inside the folder
    cp -r ./cppfitCode.tarball/* ${foldName}
    
    ## WE NEED TO CHANGE THE INPUT FILE
    fileName=$( ls ./${foldName}/*.${typeOfHadron} | xargs -n 1 basename )
    numLines=$( wc -l < ./${foldName}/*.${typeOfHadron} )
    
    ## Get the file parameters inside input,in
    sed -i "/paramsFile/c\paramsFile ${numLines} ${colSize} ${timePoints}" ${dirInput}
    sed -i "/initSeed/c\initSeed $(( RANDOM % 10000000 ))" ${dirInput}
    sed -i "/binSize/c\binSize ${binSize}" ${dirInput}
    sed -i "/fileName/c\fileName ${fileName}" ${dirInput}
    sed -i "/windowFit/c\windowFit ${inTime} ${fnTime}" ${dirInput}
    sed -i "/saveName/c\saveName bestEstimate_${resultFile}" ${dirInput}
    sed -i "/largestChiSquare/c\largestChiSquare ${maxChiSq}" ${dirInput}
    sed -i "/dimParam/c\dimParam ${dimParam}" ${dirInput}
    
    sed -i "/initGuess/c\initGuess" ${dirInput}
    sed -i "/explGuess/c\explGuess" ${dirInput}
    
    sed -i "/rescaleOrNot/c\rescaleOrNot ${resBool}" ${dirInput}
    sed -i "/rescaledFactor/c\rescaledFactor ${resFactor}" ${dirInput}
    
    for (( i = 0; i < $dimParam; i++ ))
    do
        sed -i "/initGuess/ s/$/ ${initGuess[$i]}/" ${dirInput}
        sed -i "/explGuess/ s/$/ ${initGuess[$i]}/" ${dirInput}
    done
    
    ## Change convergeFit.sh file
    sed -i "/--ntasks=/c\\\#SBATCH --ntasks=${numProc}" ./${foldName}/convergeFit.sh
    sed -i "/resultFile=/c\resultFile='${resultFile}'" ./${foldName}/convergeFit.sh
    sed -i "/Nproc=/c\Nproc=${numProc}" ./${foldName}/convergeFit.sh
    sed -i "/improvBreak=/c\improvBreak=${improvBreak}" ./${foldName}/convergeFit.sh
    sed -i "/maxIter=/c\maxIter=${maxIter}" ./${foldName}/convergeFit.sh
    sed -i "/numBootBest=/c\numBootBest=${nBootBest}" ./${foldName}/convergeFit.sh
    sed -i "/numBootCalc=/c\numBootCalc=${nBootCalc}" ./${foldName}/convergeFit.sh
    sed -i "/timeWindow=/c\timeWindow=${inTime}-${fnTime}" ./${foldName}/convergeFit.sh
    
    cd ./${foldName} && bash convergeFit.sh && cd ..

    ## Remove the unnecesary data inside the folders to avoid memory congestion
    # rm ./${foldName}/slurm*

done

