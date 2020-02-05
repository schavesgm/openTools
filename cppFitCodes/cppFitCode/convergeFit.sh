#!/bin/bash 

## ------------------------------------------------------------------------------------------------- 
##  AUTHOR : SERGIO CHAVES GARCIA-MASCARAQUE\ 
##  E-MAIL : SERGIOZTESKATE@GMAIL.COM 
##  ENERO DE 2019, SWANSEA, WALES 
## ---------------------------------------------------------------------------------

## BEFORE RUNNING THIS CODE, CHECK THAT YOU HAVE GENERATED AN 
## INPUT FILE WITH THE NAME THAT YOU WILL CLARIFY IN VARIABLE 
## 'inputFile' AND YOU HAVE CHANGED THE FUNCTION THAT YOU WOULD ## LIKE TO FIT INSIDE 'fitMain.cpp' 

## Variables we need
inputFile="input.in"            ## File containing the input definitions
Nproc=5
IFS=' '                         ## Delimiter in the string
improvBreak=0.00001
numIter=0                       ## Number of iterations until convergence
maxIter=5
resultFile='g5.cc_Gen2l_8x32.meson.g5.cc_resultsFile_10-90_15000.dat'

# Time window used in this scan
timeWindow="20-80"

## Name of the compiled files
nameComp=$( echo $resultFile | sed 's/.meson.*//' )
nameMerg=$( echo $resultFile | sed 's/.*resultsFile//' )
nameComp=${nameComp}_${nameMerg}.out

numBootCalc=20000
numBootBest=20000


## Load the modules
module purge
module load compiler/gnu/8/1.0
module load mpi/openmpi


dimParam=$( grep 'dimParam' $inputFile )    ## Dimension of the parameter space
read -ra dimParam <<< "$dimParam"
dimParam=${dimParam[1]}

## Generate the best estimate
sed -i "/numBootstrap/c\numBootstrap ${numBootBest}" ./input.in
mpic++ -o $nameComp getData.cpp
mpirun $nameComp $inputFile

# Generate the fit
mpic++ -o $nameComp fitMain.cpp 
    
declare -A oldValue
declare -A oldExplo
declare -A newValue
declare -A newExplo
 
while :
do
    
    echo "I am currently at iteration number" ${numIter}
	mpirun $nameComp $inputFile > myResult.dat
	A=$( cat myResult.dat ) && rm myResult.dat
    read -ra A <<< "$A"
    # echo "Old iter" ${A[@]}

    ## Replace the result as the initial guess in the next iteration
    countChecker=0
    for (( i = 0; i < $dimParam; i++  ))
    do  
        oldValue[$i]=${A[$countChecker]}
        oldExplo[$i]=${A[$(($countChecker+1))]}
        countChecker=$(( $countChecker + 2 ))
    done
    
    # Replace the initial guess and explore volume
    sed -i "/initGuess/c\initGuess" $inputFile
    sed -i "/explGuess/c\explGuess" $inputFile
    
    for (( i = 0; i < $dimParam; i++ ))
    do
        sed -i "/initGuess/ s/$/ ${oldValue[$i]}/" $inputFile
        sed -i "/explGuess/ s/$/ ${oldExplo[$i]}/" $inputFile
    done
    
    ## Generate a new random seed
    sed -i "/initSeed/c\initSeed $(( RANDOM % 10000000 ))" $inputFile
    
    ## Generate the next fit
	mpirun $nameComp $inputFile > myResult.dat
	A=$( cat myResult.dat ) && rm myResult.dat
    read -ra A <<< "$A"
    # echo "New iter" ${A[$((${dimParam}*2))]}
    
    ## Replace the result as the initial guess in the next iteration
    countChecker=0
    for (( i = 0; i < $dimParam; i++ ))
    do  
        newValue[$i]=${A[$countChecker]}
        newExplo[$i]=${A[$(($countChecker+1))]}
        countChecker=$(( $countChecker + 2 ))
    done
    
    # Replace the initial guess and explore volume
    sed -i "/initGuess/c\initGuess" $inputFile
    sed -i "/explGuess/c\explGuess" $inputFile
    
    for (( i = 0; i < $dimParam; i++ ))
    do
        sed -i "/initGuess/ s/$/ ${newValue[$i]}/" $inputFile
        sed -i "/explGuess/ s/$/ ${newExplo[$i]}/" $inputFile
    done
    
    ## Check the iterative condition
    declare -A valAbs
    declare -A condCheck
    
    for (( i = 0; i < $dimParam; i++ ))
    do
        valAbs[$i]=$( echo "${newValue[$i]} - ${oldValue[$i]}" | bc -l )
        valAbs[$i]=${valAbs[$i]#-}
    if (( $( bc <<< "${valAbs[$i]} <= $improvBreak" ) )); then
        condCheck[$i]=1
    else
        condCheck[$i]=0
    fi
    done
    
    checkConvergence=0
    for (( i = 0; i < $dimParam; i++ ))
    do
        checkConvergence=$(( $checkConvergence + ${condCheck[$i]} ))
    done
    
    numIter=$(( $numIter + 1 ))
    
    if (( "$checkConvergence" == "$dimParam" )); then
        echo "Convergence with eps =" $improvBreak "obtained in" \
             $numIter "iterations"
        echo "Time Window is ${timeWindow}" >> ${resultFile}
        echo "ChiSquare :" ${A[$(( 2*$dimParam ))]} >> ${resultFile}
        for (( i = 0; i < $dimParam; i++ ))
        do  
            echo "Param :" $(( $i + 1 )) "is" ${newValue[$i]} \
                 "+-" ${newExplo[$i]} >> ${resultFile}
        done
	    break
    fi
    
    if (( "$numIter" == "$maxIter" )); then
        echo "Maximum number of iterations achieved" $maxIter \
             ", I broke the program"
        echo "Time Window is ${timeWindow}" >> ${resultFile}
        echo "ChiSquare :" ${A[$(( 2*$dimParam ))]} >> ${resultFile}
        for (( i = 0; i < $dimParam; i++ ))
        do  
            echo "Param :" $(( $i + 1 )) "is" ${newValue[$i]} \
            "+-" ${newExplo[$i]} >> ${resultFile}
        done
        break
    fi
    
    ## Generate a new random seed
    sed -i "/initSeed/c\initSeed $(( RANDOM % 10000000 ))" $inputFile
	echo "Acabo de terminar una iteracion" 
    
done    # End of while true

rm $nameComp 
sed -i '/initGuess/c\initGuess 1 1' $inputFile
sed -i '/explGuess/c\explGuess 1 1' $inputFile

# Remove uneeded files when finished
rm -r ./hFiles
rm input.in
rm getData.cpp
rm fitMain.cpp


    
