#!/bin/bash 

# Properties of the data - As input arguments
checkChan=$1
if [ $checkChan == 'pseudoscalar' ]; then
    chan='g5'
elif [ $checkChan == 'vectorConc' ]; then
    chan='vec'
elif [ $checkChan == 'axialConc' ]; then
    chan='ax'
elif [ $checkChan == 'scalar' ]; then
    chan='g0'
elif [ $checkChan == 'axialZero' ]; then    
    chan='g6'
elif [ $checkChan == 'tempTensor' ]; then     # Temporal tensor
    chan='teto'
elif [ $checkChan == 'spatTensor' ]; then     # Spatial tensor
    chan='tesp'
elif [ $checkChan == 'ratioAxVec' ]; then
    chan='axvec'
elif [ $checkChan == 'ratioG0G5' ]; then
    chan='g0g5'
elif [ $checkChan == 'ratioVecG5' ]; then
    chan='vecg5'
fi

# We need to generate a cleaned file to make the analysis easier
file=$( ls ${chan}* )
fMass=mass_$file

# Now we need to do a series of commands to get rid of unwanted results

# Eliminate the Amplitude from the data
sed '/Param : 1/d' ./$file > $fMass
sed -i 's/Param : 2 is //' ./$fMass

# Eliminate the chiSquare obtained from the data
sed -i '/ChiSquare/d' ./$fMass

# Eliminate unwanted comments
sed -i '/# --/d' ./$fMass
sed -i 's/Time Window is //' ./$fMass
sed -i 's/+-//' ./$fMass

# Remove space between lines
sed -i 'N;s/\n/ /' ./$fMass

# I can remove the last point which is nan
sed -i '$d' ./$fMass

touch tempFile
# Loop to generate a new file

while read LINE
do
    tZero=$( echo $LINE | awk '{ print $1 }' | sed 's/-.*//' )
    mVlue=( $( echo $LINE | awk '{ print $2 " " $3 }' ))
    echo $tZero " " ${mVlue[@]} >> tempFile
done < $fMass

mv tempFile $fMass


