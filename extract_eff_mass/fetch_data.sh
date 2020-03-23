#!/bin/bash

NT=( 128 64 56 48 40 36 32 28 24 20 16 )
MESON=( uu us uc ss sc cc )

case $1 in
    'g5') CHANNEL='Pseudoscalar'; fold='g5';;
    'vec') CHANNEL='Vector_Spatial'; fold='vecConc';;
    'ax_plus') CHANNEL='Ax_plus_Spatial'; fold='axplusConc';;
    'ax_minus') CHANNEL='Ax_minus_Spatial'; fold='axminusConc';;
    'g0') CHANNEL='Scalar'; fold='g0';;
    *) echo 'ERROR: Channel provided does not exist!'; exit 1 ;;
esac

function gen_plot() {
# Define the variables for the plot file
chan=$1; meson=$2; nt=$3; sourc=$4; output_file=$5

echo "
set style line 1 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2

set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror

set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12

set xlabel '\$\\tau / a_\\tau\$'
set ylabel '\$M_{eff} \\cdot a_\\tau \$'
set title '${chan} -- ${meson} -- ${nt}x32 -- ${type}'

plot '${output_file}' u 1:2:3 w yerr ls 1 notitle
" > ./plot.gn
}

SCW="/home/s.989336/Analyse_Corrs/data"
ROOT=$( pwd )

for nt in ${NT[@]}
do

    echo "I am currently calculating $nt"
    for meson in ${MESON[@]}
    do
        # Create the folder that will hold the results
        mkdir -p ${CHANNEL}/${meson}/${nt}x32
        
        # Loop over sources ll and ss
        for type in ll ss
        do

            echo "I am inside $meson - $type"
            # Copy the file into the folder
            loc="$SCW/${nt}x32/concMeson.${meson}/${type}/${fold}_*/*"
            mkdir -p ${CHANNEL}/${meson}/${nt}x32/${type}
            scp scw:${loc} ${CHANNEL}/${meson}/${nt}x32/${type}

        # Then execute the script and plot
        cp extract_cosh_mass.py \
            ${CHANNEL}/${meson}/${nt}x32/${type}/

        cd ${CHANNEL}/${meson}/${nt}x32/${type}/
        name_conf=$( ls Gen2l* )
        output_file="effMass_${1}_${type}_${nt}x32.${meson}"

        python ./extract_cosh_mass.py ${name_conf} ${output_file}

        gen_plot $1 $meson $nt $type $output_file
        gnp2tex -f plot.gn -s $1_effMass_${nt}x32_${type}_${meson}.pdf

        # Remove unneeded files
        rm extract_cosh_mass.py Gen2l_* plot.gn

        cd $ROOT 
        done
    done
done
