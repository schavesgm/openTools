#!/bin/bash

NT=( 128 64 56 48 40 36 32 28 24 20 16 )
chan="tespConc"

for nt in ${NT[@]}
do

    mkdir -p ${nt}x32.cc
    for type in ll ss
    do
        # Copy the file into the folder
        locat="/home/s.989336/Fastsum-Lat/${nt}x32/concMeson.cc/${type}/${chan}_Folder/*"
        mkdir -p ${nt}x32.cc/${type}
        scp scw:${locat} ${nt}x32.cc/${type}

        # Then execute the script and plot
        cp extract_cosh_mass.py ${nt}x32.cc/${type}/
        cd ${nt}x32.cc/${type}/

        nameFile=$( ls Gen2l* )
        fileEff="effMass_${type}_${nameFile}"

        python ./extract_cosh_mass.py ${nameFile} ${fileEff}
        echo "

        set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green

        set style line 11 lc rgb '#808080' lt 1
        set border 3 back ls 11
        set tics nomirror

        set style line 12 lc rgb '#808080' lt 0 lw 1
        set grid back ls 12

        set xlabel '\$\\tau / a_\\tau\$'
        set ylabel '\$M_{eff} \\cdot a_\\tau \$'
        set title '${chan}.cc -- ${nt}x32 -- ${type}'

        plot '${fileEff}' u 1:2:3 w yerr ls 1 notitle
        " > ./plot.gn

        gnp2tex plot.gn effMass_${chan}_${nt}x32_${type}.pdf

        cd ../../
    done
done
