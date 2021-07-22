#!/bin/bash

###############################################################################
# GetResults
###############################################################################
dataset=(impactEnsemble3CTev cloud5 particularEnsemble darkSky100S volcanic2 astroTurbulence astro3DTurbulence earthquake2 isabella_velocity_goodEnsemble startingVortexGoodEnsemble seaSurfaceHeightGoodEnsemble vortexStreetGoodEnsemble2)

entries=()
pairs=()
mtTime40T=()
mtTime1T=()
pdTime1T=()
speedup=()
efficiency=()

len=`expr ${#dataset[@]} - 1`
for i in `seq 0 $len`; do
    d=${dataset[$i]}
    # Number of pairs
    pairs0=`./nohupTimeGetter.sh MT ./nohupFiles/nohupBary/nohupBaryMT.out $d 0 | grep node | cut -d' ' -f 5`
    pairs1=`./nohupTimeGetter.sh MT ./nohupFiles/nohupBary/nohupBaryMT.out $d 1 | grep node | cut -d' ' -f 5`
    pairs[$i]=`expr $pairs0 + $pairs1`
    
    # Number of entries
    entries[$i]=`./nohupTimeGetter.sh MT ./nohupFiles/nohupBary/nohupBaryMT.out $d 0 | grep node | cut -d' ' -f 1`
    
    # Execution time
    mtTime40T[$i]=`./nohupTimeGetter.sh MT ./nohupFiles/nohupBary/nohupBaryMT.out $d 0 | grep text | cut -d{ -f 2 | cut -d} -f 1 | head -1`
    mtTime1T[$i]=`./nohupTimeGetter.sh MT ./nohupFiles/nohupBary/nohupBaryMT_NT1.out $d 0 | grep text | cut -d{ -f 2 | cut -d} -f 1 | head -1`      
    pdTime1T[$i]=`./nohupTimeGetter.sh MT ./nohupFiles/nohupBary/nohupBaryMTE100_NT1.out $d 0 100 5 90 | grep text | cut -d{ -f 2 | cut -d} -f 1 | head -1`  
    
    if [ $d == "gerrisTSDefault3" -o $d == "vortexStreetGoodEnsemble2" -o $d == "seaSurfaceHeightGoodEnsemble" -o $d == "startingVortexGoodEnsemble" ]; then
        mtTimeTemp40T=`./nohupTimeGetter.sh MT ./nohupFiles/nohupBary/nohupBaryMT.out $d 1 | grep text | cut -d{ -f 2 | cut -d} -f 1 | head -1`
        mtTimeTemp1T=`./nohupTimeGetter.sh MT ./nohupFiles/nohupBary/nohupBaryMT_NT1.out $d 1 | grep text | cut -d{ -f 2 | cut -d} -f 1 | head -1`
        pdTimeTemp1T=`./nohupTimeGetter.sh MT ./nohupFiles/nohupBary/nohupBaryMTE100_NT1.out $d 1 100 5 90 | grep text | cut -d{ -f 2 | cut -d} -f 1 | head -1`
        mtTime40T[$i]=$(echo "${mtTime40T[$i]} + $mtTimeTemp40T" | bc -l)
        mtTime1T[$i]=$(echo "${mtTime1T[$i]} + $mtTimeTemp1T" | bc -l)
        pdTime1T[$i]=$(echo "${pdTime1T[$i]} + $pdTimeTemp1T" | bc -l)
    fi
    speedup[$i]=$(echo "scale=2;${mtTime1T[$i]} / (${mtTime40T[$i]})" | bc -l)
    efficiency[$i]=$(echo "scale=2;${mtTime1T[$i]} / (${mtTime40T[$i]}*20)" | bc -l)
done



###############################################################################
# Create table
###############################################################################
datasetName=("Asteroid Impact \cite{scivis2018} (3D)" "Cloud processes \cite{scivis2017} (2D)" "Viscous fingering \cite{scivis2016} (3D)" "Visualize Universe \cite{scivis2015} (3D)" "Volcanic Eruptions \cite{scivis2014} (2D)" "Ionization front \cite{scivis2008} (2D)" "Ionization front \cite{scivis2008} (3D)" "Earthquake \cite{scivis2006} (3D)" "Isabel \cite{scivis2004} (3D)" "Gerris \cite{favelier2018} (2D)" "Sea Surface Height \cite{favelier2018} (2D)" "Vortex Street \cite{favelier2018} (2D)")

echo """\begin{table} 
    \centering
    \scalebox{0.65}{
    \makebox[\linewidth]{%
    \begin{tabular}{|l||r|r||r||r|r|r|}
        \hline
        \textbf{Dataset} & \$N$ & $|\branchtree|$ & $\wasserstein{2}$ (1 c.) & $\wassersteinTree$ (1 c.) & $\wassersteinTree$ (20 c.) & Speedup \\\\
        \hline"""
        
len=`expr ${#datasetName[@]} - 1`
for i in `seq 0 $len`; do
    LC_NUMERIC="en_US.UTF-8" printf "        %s & %'d & %'d & %'0.2f & %'0.2f & %'0.2f & %'0.2f \\\\\ \n" "${datasetName[$i]}" ${entries[$i]} ${pairs[$i]} ${pdTime1T[$i]} ${mtTime1T[$i]} ${mtTime40T[$i]} ${speedup[$i]}
done

echo """        \hline
    \end{tabular}
    }
    }
    \label{tab_timeSeq}
\end{table}"""
