#!/bin/bash

eps=$1
nt=$2
index=$3

nc=1
prog=0
det=1
algo="clusteringMT.py"
progSpeedDiv=4
pt=0.0025
eps2=5
eps3=90

exec="../ttk-dev/install/bin/ttkMergeTreeClusteringCmd"

for j in `seq 0 1`
do
  ####################################################################################
  if [ $index -eq 1 -o $index -eq 0 ]; then
    for coef in ST.xml JT.xml; do
      echo "_____________ $algo startingVortexGoodEnsemble $pt $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
      $exec -i ../data/startingVortexGoodEnsemble_TA/startingVortexGoodEnsemble_$coef -P $pt -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt 

      dataSet="startingVortex_$(echo $coef | cut -d. -f1)"
      mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
    done
  fi

  ####################################################################################
  if [ $index -eq 2 -o $index -eq 0 ]; then
    coef=ST.xml
    echo "_____________ $algo isabella_velocity_goodEnsemble 0.01 $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
    $exec -i ../data/isabella_velocity_goodEnsemble_TA/isabella_velocity_goodEnsemble_$coef -P 0.01 -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

    dataSet="isabel_$(echo $coef | cut -d. -f1)"
    mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
  fi

  ####################################################################################
  if [ $index -eq 3 -o $index -eq 0 ]; then
    for coef in ST.xml JT.xml; do
      echo "_____________ $algo seaSurfaceHeightGoodEnsemble 0.01 $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
      $exec -i ../data/seaSurfaceHeightGoodEnsemble_TA/seaSurfaceHeightGoodEnsemble_$coef -P 0.01 -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

      dataSet="seaSurfaceHeight_$(echo $coef | cut -d. -f1)"
      mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
    done
  fi

  ####################################################################################
  if [ $index -eq 4 -o $index -eq 0 ]; then
    for coef in ST.xml JT.xml; do
      echo "_____________ $algo vortexStreetGoodEnsemble2 $pt $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
      $exec -i ../data/vortexStreetGoodEnsemble2_TA/vortexStreetGoodEnsemble2_$coef -P $pt -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt 

      dataSet="vortexStreet_$(echo $coef | cut -d. -f1)"
      mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
    done
  fi

  ####################################################################################
  if [ $index -eq 5 -o $index -eq 0 ]; then
    coef=ST.xml
    echo "_____________ $algo particularEnsemble $pt $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
    $exec -i ../data/particularEnsemble_TA/particularEnsemble_$coef -P $pt -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

    dataSet="viscousFingers_$(echo $coef | cut -d. -f1)"
    mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
  fi

  ####################################################################################
  if [ $index -eq 6 -o $index -eq 0 ]; then
    coef=ST.xml
    echo "_____________ $algo cloud5 0.1 $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
    $exec -i ../data/cloud5_TA/cloud5_$coef -P 0.1 -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

    dataSet="cloudProcesses_$(echo $coef | cut -d. -f1)"
    mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
  fi

  ####################################################################################
  if [ $index -eq 7 -o $index -eq 0 ]; then
    coef=ST.xml
    echo "_____________ $algo astroTurbulence $pt $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
    $exec -i ../data/astroTurbulence_TA/astroTurbulence_$coef -P $pt -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

    dataSet="ionizationFront2D_$(echo $coef | cut -d. -f1)"
    mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
  fi

  ####################################################################################
  if [ $index -eq 8 -o $index -eq 0 ]; then
    coef=ST.xml
    echo "_____________ $algo impactEnsemble3CTev $pt $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
    $exec -i ../data/impactEnsemble3CTev_TA/impactEnsemble2CTev_$coef -P $pt -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

    dataSet="asteroidImpact_$(echo $coef | cut -d. -f1)"
    mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
  fi

  ####################################################################################
  if [ $index -eq 9 -o $index -eq 0 ]; then
    coef=ST.xml
    echo "_____________ $algo darkSky100S 0.1 $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
    $exec -i ../data/darkSky100S_TA/darkSky100S_$coef -P 0.1 -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

    dataSet="darkMatter_$(echo $coef | cut -d. -f1)"
    mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
  fi

  ####################################################################################
  if [ $index -eq 10 -o $index -eq 0 ]; then
    coef=ST.xml
    echo "_____________ $algo volcanic2 $pt $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
    $exec -i ../data/volcanic2_TA/volcanic2_$coef -P $pt -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

    dataSet="volcanicEruption_$(echo $coef | cut -d. -f1)"
    mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
  fi

  ####################################################################################
  if [ $index -eq 11 -o $index -eq 0 ]; then
    coef=ST.xml
    echo "_____________ $algo astro3DTurbulence 0.01 $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
    $exec -i ../data/astro3DTurbulence_TA/astro3DTurbulence_$coef -P 0.01 -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

    dataSet="ionizationFront3D_$(echo $coef | cut -d. -f1)"
    mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
  fi

  ####################################################################################
  if [ $index -eq 12 -o $index -eq 0 ]; then
    coef=ST.xml
    echo "_____________ $algo earthquake2 $pt $eps $eps2 $eps3 $coef $prog $det $progSpeedDiv $nt _____________"
    $exec -i ../data/earthquake2_TA/earthquake2_$coef -P $pt -E1 $eps -E2 $eps2 -E3 $eps3 -K $nc -D $det -T $nt

    dataSet="earthquake_$(echo $coef | cut -d. -f1)"
    mv output_port\#4.vtu outputs/"$dataSet"_barycenter.vtu
  fi
done

exit
