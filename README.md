# Wasserstein Merge Trees

This github repository contains the exact code used for the reference below.

The module has been integrated to TTK https://github.com/topology-tool-kit/ttk and has been improved since the publication of the paper.

We recommend to use the module now present in TTK to enjoy the last features (for example, the new module is in average 3 times faster).

## Reference

If you plan to use this code to generate results for a scientific document, thanks for referencing the following publication:

"Wasserstein Distances, Geodesics and Barycenters of Merge Trees"  
Mathieu Pont, Jules Vidal, Julie Delon, Julien Tierny.  
Proc. of IEEE VIS 2021.  
IEEE Transactions on Visualization and Computer Graphics, 2021  

## Installation note

The following procedure has been tested on Ubuntu 18.04.5 LTS.

## Install the dependencies

```bash
sudo apt-get install cmake-qt-gui libboost-system-dev libpython3.6-dev libxt-dev
sudo apt-get install qt5-default qttools5-dev libqt5x11extras5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools 
sudo apt-get install python3-sklearn 
```

## Install Paraview

Extract ParaView:

```bash
tar xvJf ParaView-v5.7.0.tar.xz
```

Patch ParaView with TTK:

```bash
cd ttk-dev/paraview/patch
chmod u+x patch-paraview-5.7.0.sh
./patch-paraview-5.7.0.sh ../../../ParaView-v5.7.0/
```

Install ParaView:
(replace the 4 in "make -j4" by the number of available cores on your system)

```bash
cd ../../../ParaView-v5.7.0/
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_PYTHON_VERSION=3 -DCMAKE_INSTALL_PREFIX=../install ..
make -j4
make install
```

## Install TTK

(replace the 4 in "make -j4" by the number of available cores on your system)

```bash
cd ttk-dev/
mkdir build && cd build
paraviewPath=`pwd`/../../ParaView-v5.7.0/install/lib/cmake/paraview-5.7
cmake -DCMAKE_INSTALL_PREFIX=../install -DParaView_DIR=$paraviewPath ..
make -j4
make install
```

## Get the results

Extract the data:

```bash
tar xvJf data.tar.xz
```

### table 1

To reproduce the results of Table 1 in the paper, please enter the following commands:

```bash
cd scripts
for f in *.sh; do chmod u+x $f; done
```

Run the experiments (it will take a LONG time) and print table:
(replace N with the number of available cores on your system)

```bash
./automata5.sh N
./table1.sh
```

The barycenters computed during the benchmarks were saved in the `outputs` folder, using a default planar layout.

You can also do the computation for one dataset at a time with:

```
./runVortexStreet.sh N
./runStartingVortex.sh N
./runViscousFingering.sh N
./runIonizationFront2D.sh N
./runVolcanicEruptions.sh N
./runIonizationFront3D.sh N
./runAsteroidImpact.sh N
./runEarthquake.sh N
./runCloudProcesses.sh N
./runIsabel.sh N
./runDarkMatter.sh N
./runSeaSurfaceHeight.sh N
```

And finally the dataset ordered by size (number of members * number of persistence pairs):
- Vortex Street (45 * 23)
- Starting Vortex (12 * 124)
- Viscous fingering (15 * 118)
- Ionization front 2D (16 * 135)
- Volcanic eruptions (12 * 811)
- Ionization front 3D (16 * 763)
- Asteroid Impact (7 * 1295)
- Earthquake (12 * 1203)
- Cloud processes (12 * 1209)
- Isabel (12 * 1338)
- Dark matter (12 * 2592)
- Sea Surface Height (48 * 1787)
