#!/bin/bash

mkdir outputs

nohup ./automata4.sh 5 $1 $2 > nohupBaryMT.out
nohup ./automata4.sh 5 1 $2 > nohupBaryMT_NT1.out
nohup ./automata4.sh 100 1 $2 > nohupBaryMTE100_NT1.out

mkdir nohupFiles
mkdir nohupFiles/nohupBary

if [ -f nohupFiles/nohupBary/nohupBaryMT.out ]; then
  cat nohupBaryMT.out >> nohupFiles/nohupBary/nohupBaryMT.out
  rm nohupBaryMT.out
else
  mv nohupBaryMT.out nohupFiles/nohupBary/
fi

if [ -f nohupFiles/nohupBary/nohupBaryMT_NT1.out ]; then
  cat nohupBaryMT_NT1.out >> nohupFiles/nohupBary/nohupBaryMT_NT1.out
  rm nohupBaryMT_NT1.out
else
  mv nohupBaryMT_NT1.out nohupFiles/nohupBary/
fi

if [ -f nohupFiles/nohupBary/nohupBaryMTE100_NT1.out ]; then
  cat nohupBaryMTE100_NT1.out >> nohupFiles/nohupBary/nohupBaryMTE100_NT1.out
  rm nohupBaryMTE100_NT1.out
else
  mv nohupBaryMTE100_NT1.out nohupFiles/nohupBary/
fi
