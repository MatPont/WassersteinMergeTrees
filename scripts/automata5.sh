mkdir outputs

nohup ./automata4.sh 5 $1 > nohupBaryMT.out
nohup ./automata4.sh 5 1 > nohupBaryMT_NT1.out
nohup ./automata4.sh 100 1 > nohupBaryMTE100_NT1.out

mkdir nohupFiles
mkdir nohupFiles/nohupBary
mv nohupBaryMT.out nohupFiles/nohupBary/
mv nohupBaryMT_NT1.out nohupFiles/nohupBary/
mv nohupBaryMTE100_NT1.out nohupFiles/nohupBary/
