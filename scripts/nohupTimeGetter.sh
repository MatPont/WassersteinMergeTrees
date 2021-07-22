get_min_max_avg(){
    awk 'BEGIN {max=-999999;min=999999} {max=($1>max)?$1:max;min=($1<min)?$1:min;total += $1; count++} END { count=(count==0)?1:count; printf "\\textbf{%.2f} [%.2f ; %.2f] (%d)\n", total/count,min,max,count }' $1
}

if [ $# -lt 4 ]; then
    echo "Usage PD/MT file dataset coef [e1 e2 e3]"
    exit
fi

e1=5
e2=5
e3=90
if [ $# -gt 4 ]; then
    e1=$5
    e2=$6
    e3=$7
fi

if [ $1 == "PD" ]; then
    # Get min max avg time
    echo "--- time"
    sed -e '/'$3'.*'$4' _/,/Memory/!d' $2 | grep Processed | cut -d' ' -f 4 | get_min_max_avg

    # Get min max avg iteration
    echo "--- iteration"
    sed -e '/'$3'.*'$4' _/,/Memory/!d' $2 | grep iteration -B 1 | tr '\n' ' ' | tr '\-\-' '\n' | tr ':' ' ' | rev | cut -d' ' -f 2 | rev | awk 'NF' | get_min_max_avg  
fi

if [ $1 == "MT" ]; then
    toGet=1
    if (( $(echo "$4 == 0.5" | bc -l) )); then
        toGet=2
    fi
    coef=$4
    if [ $4 != "0.5" ]; then 
        if [ $4 -eq 10 ]; then
        coef="(0|1)"
        fi 
    fi
    sed -e '/'$3'.*'${e1}.${e2}.${e3}.$coef.*' _/,/TREES OUT/!d' $2 | grep "average" | tail -$toGet

    # Get min max avg time
    echo "--- time"
    if [ $4 != "0.5" ]; then 
        if [ $4 -eq 10 ]; then
            res0=`sed -E -e '/'$3'.*'${e1}.${e2}.${e3}.0.*' _/,/TREES OUT/!d' $2 | grep "TREES OUT" -B 3 | grep "TIME BARYCENTER" | cut -d' ' -f 4 | get_min_max_avg`
            res1=`sed -E -e '/'$3'.*'${e1}.${e2}.${e3}.1.*' _/,/TREES OUT/!d' $2 | grep "TREES OUT" -B 3 | grep "TIME BARYCENTER" | cut -d' ' -f 4 | get_min_max_avg`
            avgTime0=`echo $res0 | cut -d{ -f 2 | cut -d} -f 1` 
            avgTime1=`echo $res1 | cut -d{ -f 2 | cut -d} -f 1` 
            avgTime=`echo "$avgTime0 + $avgTime1" | bc -l`
            minTime0=`echo $res0 | cut -d[ -f 2 | cut -d\; -f 1` 
            minTime1=`echo $res1 | cut -d[ -f 2 | cut -d\; -f 1` 
            minTime=`echo "$minTime0 + $minTime1" | bc -l`
            maxTime0=`echo $res0 | cut -d\; -f 2 | cut -d] -f 1` 
            maxTime1=`echo $res1 | cut -d\; -f 2 | cut -d] -f 1` 
            maxTime=`echo "$maxTime0 + $maxTime1" | bc -l`
            echo "\\textbf{"$avgTime"} ["$minTime ";" $maxTime"]" 
        else
            #sed -E -e '/'$3'.*'$coef' _/,/TREES OUT/!d' $2 | grep "TREES OUT" -B 3 | grep "TIME BARYCENTER" | cut -d' ' -f 4
            sed -E -e '/'$3'.*'${e1}.${e2}.${e3}.$coef.*' _/,/TREES OUT/!d' $2 | grep "TREES OUT" -B 3 | grep "TIME BARYCENTER" | cut -d' ' -f 4 | get_min_max_avg
        fi 
    fi
    

    # Get min max avg iteration
    echo "--- iteration"
    sed -E -e '/'$3'.*'${e1}.${e2}.${e3}.$coef.*' _/,/TREES OUT/!d' $2 | grep "TREES OUT" -B 14 | grep Iteration | cut -d' ' -f 2 | get_min_max_avg
fi
