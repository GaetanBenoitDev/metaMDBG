#!/bin/sh

# Run this from the build subdirectory.
# Usage: benchmark.sh filename

file=$1
file2=`echo $1 | sed 's#.*/##'`
test_dir=${TEST_DIR:-./tests}
r4x8=$test_dir/rans4x8
r4x16=$test_dir/rans4x16pr
ntrials=${ntrials:-5}

awkscript='BEGIN {e1=99999;e2=0;d1=99999;d2=0} /bytes/ {if (e1 > $1) {e1 = $1} if (e2 < $1) {e2 = $1} if (d1 > $4) {d1 = $4} if (d2 < $4) {d2 = $4};s=$10} END {print e1,e2,d1,d2,s}'

echo   "Program Opts             Size     Encode Decode"
echo   "-----------------------------------------------"

# Order-0
set -- $(for i in `seq 1 $ntrials`;do
             $r4x8 -t -o0 $file 2>&1
         done | awk "$awkscript")
printf "r4x8    -o0           %10d %6.1f %6.1f\n" $5 $2 $4

set -- $(for i in `seq 1 $ntrials`;do
             $r4x16 -t -o0 $file 2>&1
        done | awk "$awkscript")
printf "r4x16   -o0           %10d %6.1f %6.1f\n" $5 $2 $4

for c in 0x0000 0x0101 0x0202 0x0404
do
    set -- $(for i in `seq 1 $ntrials`;do
             $r4x16 -t -o4 -c$c $file 2>&1
        done | awk "$awkscript")
    printf "r32x16  -o4 -c %-4s %10d %6.1f %6.1f\n" $c $5 $2 $4
done

echo

# Order-1
set -- $(for i in `seq 1 $ntrials`;do
             $r4x8 -t -o1 $file 2>&1
         done | awk "$awkscript")
printf "r4x8    -o1           %10d %6.1f %6.1f\n" $5 $2 $4

set -- $(for i in `seq 1 $ntrials`;do
             $r4x16 -t -o1 $file 2>&1
         done | awk "$awkscript")
printf "r4x16   -o1           %10d %6.1f %6.1f\n" $5 $2 $4

for c in 0x0000 0x0101 0x0202 0x0404
do
    set -- $(for i in `seq 1 $ntrials`;do
             $r4x16 -t -o5 -c$c $file 2>&1
        done | awk "$awkscript")
    printf "r32x16  -o5 -c %-4s %10d %6.1f %6.1f\n" $c $5 $2 $4
done
