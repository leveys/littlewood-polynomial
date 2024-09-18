#!/bin/bash

gcc -Wall littlewood_alg.c -lm

for ((i=1; i<=$1; i++)); do
   ./a.out $i images/lw_alg_$i.pbm $2
   convert images/lw_alg_$i.pbm images/lw_alg_$i.png
   rm images/lw_alg_$i.pbm
done