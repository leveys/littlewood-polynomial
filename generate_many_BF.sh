#!/bin/bash

export LD_LIBRARY_PATH=/usr/local/include/gsl/lib

gcc -Wall -I/usr/local/include/gsl/include -c littlewood_BF.c
gcc -L/usr/local/include/gsl/lib littlewood_BF.o -lgsl -lgslcblas -lm

for ((i=1; i<=$1; i++)); do
   ./a.out $i images/lw_BF_$i.pbm $2
   convert images/lw_BF_$i.pbm images/lw_BF_$i.png
   rm images/lw_BF_$i.pbm
done