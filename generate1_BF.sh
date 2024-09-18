#!/bin/bash

export LD_LIBRARY_PATH=/usr/local/include/gsl/lib

gcc -Wall -I/usr/local/include/gsl/include -c littlewood_BF.c
gcc -L/usr/local/include/gsl/lib littlewood_BF.o -lgsl -lgslcblas -lm

./a.out $1 images/lw_BF_$1.pbm $2
convert images/lw_BF_$1.pbm images/lw_BF_$1.png
rm images/lw_BF_$1.pbm