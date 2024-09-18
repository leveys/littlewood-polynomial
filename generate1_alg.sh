#!/bin/bash

gcc -Wall littlewood_alg.c -lm

./a.out $1 images/lw_alg_$1.pbm $2
convert images/lw_alg_$1.pbm images/lw_alg_$1.png
rm images/lw_alg_$1.pbm