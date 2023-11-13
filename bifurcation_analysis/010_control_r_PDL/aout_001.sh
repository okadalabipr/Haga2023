#!/bin/sh

rm -r 01_outgoing
mkdir 01_outgoing
gcc  -O3  pg_001_main.c  pg_001_func.c  pg_001_outg.c  -lm  -Wall  -o  a_001.out
./a_001.out
rm a_001.out
