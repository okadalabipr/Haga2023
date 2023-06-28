#!/bin/sh

cat 01_outgoing/x.txt | sort | uniq > x.txt
mv x.txt 01_outgoing/x.txt
