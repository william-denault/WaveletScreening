#!/bin/bash

parallel -j 20 --progress -a slice_definitions.txt --colsep ' ' ./make_slice.sh 
