#!/bin/bash
#!
## addGrch37coords.sh
## author: Courtney French cf458@cam.ac.uk
## updated last: 23 Apr 2019
## This script uses CrossMap to annotate family filtered files
## with grch37 coordinates


. /etc/profile.d/modules.sh                # Leave this line (enables the module command)

#Defining variables

input=$1
bed=${input/.txt/.tmp.bed}
bed2=${input/.txt/_liftover.tmp.bed}
tmp=${input/.txt/.tmp.txt}

script_path=$2
crossmap=$3
chain=$4
new_col=$5
python_script="${script_path}/addLiftoverCoords.py"

# make bed file

tail -n+2 $input | awk -F'\t' '{ OFS = "\t"; print $1, $2, $2, $1":"$2 }' > $bed

# CrossMap

python $crossmap bed $chain $bed $bed2

# add coords to original file

python $python_script --annot_file $bed2 --in_file $input --new_col $new_col > $tmp

mv $tmp $input

rm $bed
rm $bed2
