#!/bin/bash

STARTING_GPT=$1
INTERCELL_TECHNIQUE=$2
Nphot=$3
Pesc_mode=$4
scatter=$5

for ((i=$STARTING_GPT; i<36; i++));  do
  ../build/lwproj "r3D${i}" "$INTERCELL_TECHNIQUE" "$Nphot" "$Pesc_mode" "$scatter"
done
