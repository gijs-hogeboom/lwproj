#!/bin/bash

CASE="s3D1"
INTERCELL_TECHNIQUE="power"
Nphot=20
Pesc=0
scatter=0
N=10

# For 1D heating rates
for Nphot in {16..28}; do
	for Pesc in 0 1; do
		for ((i=0; i<N; i++)); do
			../build/lwproj $CASE $INTERCELL_TECHNIQUE $Nphot $Pesc $scatter
			PATH_IN="../data_output/heating_rates/HR_MC_${CASE}_${INTERCELL_TECHNIQUE}_Nphot${Nphot}_Pesc${Pesc}_scatter${scatter}.csv"
			PATH_OUT="../data_output/heating_rates/HR_MC_${CASE}_${INTERCELL_TECHNIQUE}_Nphot${Nphot}_Pesc${Pesc}_scatter${scatter}_${i}.csv"
			mv "$PATH_IN" "$PATH_OUT"
		done
	done
done


