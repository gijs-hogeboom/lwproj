#!/bin/bash

CASE="s3D1"
INTERCELL_TECHNIQUE="power"
Nphot=28
Pesc=1
scatter=0
N=10

# For 3D heating rates outputs
for Nphot in {10..16}; do
	for Pesc in 0 1; do
		for ((i=0; i<N; i++)); do
			../build/lwproj $CASE $INTERCELL_TECHNIQUE $Nphot $Pesc $scatter
			PATH_IN_1="../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot}_${INTERCELL_TECHNIQUE}_Pesc${Pesc}_scatter${scatter}.dat"
			PATH_OUT_1="../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot}_${INTERCELL_TECHNIQUE}_Pesc${Pesc}_scatter${scatter}_${i}.dat"

			PATH_IN_2="../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot}_${INTERCELL_TECHNIQUE}_Pesc${Pesc}_scatter${scatter}.dat"
			PATH_OUT_2="../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot}_${INTERCELL_TECHNIQUE}_Pesc${Pesc}_scatter${scatter}_${i}.dat"

			PATH_IN_3="../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot}_${INTERCELL_TECHNIQUE}_Pesc${Pesc}_scatter${scatter}.dat"
			PATH_OUT_3="../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot}_${INTERCELL_TECHNIQUE}_Pesc${Pesc}_scatter${scatter}_${i}.dat"
			mv "$PATH_IN_1" "$PATH_OUT_1"
			mv "$PATH_IN_2" "$PATH_OUT_2"
			mv "$PATH_IN_3" "$PATH_OUT_3"
		done
	done
done

