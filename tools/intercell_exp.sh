#!/bin/bash

Niter=10
NPhot_start=$1
NPhot_end=$2
Pesc_mode=0
scatter=0

for ((i=7; i<Niter; i++)); do
	for ((j=0; j<36; j++)) do
		for ((k=NPhot_start; k<=NPhot_end; k++)); do
			for INTERCELL_TECHNIQUE in power uniform power-gradient; do
				../build/lwproj "gpt${j}" "$INTERCELL_TECHNIQUE" "$k" "$Pesc_mode" "$scatter"
				mv "../data_output/raw_output_3D/hr_3D_atm_gpt${j}_Nphot${k}.00_${INTERCELL_TECHNIQUE}_Pesc${Pesc_mode}_scatter${scatter}.dat" "../data_output/raw_output_3D/IT_exp/hr_3D_atm_gpt${j}_Nphot${k}.00_${INTERCELL_TECHNIQUE}_Pesc${Pesc_mode}_scatter${scatter}_${i}.dat"
				mv "../data_output/raw_output_3D/flux_3D_sfc_gpt${j}_Nphot${k}.00_${INTERCELL_TECHNIQUE}_Pesc${Pesc_mode}_scatter${scatter}.dat" "../data_output/raw_output_3D/IT_exp/flux_3D_sfc_gpt${j}_Nphot${k}.00_${INTERCELL_TECHNIQUE}_Pesc${Pesc_mode}_scatter${scatter}_${i}.dat"
				mv "../data_output/raw_output_3D/flux_3D_TOA_gpt${j}_Nphot${k}.00_${INTERCELL_TECHNIQUE}_Pesc${Pesc_mode}_scatter${scatter}.dat" "../data_output/raw_output_3D/IT_exp/flux_3D_TOA_gpt${j}_Nphot${k}.00_${INTERCELL_TECHNIQUE}_Pesc${Pesc_mode}_scatter${scatter}_${i}.dat"
			done
		done
	done
done
