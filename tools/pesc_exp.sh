#!/bin/bash

Niter=12
Niter_start=$1
NPhot_start=$2
NPhot_end=32
GPT_start=$3
GPT_end=8
Pesc_start=$4
Pesc_end=1
scatter=$5

flag_firstloop_Niter=0
flag_firstloop_NPhot=0
flag_firstloop_GPT=0
flag_firstloop_Pesc=0


if [ $flag_firstloop_Niter -eq 1 ]; then
	Niter_start=0
fi
for ((i=0; i<Niter; i++)); do
	if [ $flag_firstloop_GPT -eq 1 ]; then
		GPT_start=4
	fi
	for ((j=GPT_start; j<=GPT_end; j++)) do

		if [ $flag_firstloop_NPhot -eq 1 ]; then
			NPhot_start=22
		fi
		for ((k=NPhot_start; k<=NPhot_end; k++)); do
			if [ $flag_firstloop_Pesc -eq 1 ]; then
				Pesc=0
			fi
			for ((l=Pesc_start; l<=Pesc_end; l++)); do
				../build/lwproj "r3D${j}" "power" "$k" "$l" "$scatter"
				mv "../data_output/raw_output_3D/hr_3D_atm_r3D${j}_Nphot${k}.00_power_Pesc${l}_scatter${scatter}.dat" "../data_output/raw_output_3D/pesc_exp/hr_3D_atm_r3D${j}_Nphot${k}.00_power_Pesc${l}_scatter${scatter}_${i}.dat"
				mv "../data_output/raw_output_3D/flux_3D_sfc_r3D${j}_Nphot${k}.00_power_Pesc${l}_scatter${scatter}.dat" "../data_output/raw_output_3D/pesc_exp/flux_3D_sfc_r3D${j}_Nphot${k}.00_power_Pesc${l}_scatter${scatter}_${i}.dat"
				mv "../data_output/raw_output_3D/flux_3D_TOA_r3D${j}_Nphot${k}.00_power_Pesc${l}_scatter${scatter}.dat" "../data_output/raw_output_3D/pesc_exp/flux_3D_TOA_r3D${j}_Nphot${k}.00_power_Pesc${l}_scatter${scatter}_${i}.dat"
			done
			flag_firstloop_Pesc=1
		done
		flag_firstloop_NPhot=1
	done
	flag_firstloop_GPT=1 
done
flag_firstloop_Niter=1
