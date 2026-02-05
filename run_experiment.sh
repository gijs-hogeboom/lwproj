#!/bin/bash

#Niter=10
#for ((i=0;i<Niter;i++)); do
#	for ((j=26;j<29;j++)); do
#		for ((P=0;P<2;P++)); do
#			./lwproj "s3D1" "power" "$j" "$P" "0"
#			mv "../data_output/raw_output_3D/hr_3D_atm_s3D1_Nphot${j}.00_power_Pesc${P}_scatter0.dat" "../data_output/raw_output_3D/hr_3D_atm_s3D1_Nphot${j}.00_power_Pesc${P}_scatter0_${i}.dat"
#			mv "../data_output/raw_output_3D/flux_3D_sfc_s3D1_Nphot${j}.00_power_Pesc${P}_scatter0.dat" "../data_output/raw_output_3D/flux_3D_sfc_s3D1_Nphot${j}.00_power_Pesc${P}_scatter0_${i}.dat" 
#			mv "../data_output/raw_output_3D/flux_3D_TOA_s3D1_Nphot${j}.00_power_Pesc${P}_scatter0.dat" "../data_output/raw_output_3D/flux_3D_TOA_s3D1_Nphot${j}.00_power_Pesc${P}_scatter0_${i}.dat"
#		done
#	done
#done


CASE="r3D4"
Nphot0="32.00"
Nphot1="31.54"
Niter=10
for ((i=0;i<Niter;i++)); do
	./lwproj "${CASE}" "power" "${Nphot0}" "0" "0"
	./lwproj "${CASE}" "power" "${Nphot1}" "1" "0"
	mv "../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0.dat" "../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0_${i}.dat"
	mv "../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0.dat" "../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0_${i}.dat" 
	mv "../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0.dat" "../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0_${i}.dat"

	mv "../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0.dat" "../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0_${i}.dat"
	mv "../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0.dat" "../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0_${i}.dat" 
	mv "../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0.dat" "../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0_${i}.dat"
done



CASE="r3D5"
Nphot0="32.00"
Nphot1="31.00"
Niter=10
for ((i=0;i<Niter;i++)); do
	./lwproj "${CASE}" "power" "${Nphot0}" "0" "0"
	./lwproj "${CASE}" "power" "${Nphot1}" "1" "0"
	mv "../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0.dat" "../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0_${i}.dat"
	mv "../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0.dat" "../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0_${i}.dat" 
	mv "../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0.dat" "../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot0}_power_Pesc0_scatter0_${i}.dat"

	mv "../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0.dat" "../data_output/raw_output_3D/hr_3D_atm_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0_${i}.dat"
	mv "../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0.dat" "../data_output/raw_output_3D/flux_3D_sfc_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0_${i}.dat" 
	mv "../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0.dat" "../data_output/raw_output_3D/flux_3D_TOA_${CASE}_Nphot${Nphot1}_power_Pesc1_scatter0_${i}.dat"
done

