import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d


mode = '3D'
CASE = 's3D1'
INTERCELL_TECHNIQUE = 'power'
# Nphot = 20
# Pesc = 0
scattering_enabled = 0

Niter = 10
arr_Nphot = np.arange(16, 29, 1)


def f_converge(x, a, b, c):
    return a*x**(-1/2) + c


if mode == '1D':
    path_in = '/home/gijs-hogeboom/dev/mclw/data_output/heating_rates/'

    itot = 209

    arr_std_pesc0 = np.zeros(len(arr_Nphot))
    arr_std_pesc1 = np.zeros(len(arr_Nphot))

    for j, Nphot in enumerate(arr_Nphot):
        M_Pesc0 = np.zeros((itot, Niter))
        M_Pesc1 = np.zeros((itot, Niter))
        for i in range(Niter):
            file_Pesc0 = os.path.join(path_in, f'HR_MC_{CASE}_{INTERCELL_TECHNIQUE}_Nphot{Nphot}_Pesc0_scatter{scattering_enabled}_{i}.csv')
            df_Pesc0 = pd.read_csv(file_Pesc0)
            file_Pesc1 = os.path.join(path_in, f'HR_MC_{CASE}_{INTERCELL_TECHNIQUE}_Nphot{Nphot}_Pesc1_scatter{scattering_enabled}_{i}.csv')
            df_Pesc1 = pd.read_csv(file_Pesc1)
        
            M_Pesc0[:,i] = df_Pesc0['heating_rate'].values
            M_Pesc1[:,i] = df_Pesc1['heating_rate'].values
        
        arr_std_pesc0[j] = np.mean(np.std(M_Pesc0, axis=1))
        arr_std_pesc1[j] = np.mean(np.std(M_Pesc1, axis=1))




elif mode == '3D':
    path_in = '/home/gijs-hogeboom/dev/mclw/data_output/raw_output_3D'


    # Obtaining shape of raw files
    file_metadata = os.path.join(path_in, f'hr_3D_atm_{CASE}_Nphot{arr_Nphot[0]}_{INTERCELL_TECHNIQUE}_Pesc0_scatter{scattering_enabled}_0.dat')
    with open(os.path.join(path_in, file_metadata), 'rb') as f:
        itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)

    arr_std_pesc0 = np.zeros(len(arr_Nphot))
    arr_std_pesc1 = np.zeros(len(arr_Nphot))

    for j, Nphot in enumerate(arr_Nphot):
        M_Pesc0 = np.zeros((itot, jtot, ktot, Niter))
        M_Pesc1 = np.zeros((itot, jtot, ktot, Niter))
        for i in range(Niter):
            file_Pesc0 = os.path.join(path_in, f'hr_3D_atm_{CASE}_Nphot{Nphot}_{INTERCELL_TECHNIQUE}_Pesc0_scatter{scattering_enabled}_{i}.dat')
            with open(os.path.join(path_in, file_Pesc0), 'rb') as f:
                itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)
                hr = np.fromfile(f, dtype=np.float32)
            hr_Pesc0 = hr.reshape((itot, jtot, ktot))

            file_Pesc1 = os.path.join(path_in, f'hr_3D_atm_{CASE}_Nphot{Nphot}_{INTERCELL_TECHNIQUE}_Pesc1_scatter{scattering_enabled}_{i}.dat')
            with open(os.path.join(path_in, file_Pesc1), 'rb') as f:
                itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)
                hr = np.fromfile(f, dtype=np.float32)
            hr_Pesc1 = hr.reshape((itot, jtot, ktot))
        
            M_Pesc0[:,:,:,i] = hr_Pesc0
            M_Pesc1[:,:,:,i] = hr_Pesc1
        
        arr_std_pesc0[j] = np.std(M_Pesc0, axis=3)[:, :, :].mean()
        arr_std_pesc1[j] = np.std(M_Pesc1, axis=3)[:, :, :].mean()


f_stdev_pesc0 = interp1d(arr_Nphot, arr_std_pesc0, fill_value='extrapolate')
f_stdev_pesc1 = interp1d(arr_Nphot, arr_std_pesc1, fill_value='extrapolate')
f_Nphot_pesc0 = interp1d(arr_std_pesc0, arr_Nphot, fill_value='extrapolate')
f_Nphot_pesc1 = interp1d(arr_std_pesc1, arr_Nphot, fill_value='extrapolate')


y_Nphot_pesc1_for_stdev_pesc0 = f_Nphot_pesc1(f_stdev_pesc0(arr_Nphot))

y_Nphot_pesc1_for_stdev_pesc0_extrapolated = list(y_Nphot_pesc1_for_stdev_pesc0[y_Nphot_pesc1_for_stdev_pesc0 < arr_Nphot[0]]) + [arr_Nphot[0]]
arr_Nphot_extrapolated                     = list(arr_Nphot[y_Nphot_pesc1_for_stdev_pesc0 < arr_Nphot[0]]) + [f_Nphot_pesc0(f_stdev_pesc1(arr_Nphot[0]))]
y_Nphot_pesc1_for_stdev_pesc0_interpolated = [arr_Nphot[0]] + list(y_Nphot_pesc1_for_stdev_pesc0[y_Nphot_pesc1_for_stdev_pesc0 >= arr_Nphot[0]])
arr_Nphot_interpolated                     = [f_Nphot_pesc0(f_stdev_pesc1(arr_Nphot[0]))] + list(arr_Nphot[y_Nphot_pesc1_for_stdev_pesc0 >= arr_Nphot[0]])


print(f_Nphot_pesc1(f_stdev_pesc0(32)))

plt.figure(figsize = (15,5))
plt.subplot(1,2,1)
plt.plot(arr_Nphot, arr_std_pesc0, marker='.', color='dodgerblue', label='naive')
plt.plot(arr_Nphot, arr_std_pesc1, marker='.', color='gold', label='Pesc mode')
plt.grid(which='major', alpha=0.5)
plt.grid(which='minor', alpha=0.2)
plt.minorticks_on()
plt.xlabel('2^N photons')
plt.ylabel('std (10 runs)')
plt.yscale('log')
plt.legend()
plt.title(f"{CASE} | Standard Deviation per N photons")

plt.subplot(1,2,2)
plt.plot(arr_Nphot_extrapolated, y_Nphot_pesc1_for_stdev_pesc0_extrapolated, color='black', linestyle='dashed', label='extrapolated')
plt.plot(arr_Nphot_interpolated, y_Nphot_pesc1_for_stdev_pesc0_interpolated, color='black', label='interpolated')
plt.grid(which='major', alpha=0.5)
plt.grid(which='minor', alpha=0.2)
plt.minorticks_on()
plt.xlabel('Nphot naive')
plt.ylabel('Nphot Pesc mode')
plt.legend()
plt.title(f"{CASE} | Pesc N photons equivalence")

plt.show()
plt.close()