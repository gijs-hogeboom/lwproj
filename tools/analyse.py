import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pyvista as pv


path = '../data_output/raw_output_3D'

Case = 's3D1'
InterCellTechnique = 'power'
Nphot = '28.00'
PescMode = 1
enable_scattering = 0

file_atm = f'hr_3D_atm_{Case}_Nphot{Nphot}_{InterCellTechnique}_Pesc{PescMode}_scatter{enable_scattering}.dat'
file_sfc = f'flux_3D_sfc_{Case}_Nphot{Nphot}_{InterCellTechnique}_Pesc{PescMode}_scatter{enable_scattering}.dat'
file_TOA = f'flux_3D_TOA_{Case}_Nphot{Nphot}_{InterCellTechnique}_Pesc{PescMode}_scatter{enable_scattering}.dat'

with open(os.path.join(path, file_atm), 'rb') as f:
    itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)
    hr = np.fromfile(f, dtype=np.float32)
hr = hr.reshape((itot, jtot, ktot))

with open(os.path.join(path, file_sfc), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_sfc = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100 * 100) * 86400
hr_sfc = hr_sfc.reshape((jtot, ktot))

with open(os.path.join(path, file_TOA), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_TOA = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100 * 100) * 86400
hr_TOA = hr_TOA.reshape((jtot, ktot))

img_sfc = pv.ImageData(dimensions=(ktot, jtot, 1))
img_sfc['Values'] = hr_sfc.ravel()
img_TOA = pv.ImageData(dimensions=(ktot, jtot, 1))
img_TOA['Values'] = hr_TOA.ravel()



hr = hr[:50,:,:]


vmax_atm = np.max(np.abs(hr))
vmax_sfc = np.max(np.abs(hr_sfc))
print(vmax_atm, vmax_sfc)
vmax = np.min([vmax_atm, vmax_sfc])
vmax = 12

# pl = pv.Plotter()
# pl.add_volume(np.transpose(hr, (2,1,0)), clim=[-vmax_atm, vmax_atm], cmap='seismic', opacity=[1,0,1])
# pl.add_mesh(img_sfc, clim=[-vmax_sfc, vmax_sfc], cmap='seismic')
# pl.show()

plt.figure(figsize = (7,7))
plt.imshow(hr[:,25,:], origin='lower', cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar()
plt.show()
plt.close()







# # Total energy via PP
# itot, jtot, ktot = 209, 192, 192
# Matm = np.zeros((itot, jtot, ktot))
# Msfc = np.zeros((jtot, ktot))
# MTOA = np.zeros((jtot, ktot))

# for igpt in range(36):
#     CASE = f'r3D{igpt}'
#     file_atmPP = f'../data_output/raw_output_3DPP/hr_3DPP_atm_{CASE}.dat'
#     file_sfcPP = f'../data_output/raw_output_3DPP/flux_3DPP_sfc_{CASE}.dat'
#     file_TOAPP = f'../data_output/raw_output_3DPP/flux_3DPP_TOA_{CASE}.dat'

#     with open(file_atmPP, 'rb') as f:
#         itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)
#         hr = np.fromfile(f, dtype=np.float32)
#     hr = hr.reshape((itot, jtot, ktot))

#     with open(file_sfcPP, 'rb') as f:
#         jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
#         hr_sfc = np.fromfile(f, dtype=np.float32)
#     hr_sfc = hr_sfc.reshape((jtot, ktot))/(1.225* 1004)

#     with open(file_TOAPP, 'rb') as f:
#         jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
#         hr_TOA = np.fromfile(f, dtype=np.float32)
#     hr_TOA = hr_TOA.reshape((jtot, ktot))/(1.225* 1004)

#     Matm += hr
#     Msfc += hr_sfc
#     MTOA += hr_TOA
#     print(igpt)


# img_sfcPP = pv.ImageData(dimensions=(ktot, jtot, 1))
# img_sfcPP['Values'] = Msfc.ravel()
# img_TOAPP = pv.ImageData(dimensions=(ktot, jtot, 1))
# img_TOAPP['Values'] = MTOA.ravel()

# vmax_atm = np.max(np.abs(Matm))
# vmax_sfc = np.max(np.abs(Msfc))
# vmax = np.min([vmax_atm, vmax_sfc])


# pl = pv.Plotter()
# pl.add_volume(np.transpose(Matm, (2,1,0)), clim=[-vmax_atm, vmax_atm], cmap='coolwarm', opacity=[1,0,1])
# pl.add_mesh(img_sfcPP, clim=[-vmax_sfc, vmax_sfc], cmap='coolwarm')
# pl.show()


# # plt.figure(figsize = (5,10))
# # plt.plot(Matm.mean(axis=(1,2)), np.arange(25, 25 + 50*209, 50))
# # plt.ylim([0,5000])
# # plt.grid()
# # plt.show()
# # plt.close()