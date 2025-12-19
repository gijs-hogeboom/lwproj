import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pyvista as pv

path = '/home/gijs-hogeboom/dev/mclw/data_output/raw_output_3D'

case_pesc = 1

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
    hr_sfc = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100) * 86400
hr_sfc = hr_sfc.reshape((jtot, ktot))

with open(os.path.join(path, file_TOA), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_TOA = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100) * 86400
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

pl = pv.Plotter()
pl.add_volume(np.transpose(hr, (2,1,0)), clim=[-vmax, vmax], cmap='seismic', opacity=[1,0,1])
pl.add_mesh(img_sfc, clim=[-vmax, vmax], cmap='seismic')
pl.show()

plt.figure(figsize = (7,7))
plt.imshow(hr[:,25,:], origin='lower', cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar()
plt.show()
plt.close()


# hr1D = np.mean(hr, axis=(1,2))

# plt.figure()
# plt.plot(hr1D, 25 + 50*np.arange(len(hr1D)))
# plt.show()
# plt.close()

