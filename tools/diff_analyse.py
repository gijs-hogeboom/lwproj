import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pyvista as pv

path = '/home/gijs-hogeboom/dev/mclw/data_output/raw_output_3D'


Case = 's3D1'
InterCellTechnique = 'power'
Nphot = 25
# PescMode = 0
enable_scattering = 0

file_atm1 = f'hr_3D_atm_{Case}_Nphot27.00_{InterCellTechnique}_Pesc0_scatter{enable_scattering}.dat'
file_sfc1 = f'flux_3D_sfc_{Case}_Nphot27.00_{InterCellTechnique}_Pesc0_scatter{enable_scattering}.dat'
file_TOA1 = f'flux_3D_TOA_{Case}_Nphot27.00_{InterCellTechnique}_Pesc0_scatter{enable_scattering}.dat'

file_atm2 = f'hr_3D_atm_{Case}_Nphot22.69_{InterCellTechnique}_Pesc1_scatter{enable_scattering}.dat'
file_sfc2 = f'flux_3D_sfc_{Case}_Nphot22.69_{InterCellTechnique}_Pesc1_scatter{enable_scattering}.dat'
file_TOA2 = f'flux_3D_TOA_{Case}_Nphot22.69_{InterCellTechnique}_Pesc1_scatter{enable_scattering}.dat'

with open(os.path.join(path, file_atm1), 'rb') as f:
    itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)
    hr1 = np.fromfile(f, dtype=np.float32)
hr1 = hr1.reshape((itot, jtot, ktot))

with open(os.path.join(path, file_sfc1), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_sfc1 = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100) * 86400
hr_sfc1 = hr_sfc1.reshape((jtot, ktot))

with open(os.path.join(path, file_TOA1), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_TOA1 = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100) * 86400
hr_TOA1 = hr_TOA1.reshape((jtot, ktot))


with open(os.path.join(path, file_atm2), 'rb') as f:
    itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)
    hr2 = np.fromfile(f, dtype=np.float32)
hr2 = hr2.reshape((itot, jtot, ktot))

with open(os.path.join(path, file_sfc2), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_sfc2 = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100) * 86400
hr_sfc2 = hr_sfc2.reshape((jtot, ktot))

with open(os.path.join(path, file_TOA2), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_TOA2 = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100) * 86400
hr_TOA2 = hr_TOA2.reshape((jtot, ktot))


diff_hr = hr2 - hr1
diff_hr_sfc = hr_sfc2 - hr_sfc1
diff_hr_TOA = hr_TOA2 - hr_TOA1



diff_img_sfc = pv.ImageData(dimensions=(ktot, jtot, 1))
diff_img_sfc['Values'] = diff_hr_sfc.ravel()
diff_img_TOA = pv.ImageData(dimensions=(ktot, jtot, 1))
diff_img_TOA['Values'] = diff_hr_TOA.ravel()


vmax_atm = np.max(np.abs(diff_hr))
vmax_sfc = np.max(np.abs(diff_hr_sfc))
print('Pesc - no Pesc')
vmax = np.min([vmax_atm, vmax_sfc])
vmax = 15

diff_hr = diff_hr[:50, :, :]


pl = pv.Plotter()
pl.add_volume(np.transpose(diff_hr, (2,1,0)), clim=[-vmax, vmax], cmap='seismic', opacity=[1,0,1])
pl.add_mesh(diff_img_sfc, clim=[-vmax, vmax], cmap='seismic', show_scalar_bar=False)
pl.show()


plt.figure(figsize = (7,7))
plt.imshow(diff_hr[:,25,:], origin='lower', cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar()
plt.show()
plt.close()

plt.figure(figsize = (20,5))
plt.hist(diff_hr.flatten(),bins=1000)
plt.show()
plt.close()


# hr_hist = diff_hr.flatten()
# print(np.sum(hr_hist > 0), np.sum(hr_hist < 0))

# plt.figure(figsize = (10,5))
# plt.hist(hr_hist, bins=1000)
# plt.yscale('log')
# plt.show()
# plt.close()


# hr1D = np.mean(hr, axis=(1,2))

# plt.figure()
# plt.plot(hr1D, 25 + 50*np.arange(len(hr1D)))
# plt.show()
# plt.close()

