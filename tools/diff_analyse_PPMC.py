import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pyvista as pv

path_MC = '/home/gijs-hogeboom/dev/mclw/data_output/raw_output_3D'
path_PP = '/home/gijs-hogeboom/dev/mclw/data_output/raw_output_3DPP'


Case = 'r3D21'
InterCellTechnique = 'power'
Nphot = '28.00'
PescMode = 0
enable_scattering = 0

file_atm_MC = f'hr_3D_atm_{Case}_Nphot{Nphot}_{InterCellTechnique}_Pesc{PescMode}_scatter{enable_scattering}.dat'
file_sfc_MC = f'flux_3D_sfc_{Case}_Nphot{Nphot}_{InterCellTechnique}_Pesc{PescMode}_scatter{enable_scattering}.dat'
file_TOA_MC = f'flux_3D_TOA_{Case}_Nphot{Nphot}_{InterCellTechnique}_Pesc{PescMode}_scatter{enable_scattering}.dat'

file_atm_PP = f'hr_3DPP_atm_{Case}.dat'
file_sfc_PP = f'flux_3DPP_sfc_{Case}.dat'
file_TOA_PP = f'flux_3DPP_TOA_{Case}.dat'

with open(os.path.join(path_MC, file_atm_MC), 'rb') as f:
    itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)
    hr1 = np.fromfile(f, dtype=np.float32)
hr_MC = hr1.reshape((itot, jtot, ktot))

with open(os.path.join(path_MC, file_sfc_MC), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_sfc1 = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100 * 100) * 86400
hr_sfc_MC = hr_sfc1.reshape((jtot, ktot))

with open(os.path.join(path_MC, file_TOA_MC), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_TOA1 = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100) * 86400
hr_TOA_MC = hr_TOA1.reshape((jtot, ktot))


with open(os.path.join(path_PP, file_atm_PP), 'rb') as f:
    itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)
    hr2 = np.fromfile(f, dtype=np.float32)
hr_PP = hr2.reshape((itot, jtot, ktot))

with open(os.path.join(path_PP, file_sfc_PP), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_sfc2 = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100 * 100) * 86400
hr_sfc_PP = hr_sfc2.reshape((jtot, ktot))

with open(os.path.join(path_PP, file_TOA_PP), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_TOA2 = np.fromfile(f, dtype=np.float32) / (1.225 * 1004 * 100 * 100) * 86400
hr_TOA_PP = hr_TOA2.reshape((jtot, ktot))



diff_hr = hr_PP - hr_MC
diff_hr_sfc = hr_sfc_PP - hr_sfc_MC
diff_hr_TOA = hr_TOA_PP - hr_TOA_MC


img_sfc_MC = pv.ImageData(dimensions=(ktot, jtot, 1))
img_sfc_MC['Values'] = hr_sfc_MC.ravel()
img_TOA_MC = pv.ImageData(dimensions=(ktot, jtot, 1))
img_TOA_MC['Values'] = hr_TOA_MC.ravel()

img_sfc_PP = pv.ImageData(dimensions=(ktot, jtot, 1))
img_sfc_PP['Values'] = hr_sfc_PP.ravel()
img_TOA_PP = pv.ImageData(dimensions=(ktot, jtot, 1))
img_TOA_PP['Values'] = hr_TOA_PP.ravel()

diff_img_sfc = pv.ImageData(dimensions=(ktot, jtot, 1))
diff_img_sfc['Values'] = diff_hr_sfc.ravel()
diff_img_TOA = pv.ImageData(dimensions=(ktot, jtot, 1))
diff_img_TOA['Values'] = diff_hr_TOA.ravel()


vmax_atm = np.max(np.abs(hr_PP))
vmax_sfc = np.max(np.abs(hr_sfc_PP))
print('PP - MC')
vmax = np.min([vmax_atm, vmax_sfc])


if Case[0] == 's':
    z_upper = 50
elif Case[0] == 'r':
    z_upper = 100
    
hr_MC   = hr_MC[:z_upper, :, :]
hr_PP   = hr_PP[:z_upper, :, :]
diff_hr = diff_hr[:z_upper, :, :]


pl = pv.Plotter()
pl.add_volume(np.transpose(hr_MC, (2,1,0)), clim=[-vmax_atm, vmax_atm], cmap='seismic', opacity=[1,0,1])
pl.add_mesh(img_sfc_MC, clim=[-vmax_sfc, vmax_sfc], cmap='seismic')
pl.show()

pl = pv.Plotter()
pl.add_volume(np.transpose(hr_PP, (2,1,0)), clim=[-vmax_atm, vmax_atm], cmap='seismic', opacity=[1, 0, 1])
pl.add_mesh(img_sfc_PP, clim=[-vmax_sfc, vmax_sfc], cmap='seismic')
pl.show()

pl = pv.Plotter()
pl.add_volume(np.transpose(diff_hr, (2,1,0)), clim=[-vmax, vmax], cmap='seismic', opacity=[1,0,1])
pl.add_mesh(diff_img_sfc, clim=[-vmax, vmax], cmap='seismic')
pl.show()

if Case[0] == 's':
    y_disect = 25
elif Case[0] == 'r':
    y_disect = 100

plt.figure(figsize = (20,7))
plt.subplot(1,3,1)
plt.imshow(hr_MC[:,y_disect,:], origin='lower', cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar()
plt.title('MC')

plt.subplot(1,3,2)
plt.imshow(hr_PP[:,y_disect,:], origin='lower', cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar()
plt.title('PP')

plt.subplot(1,3,3)
plt.imshow(diff_hr[:,y_disect,:], origin='lower', cmap='seismic', vmin=-vmax, vmax=vmax)
plt.colorbar()
plt.title('PP - MC')

plt.show()
plt.close()


