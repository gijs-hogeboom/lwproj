import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pyvista as pv

path = '/home/gijs-hogeboom/dev/mclw/data_output/raw_output_3D'


Case = 'r3D21'
InterCellTechnique = 'power'
Natm = Nsfc = 21
PescMode = 1

file_atm = f'hr_3D_atm_{Case}_Natm{Natm}_Nsfc{Nsfc}_{InterCellTechnique}_Pesc{PescMode}.dat'
file_sfc = f'flux_3D_sfc_{Case}_Natm{Natm}_Nsfc{Nsfc}_{InterCellTechnique}_Pesc{PescMode}.dat'

with open(os.path.join(path, file_atm), 'rb') as f:
    itot, jtot, ktot = np.fromfile(f, dtype=np.int32, count=3)
    hr = np.fromfile(f, dtype=np.float64)
hr = hr.reshape((itot, jtot, ktot))

with open(os.path.join(path, file_sfc), 'rb') as f:
    jtot, ktot = np.fromfile(f, dtype=np.int32, count=2)
    hr_sfc = np.fromfile(f, dtype=np.float64) / (1.225 * 1004 * 100 * 100) * 86400
hr_sfc = hr_sfc.reshape((jtot, ktot))

img = pv.ImageData(dimensions=(1, jtot, ktot))
img['Values'] = hr_sfc.ravel()



# hr = hr[:50,:,:]

vmax_atm = np.max(np.abs(hr))
vmax_sfc = np.max(np.abs(hr_sfc))
vmax = np.min([vmax_atm, vmax_sfc])
# vmax = 116

pl = pv.Plotter()
pl.add_volume(hr, clim=[-vmax, vmax], cmap='seismic', opacity=[1,0,1])
pl.add_mesh(img, clim=[-vmax, vmax], cmap='seismic')
pl.show()


# hr1D = np.mean(hr, axis=(1,2))

# plt.figure()
# plt.plot(hr1D, 25 + 50*np.arange(len(hr1D)))
# plt.show()
# plt.close()

