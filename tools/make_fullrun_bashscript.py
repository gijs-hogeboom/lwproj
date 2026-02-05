import numpy as np
import pandas as pd


df = pd.read_csv('../data_output/data_statistics.csv')


Nphot_pow = 35
Nphot = 2**Nphot_pow

df['phi tot mean'] = df['phi mean'] + df['Bsfc mean']*np.pi*100*100

fraction_photons = df['phi tot mean']/np.sum(df['phi tot mean'])
ls_photons = []

for frac in fraction_photons:
    ls_photons.append(f"{np.log2(frac * Nphot):.2f}")


with open(f'../tools/fullrun_Nphot{Nphot_pow:d}.sh', 'w') as f:
    for gpt, Nphot in enumerate(ls_photons):
        f.write(f'../build/lwproj "r3D{gpt}" "power" ' + Nphot + ' 1 0\n')
