import numpy as np
import json
import matplotlib.pyplot as plt
import pyvista as pv




CASE = 'r3D21'


# reading and unpacking data
with open(f"/home/gijs-hogeboom/dev/mclw/data_input/{CASE[:3]}cases.json", 'r') as f:
    dct = json.load(f)



if CASE[:3] == 's3D':

    itot, jtot, ktot = dct[f'{CASE}_case_limits']
    dx, dy = dct[f'{CASE}_cell_size_horizontal']

    arr_z  = dct[f'{CASE}_z']
    arr_zh = dct[f'{CASE}_zh']
    arr_dz = dct[f'{CASE}_dz']

    col_kext_hom   = dct[f'{CASE}_kext_open']
    col_kext_cloud = dct[f'{CASE}_kext_cloud']
    col_B_hom      = dct[f'{CASE}_Batm_open']
    col_B_cloud    = dct[f'{CASE}_Batm_cloud']

    cloud_coords_x = dct[f'{CASE}_cloud_coords_x']
    cloud_coords_y = dct[f'{CASE}_cloud_coords_y']

    # Deducing homogenized column coordinates from cloud coordinates
    hom_coords_x = np.zeros(jtot*ktot - len(cloud_coords_x), dtype=int)
    hom_coords_y = np.zeros(jtot*ktot - len(cloud_coords_y), dtype=int)
    hom_idx = 0
    for j in range(jtot):
        for k in range(ktot):

            cloudy_column = False
            for cloud_x, cloud_y in zip(cloud_coords_x, cloud_coords_y):
                if ((j == cloud_y) and (k == cloud_x)):
                    cloudy_column = True
            
            if not cloudy_column:
                hom_coords_x[hom_idx] = k
                hom_coords_y[hom_idx] = j
                hom_idx += 1


    # Initializing fields
    arr_kext = np.zeros(itot*jtot*ktot)
    arr_B    = np.zeros(itot*jtot*ktot)

    # Filling in homogenized columns
    for cloud_x, cloud_y in zip(cloud_coords_x, cloud_coords_y):
        for i in range(itot):

            j = cloud_y
            k = cloud_x
            idx = i*jtot*ktot + j*ktot + k

            arr_kext[idx] = col_kext_cloud[i]
            arr_B[idx]    = col_B_cloud[i]

    # Filling in cloudy columns
    for hom_x, hom_y in zip(hom_coords_x, hom_coords_y):
        for i in range(itot):

            j = hom_y
            k = hom_x
            idx = i*jtot*ktot + j*ktot + k

            arr_kext[idx] = col_kext_cloud[i]
            arr_B[idx]    = col_B_cloud[i]


    for i in range(itot):
        for j in range(jtot):
            for k in range(ktot):
                
                idx = i*jtot*ktot + j*ktot + k

                # Check whether current column is a cloudy one or non-cloudy one
                column_is_cloud = False
                for cloud_x, cloud_y in zip(cloud_coords_x, cloud_coords_y):
                    if ((cloud_x == k) and (cloud_y == j)):
                        column_is_cloud = True
                
                # Filling in arrays
                if column_is_cloud:
                    arr_kext[idx] = col_kext_cloud[i]
                    arr_B[idx]    = col_B_cloud[i]
                else:
                    arr_kext[idx] = col_kext_hom[i]
                    arr_B[idx]    = col_B_hom[i]


    field_kext = arr_kext.reshape((itot, jtot, ktot))
    field_B    = arr_B.reshape((itot, jtot, ktot))



    pl = pv.Plotter()
    pl.add_volume(field_kext, cmap='jet')
    pl.show()


elif CASE[:3] == 'r3D':

    arr_z  = np.array(dct[f'z'])
    arr_zh = np.array(dct[f'zh'])
    arr_dz = np.array(dct[f'dz'])

    itot = len(arr_z)
    jtot = ktot = dct['xy_size']
    
    field_atm_kext = np.array(dct[f'kext_{CASE}']).reshape((itot, jtot, ktot))
    field_atm_B = np.array(dct[f'Batm_{CASE}']).reshape((itot, jtot, ktot))
    field_sfc_B = np.array(dct[f'Bsfc_{CASE}']).reshape((jtot, ktot))

    field_atm_phi = 4*np.pi * field_atm_kext * field_atm_B * 100 * 100
    field_sfc_phi =   np.pi * 1.0 * field_sfc_B * 100 * 100

    img = pv.ImageData(dimensions=(1,jtot, ktot))
    img['values'] = field_sfc_phi.ravel()

    my_cmap = 'jet'

    pl = pv.Plotter()
    pl.add_volume(field_atm_phi, cmap=my_cmap)#, opacity=[0,0,1])
    pl.add_mesh(img, cmap=my_cmap)
    pl.show()

