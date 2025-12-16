import numpy as np
import os
import matplotlib.pyplot as plt
import json
import netCDF4 as nc



# Single cloud scenario in boundary layer atmopshere

itot = 100
jtot = 49
ktot = jtot

dx = 100
dy = 100
dz = 50

cloud_center = jtot//2
cloud_r      = jtot//4

cloud_coords_x = [ x for x in range(cloud_center - cloud_r, cloud_center + cloud_r) for y in range(cloud_center - cloud_r, cloud_center + cloud_r) ]
cloud_coords_y = [ y for x in range(cloud_center - cloud_r, cloud_center + cloud_r) for y in range(cloud_center - cloud_r, cloud_center + cloud_r) ]


# Homogeneous column
def f_kext_open(z):
    a = -0.0005
    b = 0.01
    return 10**(a*z + b)/1000

def f_B_open(z):
    return 0.1

def f_SSA_open(z):
    return 0.

def f_ASY_open(z):
    return 0.


# Cloudy column
kext_cloud = 0.6
B_cloud = 2
SSA_cloud = 0.51
ASY_cloud = 0.75

cloud_height_lower = 1000
cloud_height_upper = 1400

def f_kext_cloud(z):
    if ((z >= cloud_height_lower) and (z <= cloud_height_upper)):
        return kext_cloud
    else:
        return f_kext_open(z)


def f_B_cloud(z):
    if ((z >= cloud_height_lower) and (z <= cloud_height_upper)):
        return B_cloud
    else:
        return f_B_open(z)


def f_SSA_cloud(z):
    if ((z >= cloud_height_lower) and (z <= cloud_height_upper)):
        return SSA_cloud
    else:
        return f_SSA_open(z)


def f_ASY_cloud(z):
    if ((z >= cloud_height_lower) and (z <= cloud_height_upper)):
        return ASY_cloud
    else:
        return f_ASY_open(z)
    

arr_zh = np.arange(0, itot*dz + dz, dz)
arr_z  = np.arange(dz/2, itot*dz + dz/2, dz)
arr_dz = np.array( [arr_zh[i+1] - arr_zh[i] for i in range(len(arr_z))] )


arr_kext_open  = np.array( [f_kext_open(z) for z in arr_z] )
arr_kext_cloud = np.array( [f_kext_cloud(z) for z in arr_z] )
arr_Batm_open  = np.array( [f_B_open(z) for z in arr_z] )
arr_Batm_cloud = np.array( [f_B_cloud(z) for z in arr_z] )
arr_SSA_open   = np.array( [f_SSA_open(z) for z in arr_z] )
arr_SSA_cloud  = np.array( [f_SSA_cloud(z) for z in arr_z] )
arr_ASY_open   = np.array( [f_ASY_open(z) for z in arr_z] )
arr_ASY_cloud  = np.array( [f_ASY_cloud(z) for z in arr_z] )

Bsfc = arr_Batm_open[0]*75


print(arr_kext_cloud)
print(arr_Batm_cloud)


dct_case = {
    'case_limits':[itot, jtot, ktot],
    'cell_size_horizontal':[dx, dy],
    'z':arr_z,
    'zh':arr_zh,
    'dz':arr_dz,
    'kext_open':arr_kext_open,
    'kext_cloud':arr_kext_cloud,
    'Batm_open':arr_Batm_open,
    'Batm_cloud':arr_Batm_cloud,
    'SSA_open':arr_SSA_open,
    'SSA_cloud':arr_SSA_cloud,
    'ASY_open':arr_ASY_open,
    'ASY_cloud':arr_ASY_cloud,
    'Bsfc':Bsfc,
    'cloud_coords_x':cloud_coords_x,
    'cloud_coords_y':cloud_coords_y
}

dct_out = {}
for key, val in dct_case.items():
    if not isinstance(val, list):
        val = val.tolist()

    dct_out[f's3D1_{key}'] = val

# Dumping case dictionary locally
with open("s3Dcases.json", 'w') as f:
    json.dump(dct_out, f)


# Dumping case dictionary at lwproj's data_input location
path = '/home/gijs-hogeboom/dev/mclw/data_input'
with open(os.path.join(path,"s3Dcases.json"), 'w') as f:
    json.dump(dct_out, f)






build_r3D_cases = True
if build_r3D_cases:
    # 3D Menno GPT cases

    path_MV3D = '/home/gijs-hogeboom/dev/Python_programmes/LWproject2'

    with nc.Dataset(os.path.join(path_MV3D, 'lw_optical_properties.nc'), 'r') as df:
        with nc.Dataset(os.path.join(path_MV3D, 'grid.nc'), 'r') as dfinfo:


            tau   = np.array(df.variables['lw_tau'][:])
            Batm  = np.array(df.variables['lay_source'][:])
            Bsfc  = np.array(df.variables['sfc_source'][:])
            ssa   = np.array(df.variables['lw_ssa'][:])
            asy   = np.array(df.variables['lw_asy'][:])

            arr_z  = np.array(dfinfo.variables['lay'][:])
            arr_zh = np.array(dfinfo.variables['lev'][:])
            arr_dz = np.array(dfinfo.variables['lev'][1:] - dfinfo.variables['lev'][:-1])


    kext = np.array(tau[:,:,:,:]) / arr_dz[np.newaxis, :, np.newaxis, np.newaxis]


    dct_cases   = {'z':arr_z.tolist(), 'zh':arr_zh.tolist(), 'dz':arr_dz.tolist()}

    ls_gpts = [3, 17, 21, 30]
    xy_size = 192


    dct_cases['xy_size'] = xy_size
    for gpt in ls_gpts:
        dct_cases[f'kext_r3D{gpt}'] = kext[gpt,:,:xy_size,:xy_size].flatten().tolist()
        dct_cases[f'Batm_r3D{gpt}'] = Batm[gpt,:,:xy_size,:xy_size].flatten().tolist()
        dct_cases[f'SSA_r3D{gpt}']  = ssa[gpt,:,:xy_size,:xy_size].flatten().tolist()
        dct_cases[f'ASY_r3D{gpt}']  = asy[gpt,:,:xy_size,:xy_size].flatten().tolist()
        dct_cases[f'Bsfc_r3D{gpt}'] = Bsfc[gpt,:xy_size,:xy_size].flatten().tolist()
        print(gpt, 'done')

    with open(os.path.join(path,"r3Dcases.json"), 'w') as f:
        json.dump(dct_cases, f)






build_MVcases = False
if build_MVcases:

    path = '/home/gijs-hogeboom/dev/Python_programmes/LWproject2'
    path_out = '/home/gijs-hogeboom/dev/mclw/data_input'


    with nc.Dataset(os.path.join(path, 'lw_optical_properties.nc'), 'r') as df:
        with nc.Dataset(os.path.join(path, 'grid.nc'), 'r') as dfinfo:

            tau    = np.array(df.variables['lw_tau'][:])
            Batm   = np.array(df.variables['lay_source'][:])
            Bsfc   = np.array(df.variables['sfc_source'][:])
            Batmh  = np.array(df.variables['lev_source'][:])

            arr_z  = np.array(dfinfo.variables['lay'][:])
            arr_zh = np.array(dfinfo.variables['lev'][:])
            arr_dz = np.array(dfinfo.variables['lev'][1:] - dfinfo.variables['lev'][:-1])

            ssa    = np.array(df.variables['lw_ssa'][:])
            asy    = np.array(df.variables['lw_asy'][:])

    kext = np.array(tau[:,:,:,:]) / arr_dz[np.newaxis, :, np.newaxis, np.newaxis]

    dct_cases   = {'z':arr_z, 'zh':arr_zh, 'dz':arr_dz}



    for gpt in range(36):
        
        print(gpt)
        
        arr_kext  = kext[gpt,:,:,:].mean(axis=(1,2))
        arr_Batm  = Batm[gpt,:,:,:].mean(axis=(1,2))
        arr_SSA   = ssa[gpt,:,:,:].mean(axis=(1,2))
        arr_ASY   = asy[gpt,:,:,:].mean(axis=(1,2))

        B_sfc     = Bsfc[gpt,:,:].mean()

        arr_phi   = 4*np.pi*arr_kext*arr_Batm*arr_dz
        
        dct_cases[f'z_gpt{gpt}']     = arr_z
        dct_cases[f'zh_gpt{gpt}']    = arr_zh
        dct_cases[f'dz_gpt{gpt}']    = arr_dz
        dct_cases[f'kext_gpt{gpt}']  = arr_kext
        dct_cases[f'Batm_gpt{gpt}']  = arr_Batm
        dct_cases[f'SSA_gpt{gpt}']   = arr_SSA
        dct_cases[f'ASY_gpt{gpt}']   = arr_ASY
        dct_cases[f'Bsfc_gpt{gpt}']  = B_sfc


    dct_jsondump = {k: v.tolist() for k,v in dct_cases.items()}


    with open(os.path.join(path_out, 'gptcases.json'), 'w') as f:
        json.dump(dct_jsondump, f)


