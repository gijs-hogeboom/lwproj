#include <iostream>
#include <cmath>
#include <vector>


std::vector<float> run_MC(const std::vector<float>& arr_z,
                          const std::vector<float>& arr_zh,
                          const std::vector<float>& arr_dz,
                          const std::vector<float>& arr_kext,
                          const std::vector<float>& arr_Batm,
                          float Bsfc,
                          float dx,
                          float dy,
                          int ktot,
                          int jtot,
                          int itot,
                          const std::string& CASE,
                          const std::string& INTERCELL_TECHNIQUE,
                          const std::string& INTRACELL_TECHNIQUE,
                          int Natm,
                          int Nsfc)
{

    // Initializing MC fields

    int n_volumes = itot * jtot * ktot;
    int n_tiles   = jtot * ktot;

    std::vector<float> arr_xh(ktot + 1);
    std::vector<float> arr_yh(jtot + 1);
    std::vector<float> arr_x(ktot);
    std::vector<float> arr_y(jtot);


    for (int j = 0; j < (jtot + 1); j++)
    {
        arr_yh[j] = j*dy;
        if (j < jtot)
        {
            arr_y[j] = (j*dy + (j+1)*dy)/2;
        }
    }
    for (int k = 0; k < (ktot + 1); k++)
    {
        arr_xh[k] = k*dx;
        if (k < ktot)
        {
            arr_x[k] = (k*dx + (k+1)*dx)/2;
        }
    }

    // Initializing fields
    std::vector<float> field_atm_kext(n_volumes);
    std::vector<float> field_atm_B(n_volumes);
    std::vector<float> field_atm_phi(n_volumes);
    std::vector<float> field_sfc_B(n_tiles);
    std::vector<float> field_sfc_phi(n_tiles);
    std::vector<float> field_sfc_eps(n_tiles, 1.0f);

    std::vector<float> field_atm_absorbed(n_volumes);
    std::vector<float> field_atm_emitted(n_volumes);
    std::vector<float> field_sfc_absorbed(n_tiles);
    std::vector<float> field_sfc_emitted(n_tiles);

    // Generating fields
    for (int i = 0; i < itot; i++)
    {
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm = i * ktot * jtot + j * ktot + k;
                field_atm_kext[idx_atm] = arr_kext[i];
                field_atm_B[idx_atm] = arr_Batm[i];
                field_atm_phi[idx_atm] = 4*

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    field_sfc_B[idx_sfc] = Bsfc;
                }
            }
        }
    }










    std::vector<float> arr_heating_rates(1, 0.0f);

    return arr_heating_rates;
}