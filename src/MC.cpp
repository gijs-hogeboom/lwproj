#include <iostream>
#include <cmath>
#include <vector>


std::vector<float> run_MC(const std::vector<float>& arr_x,
                          const std::vector<float>& arr_y,
                          const std::vector<float>& arr_z,
                          const std::vector<float>& arr_xh,
                          const std::vector<float>& arr_yh,
                          const std::vector<float>& arr_zh,
                          const std::vector<float>& arr_dz,
                          const std::vector<float>& field_kext,
                          const std::vector<float>& arr_field_bounds,
                          const std::vector<float>& arr_photons_atm_pos_x,
                          const std::vector<float>& arr_photons_atm_pos_y,
                          const std::vector<float>& arr_photons_atm_pos_z,
                          const std::vector<float>& arr_photons_atm_mu,
                          const std::vector<float>& arr_photons_atm_az,
                          const std::vector<float>& arr_photons_atm_tau,
                          const std::vector<float>& arr_power_per_photon_atm,
                          const std::vector<float>& arr_photons_sfc_pos_x,
                          const std::vector<float>& arr_photons_sfc_pos_y,
                          const std::vector<float>& arr_photons_sfc_pos_z,
                          const std::vector<float>& arr_photons_sfc_mu,
                          const std::vector<float>& arr_photons_sfc_az,
                          const std::vector<float>& arr_photons_sfc_tau,
                          const std::vector<float>& arr_power_per_photon_sfc,
                          const std::string& CASE,
                          const std::string& INTERCELL_TECHNIQUE,
                          const std::string& INTRACELL_TECHNIQUE,
                          int Natm,
                          int Nsfc)
{

    std::vector<float> arr_heating_rates;
    arr_heating_rates.reserve(arr_z.size());

    for (auto el : arr_z)
    {
        arr_heating_rates.emplace_back(1.0 * el);
    }

    return arr_heating_rates;
}