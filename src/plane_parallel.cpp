#include <iostream>
#include <cmath>
#include <vector>
#include <string>


std::vector<float> run_plane_parallel(const std::vector<float>& arr_z,
                                      const std::vector<float>& arr_zh,
                                      const std::vector<float>& arr_dz,
                                      const std::vector<float>& arr_kext,
                                      const std::vector<float>& arr_Batm,
                                      const std::vector<float>& arr_Batmh,
                                      const std::string& CASE,
                                      int N_mu)
{

    std::vector<float> arr_heating_rates;
    arr_heating_rates.reserve(arr_z.size());

    for (auto el : arr_z)
    {
        arr_heating_rates.emplace_back(1.0 * el);
    }

    return arr_heating_rates;
}