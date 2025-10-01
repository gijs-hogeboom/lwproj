#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "MC.h"
#include "plane_parallel.h"


int main()
{

    // Control variables
    std::string CASE = "gpt3";

    int Nphotpow = 21;
    std::string INTERCELL_TECHNIQUE = "power";
    std::string INTRACELL_TECHNIQUE = "naive";

    int N_mu = 20;



    // Reading MVcases from json file
    std::fstream file_MVcases("../data_input/MVcases.json");
    if (!file_MVcases.is_open())
    {
        std::cerr << "Could not open MVcases.json!" << std::endl;
        return 1;
    }

    nlohmann::json j;
    file_MVcases >> j;

    std::vector<float> arr_z = j["z_" + CASE];
    std::vector<float> arr_zh = j["zh_" + CASE];
    std::vector<float> arr_dz = j["dz_" + CASE];






    // Initializing MC
    int Nphot = pow(2, Nphotpow);
    int Natm = Nphot / 2;
    int Nsfc = Nphot / 2;



    std::vector<float> arr_x = {1,2,3,4,5};
    std::vector<float> arr_y = {1,2,3,4,5};
    // std::vector<float> arr_z = {1,2,3,4,5};
    std::vector<float> arr_xh = {1,2,3,4,5};
    std::vector<float> arr_yh = {1,2,3,4,5};
    // std::vector<float> arr_zh = {1,2,3,4,5};
    std::vector<float> field_kext = {1,2,3,4,5};
    std::vector<float> arr_field_bounds= {1,2,3,4,5};
    std::vector<float> arr_photons_atm_pos_x = {1,2,3,4,5};
    std::vector<float> arr_photons_atm_pos_y = {1,2,3,4,5};
    std::vector<float> arr_photons_atm_pos_z = {1,2,3,4,5};
    std::vector<float> arr_photons_atm_mu = {1,2,3,4,5};
    std::vector<float> arr_photons_atm_az = {1,2,3,4,5};
    std::vector<float> arr_photons_atm_tau = {1,2,3,4,5};
    std::vector<float> arr_power_per_photon_atm = {1,2,3,4,5};
    std::vector<float> arr_photons_sfc_pos_x = {1,2,3,4,5};
    std::vector<float> arr_photons_sfc_pos_y = {1,2,3,4,5};
    std::vector<float> arr_photons_sfc_pos_z = {1,2,3,4,5};
    std::vector<float> arr_photons_sfc_mu = {1,2,3,4,5};
    std::vector<float> arr_photons_sfc_az = {1,2,3,4,5};
    std::vector<float> arr_photons_sfc_tau = {1,2,3,4,5};
    std::vector<float> arr_power_per_photon_sfc = {1,2,3,4,5};

    std::vector<float> arr_kext = {1,2,3,4,5};
    std::vector<float> arr_Batm = {1,2,3,4,5};
    std::vector<float> arr_Batmh = {1,2,3,4,5};


    std::vector<float> heating_rates_PP = run_plane_parallel(arr_z,
                                                            arr_zh,
                                                            arr_dz,
                                                            arr_kext,
                                                            arr_Batm,
                                                            arr_Batmh,
                                                            CASE,
                                                            N_mu);
    
    std::vector<float> heating_rates_MC = run_MC(arr_x,
                                                 arr_y,
                                                 arr_z,
                                                 arr_xh,
                                                 arr_yh,
                                                 arr_zh,
                                                 arr_dz,
                                                 field_kext,
                                                 arr_field_bounds,
                                                 arr_photons_atm_pos_x,
                                                 arr_photons_atm_pos_y,
                                                 arr_photons_atm_pos_z,
                                                 arr_photons_atm_mu,
                                                 arr_photons_atm_az,
                                                 arr_photons_atm_tau,
                                                 arr_power_per_photon_atm,
                                                 arr_photons_sfc_pos_x,
                                                 arr_photons_sfc_pos_y,
                                                 arr_photons_sfc_pos_z,
                                                 arr_photons_sfc_mu,
                                                 arr_photons_sfc_az,
                                                 arr_photons_sfc_tau,
                                                 arr_power_per_photon_sfc,
                                                 CASE,
                                                 INTERCELL_TECHNIQUE,
                                                 INTRACELL_TECHNIQUE,
                                                 Natm,
                                                 Nsfc);



    std::cout << "Heating rates:" << std::endl;
    std::cout << "PP:" << std::endl;
    for (auto el : heating_rates_PP)
    {
        std::cout << el << ',';
    }
    std::cout << std::endl;
    std::cout << "MC:" << std::endl;

    for (auto el : heating_rates_PP)
    {
        std::cout << el << ',';
    }
    std::cout << std::endl;
    return 0;
}