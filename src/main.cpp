#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "MC.h"
#include "plane_parallel.h"


int main()
{

    // Control variables
    std::string CASE = "gpt3";                  // {gpt0, gpt1, gpt3, gpt21}
    bool ENABLE_MC = true;                      // Enables Monte Carlo algorithm
    bool ENABLE_PP = false;                     // Enables plane-parallel algorithm
    float dx = 1e8;                             // [m]
    float dy = 1e8;                             // [m]

    int Nphotpow = 21;                          // Total 2^pow number of photons
    std::string INTERCELL_TECHNIQUE = "power";  // {uniform, power, (power-gradient)}
    std::string INTRACELL_TECHNIQUE = "naive";  // {naive, (margin), (edge)}

    int N_mu = 20;                              // Number of angles to calculate at each height in PP-algorithm



    // Reading MVcases from json file
    std::fstream file_MVcases("/home/gijs-hogeboom/dev/lwproj/data_input/MVcases.json");
    if (!file_MVcases.is_open())
    {
        std::cerr << "Could not open MVcases.json!" << std::endl;
        return 1;
    }

    nlohmann::json j;
    file_MVcases >> j;

    std::vector<float> arr_z = j["z_" + CASE].get<std::vector<float>>();
    std::vector<float> arr_zh = j["zh_" + CASE].get<std::vector<float>>();
    std::vector<float> arr_dz = j["dz_" + CASE].get<std::vector<float>>();

    std::vector<float> arr_kext = j["kext_" + CASE].get<std::vector<float>>();
    std::vector<float> arr_Batm = j["Batm_" + CASE].get<std::vector<float>>();
    std::vector<float> arr_Batmh = j["Batmh_" + CASE].get<std::vector<float>>();

    float Bsfc = j["Bsfc_" + CASE].get<float>();


    // General initializations
    int Nphot = pow(2, Nphotpow);
    int Natm = Nphot / 2;
    int Nsfc = Nphot / 2;

    int itot = arr_z.size();
    int jtot = 1;
    int ktot = 1;

    float xtot = ktot * dx;
    float ytot = jtot * dy;
    float ztot = arr_zh[arr_zh.size() - 1];



    if (ENABLE_MC)
    {
        std::vector<float> heating_rates_MC = run_MC(arr_z,
                                                     arr_zh,
                                                     arr_dz,
                                                     Bsfc,
                                                     dx,
                                                     dy,
                                                     ktot,
                                                     jtot,
                                                     itot,
                                                     arr_kext,
                                                     arr_Batm,
                                                     CASE,
                                                     INTERCELL_TECHNIQUE,
                                                     INTRACELL_TECHNIQUE,
                                                     Natm,
                                                     Nsfc);
    }


    if (ENABLE_PP)
    {
        std::vector<float> heating_rates_PP = run_plane_parallel(arr_z,
                                                                arr_zh,
                                                                arr_dz,
                                                                arr_kext,
                                                                arr_Batm,
                                                                arr_Batmh,
                                                                CASE,
                                                                N_mu);
    }

    

    return 0;
}