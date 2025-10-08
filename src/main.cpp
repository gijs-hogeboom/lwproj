#include <iostream>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include "gnuplot-iostream.h"

#include "MC.h"
#include "plane_parallel.h"


int main()
{
    
    std::cout << "Start of program, reading input data" << std::endl;

    // Control variables
    std::string CASE = "gpt21";                   // {gpt0, gpt1, gpt3, gpt21}
    bool ENABLE_MC = true;                        // Enables Monte Carlo algorithm
    bool ENABLE_PP = true;                        // Enables plane-parallel algorithm
    double dx = 1e8;                               // [m]
    double dy = 1e8;                               // [m]

    int Nphot_pow = 25;                            // Total 2^pow number of photons
    std::string INTERCELL_TECHNIQUE = "power";  // {uniform, power, (power-gradient)}
    std::string INTRACELL_TECHNIQUE = "naive";    // {naive, (margin), (edge)}

    int N_mu = 30;                                // Number of angles to calculate at each height in PP-algorithm



    // Reading MVcases from json file
    std::fstream file_MVcases("/home/gijs-hogeboom/dev/lwproj/data_input/MVcases.json");
    if (!file_MVcases.is_open())
    {
        std::cerr << "Could not open MVcases.json!" << std::endl;
        return 1;
    }

    nlohmann::json j;
    file_MVcases >> j;

    std::vector<double> arr_z = j["z_" + CASE].get<std::vector<double>>();
    std::vector<double> arr_zh = j["zh_" + CASE].get<std::vector<double>>();
    std::vector<double> arr_dz = j["dz_" + CASE].get<std::vector<double>>();

    std::vector<double> arr_kext = j["kext_" + CASE].get<std::vector<double>>();
    std::vector<double> arr_Batm = j["Batm_" + CASE].get<std::vector<double>>();
    std::vector<double> arr_Batmh = j["Batmh_" + CASE].get<std::vector<double>>();

    double Bsfc = j["Bsfc_" + CASE].get<double>();


    // General initializations
    int Nphot    = pow(2, Nphot_pow);
    int Natm_pow = Nphot_pow - 1; // splitting total amount in half
    int Nsfc_pow = Nphot_pow - 1;
    int Natm     = pow(2, Natm_pow);
    int Nsfc     = pow(2, Nsfc_pow);

    int itot     = arr_z.size();
    int jtot     = 1;
    int ktot     = 1;
    

    std::vector<double> heating_rates_MC(itot);
    std::vector<double> heating_rates_PP(itot);

    if (ENABLE_MC)
    {
        std::cout << "Start of MC" << std::endl;

        heating_rates_MC = run_MC(arr_z,
                                arr_zh,
                                arr_dz,
                                arr_kext,
                                arr_Batm,
                                Bsfc,
                                dx,
                                dy,
                                ktot,
                                jtot,
                                itot,
                                CASE,
                                INTERCELL_TECHNIQUE, 
                                INTRACELL_TECHNIQUE, 
                                Natm, 
                                Nsfc);
        
        // Storing output
        std::ofstream file_MCoutput("/home/gijs-hogeboom/dev/lwproj/data_output/HR_MC_" + CASE + "_Natm" + std::to_string(Natm_pow) + "_Nsfc" + std::to_string(Nsfc_pow) + ".csv");
        if (!file_MCoutput.is_open())
        {
            std::cerr << "Error: cannot open MC output file!" << std::endl;
            return 1;
        }

        file_MCoutput << "z,heating_rate\n"; // header
        for (size_t i = 0; i < itot; i++)
        {
            file_MCoutput << arr_z[i] << ',' << heating_rates_MC[i] << std::endl;
        }
    }

    if (ENABLE_PP)
    {

        std::cout << "Start of PP" << std::endl;

        heating_rates_PP = run_plane_parallel(arr_z,
                                            arr_zh,
                                            arr_dz,
                                            arr_kext,
                                            arr_Batm,
                                            arr_Batmh,
                                            CASE,
                                            Bsfc,
                                            N_mu);
        
        // Storing output
        
        std::ofstream file_PPoutput("/home/gijs-hogeboom/dev/lwproj/data_output/HR_PP_" + CASE + ".csv");
        if (!file_PPoutput.is_open())
        {
            std::cerr << "Error: cannot open PP output file!" << std::endl;
            return 1;
        }

        file_PPoutput << "z,heating_rate\n"; // header
        for (size_t i = 0; i < itot; i++)
        {
            file_PPoutput << arr_z[i] << ',' << heating_rates_PP[i] << std::endl;
        }
    }




    
    // Loading results
    std::vector<double> heating_rates_PP_in(itot, 0.);
    std::vector<double> heating_rates_MC_in(itot, 0.);
    std::vector<double> arr_z_in(itot, 0.);

    std::fstream file_MCinput("/home/gijs-hogeboom/dev/lwproj/data_output/HR_MC_" + CASE + "_Natm" + std::to_string(Natm_pow) + "_Nsfc" + std::to_string(Nsfc_pow) + ".csv");
    std::fstream file_PPinput("/home/gijs-hogeboom/dev/lwproj/data_output/HR_PP_" + CASE + ".csv");
    
    if (!file_MCinput.is_open())
    {
        std::cout << "File |HR_MC_" + CASE + "_Natm" + std::to_string(Natm_pow) + "_Nsfc" + std::to_string(Nsfc_pow) + ".csv| does not exist!" << std::endl;
    }
    else
    {
        size_t i = 0;
        std::string line;
        std::getline(file_MCinput, line); // Skipping header
        while (std::getline(file_MCinput, line))
        {
            std::stringstream ss(line);
            std::string cell;

            std::getline(ss, cell, ',');
            arr_z_in[i] = std::stod(cell);

            std::getline(ss, cell, ',');
            heating_rates_MC_in[i] = std::stod(cell);
            std::cout << heating_rates_MC_in[i] << std::endl;
            i++;
        }
        file_MCinput.close();
    }

    if (!file_PPinput.is_open())
    {
        std::cout << "File |HR_PP_" + CASE + ".csv| does not exist!" << std::endl;
    }
    else
    {
        size_t i = 0;
        std::string line;
        std::getline(file_PPinput, line); // Skipping header
        while (std::getline(file_PPinput, line))
        {
            std::stringstream ss(line);
            std::string cell;

            std::getline(ss, cell, ',');
            arr_z_in[i] = std::stod(cell);

            std::getline(ss, cell, ',');
            heating_rates_PP_in[i] = std::stod(cell);

            i++;
        }
        file_PPinput.close();
    }

    // Plotting results
    Gnuplot gp;

    std::vector<std::pair<double,double>> hr_1D, hr_3D;
    for (int i = 0; i < itot; i++)
    {
        hr_1D.emplace_back(heating_rates_PP_in[i], arr_z_in[i]);
        hr_3D.emplace_back(heating_rates_MC_in[i], arr_z_in[i]);
    }


    
    gp << "set yrange [0 : 5000]\n"
       << "plot '-' with lines title 'PP', '-' with lines title 'MC'\n";
    gp.send1d(hr_1D);
    gp.send1d(hr_3D);
    std::cout << "Done" << std::endl;
    return 0;
}