#include <iostream>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include <chrono>
#include "gnuplot-iostream.h"

#include "MC.h"
#include "plane_parallel.h"
#include "util.h"


int main()
{

    bool constexpr write_EB = false;

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;


    std::ofstream foutEB("/home/gijs-hogeboom/dev/lwproj/data_output/energy_balance_data_low.dat");
    if (!foutEB.is_open())
    {
        std::cerr << "Error, cannot open foutEB" << std::endl;
        return 1;
    }


    std::vector<int> arr_photonpows = {24};


    // Control variables
    std::string CASE = "gpt21";                   // {gpt0, gpt1, gpt3, gpt21}
    bool ENABLE_MC = false;                        // Enables Monte Carlo algorithm
    bool ENABLE_PP = true;                        // Enables plane-parallel algorithm
    double dx = 1e4;                              // [m]
    double dy = 1e4;                              // [m]

    std::string INTERCELL_TECHNIQUE = "power";    // {uniform, power, power-gradient}
    std::string INTRACELL_TECHNIQUE = "naive";    // {naive, (margin), (edge)}

    int N_mu = 100;                                // Number of angles to calculate at each height in PP-algorithm



    // Console output params
    bool print_EB_MC         = false;
    bool print_EB_PP         = true;
    bool verbose             = true;
    bool print_final_results = true;
    bool print_counter       = false;
    bool plot_results        = false;


    if (verbose)
    {
        std::cout << "Start of program, reading input data" << std::endl;
    }

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



    int run_counter = 1;

    for (auto n_phot_pow : arr_photonpows)
    {



        // General initializations
        // int dNphot    = std::max({0, n_phot_pow - 20});
        // int Niter     = pow(2, dNphot);                                // Number of experiments for the RMSE/ME metrics

        // int Nphot_pow = n_phot_pow - dNphot;                            // Total 2^pow number of photons

        int Niter     = 1;
        int Nphot_pow = n_phot_pow;

        int Nphot     = pow(2, Nphot_pow);
        int Natm_pow  = Nphot_pow - 1; // splitting total amount in half
        int Nsfc_pow  = Nphot_pow - 1;
        int Natm      = pow(2, Natm_pow);
        int Nsfc      = pow(2, Nsfc_pow);


        int itot      = arr_z.size();
        int jtot      = 1;
        int ktot      = 1;

        std::string fnameMetrics = "metrics_" + CASE + "_" + INTERCELL_TECHNIQUE + "_Natm" + std::to_string(Natm_pow) + "_Nsfc" + std::to_string(Nsfc_pow);
        std::ofstream foutMetrics("/home/gijs-hogeboom/dev/lwproj/data_output/metrics_profiles/" + fnameMetrics + ".csv");
        if (!foutMetrics.is_open())
        {
            std::cerr << "Error, cannot open foutMetrics!" << std::endl;
            return 1;
        }

        // Header
        foutMetrics << "z,RMSE,ME" << std::endl;

        std::vector<double> arr_RMSE(itot);
        std::vector<double> arr_ME(itot);

        for (int _iter = 0; _iter < Niter; _iter++)
        {
            
            std::vector<double> heating_rates_MC(itot);
            std::vector<double> heating_rates_PP(itot);

            auto MC_t1 = std::chrono::high_resolution_clock::now();
            if (ENABLE_MC)
            {
                if (verbose)
                {
                    std::cout << "Start of MC - Natm: " + std::to_string(Natm_pow) + ", Nsfc: " + std::to_string(Nsfc_pow) << std::endl;
                }

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
                                        INTERCELL_TECHNIQUE, 
                                        INTRACELL_TECHNIQUE,
                                        CASE,
                                        Natm, 
                                        Nsfc,
                                        print_EB_MC,
                                        verbose);
                
                // Storing output
                std::ofstream file_MCoutput("/home/gijs-hogeboom/dev/lwproj/data_output/heating_rates/HR_MC_" + CASE + "_" + INTERCELL_TECHNIQUE + "_" + INTRACELL_TECHNIQUE + "_Natm" + std::to_string(Natm_pow) + "_Nsfc" + std::to_string(Nsfc_pow) + ".csv");
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
            auto MC_t2 = std::chrono::high_resolution_clock::now();

            auto PP_t1 = std::chrono::high_resolution_clock::now();
            if (ENABLE_PP)
            {

                if (verbose)
                {
                    std::cout << "Start of PP" << std::endl;
                }

                heating_rates_PP = run_plane_parallel(arr_z,
                                                    arr_zh,
                                                    arr_dz,
                                                    arr_kext,
                                                    arr_Batm,
                                                    arr_Batmh,
                                                    Bsfc,
                                                    N_mu,
                                                    print_EB_PP,
                                                    verbose);
                
                // Storing output
                
                std::ofstream file_PPoutput("/home/gijs-hogeboom/dev/lwproj/data_output/heating_rates/HR_PP_" + CASE + ".csv");
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
            auto PP_t2 = std::chrono::high_resolution_clock::now();


            duration<double, std::milli> MC_time = MC_t2 - MC_t1;
            duration<double, std::milli> PP_time = PP_t2 - PP_t1;
            
            // Loading results
            std::vector<double> heating_rates_PP_in(itot, 0.);
            std::vector<double> heating_rates_MC_in(itot, 0.);
            std::vector<double> arr_z_in(itot, 0.);

            std::fstream file_MCinput("/home/gijs-hogeboom/dev/lwproj/data_output/heating_rates/HR_MC_" + CASE + "_" + INTERCELL_TECHNIQUE + "_" + INTRACELL_TECHNIQUE + "_Natm" + std::to_string(Natm_pow) + "_Nsfc" + std::to_string(Nsfc_pow) + ".csv");
            std::fstream file_PPinput("/home/gijs-hogeboom/dev/lwproj/data_output/heating_rates/HR_PP_" + CASE + ".csv");
            
            if (!file_MCinput.is_open())
            {
                std::cout << "File |HR_MC_" + CASE + "_" + INTERCELL_TECHNIQUE + "_" + INTRACELL_TECHNIQUE + "_Natm" + std::to_string(Natm_pow) + "_Nsfc" + std::to_string(Nsfc_pow) + ".csv| does not exist!" << std::endl;
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


            ////////////// Calculating metrics ///////////////////////
            double RMSE_tot, ME_tot;

            for (size_t i = 0; i < itot; i++)
            {
                RMSE_tot += (1./itot)*pow(heating_rates_MC_in[i] - heating_rates_PP_in[i], 2);
                ME_tot   += (1./itot)*(heating_rates_MC_in[i] - heating_rates_PP_in[i]);

                arr_RMSE[i] += (1./Niter)*pow(heating_rates_MC_in[i] - heating_rates_PP_in[i], 2);
                arr_ME[i]   += (1./Niter)*(heating_rates_MC_in[i] - heating_rates_PP_in[i]);

            }

            RMSE_tot = std::sqrt(RMSE_tot);

            if (print_counter)
            {
                std::cout << run_counter << " / " << (Niter * arr_photonpows.size()) << std::endl;
                run_counter++;
            }

            ////////////// Plotting results //////////////////////////
            if (plot_results)
            {
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
            }




            ////////////// Printing results //////////////////////////
            if (print_final_results)
            {
                std::cout << "====================================" << std::endl;
                std::cout << "                Done!" << std::endl;
                std::cout << "Parameters:" << std::endl;
                std::cout << "   | Natm:        " << Natm_pow << std::endl;
                std::cout << "   | Nsfc:        " << Nsfc_pow << std::endl;
                std::cout << "   | N_Mu:        " << N_mu << std::endl;
                std::cout << "Timers:" << std::endl;
                std::cout << "   | MC time:     " << MC_time.count()/1000 << std::endl;
                std::cout << "   | PP time:     " << PP_time.count()/1000 << std::endl;
                std::cout << "Metrics:" << std::endl;
                std::cout << "   | RMSE:        " << RMSE_tot << std::endl;
                std::cout << "   | ME:          " << ME_tot << std::endl;
                std::cout << "====================================" << std::endl;
            }
        }
        
        
        
        for (size_t i = 0; i < itot; i++)
        {
            arr_RMSE[i] = std::sqrt(arr_RMSE[i]);
            foutMetrics << arr_z[i] << ',' << arr_RMSE[i] << ',' << arr_ME[i] << std::endl;
        }

        foutMetrics.close();

    }


    foutEB.close();


    return 0;
}