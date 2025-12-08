#include <iostream>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include <chrono>
#include "gnuplot-iostream.h"

#include "MC.h"
#include "plane_parallel.h"
#include "util.h"


int main(int argc, char* argv[])
{

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    


    // Handling input
    std::string arg1 = "gpt21";
    std::string arg2 = "power";
    int arg3 = 21;
    bool arg4 = false;


    if (argc > 1)
    {
        arg1 = argv[1];
    }
    if (argc > 2)
    {
        arg2 = argv[2];
    }
    if (argc > 3)
    {
        arg3 = std::stoi(argv[3]);
    }
    if (argc > 4)
    {
        std::istringstream(argv[4]) >> arg4;
    }

    std::cout << "Arguments | CASE: " << arg1 << ", INTERCELL_TECHNIQUE: " << arg2 << ", Nphot_pow: " << arg3 << ", Pesc_mode: " << arg4 << std::endl;



    // Run settings
    int Nphot_pow = arg3;

    int Nphot     = pow(2, Nphot_pow);
    int Natm_pow  = Nphot_pow - 1; // splitting total amount in half
    int Nsfc_pow  = Nphot_pow - 1;
    int Natm      = pow(2, Natm_pow);
    int Nsfc      = pow(2, Nsfc_pow);

    std::string CASE = arg1;
    bool ENABLE_MC = true;                        // Enables Monte Carlo algorithm
    bool ENABLE_PP = true;                        // Enables plane-parallel algorithm

    std::string INTERCELL_TECHNIQUE = arg2;       // {uniform, power, power-gradient}
    std::string INTRACELL_TECHNIQUE = "naive";    // {naive, (margin), (edge)}

    bool Pesc_mode = arg4;

    bool enable_full_counter_matrix = false;

    std::string OUTPUT_MODE = "3D";               // {1D, 3D}



    // Console output params
    bool print_EB_MC         = true;
    bool print_EB_PP         = true;
    bool verbose             = true;
    bool print_final_results = true;
    bool plot_results        = true;


    /////////////// READING CASE //////////////////

    if (verbose)
    {
        std::cout << "Start of program, reading input data" << std::endl;
    }
    
    // Initializing variables to load
    std::vector<double> arr_z;
    std::vector<double> arr_zh;
    std::vector<double> arr_dz;
    std::vector<double> field_atm_kext;
    std::vector<double> field_atm_B;
    std::vector<double> field_sfc_B;
    double Bsfc, dx, dy;
    int itot, jtot, ktot;

    // Reading the case json file
    std::string caseID = CASE.substr(0, 3);
    std::fstream file_cases("/home/gijs-hogeboom/dev/mclw/data_input/" + caseID + "cases.json");
    
    if (!file_cases.is_open())
    {
        std::cerr << "Could not open " + caseID + "cases.json!" << std::endl;
        return 1;
    }
    
    nlohmann::json j;
    file_cases >> j;

    // File input handling
    if (caseID == "gpt") // MV cases
    {
        arr_z = j["z_" + CASE].get<std::vector<double>>();
        arr_zh = j["zh_" + CASE].get<std::vector<double>>();
        arr_dz = j["dz_" + CASE].get<std::vector<double>>();

        // jtot, ktot, dx and dy are all predetermined in this case
        itot = arr_z.size();
        jtot = 11;
        ktot = 11;
        int n_volumes = itot * jtot * ktot;
        int n_tiles = jtot * ktot;

        dx = 100;
        dy = 100;

        // Generating fields from 1D kext and Batm data
        std::vector<double> arr_kext = j["kext_" + CASE].get<std::vector<double>>();
        std::vector<double> arr_Batm = j["Batm_" + CASE].get<std::vector<double>>();
        Bsfc = j["Bsfc_" + CASE].get<double>();

        field_atm_kext.resize(n_volumes);
        field_atm_B.resize(n_volumes);
        field_sfc_B.resize(n_tiles);

        for (int i = 0; i < itot; i++)
        {
            for (int j = 0; j < jtot; j++)
            {
                for (int k = 0; k < ktot; k++)
                {
                    int idx = i*jtot*ktot + j*ktot + k;
                    field_atm_kext[idx] = arr_kext[i];
                    field_atm_B[idx] = arr_Batm[i];
                    if (i == 0) field_sfc_B[idx] = Bsfc;
                }
            }
        }


    }
    else if (caseID == "s3D") // simple 3D cases
    {
        arr_z = j[CASE + "_z"].get<std::vector<double>>();
        arr_zh = j[CASE + "_zh"].get<std::vector<double>>();
        arr_dz = j[CASE + "_dz"].get<std::vector<double>>();

        std::vector<int> case_limits = j[CASE + "_case_limits"].get<std::vector<int>>();
        itot = case_limits[0];
        jtot = case_limits[1];
        ktot = case_limits[2];
        int n_volumes = itot * jtot * ktot;
        int n_tiles = jtot * ktot;

        std::vector<int> cell_size_horizontal = j[CASE + "_cell_size_horizontal"].get<std::vector<int>>();
        dx = cell_size_horizontal[0];
        dy = cell_size_horizontal[1];

        Bsfc = j[CASE + "_Bsfc"].get<double>();

        // Generating kext and Batm arrays
        std::vector<double> col_kext_open  = j[CASE + "_kext_open"].get<std::vector<double>>();
        std::vector<double> col_kext_cloud = j[CASE + "_kext_cloud"].get<std::vector<double>>();
        std::vector<double> col_Batm_open  = j[CASE + "_Batm_open"].get<std::vector<double>>();
        std::vector<double> col_Batm_cloud = j[CASE + "_Batm_cloud"].get<std::vector<double>>();

        std::vector<int> cloud_coords_x = j[CASE + "_cloud_coords_x"].get<std::vector<int>>();
        std::vector<int> cloud_coords_y = j[CASE + "_cloud_coords_y"].get<std::vector<int>>();
        int Ncloud_coords = cloud_coords_x.size();
        int Nopen_coords  = jtot*ktot - Ncloud_coords;


        std::vector<int> open_coords_x(Nopen_coords);
        std::vector<int> open_coords_y(Nopen_coords);
        int open_idx = 0;
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                // going through each cloud coord pair
                bool is_cloud_column = false;
                for (int idx_cloud = 0; idx_cloud < Ncloud_coords; idx_cloud++)
                {
                    int cloud_x = cloud_coords_x[idx_cloud];
                    int cloud_y = cloud_coords_y[idx_cloud];

                    if ((j == cloud_y) && (k == cloud_x)) {is_cloud_column = true; break; }
                }

                if (!is_cloud_column)
                {
                    open_coords_x[open_idx] = k;
                    open_coords_y[open_idx] = j;
                    open_idx += 1;
                }
            }
        }

        field_sfc_B.resize(n_tiles);
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_sfc = j*ktot + k;
                field_sfc_B[idx_sfc] = Bsfc;
            }
        }

        // Filling in the kext and batm arrays
        field_atm_kext.resize(n_volumes);
        field_atm_B.resize(n_volumes);
        // Open columns
        for (int idx_open = 0; idx_open < Nopen_coords; idx_open++)
        {
            int k = open_coords_x[idx_open];
            int j = open_coords_y[idx_open];
            for (int i = 0; i < itot; i++)
            {
                int idx = i*jtot*ktot + j*ktot + k;

                field_atm_kext[idx] = col_kext_open[i];
                field_atm_B[idx]    = col_Batm_open[i];
            }
        }
        // Cloudy columns
        for (int idx_cloud = 0; idx_cloud < Ncloud_coords; idx_cloud++)
        {
            int k = cloud_coords_x[idx_cloud];
            int j = cloud_coords_y[idx_cloud];
            for (int i = 0; i < itot; i++)
            {
                int idx = i*jtot*ktot + j*ktot + k;

                field_atm_kext[idx] = col_kext_cloud[i];
                field_atm_B[idx]    = col_Batm_cloud[i];
            }
        }

    }
    else if (caseID == "r3D") // "Real" 3D cases (from the gpoints data)
    {
        arr_z = j["z"].get<std::vector<double>>();
        arr_zh = j["zh"].get<std::vector<double>>();
        arr_dz = j["dz"].get<std::vector<double>>();

        int xy_size = j["xy_size"].get<int>();
        itot = arr_z.size();
        jtot = xy_size;
        ktot = xy_size;

        dx = 100;
        dy = 100;

        field_atm_B    = j["Batm_" + CASE].get<std::vector<double>>();
        field_atm_kext = j["kext_" + CASE].get<std::vector<double>>();
        field_sfc_B    = j["Bsfc_" + CASE].get<std::vector<double>>();
        

    }

    

        
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
                                field_atm_kext,
                                field_atm_B,
                                field_sfc_B,
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
                                verbose,
                                enable_full_counter_matrix,
                                Pesc_mode,
                                OUTPUT_MODE);
        
        // Storing output
        std::ofstream file_MCoutput("/home/gijs-hogeboom/dev/mclw/data_output/heating_rates/HR_MC_" + CASE + "_" + INTERCELL_TECHNIQUE + "_" + INTRACELL_TECHNIQUE + "_Natm" + std::to_string(Natm_pow) + "_Nsfc" + std::to_string(Nsfc_pow) + ".csv");
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
                                            field_atm_kext,
                                            field_atm_B,
                                            Bsfc,
                                            100,
                                            print_EB_PP,
                                            verbose);
        
        // Storing output
        
        std::ofstream file_PPoutput("/home/gijs-hogeboom/dev/mclw/data_output/heating_rates/HR_PP_" + CASE + ".csv");
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

    std::fstream file_MCinput("/home/gijs-hogeboom/dev/mclw/data_output/heating_rates/HR_MC_" + CASE + "_" + INTERCELL_TECHNIQUE + "_" + INTRACELL_TECHNIQUE + "_Natm" + std::to_string(Natm_pow) + "_Nsfc" + std::to_string(Nsfc_pow) + ".csv");
    std::fstream file_PPinput("/home/gijs-hogeboom/dev/mclw/data_output/heating_rates/HR_PP_" + CASE + ".csv");
    
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
    }

    RMSE_tot = std::sqrt(RMSE_tot);



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
        std::cout << "Timers:" << std::endl;
        std::cout << "   | MC time:     " << MC_time.count()/1000 << std::endl;
        std::cout << "   | PP time:     " << PP_time.count()/1000 << std::endl;
        std::cout << "Metrics:" << std::endl;
        std::cout << "   | RMSE:        " << RMSE_tot << std::endl;
        std::cout << "   | ME:          " << ME_tot << std::endl;
        std::cout << "====================================" << std::endl;
    }


    return 0;
}