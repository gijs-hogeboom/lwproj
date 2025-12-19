#include <iostream>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include <chrono>
#include <netcdf>
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

    using namespace netCDF;
    using namespace netCDF::exceptions;
    


    // Handling input
    std::string arg1 = "gpt21";
    std::string arg2 = "power";
    float arg3 = 20.;
    bool arg4 = false;
    bool arg5 = false;

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
        arg3 = std::stof(argv[3]);
    }
    if (argc > 4)
    {
        std::istringstream(argv[4]) >> arg4;
    }
    if (argc > 5)
    {
        std::istringstream(argv[5]) >> arg5;
    }

    std::cout << "Arguments | CASE: " << arg1 << ", INTERCELL_TECHNIQUE: " << arg2 << ", Nphot_pow: " << arg3 << ", Pesc_mode: " << arg4 << ", enable_scattering: " << arg5 << std::endl;



    // Run settings
    float Nphot_pow = arg3;
    long int Nphot = (long int) pow(2, Nphot_pow);


    std::string CASE = arg1;
    bool ENABLE_MC = true;                        // Enables Monte Carlo algorithm
    bool ENABLE_PP = true;                        // Enables plane-parallel algorithm

    std::string INTERCELL_TECHNIQUE = arg2;       // {uniform, power, power-gradient}

    bool Pesc_mode = arg4;
    bool enable_scattering = arg5;

    // Output parameters
    bool enable_full_counter_matrix = false;
    bool OUTPUT_3D = true;



    // Console output params
    bool print_EB_MC         = false;
    bool print_EB_PP         = false;
    bool verbose             = true;
    bool print_final_results = true;
    bool plot_results        = false;


    /////////////// READING CASE //////////////////

    if (verbose)
    {
        std::cout << "Start of program, reading input data" << std::endl;
    }

    std::string hr_filename_MC = "../data_output/heating_rates/HR_MC_" + CASE + "_" + INTERCELL_TECHNIQUE + 
                                    "_Nphot" + std::to_string(Nphot_pow) + "_Pesc" + std::to_string(Pesc_mode) + "_scatter" + std::to_string(enable_scattering) + ".csv";
    std::string hr_filename_PP = "../data_output/heating_rates/HR_PP_" + CASE + ".csv";

    // Initializing variables to load
    std::vector<float> arr_z;
    std::vector<float> arr_zh;
    std::vector<float> arr_dz;
    std::vector<float> field_atm_kext;
    std::vector<float> field_atm_B;
    std::vector<float> field_sfc_B;
    std::vector<float> field_atm_SSA;
    std::vector<float> field_atm_ASY;
    float Bsfc, dx, dy;
    int itot, jtot, ktot;

    // Reading the case json file
    std::string caseID = CASE.substr(0, 3);


    // File input handling
    if (caseID == "gpt") // MV cases
    {
        // Opening MV cases .json
        std::fstream file_cases("../data_input/" + caseID + "cases.json");
    
        if (!file_cases.is_open())
        {
            std::cerr << "Could not open " + caseID + "cases.json!" << std::endl;
            return 1;
        }
        
        nlohmann::json j;
        file_cases >> j;

        arr_z = j["z_" + CASE].get<std::vector<float>>();
        arr_zh = j["zh_" + CASE].get<std::vector<float>>();
        arr_dz = j["dz_" + CASE].get<std::vector<float>>();

        // jtot, ktot, dx and dy are all predetermined in this case
        itot = arr_z.size();
        jtot = 1;
        ktot = 1;
        int n_volumes = itot * jtot * ktot;
        int n_tiles = jtot * ktot;

        dx = 100;
        dy = 100;

        // Generating fields from 1D kext and Batm data
        std::vector<float> arr_kext = j["kext_" + CASE].get<std::vector<float>>();
        std::vector<float> arr_Batm = j["Batm_" + CASE].get<std::vector<float>>();
        std::vector<float> arr_SSA = j["SSA_" + CASE].get<std::vector<float>>();
        std::vector<float> arr_ASY = j["ASY_" + CASE].get<std::vector<float>>();
        Bsfc = j["Bsfc_" + CASE].get<float>();

        field_atm_kext.resize(n_volumes);
        field_atm_B.resize(n_volumes);
        field_atm_SSA.resize(n_volumes);
        field_atm_ASY.resize(n_volumes);
        field_sfc_B.resize(n_tiles);

        for (int i = 0; i < itot; i++)
        {
            for (int j = 0; j < jtot; j++)
            {
                for (int k = 0; k < ktot; k++)
                {
                    int idx = i*jtot*ktot + j*ktot + k;
                    field_atm_kext[idx] = arr_kext[i];
                    field_atm_B[idx]    = arr_Batm[i];
                    field_atm_SSA[idx]  = arr_SSA[i]; // Temporary, due to lack of data
                    field_atm_ASY[idx]  = arr_ASY[i]; // Temporary, due to lack of data
                    if (i == 0) field_sfc_B[idx] = Bsfc;
                }
            }
        }


    }
    else if (caseID == "s3D") // simple 3D cases
    {
        // Opening s3Dcases.json
        std::fstream file_cases("../mclw/data_input/" + caseID + "cases.json");
    
        if (!file_cases.is_open())
        {
            std::cerr << "Could not open " + caseID + "cases.json!" << std::endl;
            return 1;
        }
        
        nlohmann::json j;
        file_cases >> j;
        
        std::string caseNr = CASE.substr(3, CASE.length());

        arr_z = j[caseNr]["z"].get<std::vector<float>>();
        arr_zh = j[caseNr]["zh"].get<std::vector<float>>();
        arr_dz = j[caseNr]["dz"].get<std::vector<float>>();


        std::vector<int> case_limits = j[caseNr]["case_limits"].get<std::vector<int>>();
        itot = case_limits[0];
        jtot = case_limits[1];
        ktot = case_limits[2];
        int n_volumes = itot * jtot * ktot;
        int n_tiles = jtot * ktot;

        std::vector<int> cell_size_horizontal = j[caseNr]["cell_size_horizontal"].get<std::vector<int>>();
        dx = cell_size_horizontal[0];
        dy = cell_size_horizontal[1];

        Bsfc = j[caseNr]["Bsfc"].get<float>();

        // Generating kext and Batm arrays
        std::vector<float> col_kext_open  = j[caseNr]["kext_open"].get<std::vector<float>>();
        std::vector<float> col_kext_cloud = j[caseNr]["kext_cloud"].get<std::vector<float>>();
        std::vector<float> col_Batm_open  = j[caseNr]["Batm_open"].get<std::vector<float>>();
        std::vector<float> col_Batm_cloud = j[caseNr]["Batm_cloud"].get<std::vector<float>>();
        std::vector<float> col_SSA_open   = j[caseNr]["SSA_open"].get<std::vector<float>>();
        std::vector<float> col_SSA_cloud  = j[caseNr]["SSA_cloud"].get<std::vector<float>>();
        std::vector<float> col_ASY_open   = j[caseNr]["ASY_open"].get<std::vector<float>>();
        std::vector<float> col_ASY_cloud  = j[caseNr]["ASY_cloud"].get<std::vector<float>>();

        std::vector<int> cloud_coords_x = j[caseNr]["cloud_coords_x"].get<std::vector<int>>();
        std::vector<int> cloud_coords_y = j[caseNr]["cloud_coords_y"].get<std::vector<int>>();
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

        // Filling in the data arrays
        field_atm_kext.resize(n_volumes);
        field_atm_B.resize(n_volumes);
        field_atm_SSA.resize(n_volumes);
        field_atm_ASY.resize(n_volumes);
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
                field_atm_SSA[idx]  = col_SSA_open[i];
                field_atm_ASY[idx]  = col_ASY_open[i];
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
                field_atm_SSA[idx]  = col_SSA_cloud[i];
                field_atm_ASY[idx]  = col_ASY_cloud[i];
            }
        }


    }
    else if (caseID == "r3D") // "Real" 3D cases (from the gpoints data)
    {
        // Opening raw 3D lw optics + grid info .nc files
        NcFile nc_optics("../data_input/lw_optical_properties.nc", NcFile::read); 
        NcFile nc_gridinfo("../data_input/grid.nc", NcFile::read);

        // Loading vars
        NcVar nc_tau = nc_optics.getVar("lw_tau");
        NcVar nc_Batm = nc_optics.getVar("lay_source");
        NcVar nc_Bsfc = nc_optics.getVar("sfc_source");
        NcVar nc_SSA = nc_optics.getVar("lw_ssa");
        NcVar nc_ASY = nc_optics.getVar("lw_asy");

        NcVar nc_arr_z  = nc_gridinfo.getVar("lay");
        NcVar nc_arr_zh = nc_gridinfo.getVar("lev");


        // Obtaining information
        auto dims = nc_tau.getDims();
        long unsigned int nGpts = dims[0].getSize();
        long unsigned int itot_local = dims[1].getSize();
        long unsigned int jtot_local = dims[2].getSize();
        long unsigned int ktot_local = dims[3].getSize();
        int n_volumes = itot_local * jtot_local * ktot_local;
        int n_tiles = jtot_local * ktot_local;

        std::vector<double> darr_z;
        std::vector<double> darr_zh;
        std::vector<double> darr_dz;
        std::vector<double> dfield_atm_kext;
        std::vector<double> dfield_atm_B;
        std::vector<double> dfield_sfc_B;
        std::vector<double> dfield_atm_SSA;
        std::vector<double> dfield_atm_ASY;

        // Resizing vectors to fit data
        dfield_atm_kext.resize(n_volumes);
        dfield_atm_B.resize(n_volumes);
        dfield_atm_SSA.resize(n_volumes);
        dfield_atm_ASY.resize(n_volumes);
        dfield_sfc_B.resize(n_tiles);

        darr_z.resize(itot_local);
        darr_zh.resize(itot_local + 1);
        darr_dz.resize(itot_local);


        // Filling in data vectors
        nc_arr_z.getVar(darr_z.data());
        nc_arr_zh.getVar(darr_zh.data());

        for (int i = 0; i < itot_local; i++)
        {
            darr_dz[i] = darr_zh[i+1] - darr_zh[i];
        }

        long unsigned int chosen_gpt = std::stoi(CASE.substr(3, CASE.length()));
        nc_tau.getVar(  {chosen_gpt, 0, 0, 0}, {1, itot_local, jtot_local, ktot_local}, dfield_atm_kext.data());
        nc_Batm.getVar( {chosen_gpt, 0, 0, 0}, {1, itot_local, jtot_local, ktot_local}, dfield_atm_B.data());
        nc_SSA.getVar(  {chosen_gpt, 0, 0, 0}, {1, itot_local, jtot_local, ktot_local}, dfield_atm_SSA.data());
        nc_ASY.getVar(  {chosen_gpt, 0, 0, 0}, {1, itot_local, jtot_local, ktot_local}, dfield_atm_ASY.data());
        nc_Bsfc.getVar( {chosen_gpt, 0, 0},    {1, jtot_local, ktot_local},             dfield_sfc_B.data());
        
        dx = 100.;
        dy = 100.;     
        
        // casting doubles to floats, filling in actual values
        itot = (int) itot_local;
        jtot = (int) jtot_local;
        ktot = (int) ktot_local;

        field_atm_kext.resize(n_volumes);
        field_atm_B.resize(n_volumes);
        field_atm_SSA.resize(n_volumes);
        field_atm_ASY.resize(n_volumes);
        field_sfc_B.resize(n_tiles);

        arr_z.resize(itot_local);
        arr_zh.resize(itot_local + 1);
        arr_dz.resize(itot_local);

        for (int i = 0; i < itot; i++)
        {
            arr_z[i] = (float) darr_z[i];
            arr_zh[i] = (float) darr_zh[i];
            arr_dz[i] = (float) darr_dz[i];
            for (int j = 0; j < jtot; j++)
            {
                for (int k = 0; k < ktot; k++)
                {
                    int idx = i * n_tiles + j*ktot + k;

                    field_atm_kext[idx] = (float) (dfield_atm_kext[idx] / arr_dz[i]);
                    field_atm_B[idx] = (float) dfield_atm_B[idx];
                    field_atm_SSA[idx] = (float) dfield_atm_SSA[idx];
                    field_atm_ASY[idx] = (float) dfield_atm_ASY[idx];

                    if (i == 0)
                    {
                        field_sfc_B[idx] = (float) dfield_sfc_B[idx];
                    }
                }
            }
        }
        arr_zh[itot] = (float) darr_zh[itot]; // filling in final upper value, not captured in the loop
    }

    

        
    std::vector<float> heating_rates_MC(itot);
    std::vector<float> heating_rates_PP(itot);

    auto MC_t1 = std::chrono::high_resolution_clock::now();
    if (ENABLE_MC)
    {
        if (verbose)
        {
            std::cout << "Start of MC - Nphot: " + std::to_string(Nphot_pow) << std::endl;
        }

        heating_rates_MC = run_MC(arr_z,
                                arr_zh,
                                arr_dz,
                                field_atm_kext,
                                field_atm_B,
                                field_sfc_B,
                                field_atm_SSA,
                                field_atm_ASY,
                                dx,
                                dy,
                                ktot,
                                jtot,
                                itot,
                                INTERCELL_TECHNIQUE,
                                CASE,
                                Nphot,
                                print_EB_MC,
                                verbose,
                                enable_full_counter_matrix,
                                Pesc_mode,
                                OUTPUT_3D,
                                enable_scattering);
        
        // Storing output
        std::ofstream file_MCoutput(hr_filename_MC);
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
                                              field_sfc_B,
                                              CASE,
                                              dx,
                                              dy,
                                              jtot,
                                              ktot,
                                              print_EB_PP,
                                              verbose,
                                              OUTPUT_3D);
        
        // Storing output
        
        std::ofstream file_PPoutput(hr_filename_PP);
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
    std::vector<float> heating_rates_PP_in(itot, 0.);
    std::vector<float> heating_rates_MC_in(itot, 0.);
    std::vector<float> arr_z_in(itot, 0.);

    std::fstream file_MCinput("../data_output/heating_rates/HR_MC_" + CASE + "_" + INTERCELL_TECHNIQUE + 
                              "_Nphot" + std::to_string(Nphot_pow) + "_Pesc" + std::to_string(Pesc_mode) + "_scatter" + std::to_string(enable_scattering) + ".csv");
    std::fstream file_PPinput("../data_output/heating_rates/HR_PP_" + CASE + ".csv");
    
    if (!file_MCinput.is_open())
    {
        std::cout << "File |HR_MC_" + CASE + "_" + INTERCELL_TECHNIQUE + "_Nphot" + std::to_string(Nphot_pow) + "_Pesc" + 
                      std::to_string(Pesc_mode) + "_scatter" + std::to_string(enable_scattering) + ".csv| does not exist!" << std::endl;
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
    float RMSE_tot, ME_tot;

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

        std::vector<std::pair<float,float>> hr_1D, hr_3D;
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
        std::cout << "   | Case:        " << CASE << std::endl;
        std::cout << "   | IT sampling: " << INTERCELL_TECHNIQUE << std::endl;
        std::cout << "   | Nphot:       " << Nphot_pow << std::endl;
        std::cout << "   | Pesc_mode:   " << Pesc_mode << std::endl;
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