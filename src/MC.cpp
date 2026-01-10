#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>
#include <fstream>
#include <limits>
#include <string>
#include <numeric>
#include <algorithm>
#include <omp.h>
#include <chrono>

#include "util.h"
#include "photprop.h"




std::vector<float> run_MC(const std::vector<float>& arr_z,
                          const std::vector<float>& arr_zh,
                          const std::vector<float>& arr_dz,
                          const std::vector<float>& field_atm_kext,
                          const std::vector<float>& field_atm_B,
                          const std::vector<float>& field_sfc_B,
                          const std::vector<float>& field_atm_SSA,
                          const std::vector<float>& field_atm_ASY,
                          const float dx,
                          const float dy,
                          const int ktot,
                          const int jtot,
                          const int itot,
                          const std::string& INTERCELL_TECHNIQUE,
                          const std::string& CASE,
                          const long int Nphot,
                          const bool print_EB,
                          const bool verbose,
                          const bool enable_full_counter_matrix,
                          const bool Pesc_mode,
                          const bool OUTPUT_3D,
                          const bool enable_scattering)
{

    ///////////////////// INITIALIZING DOMAIN ///////////////////////
    if (verbose)
    {
        std::cout << "  MC: Initializing domain" << std::endl; 
    }
 
    // Randomization setup
    FastRNG rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    // Initializing domain skeleton
    int n_volumes = itot * jtot * ktot;
    int n_tiles   = jtot * ktot;

    float x_max = ktot * dx;
    float y_max = jtot * dy;
    float z_max = arr_zh[arr_zh.size() - 1];

    std::vector<float> arr_xh(ktot + 1);
    std::vector<float> arr_yh(jtot + 1);
    std::vector<float> arr_x(ktot);
    std::vector<float> arr_y(jtot);


    // Generating arr_x(h) and arr_y(h)
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


    // Initializing optical property fields
    std::vector<float> field_atm_phi(n_volumes);
    std::vector<float> field_atm_netto_power(n_volumes);

    std::vector<float> field_sfc_phi(n_tiles);
    std::vector<float> field_sfc_eps(n_tiles, 1.0);
    std::vector<float> field_sfc_netto_power(n_tiles);

    std::vector<float> field_TOA_netto_power(n_tiles);


    // Generating fields
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < itot; i++)
    {
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm = i * ktot * jtot + j * ktot + k;
                float current_kext = field_atm_kext[idx_atm];
                float current_Batm = field_atm_B[idx_atm];
                field_atm_phi[idx_atm] = 4*cfloat::PI * current_kext * current_Batm * dx * dy * arr_dz[i];

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    field_sfc_phi[idx_sfc] = cfloat::PI * field_sfc_eps[idx_sfc] * field_sfc_B[idx_sfc] * dx * dy;
                }
            }
        }
    }


    // Loading Pesc-kext curves

    if (Pesc_mode)
    {
        if (verbose)
        {
            std::cout << "  MC: Attributing Pesc to cells" << std::endl;
        }
        for (int i = 0; i < itot; i++)
        {

            float dz = arr_dz[i];
            std::string Pesc_path_name = f_Pesccurve_name(dx, dy, dz);  

            std::fstream Pesc_curve(Pesc_path_name);

            if (!Pesc_curve.is_open())
            {
                std::cout << "WARNING! Pesc curve '" << Pesc_path_name << "' not found!" << std::endl;
                std::vector<float> exit_vec(1);
                return exit_vec;
            }
            
            size_t N_points = count_lines(Pesc_curve) - 1;
            Pesc_curve.clear();
            Pesc_curve.seekg(0, std::ios::beg);
            

            std::vector<float> arr_kext(N_points);
            std::vector<float> arr_Pesc(N_points);

            std::string line;
            int idx = 0;

            std::getline(Pesc_curve, line); // Skipping the header

            while (std::getline(Pesc_curve, line))
            {
                std::stringstream ss(line);
                std::string cell;
                std::getline(ss, cell, ','); // index col
                std::getline(ss, cell, ','); // kext value
                float kext = std::stod(cell);
                if ((kext >= 1e-15) && (kext <= 1e5))
                {
                    arr_kext[idx] = kext;
                }
                std::getline(ss, cell, ','); // Pesc value
                float Pesc = std::stod(cell);
                if ((kext >= 1e-15) && (kext <= 1e5))
                {
                    arr_Pesc[idx] = Pesc;
                }
                idx++;
            }
    
            LinearInterpolator_float f_Pesc(arr_kext, arr_Pesc);

            for (int j = 0; j < jtot; j++)
            {
                for (int k = 0; k < ktot; k++)
                {
                    int idx = i*jtot*ktot + j*ktot + k;

                    float kext = field_atm_kext[idx];

                    float Pesc = f_Pesc(kext);

                    float emitted_power = Pesc * field_atm_phi[idx];

                    field_atm_phi[idx] = emitted_power;
                }
            }
        }
    }

    // Partitioning photons between Atm and Sfc
    float phi_atm_total = std::accumulate(field_atm_phi.begin(), field_atm_phi.end(), 0.0);
    float phi_sfc_total = std::accumulate(field_sfc_phi.begin(), field_sfc_phi.end(), 0.0);
    float phi_ratio_atm_sfc = phi_atm_total / (phi_atm_total + phi_sfc_total);

    long int Natm = (long int)  (phi_ratio_atm_sfc * Nphot);  
    long int Nsfc = Nphot - Natm;

    // int Natm = Nphot/2;
    // int Nsfc = Nphot/2;


    if (verbose)
    {
        std::cout << "  MC: Ratio Natm / Nsfc: " << phi_ratio_atm_sfc * 100.0 << '%' << std::endl;
    }





    /////////////////////////// ALIAS TABLES /////////////////////////// 
    
    if (verbose)
    {
        std::cout << "  MC: Preparing photon sampling Alias Tables" << std::endl;
    }
    
    
    // Generating AliasTabel: this will act as the PDF for where photons are
    // sampled within the domain
    std::vector<float> aliastable_weights_atm;
    std::vector<float> aliastable_weights_sfc;

    
    if (INTERCELL_TECHNIQUE == "uniform")
    {
        aliastable_weights_atm = std::vector<float>(n_volumes, 1.0);
        aliastable_weights_sfc = std::vector<float>(n_tiles, 1.0);
    } else
    if (INTERCELL_TECHNIQUE == "power")
    {
        aliastable_weights_atm = field_atm_phi;
        aliastable_weights_sfc = field_sfc_phi;

    } else 
    if (INTERCELL_TECHNIQUE == "power-gradient")
    {
        // Creating the power-gradient field (atmosphere)
        std::vector<float> field_atm_phi_gradient(n_volumes);
        for (size_t i = 0; i < itot; i++)
        {
            for (size_t j = 0; j < jtot; j++)
            {
                for (size_t k = 0; k < ktot; k++)
                {

                    // Creating the supposed indexes for accessing the values for calculating the gradient
                    size_t idx_atm    = i*jtot*ktot + j*ktot + k;
                    size_t idx_z_up   = (i+1)*jtot*ktot + j*ktot + k;
                    size_t idx_z_down = (i-1)*jtot*ktot + j*ktot + k;
                    size_t idx_y_pos  = i*jtot*ktot + (j+1)*ktot + k;
                    size_t idx_y_neg  = i*jtot*ktot + (j-1)*ktot + k;
                    size_t idx_x_pos  = i*jtot*ktot + j*ktot + (k+1);
                    size_t idx_x_neg  = i*jtot*ktot + j*ktot + (k-1);
                    
                    float dz_tot, dx_tot, dy_tot;
                    float phi_z_up, phi_z_down, phi_y_pos, phi_y_neg, phi_x_pos, phi_x_neg;

                    dx_tot = 2*dx;
                    dy_tot = 2*dy;

                    // Fetching values
                    // z-direction
                    if (i == (itot-1))
                    {
                        dz_tot     = 0.5*arr_dz[i-1] + 1.5*arr_dz[i];
                        phi_z_up   = 0.;
                        phi_z_down = field_atm_phi[idx_z_down];
                    } else
                    if (i == 0)
                    {
                        dz_tot     = 0.5*arr_dz[i];
                        size_t idx_sfc = j*ktot + k;
                        phi_z_up   = field_atm_phi[idx_z_up];
                        phi_z_down = field_atm_phi[idx_atm]; // ignoring surface
                    }
                    else
                    {
                        dz_tot     = 0.5*arr_dz[i-1] + arr_dz[i] + 0.5*arr_dz[i+1];
                        phi_z_up   = field_atm_phi[idx_z_up];
                        phi_z_down = field_atm_phi[idx_z_down];
                    }

                    // y-direction
                    if (j == (jtot-1))
                    {
                        idx_y_pos = i*jtot*ktot + k;
                    }
                    if (j == 0)
                    {
                        idx_y_neg = i*jtot*ktot + (jtot-1)*ktot + k;
                    }

                    phi_y_pos = field_atm_phi[idx_y_pos];
                    phi_y_neg = field_atm_phi[idx_y_neg];


                    // x-direction
                    if (k == (ktot-1))
                    {
                        idx_x_pos = i*jtot*ktot + j*ktot;
                    }
                    if (k == 0)
                    {
                        idx_x_neg = i*jtot*ktot + j*ktot + (ktot-1);
                    }

                    phi_x_pos = field_atm_phi[idx_x_pos];
                    phi_x_neg = field_atm_phi[idx_x_neg];
                    


                    // Calculating the gradient
                    float term_x = (phi_x_pos - phi_x_neg)/dx_tot;
                    float term_y = (phi_y_pos - phi_y_neg)/dy_tot;
                    float term_z = (phi_z_up - phi_z_down)/dz_tot;

                    float gradient_magnitude = sqrt(term_x*term_x + term_y*term_y + term_z*term_z);

                    // Inserting value
                    field_atm_phi_gradient[idx_atm] = gradient_magnitude;

                }
            }
        }

        // Surface (only gradient in z direction is considered important)
        std::vector<float> field_sfc_phi_gradient(n_tiles);
        for (size_t j = 0; j < jtot; j++)
        {
            for (size_t k = 0; k < ktot; k++)
            {
                size_t idx = j*ktot + k;
                field_sfc_phi_gradient[idx] = std::abs((field_sfc_phi[idx] - field_atm_phi[idx])/arr_z[0]);
            }
        }


        // Atmosphere
        // Generating weighted choice
        aliastable_weights_atm = field_atm_phi_gradient;
        aliastable_weights_sfc = field_sfc_phi_gradient;

    }



    AliasTable_float AliasTable_atm(aliastable_weights_atm);
    AliasTable_float AliasTable_sfc(aliastable_weights_sfc);   


    
    /////////////////////// PHOTON RELEASE /////////////////////// 
    if (verbose) 
    {
        std::cout << "  MC: Photon release - atm" << std::endl;
    }


    // Photons from the atmosphere
    photon_propagation(AliasTable_atm,
                       field_atm_kext,
                       field_sfc_eps,
                       field_atm_SSA,
                       field_atm_ASY,
                       arr_xh,
                       arr_yh,
                       arr_zh,
                       arr_x,
                       arr_y,
                       arr_z,
                       arr_dz,
                       field_atm_phi,
                       field_atm_netto_power,
                       field_sfc_netto_power,
                       field_TOA_netto_power,
                       Natm,
                       0,
                       INTERCELL_TECHNIQUE,
                       Pesc_mode,
                       enable_scattering);

    if (verbose)
    {
        std::cout << "  MC: Photon release - sfc" << std::endl;
    }
    // Photons from the surface
    photon_propagation(AliasTable_sfc,
                       field_atm_kext,
                       field_sfc_eps,
                       field_atm_SSA,
                       field_atm_ASY,
                       arr_xh,
                       arr_yh,
                       arr_zh,
                       arr_x,
                       arr_y,
                       arr_z,
                       arr_dz,
                       field_sfc_phi,
                       field_atm_netto_power,
                       field_sfc_netto_power,
                       field_TOA_netto_power,
                       Nsfc,
                       1,
                       INTERCELL_TECHNIQUE,
                       Pesc_mode,
                       enable_scattering);

            
           
    ///////////////// CALCULATING HEATING RATES /////////////////

    if (verbose)
    {
        std::cout << "  MC: Calculating heating rates" << std::endl;
    }


    std::vector<float> field_atm_heating_rates(n_volumes);
    std::vector<float> field_sfc_heating_rates(n_tiles);
    for (int i = 0; i < itot; i++)
    {
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm =  i * ktot * jtot + j * ktot + k;
                float dz = arr_dz[i];
                field_atm_heating_rates[idx_atm] = field_atm_netto_power[idx_atm] / (cfloat::RHO * cfloat::CP * dx * dy * dz) * 86400;

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    field_sfc_heating_rates[idx_sfc] = field_sfc_netto_power[idx_sfc] / (cfloat::RHO * cfloat::CP * dx * dy) * 86400;
                }
            }
        }
    }


    // Output
    if (OUTPUT_3D)
    {
        std::ostringstream oss_Nphot;
        oss_Nphot << std::fixed << std::setprecision(2) << (float) std::log2(Nphot);
        std::string atm_output_name = "hr_3D_atm_"   + CASE + "_Nphot" + oss_Nphot.str() + "_" + INTERCELL_TECHNIQUE + "_Pesc" + std::to_string(Pesc_mode) + "_scatter" + std::to_string(enable_scattering) + ".dat";
        std::string sfc_output_name = "flux_3D_sfc_" + CASE + "_Nphot" + oss_Nphot.str() + "_" + INTERCELL_TECHNIQUE + "_Pesc" + std::to_string(Pesc_mode) + "_scatter" + std::to_string(enable_scattering) + ".dat";
        std::string TOA_output_name = "flux_3D_TOA_" + CASE + "_Nphot" + oss_Nphot.str() + "_" + INTERCELL_TECHNIQUE + "_Pesc" + std::to_string(Pesc_mode) + "_scatter" + std::to_string(enable_scattering) + ".dat";
        std::ofstream atm_output("../data_output/raw_output_3D/" + atm_output_name, std::ios::binary);
        std::ofstream sfc_output("../data_output/raw_output_3D/" + sfc_output_name, std::ios::binary);
        std::ofstream TOA_output("../data_output/raw_output_3D/" + TOA_output_name, std::ios::binary);
        int atm_dims[3] = {itot, jtot, ktot};
        int sfc_dims[2] = {jtot, ktot};

        
        atm_output.write(reinterpret_cast<char*>(atm_dims), sizeof(atm_dims));
        atm_output.write(reinterpret_cast<char*>(field_atm_heating_rates.data()), sizeof(float)*n_volumes);
        sfc_output.write(reinterpret_cast<char*>(sfc_dims), sizeof(sfc_dims));
        sfc_output.write(reinterpret_cast<char*>(field_sfc_netto_power.data()), sizeof(float)*n_tiles);
        TOA_output.write(reinterpret_cast<char*>(sfc_dims), sizeof(sfc_dims));
        TOA_output.write(reinterpret_cast<char*>(field_TOA_netto_power.data()), sizeof(float)*n_tiles);
        atm_output.close();
        sfc_output.close();
        TOA_output.close();
    }


    // Making 1D for comparing against plane-parallel
    std::vector<float> arr_atm_heating_rates_1D(itot);

    // Averaging the horizontal directions
    for (int i = 0; i < itot; i++)
    {
        std::vector<float> horizontal_values(jtot*ktot);
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm =  i * ktot * jtot + j * ktot + k;
                int idx_horizontal = j*ktot + k;
                horizontal_values[idx_horizontal] = field_atm_heating_rates[idx_atm];
            }
        }
        arr_atm_heating_rates_1D[i] = std::accumulate(horizontal_values.begin(), horizontal_values.end(), 0.0) / (jtot*ktot);

    }

    // MC energy balance
    float sum_of_net_phi_atm = std::accumulate(field_atm_netto_power.begin(), field_atm_netto_power.end(), 0.0);
    float sum_of_net_phi_sfc = std::accumulate(field_sfc_netto_power.begin(), field_sfc_netto_power.end(), 0.0);
    float sum_of_net_phi_TOA = std::accumulate(field_TOA_netto_power.begin(), field_TOA_netto_power.end(), 0.0);
    float sum_of_sums_net_phi = sum_of_net_phi_atm + sum_of_net_phi_sfc + sum_of_net_phi_TOA;


    if (print_EB)
    {
        std::cout << "+++ MC ENERGY BALANCE ++++++++++++++" << std::endl;
        std::cout << "+ -- net ---------------------------" << std::endl;
        std::cout << "+ sfc net:         " << sum_of_net_phi_sfc << std::endl;
        std::cout << "+ atm net:         " << sum_of_net_phi_atm << std::endl;
        std::cout << "+ TOA net:         " << sum_of_net_phi_TOA << std::endl;
        std::cout << "+ -- sum ---------------------------" << std::endl;
        std::cout << "+ netto:           " << sum_of_sums_net_phi << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    }


    

    return arr_atm_heating_rates_1D;
}