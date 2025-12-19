#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <boost/math/quadrature/gauss.hpp>
#include <chrono>
#include <ctime>

#include "util.h"


using namespace boost::math::quadrature;




std::vector<float> run_plane_parallel(const std::vector<float>& arr_z,
                                      const std::vector<float>& arr_zh,
                                      const std::vector<float>& arr_dz,
                                      const std::vector<float>& field_atm_kext,
                                      const std::vector<float>& field_atm_B,
                                      const std::vector<float>& field_sfc_B,
                                      const std::string& CASE,
                                      const float dx,
                                      const float dy,
                                      const int jtot,
                                      const int ktot,
                                      const bool print_EB,
                                      const bool verbose,
                                      const bool OUTPUT_3D)
{

    // Initializing domain
    const int N_mu      = 100;

    const int itot      = arr_z.size();
    const int itoth     = arr_zh.size();
    const int n_volumes = field_atm_kext.size();
    const int n_tiles   = field_sfc_B.size();

    std::vector<float> field_atm_heating_rates(n_volumes);
    std::vector<float> field_sfc_fluxes(n_tiles);
    std::vector<float> field_TOA_fluxes(n_tiles);

    if (verbose)
    {
        std::time_t starting_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::cout << "  PP: Starting with " + std::to_string(n_tiles) + " columns. Start time: " + std::ctime(&starting_time) << std::endl;
    }

    #pragma omp parallel
    {

        // Going through each tile in the field
        #pragma omp for
        for (int jk = 0; jk < n_tiles; jk++)
        {

            // Handling the fields into 1D
            std::vector<float> arr_kext(itot);
            std::vector<float> arr_Batm(itot);

            for (int i = 0; i < itot; i++)
            {
                int idx_atm = i*n_tiles + jk;
                arr_kext[i] = field_atm_kext[idx_atm];
                arr_Batm[i] = field_atm_B[idx_atm];
            }


            // Generating tauh (cumulative sum of kext*dz from TOA to bottom)
            std::vector<float> arr_tauh(itoth);
            arr_tauh[itoth-1] = 0.;
            for (int i = (itoth-2); i >= 0; i--)
            {
                arr_tauh[i] = arr_tauh[i+1] + (arr_dz[i]*arr_kext[i]);
            }


            float tau_at_sfc = arr_tauh[0];
            float tau_at_TOA = arr_tauh[itoth-1];

            float I_at_sfc = field_sfc_B[jk];
            float I_at_TOA = 0.;
            

            // Angles
            std::vector<float> arr_mu = linspace(0.01, 1.0, N_mu);
            int jtot = N_mu;


            // Calculating I(tau, mu) (p.32 of K. N. Liou, 2002)
            std::vector<float> M_I_uph(itoth*jtot);
            std::vector<float> M_I_downh(itoth*jtot);

            for (int ih = 0; ih < itoth; ih++)
            {
                // Loading tau
                float tau = arr_tauh[ih];

                for (int j = 0; j < jtot; j++)
                {
                    // Loading mu
                    float mu = arr_mu[j];

                    // Calculating extinction term
                    float extinction_term_uph   = I_at_sfc * std::exp( -(tau_at_sfc - tau)/mu );
                    float extinction_term_downh = I_at_TOA * std::exp( -(tau - tau_at_TOA)/mu );

                    // Calculating emission term
                    float emission_term_uph = 0.;
                    float emission_term_downh = 0.;
                    for (int k = 0; k < ih; k++)
                    {
                        float temp = (arr_Batm[k] * (std::exp(-(arr_tauh[k+1] - tau)/mu) - std::exp(-(arr_tauh[k] - tau)/mu))); 
                        emission_term_uph += temp;
                    }
                    for (int k = (itot-1); k > (ih - 1); k--)
                    {
                        float temp = (arr_Batm[k] * (std::exp(-(tau - arr_tauh[k])/mu) - std::exp(-(tau - arr_tauh[k+1])/mu)));
                        emission_term_downh += temp;
                    }    

                    // Combining extinction + emission term
                    int idx        = ih*jtot + j;
                    M_I_uph[idx]   = extinction_term_uph + emission_term_uph;
                    M_I_downh[idx] = extinction_term_downh + emission_term_downh;
                }
            }


            // Calculating F(z)
            std::vector<float> arr_F_uph(itoth);
            std::vector<float> arr_F_downh(itoth);

            for (size_t ih = 0; ih < itoth; ih++)
            {
                // Calculating and storing I(tau, mu) * mu for each angle per level
                std::vector<float> arr_Imu_uph(jtot);
                std::vector<float> arr_Imu_downh(jtot);
                for (size_t j = 0; j < jtot; j++)
                {
                    float mu = arr_mu[j];

                    int idx = ih*jtot + j;

                    arr_Imu_uph[j]   = M_I_uph[idx]*mu;
                    arr_Imu_downh[j] = M_I_downh[idx]*mu;
                }
                arr_F_uph[ih]   = 2.*cfloat::PI*trapezoid(arr_mu, arr_Imu_uph);
                arr_F_downh[ih] = 2.*cfloat::PI*trapezoid(arr_mu, arr_Imu_downh);
            }

            // Calculating net flux at each cell
            std::vector<float> arr_F_net(itot);
            for (size_t i = 0; i < itot; i++)
            {
                arr_F_net[i] = arr_F_uph[i] + arr_F_downh[i+1] - arr_F_uph[i+1] - arr_F_downh[i];
            }


            // Calculating heating rates at each cell
            // if (verbose)
            // {
            //     std::cout << "  PP: Calculating heating rates" << std::endl;
            // }
            
            std::vector<float> arr_heating_rates(itot);
            for (size_t i = 0; i < itot; i++)
            {
                arr_heating_rates[i] = 1 / (cfloat::RHO * cfloat::CP * arr_dz[i]) * arr_F_net[i] * 86400.;
            }

            // Filling up the local 3D fields
            for (int i = 0; i < itot; i++)
            {
                int idx_atm = i*n_tiles + jk;
                field_atm_heating_rates[idx_atm] = arr_heating_rates[i];
            }
            field_sfc_fluxes[jk] = (arr_F_downh[0] - arr_F_uph[0])*dx*dy;
            field_TOA_fluxes[jk] = (arr_F_uph[itoth-1] - arr_F_downh[itoth-1])*dx*dy;
        }
    }



    if (OUTPUT_3D)
    {
        std::string atm_output_name = "hr_3DPP_atm_"   + CASE + ".dat";
        std::string sfc_output_name = "flux_3DPP_sfc_" + CASE + ".dat";
        std::string TOA_output_name = "flux_3DPP_TOA_" + CASE + ".dat";
        std::ofstream atm_output("/home/gijs-hogeboom/dev/mclw/data_output/raw_output_3DPP/" + atm_output_name, std::ios::binary);
        std::ofstream sfc_output("/home/gijs-hogeboom/dev/mclw/data_output/raw_output_3DPP/" + sfc_output_name, std::ios::binary);
        std::ofstream TOA_output("/home/gijs-hogeboom/dev/mclw/data_output/raw_output_3DPP/" + TOA_output_name, std::ios::binary);
        int atm_dims[3] = {itot, jtot, ktot};
        int sfc_dims[2] = {jtot, ktot};

        
        atm_output.write(reinterpret_cast<char*>(atm_dims), sizeof(atm_dims));
        atm_output.write(reinterpret_cast<char*>(field_atm_heating_rates.data()), sizeof(float)*n_volumes);
        sfc_output.write(reinterpret_cast<char*>(sfc_dims), sizeof(sfc_dims));
        sfc_output.write(reinterpret_cast<char*>(field_sfc_fluxes.data()), sizeof(float)*n_tiles);
        TOA_output.write(reinterpret_cast<char*>(sfc_dims), sizeof(sfc_dims));
        TOA_output.write(reinterpret_cast<char*>(field_TOA_fluxes.data()), sizeof(float)*n_tiles);
        atm_output.close();
        sfc_output.close();
        TOA_output.close();
    }


    // Averaging fields to 1D for quicklook
    std::vector<float> arr_atm_heating_rates(itot);
    for (int i = 0; i < itot; i++)
    {
        std::vector<float> temp_heating_rates(n_tiles);
        for (int jk = 0; jk < n_tiles; jk++)
        {
            int idx_atm = i*n_tiles + jk;
            temp_heating_rates[jk] = field_atm_heating_rates[idx_atm];
        }

        arr_atm_heating_rates[i] = std::accumulate(temp_heating_rates.begin(), temp_heating_rates.end(), 0.0) / n_tiles;
    }

    


    // PP energy balance
    // sfc_source          + TOA_source                       = sfc_sink                 + TOA_sink                 + atm_netto
    // pi*I_at_sfc         + 0                                = F_downh[0]               + F_uph[-1]                + F_net.sum()

    // float sfc_source = cfloat::PI * I_at_sfc;
    // float sfc_sink   = arr_F_downh[0];
    // float atm_netto  = kahan_sum(arr_F_net);
    // float TOA_source = 0.;
    // float TOA_sink   = arr_F_uph[itoth-1];


    // float netto_phi = sfc_source + TOA_source - sfc_sink - TOA_sink - atm_netto;
    // float netto_phi_percentage = netto_phi/(sfc_source + TOA_source) * 100.;


    // if (print_EB)
    // {
    //     std::cout << "+++ PP ENERGY BALANCE ++++++++++++++" << std::endl;
    //     std::cout << "+ -- source ------------------------" << std::endl;
    //     std::cout << "+ sfc source:      " << sfc_source << std::endl;
    //     std::cout << "+ atm source:      -" << std::endl;
    //     std::cout << "+ TOA source:      " << TOA_source << std::endl;
    //     std::cout << "+ -- sinks -------------------------" << std::endl;
    //     std::cout << "+ sfc sink:        " << sfc_sink << std::endl;
    //     std::cout << "+ atm sink:        -" << std::endl;
    //     std::cout << "+ TOA sink:        " << TOA_sink << std::endl;
    //     std::cout << "+ -- net ---------------------------" << std::endl;
    //     std::cout << "+ sfc net:         " << sfc_sink - sfc_source << std::endl;
    //     std::cout << "+ atm net:         " << atm_netto << std::endl;
    //     std::cout << "+ TOA net:         " << TOA_sink - TOA_source << std::endl;
    //     std::cout << "+ -- sums --------------------------" << std::endl;
    //     std::cout << "+ sources:          unknown" << std::endl;
    //     std::cout << "+ sinks:            unknown" << std::endl;
    //     std::cout << "+ sources - sinks: " << netto_phi << std::endl;
    //     std::cout << "+ as percentage:   " << netto_phi_percentage << std::endl;
    //     std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    // }


    return arr_atm_heating_rates;
}