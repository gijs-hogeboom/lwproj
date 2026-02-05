#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <omp.h>
#include <atomic>

#include "util.h"
#include "precision.h"



std::vector<Real> run_plane_parallel(const std::vector<Real>& arr_z,
                                     const std::vector<Real>& arr_zh,
                                     const std::vector<Real>& arr_dz,
                                     const std::vector<Real>& field_atm_kext,
                                     const std::vector<Real>& field_atm_B,
                                     const std::vector<Real>& field_sfc_B,
                                     const std::string& CASE,
                                     const Real dx,
                                     const Real dy,
                                     const int jtot,
                                     const int ktot,
                                     const bool print_EB,
                                     const bool verbose,
                                     const bool OUTPUT_3D)
{

    // Initializing domain
    const int N_mu      = 20;

    const int itot      = arr_z.size();
    const int itoth     = arr_zh.size();
    const int n_volumes = field_atm_kext.size();
    const int n_tiles   = field_sfc_B.size();

    const int Column_batch = 1;

    std::vector<Real> field_atm_heating_rates(n_volumes);
    std::vector<Real> field_sfc_fluxes(n_tiles);
    std::vector<Real> field_TOA_fluxes(n_tiles);

    // Initializing progress bar
    if (verbose)
    {
        std::time_t starting_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::cout << "  PP: Starting with " + std::to_string(n_tiles) + " columns. Start time: " + std::ctime(&starting_time) << std::endl;
    }
    const int num_progress_bar = 100;
    std::vector<std::int64_t> progress_bar(num_progress_bar);

    for (int i = 0; i < num_progress_bar; i++)
    {
        progress_bar[i] = (n_tiles * (i + 1)) / num_progress_bar;
    }
    int idx_progress = 0;

    std::cout << '|';
    for (int _ = 0; _ < num_progress_bar; _++) std::cout << '.' << std::flush;
    std::cout << '|' << '\b';
    for (int _ = 0; _ < num_progress_bar; _++) std::cout << '\b' << std::flush;

    std::atomic<std::int64_t> columns_done{0};



    #pragma omp parallel
    {

        int column_counter = 0;

        // Going through each tile in the field
        #pragma omp for schedule(dynamic, Column_batch)
        for (int jk = 0; jk < n_tiles; jk++)
        {

            // Handling the fields into 1D
            std::vector<Real> arr_kext(itot);
            std::vector<Real> arr_Batm(itot);

            for (int i = 0; i < itot; i++)
            {
                int idx_atm = i*n_tiles + jk;
                arr_kext[i] = field_atm_kext[idx_atm];
                arr_Batm[i] = field_atm_B[idx_atm];
            }


            // Generating tauh (cumulative sum of kext*dz from TOA to bottom)
            std::vector<Real> arr_tauh(itoth);
            arr_tauh[itoth-1] = 0.;
            for (int i = (itoth-2); i >= 0; i--)
            {
                arr_tauh[i] = arr_tauh[i+1] + (arr_dz[i]*arr_kext[i]);
            }


            Real tau_at_sfc = arr_tauh[0];
            Real tau_at_TOA = arr_tauh[itoth-1];

            Real I_at_sfc = field_sfc_B[jk];
            Real I_at_TOA = 0.;
            

            // Angles
            std::vector<Real> arr_mu = linspace(0.01, 1.0, N_mu);
            int jtot = N_mu;


            // Calculating I(tau, mu) (p.32 of K. N. Liou, 2002)
            std::vector<Real> M_I_uph(itoth*jtot);
            std::vector<Real> M_I_downh(itoth*jtot);

            for (int ih = 0; ih < itoth; ih++)
            {
                // Loading tau
                Real tau = arr_tauh[ih];

                for (int j = 0; j < jtot; j++)
                {
                    // Loading mu
                    Real mu = arr_mu[j];

                    // Calculating extinction term
                    Real extinction_term_uph   = I_at_sfc * std::exp( -(tau_at_sfc - tau)/mu );
                    Real extinction_term_downh = I_at_TOA * std::exp( -(tau - tau_at_TOA)/mu );

                    // Calculating emission term
                    Real emission_term_uph = 0.;
                    Real emission_term_downh = 0.;
                    for (int k = 0; k < ih; k++)
                    {
                        Real temp = (arr_Batm[k] * (std::exp(-(arr_tauh[k+1] - tau)/mu) - std::exp(-(arr_tauh[k] - tau)/mu))); 
                        emission_term_uph += temp;
                    }
                    for (int k = (itot-1); k > (ih - 1); k--)
                    {
                        Real temp = (arr_Batm[k] * (std::exp(-(tau - arr_tauh[k])/mu) - std::exp(-(tau - arr_tauh[k+1])/mu)));
                        emission_term_downh += temp;
                    }    

                    // Combining extinction + emission term
                    int idx        = ih*jtot + j;
                    M_I_uph[idx]   = extinction_term_uph + emission_term_uph;
                    M_I_downh[idx] = extinction_term_downh + emission_term_downh;
                }
            }


            // Calculating F(z)
            std::vector<Real> arr_F_uph(itoth);
            std::vector<Real> arr_F_downh(itoth);

            for (size_t ih = 0; ih < itoth; ih++)
            {
                // Calculating and storing I(tau, mu) * mu for each angle per level
                std::vector<Real> arr_Imu_uph(jtot);
                std::vector<Real> arr_Imu_downh(jtot);
                for (size_t j = 0; j < jtot; j++)
                {
                    Real mu = arr_mu[j];

                    int idx = ih*jtot + j;

                    arr_Imu_uph[j]   = M_I_uph[idx]*mu;
                    arr_Imu_downh[j] = M_I_downh[idx]*mu;
                }
                arr_F_uph[ih]   = static_cast<Real>(2.)*constants::PI*trapezoid(arr_mu, arr_Imu_uph);
                arr_F_downh[ih] = static_cast<Real>(2.)*constants::PI*trapezoid(arr_mu, arr_Imu_downh);
            }

            // Calculating net flux at each cell
            std::vector<Real> arr_F_net(itot);
            for (size_t i = 0; i < itot; i++)
            {
                arr_F_net[i] = arr_F_uph[i] + arr_F_downh[i+1] - arr_F_uph[i+1] - arr_F_downh[i];
            }
            
            std::vector<double> arr_heating_rates(itot);
            for (size_t i = 0; i < itot; i++)
            {
                arr_heating_rates[i] = 1 / (constants::RHO * constants::CP * arr_dz[i]) * arr_F_net[i] * static_cast<Real>(86400.);
            }

            // Filling up the local 3D fields
            for (int i = 0; i < itot; i++)
            {
                int idx_atm = i*n_tiles + jk;
                field_atm_heating_rates[idx_atm] = arr_heating_rates[i];
            }
            field_sfc_fluxes[jk] = arr_F_downh[0] - arr_F_uph[0];
            field_TOA_fluxes[jk] = arr_F_uph[itoth-1] - arr_F_downh[itoth-1];


            column_counter++;

            if (column_counter == Column_batch)
            {
                column_counter = 0;
                
                columns_done.fetch_add(Column_batch, std::memory_order_relaxed);

                if (omp_get_thread_num() == 0)
                {
                    auto done = columns_done.load(std::memory_order_relaxed);

                    #pragma omp critical
                    {
                        while (idx_progress < num_progress_bar && done >= progress_bar[idx_progress])
                        {
                            std::cout << '#' << std::flush;
                            idx_progress++;
                        }
                    }
                }
            }
        }
    }
    std::cout << '\n';



    if (OUTPUT_3D)
    {
        std::string atm_output_name = "hr_3DPP_atm_"   + CASE + ".dat";
        std::string sfc_output_name = "flux_3DPP_sfc_" + CASE + ".dat";
        std::string TOA_output_name = "flux_3DPP_TOA_" + CASE + ".dat";
        std::ofstream atm_output("../data_output/raw_output_3DPP/" + atm_output_name, std::ios::binary);
        std::ofstream sfc_output("../data_output/raw_output_3DPP/" + sfc_output_name, std::ios::binary);
        std::ofstream TOA_output("../data_output/raw_output_3DPP/" + TOA_output_name, std::ios::binary);
        int atm_dims[3] = {itot, jtot, ktot};
        int sfc_dims[2] = {jtot, ktot};

        std::vector<float> field_atm_heating_rates_f(n_volumes);
        std::vector<float> field_sfc_fluxes_f(n_tiles);
        std::vector<float> field_TOA_fluxes_f(n_tiles);

        for (int i = 0; i < itot; i++)
        {
            for (int j = 0; j < jtot; j++)
            {
                for (int k = 0; k < ktot; k++)
                {
                    long int idx = i*jtot*ktot + j*ktot + k;
                    field_atm_heating_rates_f[idx] = static_cast<float>(field_atm_heating_rates[idx]);

                    if (i == 0)
                    {
                        field_sfc_fluxes_f[idx] = static_cast<float>(field_sfc_fluxes[idx]);
                        field_TOA_fluxes_f[idx] = static_cast<float>(field_TOA_fluxes[idx]);
                    }
                }
            }
        }

        atm_output.write(reinterpret_cast<char*>(atm_dims), sizeof(atm_dims));
        atm_output.write(reinterpret_cast<char*>(field_atm_heating_rates_f.data()), sizeof(float)*n_volumes);
        sfc_output.write(reinterpret_cast<char*>(sfc_dims), sizeof(sfc_dims));
        sfc_output.write(reinterpret_cast<char*>(field_sfc_fluxes_f.data()), sizeof(float)*n_tiles);
        TOA_output.write(reinterpret_cast<char*>(sfc_dims), sizeof(sfc_dims));
        TOA_output.write(reinterpret_cast<char*>(field_TOA_fluxes_f.data()), sizeof(float)*n_tiles);
        atm_output.close();
        sfc_output.close();
        TOA_output.close();
    }


    // Averaging fields to 1D for quicklook
    std::vector<Real> arr_atm_heating_rates(itot);
    for (int i = 0; i < itot; i++)
    {
        std::vector<Real> temp_heating_rates(n_tiles);
        for (int jk = 0; jk < n_tiles; jk++)
        {
            int idx_atm = i*n_tiles + jk;
            temp_heating_rates[jk] = field_atm_heating_rates[idx_atm];
        }

        arr_atm_heating_rates[i] = std::accumulate(temp_heating_rates.begin(), temp_heating_rates.end(), static_cast<Real>(0.)) / n_tiles;
    }

    


    // PP energy balance
    // sfc_source          + TOA_source                       = sfc_sink                 + TOA_sink                 + atm_netto
    // pi*I_at_sfc         + 0                                = F_downh[0]               + F_uph[-1]                + F_net.sum()

    // double sfc_source = cdouble::PI * I_at_sfc;
    // double sfc_sink   = arr_F_downh[0];
    // double atm_netto  = kahan_sum(arr_F_net);
    // double TOA_source = 0.;
    // double TOA_sink   = arr_F_uph[itoth-1];


    // double netto_phi = sfc_source + TOA_source - sfc_sink - TOA_sink - atm_netto;
    // double netto_phi_percentage = netto_phi/(sfc_source + TOA_source) * 100.;


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