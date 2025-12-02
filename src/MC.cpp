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



void photon_propagation(const AliasTable_double& aliastable,
                        FastRNG rng,
                        const std::vector<double>& field_kext,
                        const std::vector<double>& arr_xh,
                        const std::vector<double>& arr_yh,
                        const std::vector<double>& arr_zh,
                        const std::vector<double>& arr_x,
                        const std::vector<double>& arr_y,
                        const std::vector<double>& arr_z,
                        const std::vector<double>& arr_dz,
                        std::vector<double>& field_phi,
                        std::vector<double>& field_atm_net_phi,
                        std::vector<double>& field_sfc_net_phi,
                        std::vector<double>& field_TOA_net_phi,
                        const int N,
                        const int domain_section,
                        const std::string& INTERCELL_TECHNIQUE,
                        const bool Pesc_mode)
{

    const int itot = arr_z.size();
    const int jtot = arr_y.size();
    const int ktot = arr_x.size();
    const double x_max = arr_xh[ktot];
    const double y_max = arr_yh[jtot];
    const double z_max = arr_zh[itot];
    const double cell_dx = arr_xh[1] - arr_xh[0];
    const double cell_dy = arr_yh[1] - arr_yh[0];
    
    const double eps = 2e-5;
    const int jktot = jtot*ktot;

    double photon_power;

    if (INTERCELL_TECHNIQUE == "power")
    {
        photon_power = std::accumulate(field_phi.begin(), field_phi.end(), 0.0) / N;
    }
    
    int idx_1_counter = 0;

    u_int64_t photon_not_tracked_counter = 0;

    int idx_photon = 0;
    while (idx_photon < N)
    {
        // Tracking wheter each photon is accounted for inside the photon counter matrix
        bool photon_is_tracked = false;

        // Tracking whether photon leaves the cell
        bool out_of_cell = false;
        if (domain_section == 1) out_of_cell = true; // Surface photons are always "out" of the surface

        
        // Sampling location within domain, determining photon power
        int idx_flat = aliastable.sample(rng);
        int idx_original = idx_flat; // storing starting position

        if (!(INTERCELL_TECHNIQUE == "power"))
        {
            double sample_weight = aliastable.weights[idx_flat];
            photon_power = field_phi[idx_flat] / (sample_weight * N); // either field_atm_phi or field_sfc_phi
        }


        // Initializing position/direction/optical thickness
        int idx_z, idx_y, idx_x;
        double x, y, z, mu, az, tau;

        if (domain_section == 0)
        {
            // Atmosphere
            idx_z  = idx_flat / jktot;
            int idx_2D = idx_flat % jktot;
            idx_y  = idx_2D / ktot;
            idx_x  = idx_2D % ktot;

            x = (idx_x + rng.uniform()) * cell_dx;
            y = (idx_y + rng.uniform()) * cell_dy;
            z = arr_zh[idx_z] + rng.uniform()*arr_dz[idx_z];

            mu = rng.uniform()*2 - 1;
            az = rng.uniform()*2*cdouble::PI;

            tau = -std::log(rng.uniform());
        }
        else if (domain_section == 1)
        {
            // Surface
            idx_z = 0;
            idx_y = idx_flat / ktot;
            idx_x = idx_flat % ktot;

            x = (idx_x + rng.uniform()) * cell_dx;
            y = (idx_y + rng.uniform()) * cell_dy;
            z = 0.;

            mu = std::sqrt(rng.uniform());
            az = rng.uniform()*2*cdouble::PI;

            tau = -std::log(rng.uniform());
        }

        // Calculating cartesian direction vector
        double s = std::sqrt(1 - mu*mu);
        double dx = s*std::cos(az);
        double dy = s*std::sin(az);
        double dz = mu;
        
        int counter = 0;
        
        // Starting propegation...
        while (tau > 1e-10)
        {
            

            // field boundary detection in x direction - wrapping
            bool at_far_wall_x     = (std::abs(x - x_max) < eps);
            bool going_forwards_x  = (dx >= 0.);
            if (at_far_wall_x && going_forwards_x) 
            { 
                x = 0.;
                idx_x = 0;
            }
            bool at_near_wall_x    = (std::abs(x) < eps);
            bool going_backwards_x = (dx < 0.);
            if (at_near_wall_x && going_backwards_x) 
            { 
                x = x_max; 
                idx_x = ktot - 1;
            }

            // field boundary detection in y direction - wrapping
            bool at_far_wall_y     = (std::abs(y - y_max) < eps);
            bool going_forwards_y  = (dy >= 0.);
            if (at_far_wall_y && going_forwards_y) 
            { 
                y = 0.; 
                idx_y = 0;
            }
            bool at_near_wall_y    = (std::abs(y) < eps);
            bool going_backwards_y = (dy < 0.);
            if (at_near_wall_y && going_backwards_y) 
            { 
                y = y_max; 
                idx_y = jtot - 1;
            }


            // field boundary detection in z direction - loss through TOA or absorbtion by surface
            bool at_TOA            = (std::abs(z - z_max) < eps);
            bool going_up          = (dz >= 0.);
            if (at_TOA && going_up)
            {
                tau = 0.;
                int idx_tile = idx_y * ktot + idx_x;

                field_TOA_net_phi[idx_tile] += photon_power;
                if (domain_section == 0)
                {
                    field_atm_net_phi[idx_original] -= photon_power;
                } else if (domain_section == 1)
                {
                    field_sfc_net_phi[idx_original] -= photon_power;
                }

                photon_is_tracked = true;
                break;
            }
            bool at_surface        = (std::abs(z) < eps);
            bool going_down        = (dz < 0.);
            if (at_surface && going_down)
            {
                tau = 0.;
                int idx_tile = idx_y * ktot + idx_x;
                
                field_sfc_net_phi[idx_tile] += photon_power;
                if (domain_section == 0)
                {
                    field_atm_net_phi[idx_original] -= photon_power;
                } else if (domain_section == 1)
                {
                    field_sfc_net_phi[idx_original] -= photon_power;
                }

                photon_is_tracked = true;
                break;
            }

            // Updating position
            idx_flat = idx_z*jktot + idx_y*ktot + idx_x;
            
            // Loading kext
            double current_kext = field_kext[idx_flat];



            // Scanning collision with cell boundaries
            double time_x, time_y, time_z;
            double dnx, dny, dnz;

            if (dx >= 0.) // x
            {
                dnx = arr_xh[idx_x + 1] - x;
            } else {
                dnx = arr_xh[idx_x] - x;
            }
            time_x = dnx/dx;

            if (dy >= 0.) // y
            {
                dny = arr_yh[idx_y + 1] - y;
            } else {
                dny = arr_yh[idx_y] - y;
            }
            time_y = dny/dy;

            if (dz >= 0.) // z
            {
                dnz = arr_zh[idx_z + 1] - z;
            } else {
                dnz = arr_zh[idx_z] - z;
            }
            time_z = dnz/dz;



            // Determinig the scaling factor based on which cell is hit (i.e., which direction takes the least amount of time)
            // Additionally, updating photon index position for next iteration (if photon extincts within the cell, the idx will not be used anyways)
            double dist_x, dist_y, dist_z;

            bool hit_x_wall = ((time_x <= time_y) && (time_x <= time_z));
            bool hit_y_wall = ((time_y <= time_x) && (time_y <= time_z));
            bool hit_z_wall = ((time_z <= time_x) && (time_z <= time_y));

            if (hit_x_wall)
            {
                dist_x = dnx;
                dist_y = time_x * dy;
                dist_z = time_x * dz;
                if (going_forwards_x) {idx_x += 1;} else {idx_x -= 1;}
            }
            else if (hit_y_wall)
            {
                dist_x = time_y * dx;
                dist_y = dny;
                dist_z = time_y * dz;
                if (going_forwards_y) {idx_y += 1;} else {idx_y -= 1;}
            }
            else if (hit_z_wall)
            {
                dist_x = time_z * dx;
                dist_y = time_z * dy;
                dist_z = dnz;
                if (going_up) {idx_z += 1;} else {idx_z -= 1;}
            }


            // Calculating distance traveled
            double ds = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
            double max_s = tau/current_kext;
            double tau_absorbed = current_kext*ds;


            if (ds < max_s)
            {
                tau -= tau_absorbed;
                x += dist_x;
                y += dist_y;
                z += dist_z;

                out_of_cell = true;
            }
            else
            {
                tau = 0.;
                double fs = max_s / ds;
                x += dist_x*fs;
                y += dist_y*fs;
                z += dist_z*fs;

                if (out_of_cell)
                {
                    field_atm_net_phi[idx_flat] += photon_power;
                    if (domain_section == 0)
                    {
                        field_atm_net_phi[idx_original] -= photon_power;
                    } else if (domain_section == 1)
                    {
                        field_sfc_net_phi[idx_original] -= photon_power;
                    }
                }

                photon_is_tracked = true;
            }
            
            counter++;
        }

        // Photon tracking metric
        if (!photon_is_tracked)
        {
            photon_not_tracked_counter++;
        }

        // Counting to the next photon
        if (Pesc_mode)
        {
            // For Pesc mode, only count to the next photon if the photon actually crossed cell boundaries
            if (out_of_cell) idx_photon++;
        }
        else
        {
            idx_photon++;
        }
    }

    if (photon_not_tracked_counter > 0)
    {
        std::cout << "!!!!!!!!!! WARNING !!!!!!!!!    " << photon_not_tracked_counter << " photons are not tracked!! " << std::endl;
    }


    std::cout << domain_section << ',' << idx_1_counter << std::endl;
}














std::vector<double> run_MC(const std::vector<double>& arr_z,
                          const std::vector<double>& arr_zh,
                          const std::vector<double>& arr_dz,
                          const std::vector<double>& field_atm_kext,
                          const std::vector<double>& field_atm_B,
                          const double Bsfc,
                          const double dx,
                          const double dy,
                          const int ktot,
                          const int jtot,
                          const int itot,
                          const std::string& INTERCELL_TECHNIQUE,
                          const std::string& INTRACELL_TECHNIQUE,
                          const std::string& CASE,
                          const int Natm,
                          const int Nsfc,
                          const bool print_EB,
                          const bool verbose,
                          const bool enable_full_counter_matrix,
                          const bool Pesc_mode)
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

    double x_max = ktot * dx;
    double y_max = jtot * dy;
    double z_max = arr_zh[arr_zh.size() - 1];

    std::vector<double> arr_xh(ktot + 1);
    std::vector<double> arr_yh(jtot + 1);
    std::vector<double> arr_x(ktot);
    std::vector<double> arr_y(jtot);


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
    std::vector<double> field_atm_phi(n_volumes);
    std::vector<double> field_atm_netto_power(n_volumes);

    std::vector<double> field_sfc_B(n_tiles);
    std::vector<double> field_sfc_phi(n_tiles);
    std::vector<double> field_sfc_eps(n_tiles, 1.0);
    std::vector<double> field_sfc_netto_power(n_tiles);

    std::vector<double> field_TOA_netto_power(n_tiles);


    // Generating fields
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < itot; i++)
    {
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm = i * ktot * jtot + j * ktot + k;
                double current_kext = field_atm_kext[idx_atm];
                double current_Batm = field_atm_B[idx_atm];
                field_atm_phi[idx_atm] = 4*cdouble::PI * current_kext * current_Batm * dx * dy * arr_dz[i];

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    field_sfc_B[idx_sfc] = Bsfc;
                    field_sfc_phi[idx_sfc] = cdouble::PI * field_sfc_eps[idx_sfc] * Bsfc * dx * dy;
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

            double dz = arr_dz[i];
            std::string Pesc_path_name = f_Pesccurve_name(dx, dy, dz);  

            std::fstream Pesc_curve(Pesc_path_name);

            if (!Pesc_curve.is_open())
            {
                std::cout << "WARNING! Pesc curve '" << Pesc_path_name << "' not found!" << std::endl;
                std::vector<double> exit_vec(1);
                return exit_vec;
            }
            
            size_t N_points = count_lines(Pesc_curve) - 1;
            Pesc_curve.clear();
            Pesc_curve.seekg(0, std::ios::beg);
            

            std::vector<double> arr_kext(N_points);
            std::vector<double> arr_Pesc(N_points);

            std::string line;
            int idx = 0;

            std::getline(Pesc_curve, line); // Skipping the header

            while (std::getline(Pesc_curve, line))
            {
                std::stringstream ss(line);
                std::string cell;
                std::getline(ss, cell, ','); // index col
                std::getline(ss, cell, ','); // kext value
                double kext = std::stod(cell);
                if ((kext >= 1e-15) && (kext <= 1e5))
                {
                    arr_kext[idx] = kext;
                }
                std::getline(ss, cell, ','); // Pesc value
                double Pesc = std::stod(cell);
                if ((kext >= 1e-15) && (kext <= 1e5))
                {
                    arr_Pesc[idx] = Pesc;
                }
                idx++;
            }
    
            LinearInterpolator_double f_Pesc(arr_kext, arr_Pesc);

            for (int j = 0; j < jtot; j++)
            {
                for (int k = 0; k < ktot; k++)
                {
                    int idx = i*jtot*ktot + j*ktot + k;

                    double kext = field_atm_kext[idx];

                    double Pesc = f_Pesc(kext);

                    double emitted_power = Pesc * field_atm_phi[idx];

                    field_atm_phi[idx] = emitted_power;
                }
            }
        }
    }








    /////////////////////////// SAMPLING /////////////////////////// 
    
    if (verbose)
    {
        std::cout << "  MC: Preparing photon sampling Alias Tables" << std::endl;
    }
    
    
    // Generating AliasTabel: this will act as the PDF for where photons are
    // sampled within the domain
    std::vector<double> aliastable_weights_atm;
    std::vector<double> aliastable_weights_sfc;

    
    if (INTERCELL_TECHNIQUE == "uniform")
    {
        aliastable_weights_atm = std::vector<double>(n_volumes, 1.0);
        aliastable_weights_sfc = std::vector<double>(n_tiles, 1.0);
    } else
    if (INTERCELL_TECHNIQUE == "power")
    {
        aliastable_weights_atm = field_atm_phi;
        aliastable_weights_sfc = field_sfc_phi;

    } else 
    if (INTERCELL_TECHNIQUE == "power-gradient")
    {
        // Creating the power-gradient field (atmosphere)
        std::vector<double> field_atm_phi_gradient(n_volumes);
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
        std::vector<double> field_sfc_phi_gradient(n_tiles);
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



    AliasTable_double AliasTable_atm(aliastable_weights_atm);
    AliasTable_double AliasTable_sfc(aliastable_weights_sfc);   


    
    /////////////////////// PHOTON RELEASE /////////////////////// 
    if (verbose) 
    {
        std::cout << "  MC: Photon release - atm" << std::endl;
    }


    // Photons from the atmosphere
    photon_propagation(AliasTable_atm,
                       rng,
                       field_atm_kext,
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
                       Pesc_mode);
    
    if (verbose)
    {
        std::cout << "  MC: Photon release - sfc" << std::endl;
    }
    // Photons from the surface
    photon_propagation(AliasTable_sfc,
                       rng,
                       field_atm_kext,
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
                       Pesc_mode);

            
           
    ///////////////// CALCULATING HEATING RATES /////////////////

    if (verbose)
    {
        std::cout << "  MC: Calculating heating rates" << std::endl;
    }


    std::vector<double> field_atm_heating_rates(n_volumes);
    std::vector<double> field_sfc_heating_rates(n_tiles);
    for (int i = 0; i < itot; i++)
    {
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm =  i * ktot * jtot + j * ktot + k;
                double dz = arr_dz[i];
                field_atm_heating_rates[idx_atm] = field_atm_netto_power[idx_atm] / (cdouble::RHO * cdouble::CP * dx * dy * dz) * 86400;

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    field_sfc_heating_rates[idx_sfc] = field_sfc_netto_power[idx_sfc] / (cdouble::RHO * cdouble::CP * dx * dy) * 86400;
                }
            }
        }
    }


    std::vector<double> arr_atm_heating_rates_1D(itot);
    // Averaging the horizontal directions
    for (int i = 0; i < itot; i++)
    {
        std::vector<double> horizontal_values(jtot*ktot);
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
    // sfc_source          + atm_source          + TOA_source = sfc_sink                 + atm_sink                 + TOA_sink
    // field_sfc_phi.sum() + field_atm_phi.sum() + 0          = field_sfc_absorbed.sum() + field_atm_absorbed.sum() + TOA_absorbed

    // double sfc_source = kahan_sum(field_sfc_emitted) / (dx*dy);
    // double sfc_sink   = kahan_sum(field_sfc_absorbed) / (dx*dy);
    // double atm_source = kahan_sum(field_atm_emitted) / (dx*dy);
    // double atm_sink   = kahan_sum(field_atm_absorbed) / (dx*dy);
    // double TOA_source = 0.;
    // double TOA_sink   = TOA_absorbed / (dx*dy);

    // double netto_phi = sfc_source + atm_source + TOA_source - sfc_sink - atm_sink - TOA_sink;
    // double netto_phi_percentage = netto_phi/(sfc_source + atm_source + TOA_source) * 100.;


    // if (print_EB)
    // {
    //     std::cout << "+++ MC ENERGY BALANCE ++++++++++++++" << std::endl;
    //     std::cout << "+ -- source ------------------------" << std::endl;
    //     std::cout << "+ sfc source:      " << sfc_source << std::endl;
    //     std::cout << "+ atm source:      " << atm_source << std::endl;
    //     std::cout << "+ TOA source:      " << TOA_source << std::endl;
    //     std::cout << "+ -- sinks -------------------------" << std::endl;
    //     std::cout << "+ sfc sink:        " << sfc_sink << std::endl;
    //     std::cout << "+ atm sink:        " << atm_sink << std::endl;
    //     std::cout << "+ TOA sink:        " << TOA_sink << std::endl;
    //     std::cout << "+ -- net ---------------------------" << std::endl;
    //     std::cout << "+ sfc net:         " << sfc_sink - sfc_source << std::endl;
    //     std::cout << "+ atm net:         " << atm_sink - atm_source << std::endl;
    //     std::cout << "+ TOA net:         " << TOA_sink - TOA_source << std::endl;
    //     std::cout << "+ -- sums --------------------------" << std::endl;
    //     std::cout << "+ sources:         " << sfc_source + atm_source + TOA_source << std::endl;
    //     std::cout << "+ sinks:           " << sfc_sink + atm_sink + TOA_sink << std::endl;
    //     std::cout << "+ sources - sinks: " << netto_phi << std::endl;
    //     std::cout << "+ as percentage:   " << netto_phi_percentage << std::endl;
    //     std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    // }


    

    return arr_atm_heating_rates_1D;
}