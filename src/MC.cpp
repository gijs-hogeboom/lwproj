#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>
#include <fstream>
#include <limits>
#include <string>

#include "util.h"


void photon_propagation(const std::vector<int>& arr_photons_pos_idx,
                        const std::vector<float>& arr_photons_pos_x,
                        const std::vector<float>& arr_photons_pos_y,
                        const std::vector<float>& arr_photons_pos_z,
                        const std::vector<float>& arr_photons_mu,
                        const std::vector<float>& arr_photons_az,
                        const std::vector<float>& arr_photons_tau,
                        const std::vector<float>& arr_photons_phi,
                        const std::vector<float>& field_kext,
                        const std::vector<float>& arr_xh,
                        const std::vector<float>& arr_yh,
                        const std::vector<float>& arr_zh,
                        const std::vector<float>& arr_x,
                        const std::vector<float>& arr_y,
                        const std::vector<float>& arr_z,
                        std::vector<float>& field_atm_absorbed,
                        std::vector<float>& field_sfc_absorbed,
                        float x_max,
                        float y_max,
                        float z_max,
                        int itot,
                        int jtot,
                        int ktot,
                        float dx,
                        float dy,
                        int N,
                        int domain_section)
{
    float eps = 1e-3;

    constexpr bool debug = false;
    std::ofstream fout; // debugging output
    
    if constexpr (debug) 
    {
        fout.open("debug" + std::to_string(domain_section) + ".dat");
        fout << "state,section,photon,steps,x,y,z,dx,dy,dz,tau,kext" << std::endl;
    }

    for (int idx_photon = 0; idx_photon < N; idx_photon++)
    {
        
        // Loading initial photon variables
        int idx_flat = arr_photons_pos_idx[idx_photon];
        float x      = arr_photons_pos_x[idx_photon];
        float y      = arr_photons_pos_y[idx_photon];
        float z      = arr_photons_pos_z[idx_photon];
        float mu     = arr_photons_mu[idx_photon];
        float az     = arr_photons_az[idx_photon];
        float tau    = arr_photons_tau[idx_photon];

        // Calculating cartesian direction vector
        float s = sqrt(1 - mu*mu);
        float dx = s*cos(az);
        float dy = s*sin(az);
        float dz = mu;

        // Finding photon position as idx of the field
        // TODO - does not work exactly on cell boundaries, since it does not take direction into account
        int idx_x = 0, idx_y = 0, idx_z = 0;
        for (int xi = 0; xi < ktot; xi++)
        {
            if ((arr_xh[xi] <= x) && (x < arr_xh[xi + 1]))
            {
                idx_x = xi;
                break;
            }
        }
        for (int yi = 0; yi < jtot; yi++)
        {
            if ((arr_yh[yi] <= y) && (y < arr_yh[yi + 1]))
            {
                idx_y = yi;
                break;
            }
        }
        for (int zi = 0; zi < itot; zi++)
        {
            if ((arr_zh[zi] <= z) && (z < arr_zh[zi + 1]))
            {
                idx_z = zi;
                break;
            }
        }


        int counter = 0;

        // Writing header for debug output
        if constexpr (debug)
        {
            fout << "START" << ',' << domain_section << ',' << idx_photon << ',' << counter << ',' << x << ',' << y << ',' << z << ',' << dx << ',' << dy << ',' << dz << ',' << tau << ',' << 9999. << std::endl;
        }

        

        // Starting propegation...
        while (tau > 1e-10)
        {

            // field boundary detection in x direction - wrapping
            bool at_far_wall_x     = (abs(x - x_max) < eps);
            bool going_forwards_x  = (dx >= 0.);
            if (at_far_wall_x && going_forwards_x) 
            { 
                x = 0.;
                idx_x = 0;
            }
            bool at_near_wall_x    = (abs(x) < eps);
            bool going_backwards_x = (dx < 0.);
            if (at_near_wall_x && going_backwards_x) 
            { 
                x = x_max; 
                idx_x = ktot - 1;
            }

            // field boundary detection in y direction - wrapping
            bool at_far_wall_y     = (abs(y - y_max) < eps);
            bool going_forwards_y  = (dy >= 0.);
            if (at_far_wall_y && going_forwards_y) 
            { 
                y = 0.; 
                idx_y = 0;
            }
            bool at_near_wall_y    = (abs(y) < eps);
            bool going_backwards_y = (dy < 0.);
            if (at_near_wall_y && going_backwards_y) 
            { 
                y = y_max; 
                idx_y = jtot - 1;
            }

            // field boundary detection in z direction - loss through TOA or absorbtion by surface
            bool at_TOA            = (abs(z - z_max) < eps);
            bool going_up          = (dz >= 0.);
            if (at_TOA && going_up)
            {
                tau = 0.;
                z = std::nextafter(z_max, std::numeric_limits<float>::infinity()); // transport just above TOA
                // ptr_energy_lost_at_TOA[0] += arr_photons_phi[idx_photon];
                break;
            }
            bool at_surface        = (abs(z) < eps);
            bool going_down        = (dz < 0.);
            if (at_surface && going_down)
            {
                tau = 0.;
                int idx_sfc = idx_y * ktot + idx_x;
                field_sfc_absorbed[idx_sfc] += arr_photons_phi[idx_photon];
                break;
            }

        
            // retrieving kext from flattened grid
            int idx_flat = idx_z * jtot * ktot + idx_y * ktot + idx_x;
            float current_kext = field_kext[idx_flat];



            // Scanning collision with cell boundaries
            float time_x, time_y, time_z;
            float dn = 0.;
            float f = 1.;

            if (dx >= 0.) // x
            {
                dn = arr_xh[idx_x + 1] - x;
            } else {
                dn = arr_xh[idx_x] - x;
            }
            time_x = dn/dx;

            if (dy >= 0.) // y
            {
                dn = arr_yh[idx_y + 1] - y;
            } else {
                dn = arr_yh[idx_y] - y;
            }
            time_y = dn/dy;
            if (dz >= 0.) // z
            {
                dn = arr_zh[idx_z + 1] - z;
            } else {
                dn = arr_zh[idx_z] - z;
            }
            time_z = dn/dz;


            // Determinig the scaling factor based on which cell is hit (i.e., which direction takes the least amount of time)
            // Additionally, updating photon index position for next iteration (if photon extincts within the cell, the idx will not be used anyways)
            if ((time_x <= time_y) && (time_x <= time_z))
            {
                f = time_x;
                if (going_forwards_x) {idx_x += 1;} else {idx_x -= 1;}
            }
            else if ((time_y <= time_x) && (time_y <= time_z))
            {
                f = time_y;
                if (going_forwards_y) {idx_y += 1;} else {idx_y -= 1;}
            }
            else if ((time_z <= time_x) && (time_z <= time_y))
            {
                f = time_z;
                if (going_up) {idx_z += 1;} else {idx_z -= 1;}
            }


            // Actual direction vectors
            float dist_x = f * dx;
            float dist_y = f * dy;
            float dist_z = f * dz;
            
            float ds = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
            float max_s = tau/current_kext;
            float tau_absorbed = current_kext*ds;

            if (ds < max_s)
            {
                tau -= tau_absorbed;
                x += dist_x;
                y += dist_y;
                z += dist_z;

            }
            else
            {
                tau = 0.;
                float fs = max_s / ds;
                x += dist_x*fs;
                y += dist_y*fs;
                z += dist_z*fs;

                field_atm_absorbed[idx_flat] += arr_photons_phi[idx_photon];
            }
            counter++;




            if constexpr (debug)
            {
                if (counter > 999980)
                {
                    fout << "..." << ',' << domain_section << ',' << idx_photon << ',' << counter << ',' << x << ',' << y << ',' << z << ',' << dx << ',' << dy << ',' << dz << ',' << tau << ',' << current_kext << std::endl;
                }
                if (counter > 1000000) 
                { 
                    std::cout << "Warning: counter > 1.000.000..." << std::endl;
                    fout << "END" << ',' << domain_section << ',' << idx_photon << ',' << counter << ',' << x << ',' << y << ',' << z << ',' << dx << ',' << dy << ',' << dz << ',' << tau << ',' << current_kext << std::endl;
                    exit(-1);
                }
            }            
        }

        


        if constexpr (debug)
        {
            fout << "END" << ',' << domain_section << ',' << idx_photon << ',' << counter << ',' << x << ',' << y << ',' << z << ',' << dx << ',' << dy << ',' << dz << ',' << tau << ',' << 9999. << std::endl;
        }
        
    }
}


std::vector<float> run_MC(const std::vector<float>& arr_z,
                          const std::vector<float>& arr_zh,
                          const std::vector<float>& arr_dz,
                          const std::vector<float>& arr_kext,
                          const std::vector<float>& arr_Batm,
                          float Bsfc,
                          float dx,
                          float dy,
                          int ktot,
                          int jtot,
                          int itot,
                          const std::string& CASE,
                          const std::string& INTERCELL_TECHNIQUE,
                          const std::string& INTRACELL_TECHNIQUE,
                          int Natm,
                          int Nsfc)
{

    ///////////////////// INITIALIZING DOMAIN /////////////////////// 
    std::cout << "  MC: Initializing domain" << std::endl;

    // Randomization setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> random_float(0.0f, 1.0f);

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
    std::vector<float> field_atm_kext(n_volumes);
    std::vector<float> field_atm_B(n_volumes);
    std::vector<float> field_atm_phi(n_volumes);
    std::vector<float> field_sfc_B(n_tiles);
    std::vector<float> field_sfc_phi(n_tiles);
    std::vector<float> field_sfc_eps(n_tiles, 1.0f);

    std::vector<float> field_atm_absorbed(n_volumes);
    std::vector<float> field_atm_emitted(n_volumes);
    std::vector<float> field_sfc_absorbed(n_tiles);
    std::vector<float> field_sfc_emitted(n_tiles);

    // Generating fields
    for (int i = 0; i < itot; i++)
    {
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm = i * ktot * jtot + j * ktot + k;
                float current_kext = arr_kext[i];
                float current_Batm = arr_Batm[i];

                field_atm_kext[idx_atm] = current_kext;
                field_atm_B[idx_atm] = current_Batm;
                field_atm_phi[idx_atm] = 4*cf::PI * current_kext * current_Batm * dx * dy * arr_dz[i];
                
                field_atm_emitted[idx_atm] = field_atm_phi[idx_atm]; // basically copying the power array

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    field_sfc_B[idx_sfc] = Bsfc;
                    field_sfc_phi[idx_sfc] = cf::PI * field_sfc_eps[idx_sfc] * Bsfc * dx * dy;

                    field_sfc_emitted[idx_sfc] = field_sfc_phi[idx_sfc]; // idem
                }
            }
        }
    }

    /////////////////////////// SAMPLING /////////////////////////// 
    std::cout << "  MC: Sampling" << std::endl;

    
    // Initializing photon arrays
    std::vector<int>   arr_photons_atm_pos_idx(Natm);
    std::vector<float> arr_photons_atm_pos_x(Natm);
    std::vector<float> arr_photons_atm_pos_y(Natm);
    std::vector<float> arr_photons_atm_pos_z(Natm);
    std::vector<float> arr_photons_atm_mu(Natm);
    std::vector<float> arr_photons_atm_az(Natm);
    std::vector<float> arr_photons_atm_tau(Natm);
    std::vector<float> arr_photons_atm_phi(Natm);
    std::vector<int>   field_atm_photons_per_gridcell(n_volumes, 0);

    std::vector<int>   arr_photons_sfc_pos_idx(Nsfc);
    std::vector<float> arr_photons_sfc_pos_x(Nsfc);
    std::vector<float> arr_photons_sfc_pos_y(Nsfc);
    std::vector<float> arr_photons_sfc_pos_z(Nsfc, 0.0f);
    std::vector<float> arr_photons_sfc_mu(Nsfc);
    std::vector<float> arr_photons_sfc_az(Nsfc);
    std::vector<float> arr_photons_sfc_tau(Nsfc);
    std::vector<float> arr_photons_sfc_phi(Nsfc);
    std::vector<int>   field_sfc_photons_per_gridcell(n_tiles, 0);
    

    // Sampling photons
    if (INTERCELL_TECHNIQUE == "uniform")
    {
        if (INTRACELL_TECHNIQUE == "naive")
        {

            std::uniform_int_distribution<> random_atm_idx(0, (n_volumes - 1));
            std::uniform_int_distribution<> random_sfc_idx(0, (n_tiles - 1));

            // Atmosphere
            for (int idx_photon = 0; idx_photon < Natm; idx_photon++)
            {
                // Position
                int idx_atm = random_atm_idx(gen);
                arr_photons_atm_pos_idx[idx_photon] = idx_atm;

                int jktot = ktot * jtot;
                int idx_atm_z  = idx_atm / jktot;
                int idx_atm_2D = idx_atm % jktot;
                int idx_atm_y  = idx_atm_2D / ktot;
                int idx_atm_x  = idx_atm_2D % ktot;

                float random_shift_x = random_float(gen) * dx;
                float random_shift_y = random_float(gen) * dy;
                float random_shift_z = random_float(gen) * arr_dz[idx_atm_z];

                arr_photons_atm_pos_x[idx_photon] = idx_atm_x * dx + random_shift_x;
                arr_photons_atm_pos_y[idx_photon] = idx_atm_y * dy + random_shift_y;
                arr_photons_atm_pos_z[idx_photon] = arr_zh[idx_atm_z] + random_shift_z;

                // Angles and optical thickness
                arr_photons_atm_mu[idx_photon]  = 2*random_float(gen) - 1;
                arr_photons_atm_az[idx_photon]  = 2*cf::PI*random_float(gen);
                arr_photons_atm_tau[idx_photon] = -logf(random_float(gen));

                // Countint photons per gridcell
                field_atm_photons_per_gridcell[idx_atm] += 1;
            }

            // Surface
            for (int idx_photon = 0; idx_photon < Nsfc; idx_photon++)
            {
                // Position
                int idx_sfc = random_sfc_idx(gen);
                arr_photons_sfc_pos_idx[idx_photon] = idx_sfc;

                int idx_sfc_y  = idx_sfc / ktot;
                int idx_sfc_x  = idx_sfc % ktot;

                float random_shift_x = random_float(gen) * dx;
                float random_shift_y = random_float(gen) * dy;

                arr_photons_sfc_pos_x[idx_photon] = idx_sfc_x * dx + random_shift_x;
                arr_photons_sfc_pos_y[idx_photon] = idx_sfc_y * dy + random_shift_y;

                // Angles and optical thickness
                arr_photons_sfc_mu[idx_photon]  = std::sqrt(random_float(gen));
                arr_photons_sfc_az[idx_photon]  = 2*cf::PI*random_float(gen);
                arr_photons_sfc_tau[idx_photon] = -logf(random_float(gen));

                // Countint photons per gridcell
                field_sfc_photons_per_gridcell[idx_sfc] += 1;

            }


            // Determining carrying power of each photon
            for (int idx_photon = 0; idx_photon < Natm; idx_photon++)
            {
                int idx_atm = arr_photons_atm_pos_idx[idx_photon];
                arr_photons_atm_phi[idx_photon] = field_atm_phi[idx_atm] / field_atm_photons_per_gridcell[idx_atm];
            }
            for (int idx_photon = 0; idx_photon < Nsfc; idx_photon++)
            {
                int idx_sfc = arr_photons_sfc_pos_idx[idx_photon];
                arr_photons_sfc_phi[idx_photon] = field_sfc_phi[idx_sfc] / field_sfc_photons_per_gridcell[idx_sfc];
            }


        } else
        if (INTRACELL_TECHNIQUE == "margin")
        {
            ;
        } else
        if (INTRACELL_TECHNIQUE == "edge")
        {
            ;
        }
    } else
    if (INTERCELL_TECHNIQUE == "power")
    {
        if (INTRACELL_TECHNIQUE == "naive")
        {
            ;
        } else
        if (INTRACELL_TECHNIQUE == "margin")
        {
            ;
        } else
        if (INTRACELL_TECHNIQUE == "edge")
        {
            ;
        }
    } else 
    if (INTERCELL_TECHNIQUE == "power-gradient")
    {
        if (INTRACELL_TECHNIQUE == "naive")
        {
            ;
        } else
        if (INTRACELL_TECHNIQUE == "margin")
        {
            ;
        } else
        if (INTRACELL_TECHNIQUE == "edge")
        {
            ;
        }
    }



    /////////////////////// PHOTON RELEASE /////////////////////// 
    std::cout << "  MC: Photon release - atm" << std::endl;


    // Photons from the atmosphere
    photon_propagation(arr_photons_atm_pos_idx,
                       arr_photons_atm_pos_x,
                       arr_photons_atm_pos_y,
                       arr_photons_atm_pos_z,
                       arr_photons_atm_mu,
                       arr_photons_atm_az,
                       arr_photons_atm_tau,
                       arr_photons_atm_phi,
                       field_atm_kext,
                       arr_xh,
                       arr_yh,
                       arr_zh,
                       arr_x,
                       arr_y,
                       arr_z,
                       field_atm_absorbed,
                       field_sfc_absorbed,
                       x_max,
                       y_max,
                       z_max,
                       itot,
                       jtot,
                       ktot,
                       dx,
                       dy,
                       Natm,
                       0);
    
    
    std::cout << "  MC: Photon release - sfc" << std::endl;

    // Photons from the surface
    photon_propagation(arr_photons_sfc_pos_idx,
                       arr_photons_sfc_pos_x,
                       arr_photons_sfc_pos_y,
                       arr_photons_sfc_pos_z,
                       arr_photons_sfc_mu,
                       arr_photons_sfc_az,
                       arr_photons_sfc_tau,
                       arr_photons_sfc_phi,
                       field_atm_kext,
                       arr_xh,
                       arr_yh,
                       arr_zh,
                       arr_x,
                       arr_y,
                       arr_z,
                       field_atm_absorbed,
                       field_sfc_absorbed,
                       x_max,
                       y_max,
                       z_max,
                       itot,
                       jtot,
                       ktot,
                       dx,
                       dy,
                       Nsfc,
                       1);


    ///////////////// CALCULATING HEATING RATES /////////////////

    std::cout << "  MC: Calculating heating rates" << std::endl;

    std::vector<float> field_atm_heating_rates(n_volumes);
    std::vector<float> field_sfc_heating_rates(n_tiles);

    for (int i = 0; i < itot; i++)
    {
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm =  i * ktot * jtot + j * ktot + k;
                float atm_net_phi = field_atm_absorbed[idx_atm] - field_atm_emitted[idx_atm];
                float dz = arr_dz[i];
                field_atm_heating_rates[idx_atm] = atm_net_phi / (cf::RHO * cf::CP * dx * dy * dz) * 86400;

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    float sfc_net_phi = field_sfc_absorbed[idx_sfc] - field_sfc_emitted[idx_sfc];
                    field_sfc_heating_rates[idx_sfc] = sfc_net_phi / (cf::RHO * cf::CP * dx * dy) * 86400;
                }
            }
        }
    }



    


    return field_atm_heating_rates;
}