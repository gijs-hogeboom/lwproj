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
#include <atomic>

#include "util.h"
#include "precision.h"

    

void photon_propagation(const AliasTable& aliastable,
                        const std::vector<Real>& field_kext,
                        const std::vector<Real>& field_sfc_eps,
                        const std::vector<Real>& field_SSA,
                        const std::vector<Real>& field_ASY,
                        const std::vector<Real>& arr_xh,
                        const std::vector<Real>& arr_yh,
                        const std::vector<Real>& arr_zh,
                        const std::vector<Real>& arr_x,
                        const std::vector<Real>& arr_y,
                        const std::vector<Real>& arr_z,
                        const std::vector<Real>& arr_dz,
                        const std::vector<Real>& field_phi,
                        std::vector<Real>& field_atm_net_phi,
                        std::vector<Real>& field_sfc_net_phi,
                        std::vector<Real>& field_TOA_net_phi,
                        const long int N,
                        const int domain_section,
                        const std::string& INTERCELL_TECHNIQUE,
                        const bool Pesc_mode,
                        const bool enable_scattering)
{

    // Initializing parameters
    const int itot        = arr_z.size();
    const int jtot        = arr_y.size();
    const int ktot        = arr_x.size();
    const int n_volumes   = itot*jtot*ktot;
    const int n_tiles     = jtot*ktot;

    const Real x_max      = arr_xh[ktot];
    const Real y_max      = arr_yh[jtot];
    const Real z_max      = arr_zh[itot];
    const Real cell_dx    = arr_xh[1] - arr_xh[0];
    const Real cell_dy    = arr_yh[1] - arr_yh[0];
    
    const Real eps        = 1e-3;
    const Real tiny_eps   = 1e-12;

    const int iter_max    = 1e6;
    const Real dx_min     = tiny_eps;
    const Real dy_min     = tiny_eps;
    const Real dz_min     = std::sin(std::atan(z_max / static_cast<Real>(iter_max * cell_dx)));

    const Real w_crit     = 0.5;
    const Real Nreal      = static_cast<Real>(N);
    const int jktot       = n_tiles;

    const Real rng_smallest_val = 5.96046e-08;
    const Real tau_maxval       = 20.;

    const int Nphot_batch  = 1024;


    // Initializing progress bar
    const int num_progress_bar = 100;
    std::vector<std::int64_t> progress_bar(num_progress_bar);

    for (int i = 0; i < num_progress_bar; i++)
    {
        progress_bar[i] = (N * (i + 1)) / num_progress_bar;
    }
    int idx_progress = 0;

    std::cout << '|';
    for (int _ = 0; _ < num_progress_bar; _++) std::cout << '.' << std::flush;
    std::cout << '|' << '\b';
    for (int _ = 0; _ < num_progress_bar; _++) std::cout << '\b' << std::flush;

    std::atomic<std::int64_t> photons_done{0};

    

    // Parallelization
    #pragma omp parallel
    {
        // Initializing local rng and net phi vectors
        int tid = omp_get_thread_num();
        FastRNG rng_local(std::chrono::high_resolution_clock::now().time_since_epoch().count() + tid);

        std::vector<Real> field_atm_net_phi_local(n_volumes);
        std::vector<Real> field_sfc_net_phi_local(n_tiles);
        std::vector<Real> field_TOA_net_phi_local(n_tiles);

        // Initializing photon power
        Real photon_power;

        if (INTERCELL_TECHNIQUE == "power")
        {
            photon_power = std::accumulate(field_phi.begin(), field_phi.end(), 0.0) / N;
        }

        // Keeping track of current photon batch
        int photon_counter = 0;



        #pragma omp for schedule(dynamic, Nphot_batch)
        for (long int idx_photon = 0; idx_photon < N; idx_photon++)
        {

            // Tracking whether photon leaves the cell
            bool out_of_cell = false;
            if (domain_section == 1) out_of_cell = true; // Surface photons are always "out" of their cells
            
            // Sampling location within domain, determining photon power
            int idx_flat = aliastable.sample(rng_local);
            int idx_original = idx_flat; // storing starting position


            // Initializing position/direction/optical thickness
            int idx_z, idx_y, idx_x;
            Real x, y, z, mu, az, tau;

            // Initializing scattering-related properties
            Real absorbed_photon_power, current_ssa, current_asy;

            // Keeping track of cycles per photon
            std::uint64_t iter_counter = 0;




            ////////////////////////////////// SAMPLING PHOTON //////////////////////////////////
            if (domain_section == 0)
            {
                // Atmosphere
                idx_z  = idx_flat / jktot;
                int idx_2D = idx_flat % jktot;
                idx_y  = idx_2D / ktot;
                idx_x  = idx_2D % ktot;

                x = (idx_x + rng_local.uniform()) * cell_dx;
                y = (idx_y + rng_local.uniform()) * cell_dy;
                z = arr_zh[idx_z] + rng_local.uniform()*arr_dz[idx_z];

                mu = rng_local.uniform()*static_cast<Real>(2.) - static_cast<Real>(1.);
                az = rng_local.uniform()*static_cast<Real>(2.)*constants::PI;

                tau = std::min(-std::log(rng_local.uniform()), tau_maxval);
            }
            else if (domain_section == 1)
            {
                // Surface
                idx_z = 0;
                idx_y = idx_flat / ktot;
                idx_x = idx_flat % ktot;

                x = (idx_x + rng_local.uniform()) * cell_dx;
                y = (idx_y + rng_local.uniform()) * cell_dy;
                z = 0.;

                mu = std::sqrt(rng_local.uniform());
                az = rng_local.uniform()*static_cast<Real>(2.)*constants::PI;

                tau = std::min(-std::log(rng_local.uniform()), tau_maxval);
            }

            // Calculating cartesian direction vector
            Real s = std::sqrt(static_cast<Real>(1.) - mu*mu);
            Real dx = s*std::cos(az);
            Real dy = s*std::sin(az);
            Real dz = mu;

            dx = std::copysign(std::max(std::abs(dx), dx_min), dx);
            dy = std::copysign(std::max(std::abs(dy), dy_min), dy);
            dz = std::copysign(std::max(std::abs(dz), dz_min), dz);

            Real w = 1.;

            // Initializing working variables
            Real ds, max_s, tau_absorbed, fs, current_kext;
            Real dist_x, dist_y, dist_z;
            

            // Calculating photon power for uniform and power-gradient method
            if (!(INTERCELL_TECHNIQUE == "power"))
            {
                Real sample_weight = aliastable.weights[idx_flat];
                photon_power = field_phi[idx_flat] / (sample_weight * Nreal);
            }

            




            ////////////////////////////////// START OF PROPAGATION //////////////////////////////////
            while (w > static_cast<Real>(0.))
            {
                while (tau > static_cast<Real>(0.))
                {
                    iter_counter++;
                    if (iter_counter >= 1000000000)
                    {
                        // LOGvars(idx_photon);
                        // LOGvars(idx_x, idx_y, idx_z);
                        // LOGvars(x, y, z);
                        // LOGvars(dx, dy, dz);
                        // LOGvars(tau, w, idx_flat, current_kext*1e6, ds, max_s, tau_absorbed, fs, photon_power);
                        // LOGvars(dist_x, dist_y, dist_z);
                        w = 0.;
                        break;
                    }
                    // field boundary detection in x direction - wrapping
                    bool at_far_wall_x     = (std::abs(x - x_max) < eps);
                    bool going_forwards_x  = (dx >= static_cast<Real>(0.));
                    if (at_far_wall_x && going_forwards_x) 
                    { 
                        x = 0.;
                        idx_x = 0;
                    }
                    bool at_near_wall_x    = (std::abs(x) < eps);
                    bool going_backwards_x = (dx < static_cast<Real>(0.));
                    if (at_near_wall_x && going_backwards_x) 
                    { 
                        x = x_max; 
                        idx_x = ktot - 1;
                    }

                    // field boundary detection in y direction - wrapping
                    bool at_far_wall_y     = (std::abs(y - y_max) < eps);
                    bool going_forwards_y  = (dy >= static_cast<Real>(0.));
                    if (at_far_wall_y && going_forwards_y) 
                    { 
                        y = 0.; 
                        idx_y = 0;
                    }
                    bool at_near_wall_y    = (std::abs(y) < eps);
                    bool going_backwards_y = (dy < static_cast<Real>(0.));
                    if (at_near_wall_y && going_backwards_y) 
                    { 
                        y = y_max; 
                        idx_y = jtot - 1;
                    }


                    // field boundary detection in z direction - loss through TOA or absorbtion by surface
                    bool at_TOA            = (std::abs(z - z_max) < eps);
                    bool going_up          = (dz >= static_cast<Real>(0.));
                    if (at_TOA && going_up)
                    {
                        tau = 0.;
                        int idx_tile = idx_y * ktot + idx_x;
                        
                        Real absorbed_photon_power = w*photon_power;
                        w = 0.;

                        field_TOA_net_phi_local[idx_tile] += absorbed_photon_power;
                        if (domain_section == 0)
                        {
                            field_atm_net_phi_local[idx_original] -= absorbed_photon_power;
                        } else if (domain_section == 1)
                        {
                            field_sfc_net_phi_local[idx_original] -= absorbed_photon_power;
                        }
                        break;
                    }
                    bool at_surface        = (std::abs(z) < eps);
                    bool going_down        = (dz < static_cast<Real>(0.));
                    if (at_surface && going_down)
                    {
                        tau = 0.;
                        int idx_tile = idx_y * ktot + idx_x;

                        if (enable_scattering)
                        {
                            current_ssa = static_cast<Real>(1.) - field_sfc_eps[idx_tile]; // surface reflectivity is 1 - emissivity
                        }
                        else
                        {
                            current_ssa = 0.;
                        }
                        
                        absorbed_photon_power = (static_cast<Real>(1.) - current_ssa)*w*photon_power;
                        w *= current_ssa;

                        field_sfc_net_phi_local[idx_tile] += absorbed_photon_power;
                        if (domain_section == 0)
                        {
                            field_atm_net_phi_local[idx_original] -= absorbed_photon_power;
                        } else if (domain_section == 1)
                        {
                            field_sfc_net_phi_local[idx_original] -= absorbed_photon_power;
                        }
                        
                        if (enable_scattering)
                        {
                            if (w < w_crit)
                            {
                                Real rhow = rng_local.uniform();
                                if (rhow >= w)
                                {
                                    w = 0.;
                                    break;
                                } else {
                                    w = 1.;
                                    idx_z = 0;

                                    mu = std::sqrt(rng_local.uniform());
                                    az = rng_local.uniform()*static_cast<Real>(2.)*constants::PI;

                                    s = std::sqrt(static_cast<Real>(1.) - mu*mu);
                                    dx = s*std::cos(az);
                                    dy = s*std::sin(az);
                                    dz = mu;

                                    dx = std::copysign(std::max(std::abs(dx), dx_min), dx);
                                    dy = std::copysign(std::max(std::abs(dy), dy_min), dy);
                                    dz = std::copysign(std::max(std::abs(dz), dz_min), dz);

                                    tau = std::min(-std::log(rng_local.uniform()), tau_maxval);
                                    continue;
                                }
                            }
                            else
                            {
                                idx_z = 0;

                                mu = std::sqrt(rng_local.uniform());
                                az = rng_local.uniform()*static_cast<Real>(2.)*constants::PI;
                                
                                s = std::sqrt(static_cast<Real>(1.) - mu*mu);
                                dx = s*std::cos(az);
                                dy = s*std::sin(az);
                                dz = mu;

                                dx = std::copysign(std::max(std::abs(dx), dx_min), dx);
                                dy = std::copysign(std::max(std::abs(dy), dy_min), dy);
                                dz = std::copysign(std::max(std::abs(dz), dz_min), dz);

                                tau = std::min(-std::log(rng_local.uniform()), tau_maxval);
                                continue;
                            }
                        }
                        break;
                    }
                    


                    // Updating position
                    idx_flat = idx_z*jktot + idx_y*ktot + idx_x;

                    // Loading kext
                    Real current_kext = field_kext[idx_flat];

                    // Scanning collision with cell boundaries
                    Real time_x, time_y, time_z;
                    Real dnx, dny, dnz;

                    if (dx >= static_cast<Real>(0.)) // x
                    {
                        dnx = arr_xh[idx_x + 1] - x;
                    } else {
                        dnx = arr_xh[idx_x] - x;
                    }
                    time_x = dnx/dx;

                    if (dy >= static_cast<Real>(0.)) // y
                    {
                        dny = arr_yh[idx_y + 1] - y;
                    } else {
                        dny = arr_yh[idx_y] - y;
                    }
                    time_y = dny/dy;

                    if (dz >= static_cast<Real>(0.)) // z
                    {
                        dnz = arr_zh[idx_z + 1] - z;
                    } else {
                        dnz = arr_zh[idx_z] - z;
                    }
                    time_z = dnz/dz;


                    // Determinig the scaling factor based on which cell is hit (i.e., which direction takes the least amount of time)
                    bool hit_x_wall = ((time_x <= time_y) && (time_x <= time_z));
                    bool hit_y_wall = ((time_y <= time_x) && (time_y <= time_z));
                    bool hit_z_wall = ((time_z <= time_x) && (time_z <= time_y));

                    if (hit_x_wall)
                    {
                        dist_x = dnx;
                        dist_y = time_x * dy;
                        dist_z = time_x * dz;
                    }
                    else if (hit_y_wall)
                    {
                        dist_x = time_y * dx;
                        dist_y = dny;
                        dist_z = time_y * dz;
                    }
                    else if (hit_z_wall)
                    {
                        dist_x = time_z * dx;
                        dist_y = time_z * dy;
                        dist_z = dnz;
                    }


                    // Calculating distance traveled
                    Real ds           = std::sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
                    Real max_s        = tau/current_kext;
                    Real tau_absorbed = current_kext*ds;


                    if (Pesc_mode && !out_of_cell)
                    {
                        // Updating distance; traveling to the cell-border
                        x += dist_x;
                        y += dist_y;
                        z += dist_z;
                        
                        // Updating cell index
                        if (hit_x_wall)      {if (going_forwards_x) {idx_x += 1;} else {idx_x -= 1;}}
                        else if (hit_y_wall) {if (going_forwards_y) {idx_y += 1;} else {idx_y -= 1;}}
                        else if (hit_z_wall) {if (going_up)         {idx_z += 1;} else {idx_z -= 1;}}

                        out_of_cell = true;
                    }
                    else
                    {
                        if (ds < max_s)
                        {
                            // Photon transfers through the whole cell
                            tau -= tau_absorbed;
                            x += dist_x;
                            y += dist_y;
                            z += dist_z;

                            // Updating cell index
                            if (hit_x_wall)      {if (going_forwards_x) {idx_x += 1;} else {idx_x -= 1;}}
                            else if (hit_y_wall) {if (going_forwards_y) {idx_y += 1;} else {idx_y -= 1;}}
                            else if (hit_z_wall) {if (going_up)         {idx_z += 1;} else {idx_z -= 1;}}
                            
                            out_of_cell = true;
                        }
                        else
                        {
                            // Cell travels only partly through the cell and is absorbed/scattered
                            tau = 0.;

                            // Traveling fraction of a distance
                            Real fs = max_s / ds;
                            x += dist_x*fs;
                            y += dist_y*fs;
                            z += dist_z*fs;

                            // Retrieving scattering properties to update photon weight
                            if (enable_scattering)
                            {
                                current_ssa = field_SSA[idx_flat];
                                current_asy = field_ASY[idx_flat];
                            }
                            else
                            {
                                current_ssa = 0.;
                                current_asy = 0.;
                            }
                            
                            // Updating photon weight and carrying power
                            absorbed_photon_power = (static_cast<Real>(1.) - current_ssa)*w*photon_power;
                            w *= current_ssa;

                            // Updating absorbed/emitted power field
                            if (out_of_cell)
                            {
                                field_atm_net_phi_local[idx_flat] += absorbed_photon_power;
                                if (domain_section == 0)
                                {
                                    field_atm_net_phi_local[idx_original] -= absorbed_photon_power;
                                } else if (domain_section == 1)
                                {
                                    field_sfc_net_phi_local[idx_original] -= absorbed_photon_power;
                                }
                            }

                            // Scattering event
                            if (enable_scattering)
                            {
                                if (w < w_crit)
                                {
                                    Real rhow = rng_local.uniform();
                                    if (rhow >= w)
                                    {
                                        w = 0.;
                                        break;
                                    } else {
                                        w = 1.;
                                        Vec3 vec_new = generate_angle_HG(dx, dy, dz, current_asy, rng_local);
                                        dx = vec_new.x;
                                        dy = vec_new.y;
                                        dz = vec_new.z;

                                        dx = std::copysign(std::max(std::abs(dx), dx_min), dx);
                                        dy = std::copysign(std::max(std::abs(dy), dy_min), dy);
                                        dz = std::copysign(std::max(std::abs(dz), dz_min), dz);

                                        tau = std::min(-std::log(rng_local.uniform()), tau_maxval);
                                    }
                                }
                                else
                                {
                                    Vec3 vec_new = generate_angle_HG(dx, dy, dz, current_asy, rng_local);
                                    dx = vec_new.x;
                                    dy = vec_new.y;
                                    dz = vec_new.z;

                                    dx = std::copysign(std::max(std::abs(dx), dx_min), dx);
                                    dy = std::copysign(std::max(std::abs(dy), dy_min), dy);
                                    dz = std::copysign(std::max(std::abs(dz), dz_min), dz);

                                    tau = std::min(-std::log(rng_local.uniform()), tau_maxval);
                                }
                            }
                        }
                    }
                }
                if (tau <= static_cast<Real>(0.))
                {
                    if (w > static_cast<Real>(0.))
                    {
                        break; // caused by numerical error in 'tau -= tau_absorbed (= 0)'
                    }
                }
            }

            photon_counter++;

            if (photon_counter == Nphot_batch)
            {
                photon_counter = 0;
                
                photons_done.fetch_add(Nphot_batch, std::memory_order_relaxed);

                if (omp_get_thread_num() == 0)
                {
                    auto done = photons_done.load(std::memory_order_relaxed);

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

        // Adding remaining photons
        if (photon_counter > 0)
        {
            photons_done.fetch_add(photon_counter, std::memory_order_relaxed);
        }   

        // Combining the results of each thread
        #pragma omp critical
        {
            for (int i = 0; i < itot; i++)
            {
                for (int j = 0; j < jtot; j++)
                {
                    for (int k = 0; k < ktot; k++)
                    {
                        int idx = i*jktot + j*ktot + k;
                        field_atm_net_phi[idx] += field_atm_net_phi_local[idx];
                        if (i == 0)
                        {
                            int idx_tile = j*ktot + k;
                            field_sfc_net_phi[idx_tile] += field_sfc_net_phi_local[idx_tile];
                            field_TOA_net_phi[idx_tile] += field_TOA_net_phi_local[idx_tile];
                        }
                    }
                }
            }
        }
    }

    std::cout << "#|\n";

}


