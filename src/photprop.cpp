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

    

void photon_propagation(const AliasTable_float& aliastable,
                        const std::vector<float>& field_kext,
                        const std::vector<float>& field_sfc_eps,
                        const std::vector<float>& field_SSA,
                        const std::vector<float>& field_ASY,
                        const std::vector<float>& arr_xh,
                        const std::vector<float>& arr_yh,
                        const std::vector<float>& arr_zh,
                        const std::vector<float>& arr_x,
                        const std::vector<float>& arr_y,
                        const std::vector<float>& arr_z,
                        const std::vector<float>& arr_dz,
                        const std::vector<float>& field_phi,
                        std::vector<float>& field_atm_net_phi,
                        std::vector<float>& field_sfc_net_phi,
                        std::vector<float>& field_TOA_net_phi,
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

    const float x_max     = arr_xh[ktot];
    const float y_max     = arr_yh[jtot];
    const float z_max     = arr_zh[itot];
    const float cell_dx   = arr_xh[1] - arr_xh[0];
    const float cell_dy   = arr_yh[1] - arr_yh[0];
    
    const float eps       = 1e-3f;
    const float tiny_eps  = 1e-12f;
    const float w_crit    = 0.5f;
    const float Nfloat    = (float) N;
    const int jktot       = n_tiles;

    const float rng_smallest_val = 5.96046e-08f;
    const float mu_minval_atm    = -1.f + rng_smallest_val;
    const float mu_minval_sfc    = rng_smallest_val;
    const float tau_maxval       = 20.f;

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

        std::vector<float> field_atm_net_phi_local(n_volumes);
        std::vector<float> field_sfc_net_phi_local(n_tiles);
        std::vector<float> field_TOA_net_phi_local(n_tiles);

        // Initializing photon power
        float photon_power;

        if (INTERCELL_TECHNIQUE == "power")
        {
            photon_power = std::accumulate(field_phi.begin(), field_phi.end(), 0.0f) / N;
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


            // Calculating photon power for uniform and power-gradient method
            if (!(INTERCELL_TECHNIQUE == "power"))
            {
                float sample_weight = aliastable.weights[idx_flat];
                photon_power = field_phi[idx_flat] / (sample_weight * Nfloat);
            }


            // Initializing position/direction/optical thickness
            int idx_z, idx_y, idx_x;
            float x, y, z, mu, az, tau;

            // Initializing scattering-related properties
            float absorbed_photon_power, current_ssa, current_asy;

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

                mu = rng_local.uniform()*2.f - 1.f, mu_minval_atm;
                az = rng_local.uniform()*2.f*cfloat::PI;

                tau = fminf32(-logf(rng_local.uniform()), tau_maxval);
            }
            else if (domain_section == 1)
            {
                // Surface
                idx_z = 0;
                idx_y = idx_flat / ktot;
                idx_x = idx_flat % ktot;

                x = (idx_x + rng_local.uniform()) * cell_dx;
                y = (idx_y + rng_local.uniform()) * cell_dy;
                z = 0.f;

                mu = sqrtf(rng_local.uniform());
                az = rng_local.uniform()*2.f*cfloat::PI;

                tau = fminf32(-logf(rng_local.uniform()), tau_maxval);
            }

            // Calculating cartesian direction vector
            float s = sqrtf(1.f - mu*mu);
            float dx = s*cosf(az);
            float dy = s*sinf(az);
            float dz = mu;

            dx = copysignf(fmaxf(fabsf(dx), tiny_eps), dx);
            dy = copysignf(fmaxf(fabsf(dy), tiny_eps), dy);
            dz = copysignf(fmaxf(fabsf(dz), tiny_eps), dz);

            float w = 1.f;

            // Initializing working variables
            float ds, max_s, tau_absorbed, fs, current_kext;
            float dist_x, dist_y, dist_z;
            





            ////////////////////////////////// START OF PROPAGATION //////////////////////////////////
            while (w > 0.f)
            {
                while (tau > 0.f)
                {
                    iter_counter++;
                    if (iter_counter >= 1000000000)
                    {
                        LOGvars(idx_photon);
                        LOGvars(idx_x, idx_y, idx_z);
                        LOGvars(x, y, z);
                        LOGvars(dx, dy, dz);
                        LOGvars(tau, w, idx_flat, current_kext, ds, max_s, tau_absorbed, fs);
                        LOGvars(dist_x, dist_y, dist_z);
                    }
                    // field boundary detection in x direction - wrapping
                    bool at_far_wall_x     = (fabsf(x - x_max) < eps);
                    bool going_forwards_x  = (dx >= 0.f);
                    if (at_far_wall_x && going_forwards_x) 
                    { 
                        x = 0.f;
                        idx_x = 0;
                    }
                    bool at_near_wall_x    = (fabsf(x) < eps);
                    bool going_backwards_x = (dx < 0.f);
                    if (at_near_wall_x && going_backwards_x) 
                    { 
                        x = x_max; 
                        idx_x = ktot - 1;
                    }

                    // field boundary detection in y direction - wrapping
                    bool at_far_wall_y     = (fabsf(y - y_max) < eps);
                    bool going_forwards_y  = (dy >= 0.f);
                    if (at_far_wall_y && going_forwards_y) 
                    { 
                        y = 0.f; 
                        idx_y = 0;
                    }
                    bool at_near_wall_y    = (fabsf(y) < eps);
                    bool going_backwards_y = (dy < 0.f);
                    if (at_near_wall_y && going_backwards_y) 
                    { 
                        y = y_max; 
                        idx_y = jtot - 1;
                    }


                    // field boundary detection in z direction - loss through TOA or absorbtion by surface
                    bool at_TOA            = (fabsf(z - z_max) < eps);
                    bool going_up          = (dz >= 0.f);
                    if (at_TOA && going_up)
                    {
                        tau = 0.f;
                        int idx_tile = idx_y * ktot + idx_x;
                        
                        float absorbed_photon_power = w*photon_power;
                        w = 0.f;

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
                    bool at_surface        = (fabsf(z) < eps);
                    bool going_down        = (dz < 0.f);
                    if (at_surface && going_down)
                    {
                        tau = 0.f;
                        int idx_tile = idx_y * ktot + idx_x;

                        if (enable_scattering)
                        {
                            current_ssa = 1.f - field_sfc_eps[idx_tile]; // surface reflectivity is 1 - emissivity
                        }
                        else
                        {
                            current_ssa = 0.f;
                        }
                        
                        absorbed_photon_power = (1.f - current_ssa)*w*photon_power;
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
                                float rhow = rng_local.uniform();
                                if (rhow >= w)
                                {
                                    w = 0.f;
                                    break;
                                } else {
                                    w = 1.f;
                                    idx_z = 0;

                                    mu = sqrtf(rng_local.uniform());
                                    az = rng_local.uniform()*2.f*cfloat::PI;

                                    s = sqrtf(1.f - mu*mu);
                                    dx = s*cosf(az);
                                    dy = s*sinf(az);
                                    dz = mu;

                                    dx = copysignf(fmaxf(fabsf(dx), tiny_eps), dx);
                                    dy = copysignf(fmaxf(fabsf(dy), tiny_eps), dy);
                                    dz = copysignf(fmaxf(fabsf(dz), tiny_eps), dz);

                                    tau = fminf32(-logf(rng_local.uniform()), tau_maxval);
                                    continue;
                                }
                            }
                            else
                            {
                                idx_z = 0;

                                mu = sqrtf(rng_local.uniform());
                                az = rng_local.uniform()*2.f*cfloat::PI;
                                
                                s = sqrtf(1.f - mu*mu);
                                dx = s*cosf(az);
                                dy = s*sinf(az);
                                dz = mu;

                                dx = copysignf(fmaxf(fabsf(dx), tiny_eps), dx);
                                dy = copysignf(fmaxf(fabsf(dy), tiny_eps), dy);
                                dz = copysignf(fmaxf(fabsf(dz), tiny_eps), dz);

                                tau = fminf32(-logf(rng_local.uniform()), tau_maxval);
                                continue;
                            }
                        }
                        break;
                    }
                    
                    // Updating position
                    idx_flat = idx_z*jktot + idx_y*ktot + idx_x;

                    // Loading kext
                    float current_kext = field_kext[idx_flat];


                    if ((idx_flat < 0) | (idx_flat >= n_volumes)) {LOGvars(idx_photon); LOGvars(iter_counter, idx_x, idx_y, idx_z, tau, w, current_kext, dx, dy, dz, idx_flat);}


                    // Scanning collision with cell boundaries
                    float time_x, time_y, time_z;
                    float dnx, dny, dnz;

                    if (dx >= 0.f) // x
                    {
                        dnx = arr_xh[idx_x + 1] - x;
                    } else {
                        dnx = arr_xh[idx_x] - x;
                    }
                    time_x = dnx/dx;

                    if (dy >= 0.f) // y
                    {
                        dny = arr_yh[idx_y + 1] - y;
                    } else {
                        dny = arr_yh[idx_y] - y;
                    }
                    time_y = dny/dy;

                    if (dz >= 0.f) // z
                    {
                        dnz = arr_zh[idx_z + 1] - z;
                    } else {
                        dnz = arr_zh[idx_z] - z;
                    }
                    time_z = dnz/dz;


                    // Determinig the scaling factor based on which cell is hit (i.e., which direction takes the least amount of time)
                    // Additionally, updating photon index position for next iteration (if photon extincts within the cell, the idx will not be used anyways)

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
                    float ds           = sqrtf(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
                    float max_s        = tau/current_kext;
                    float tau_absorbed = current_kext*ds;


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
                            tau = 0.f;

                            // Traveling fraction of a distance
                            float fs = max_s / ds;
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
                                current_ssa = 0.f;
                                current_asy = 0.f;
                            }
                            
                            // Updating photon weight and carrying power
                            absorbed_photon_power = (1 - current_ssa)*w*photon_power;
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
                                    float rhow = rng_local.uniform();
                                    if (rhow >= w)
                                    {
                                        w = 0.f;
                                        break;
                                    } else {
                                        w = 1.f;
                                        Vec3float vec_new = generate_angle_HG_float(dx, dy, dz, current_asy, rng_local);
                                        dx = vec_new.x;
                                        dy = vec_new.y;
                                        dz = vec_new.z;

                                        dx = copysignf(fmaxf(fabsf(dx), tiny_eps), dx);
                                        dy = copysignf(fmaxf(fabsf(dy), tiny_eps), dy);
                                        dz = copysignf(fmaxf(fabsf(dz), tiny_eps), dz);

                                        tau = fminf32(-logf(rng_local.uniform()), tau_maxval);
                                    }
                                }
                                else
                                {
                                    Vec3float vec_new = generate_angle_HG_float(dx, dy, dz, current_asy, rng_local);
                                    dx = vec_new.x;
                                    dy = vec_new.y;
                                    dz = vec_new.z;

                                    dx = copysignf(fmaxf(fabsf(dx), tiny_eps), dx);
                                    dy = copysignf(fmaxf(fabsf(dy), tiny_eps), dy);
                                    dz = copysignf(fmaxf(fabsf(dz), tiny_eps), dz);

                                    tau = fminf32(-logf(rng_local.uniform()), tau_maxval);
                                }
                            }
                        }
                    }
                }
                if (tau <= 0.f)
                {
                    if (w > 0.f)
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

