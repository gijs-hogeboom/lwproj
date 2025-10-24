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




void photon_propagation(const std::vector<int>& arr_photons_pos_idx,
                        const std::vector<double>& arr_photons_pos_x,
                        const std::vector<double>& arr_photons_pos_y,
                        const std::vector<double>& arr_photons_pos_z,
                        const std::vector<double>& arr_photons_mu,
                        const std::vector<double>& arr_photons_az,
                        const std::vector<double>& arr_photons_tau,
                        const std::vector<double>& arr_photons_phi,
                        const std::vector<double>& arr_photons_phi_db,
                        const std::vector<double>& field_kext,
                        const std::vector<double>& arr_xh,
                        const std::vector<double>& arr_yh,
                        const std::vector<double>& arr_zh,
                        const std::vector<double>& arr_x,
                        const std::vector<double>& arr_y,
                        const std::vector<double>& arr_z,
                        std::vector<double>& field_atm_absorbed,
                        std::vector<double>& field_sfc_absorbed,
                        std::vector<double>& field_atm_absorbed_db,
                        std::vector<double>& field_sfc_absorbed_db,
                        double& TOA_absorbed,
                        double& TOA_absorbed_db,
                        double& tot_error,
                        std::vector<double>& arr_z34_phi,
                        double x_max,
                        double y_max,
                        double z_max,
                        int itot,
                        int jtot,
                        int ktot,
                        double dx,
                        double dy,
                        int N,
                        int domain_section,
                        std::vector<int>& M_photon_origin_counter)
{
    const double eps = 2e-5;
    const int jktot = jtot*ktot;
    
    const int Mstride = itot*jtot*ktot + jtot*ktot + 1; // (n_volumes + n_tiles + TOA)
    int sfc_offset = 0;
    if (domain_section == 0)
    {
        sfc_offset = jktot; // n_tiles;
    }

    

    u_int64_t photon_not_tracked_counter = 0;
    u_int64_t photon_is_tracked_counter = 0;

    // #pragma omp for schedule(dynamic)
    for (int idx_photon = 0; idx_photon < N; idx_photon++)
    {

        bool photon_is_tracked = false;
        
        // Loading initial photon variables
        int idx_flat = arr_photons_pos_idx[idx_photon];
        int idx_original = 0 + idx_flat; // storing starting position
        double x      = arr_photons_pos_x[idx_photon];
        double y      = arr_photons_pos_y[idx_photon];
        double z      = arr_photons_pos_z[idx_photon];
        double mu     = arr_photons_mu[idx_photon];
        double az     = arr_photons_az[idx_photon];
        double tau    = arr_photons_tau[idx_photon];

        // Calculating cartesian direction vector
        double s = std::sqrt(1 - mu*mu);
        double dx = s*std::cos(az);
        double dy = s*std::sin(az);
        double dz = mu;


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

        double ds_previous = 0;
        double x_prev = 0;
        double y_prev = 0;
        double distx_previous = 0;
        double disty_previous = 0;
        double dnx_prev = 0;
        double dny_prev = 0;
        
        int counter = 0;
        
        // Starting propegation...
        while (tau > 1e-10)
        {



            // if (x >= x_max) {std::cout << std::setprecision(10) << "x:\t" << x << '\t' << x_max << '\t' << x - x_max << '\t' << eps << '\t' << distx_previous << '\t' << x_prev << std::endl;}
            // if (y >= y_max) {std::cout << std::setprecision(10) << "y:\t" << y << '\t' << y_max << '\t' << y - y_max << '\t' << eps << '\t' << disty_previous << '\t' << y_prev << std::endl;}
            // if ((x <= eps) & (x != 0.)) {std::cout << std::setprecision(16) << "x:\t" << x << '\t' << eps << '\t' << distx_previous << '\t' << x_prev << '\t' << dnx_prev << '\t' << dx << std::endl;}
            // if ((y <= eps) & (y != 0.)) {std::cout << std::setprecision(16) << "y:\t" << y << '\t' << eps << '\t' << disty_previous << '\t' << y_prev << '\t' << dny_prev << '\t' << dy << std::endl;}


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

            // if (idx_x == -1) {std::cout << "x\t" << idx_x << '\t' << x << '\t' << eps << '\t' << ds_previous << std::endl;}
            // if (idx_y == -1) {std::cout << "y\t" << idx_y << '\t' << y << '\t' << eps << '\t' << ds_previous << std::endl;}
            // if (idx_x == 1) {std::cout << "x\t" << idx_x << '\t' << x << '\t' << eps << '\t' << ds_previous << std::endl;}
            // if (idx_y == 1) {std::cout << "y\t" << idx_y << '\t' << y << '\t' << eps << '\t' << ds_previous << std::endl;}


            // field boundary detection in z direction - loss through TOA or absorbtion by surface
            bool at_TOA            = (std::abs(z - z_max) < eps);
            bool going_up          = (dz >= 0.);
            if (at_TOA && going_up)
            {
                tau = 0.;

                // #pragma omp atomic
                TOA_absorbed += arr_photons_phi[idx_photon];
                // TOA_absorbed_db += static_cast<double>(arr_photons_phi[idx_photon]);
                TOA_absorbed_db += arr_photons_phi_db[idx_photon];

                int idx_M = itot + jktot + Mstride*(sfc_offset + idx_original);
                M_photon_origin_counter[idx_M]++;
                tot_error += estimate_max_numerical_double_error(TOA_absorbed);
                photon_is_tracked = true;
                break;
            }
            bool at_surface        = (std::abs(z) < eps);
            bool going_down        = (dz < 0.);
            if (at_surface && going_down)
            {
                tau = 0.;
                int idx_sfc = idx_y * ktot + idx_x;

                // #pragma omp atomic
                field_sfc_absorbed[idx_sfc] += arr_photons_phi[idx_photon];
                // field_sfc_absorbed_db[idx_sfc] += static_cast<double>(arr_photons_phi[idx_photon]);
                field_sfc_absorbed_db[idx_sfc] += arr_photons_phi_db[idx_photon];

                int idx_M = idx_sfc + Mstride*(sfc_offset + idx_original);
                // std::cout << domain_section << ',' << "absorbed by sfc" << ',' << idx_original << ',' << idx_sfc << ',' << idx_M << std::endl;
                M_photon_origin_counter[idx_M]++;
                tot_error += estimate_max_numerical_double_error(field_sfc_absorbed[idx_sfc]);
                photon_is_tracked = true;
                break;
            }

        
            // retrieving kext from flattened grid
            int idx_flat = idx_z * jktot + idx_y * ktot + idx_x;
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
            dnx_prev = dnx;

            if (dy >= 0.) // y
            {
                dny = arr_yh[idx_y + 1] - y;
            } else {
                dny = arr_yh[idx_y] - y;
            }
            time_y = dny/dy;
            dny_prev = dny;

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

            
            double ds = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
            double max_s = tau/current_kext;
            double tau_absorbed = current_kext*ds;

            ds_previous = ds;
            distx_previous = dist_x;
            disty_previous = dist_y;

            x_prev = x;
            y_prev = y;
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
                double fs = max_s / ds;
                x += dist_x*fs;
                y += dist_y*fs;
                z += dist_z*fs;

                // #pragma omp atomic
                field_atm_absorbed[idx_flat] += arr_photons_phi[idx_photon];
                // field_atm_absorbed_db[idx_flat] += static_cast<double>(arr_photons_phi[idx_photon]);
                field_atm_absorbed_db[idx_flat] += arr_photons_phi_db[idx_photon];

                int idx_M = jktot + idx_flat + Mstride*(sfc_offset + idx_original);
                M_photon_origin_counter[idx_M]++;

                photon_is_tracked = true;

                // std::cout << domain_section << ',' << "absorbed by atm" << ',' << idx_original << ',' << idx_flat << ',' << idx_M << ',' <<  sfc_offset << std::endl;

                tot_error += estimate_max_numerical_double_error(field_atm_absorbed[idx_flat]);

                if (idx_flat == 34) {arr_z34_phi.push_back(static_cast<double>(arr_photons_phi[idx_photon]));}
            }
            counter++;
        }

        if (photon_is_tracked == true)
        {
            photon_is_tracked_counter++;
        } else {
            photon_not_tracked_counter++;
        }
    }

    if (photon_not_tracked_counter > 0)
    {
        std::cout << "!!!!!!!!!! WARNING !!!!!!!!!    " << photon_not_tracked_counter << " photons are not tracked!! " << std::endl;
    }
}


std::vector<double> run_MC(const std::vector<double>& arr_z,
                          const std::vector<double>& arr_zh,
                          const std::vector<double>& arr_dz,
                          const std::vector<double>& arr_kext,
                          const std::vector<double>& arr_Batm,
                          double Bsfc,
                          double dx,
                          double dy,
                          int ktot,
                          int jtot,
                          int itot,
                          const std::string& INTERCELL_TECHNIQUE,
                          const std::string& INTRACELL_TECHNIQUE,
                          const std::string& CASE,
                          int Natm,
                          int Nsfc,
                          bool print_EB,
                          bool verbose)
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
    std::vector<double> field_atm_kext(n_volumes);
    std::vector<double> field_atm_B(n_volumes);
    std::vector<double> field_atm_phi(n_volumes);
    std::vector<double> field_sfc_B(n_tiles);
    std::vector<double> field_sfc_phi(n_tiles);
    std::vector<double> field_sfc_eps(n_tiles, 1.0);

    std::vector<double> field_atm_absorbed(n_volumes);
    std::vector<double> field_atm_emitted(n_volumes);
    std::vector<double> field_sfc_absorbed(n_tiles);
    std::vector<double> field_sfc_emitted(n_tiles);

    std::vector<double> field_atm_absorbed_db(n_volumes);
    std::vector<double> field_sfc_absorbed_db(n_tiles);


    // Generating fields
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < itot; i++)
    {
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm = i * ktot * jtot + j * ktot + k;
                double current_kext = arr_kext[i];
                double current_Batm = arr_Batm[i];

                field_atm_kext[idx_atm] = current_kext;
                field_atm_B[idx_atm] = current_Batm;
                field_atm_phi[idx_atm] = 4*cdouble::PI * current_kext * current_Batm * dx * dy * arr_dz[i];
                
                field_atm_emitted[idx_atm] = field_atm_phi[idx_atm]; // basically copying the power array

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    field_sfc_B[idx_sfc] = Bsfc;
                    field_sfc_phi[idx_sfc] = cdouble::PI * field_sfc_eps[idx_sfc] * Bsfc * dx * dy;

                    field_sfc_emitted[idx_sfc] = field_sfc_phi[idx_sfc]; // idem
                }
            }
        }
    }

    /////////////////////////// SAMPLING /////////////////////////// 
    
    if (verbose)
    {
        std::cout << "  MC: Sampling" << std::endl;
    }
    
    // Initializing photon arrays
    std::vector<int>   arr_photons_atm_pos_idx(Natm);
    std::vector<double> arr_photons_atm_pos_x(Natm);
    std::vector<double> arr_photons_atm_pos_y(Natm);
    std::vector<double> arr_photons_atm_pos_z(Natm);
    std::vector<double> arr_photons_atm_mu(Natm);
    std::vector<double> arr_photons_atm_az(Natm);
    std::vector<double> arr_photons_atm_tau(Natm);
    std::vector<double> arr_photons_atm_phi(Natm);
    std::vector<double> arr_photons_atm_phi_db(Natm);
    std::vector<double> field_atm_phi_per_photon(n_volumes);
    std::vector<int>   field_atm_photons_per_gridcell(n_volumes, 0);

    std::vector<int>   arr_photons_sfc_pos_idx(Nsfc);
    std::vector<double> arr_photons_sfc_pos_x(Nsfc);
    std::vector<double> arr_photons_sfc_pos_y(Nsfc);
    std::vector<double> arr_photons_sfc_pos_z(Nsfc, 0.0);
    std::vector<double> arr_photons_sfc_mu(Nsfc);
    std::vector<double> arr_photons_sfc_az(Nsfc);
    std::vector<double> arr_photons_sfc_tau(Nsfc);
    std::vector<double> arr_photons_sfc_phi(Nsfc);
    std::vector<double> arr_photons_sfc_phi_db(Nsfc);
    std::vector<double> field_sfc_phi_per_photon(n_tiles);
    std::vector<int>   field_sfc_photons_per_gridcell(n_tiles, 0);

    double TOA_absorbed = 0.0; // This variable keeps track of photons lost to TOA
    double TOA_absorbed_db = 0.0;

    // Sampling photons
    if (INTERCELL_TECHNIQUE == "uniform")
    {
        if (INTRACELL_TECHNIQUE == "naive")
        {

            // Atmosphere
            for (size_t idx_photon = 0; idx_photon < Natm; idx_photon++)
            {
                // Position
                int idx_atm = rng.uniform_int(n_volumes - 1);
                arr_photons_atm_pos_idx[idx_photon] = idx_atm;

                int jktot = ktot * jtot;
                int idx_atm_z  = idx_atm / jktot;
                int idx_atm_2D = idx_atm % jktot;
                int idx_atm_y  = idx_atm_2D / ktot;
                int idx_atm_x  = idx_atm_2D % ktot;

                double random_shift_x = rng.uniform() * dx;
                double random_shift_y = rng.uniform() * dy;
                double random_shift_z = rng.uniform() * arr_dz[idx_atm_z];

                arr_photons_atm_pos_x[idx_photon] = idx_atm_x * dx + random_shift_x;
                arr_photons_atm_pos_y[idx_photon] = idx_atm_y * dy + random_shift_y;
                arr_photons_atm_pos_z[idx_photon] = arr_zh[idx_atm_z] + random_shift_z;

                // Angles and optical thickness
                arr_photons_atm_mu[idx_photon]  = 2*rng.uniform() - 1;
                arr_photons_atm_az[idx_photon]  = 2*cdouble::PI*rng.uniform();
                arr_photons_atm_tau[idx_photon] = -std::log(rng.uniform());

                // Countint photons per gridcell
                field_atm_photons_per_gridcell[idx_atm] += 1;

            }


            // Surface
            for (size_t idx_photon = 0; idx_photon < Nsfc; idx_photon++)
            {
                // Position
                int idx_sfc = rng.uniform_int(n_tiles - 1);
                arr_photons_sfc_pos_idx[idx_photon] = idx_sfc;

                int idx_sfc_y  = idx_sfc / ktot;
                int idx_sfc_x  = idx_sfc % ktot;

                double random_shift_x = rng.uniform() * dx;
                double random_shift_y = rng.uniform() * dy;

                arr_photons_sfc_pos_x[idx_photon] = idx_sfc_x * dx + random_shift_x;
                arr_photons_sfc_pos_y[idx_photon] = idx_sfc_y * dy + random_shift_y;

                // Angles and optical thickness
                arr_photons_sfc_mu[idx_photon]  = std::sqrt(rng.uniform());
                arr_photons_sfc_az[idx_photon]  = 2*cdouble::PI*rng.uniform();
                arr_photons_sfc_tau[idx_photon] = -std::log(rng.uniform());

                // Countint photons per gridcell
                field_sfc_photons_per_gridcell[idx_sfc] += 1;
                

            }

            // Calculating the power per photon for each cell in the field
            for (int idx_atm = 0; idx_atm < n_volumes; idx_atm++)
            {
                field_atm_phi_per_photon[idx_atm] = field_atm_phi[idx_atm] / field_atm_photons_per_gridcell[idx_atm];
            }
            for (int idx_sfc = 0; idx_sfc < n_tiles; idx_sfc++)
            {
                field_sfc_phi_per_photon[idx_sfc] = field_sfc_phi[idx_sfc] / field_sfc_photons_per_gridcell[idx_sfc];
            }


            // Determining carrying power of each photon
            for (size_t idx_photon = 0; idx_photon < Natm; idx_photon++)
            {
                int idx_atm = arr_photons_atm_pos_idx[idx_photon];
                arr_photons_atm_phi[idx_photon] = field_atm_phi_per_photon[idx_atm];
            }
            for (size_t idx_photon = 0; idx_photon < Nsfc; idx_photon++)
            {
                int idx_sfc = arr_photons_sfc_pos_idx[idx_photon];
                arr_photons_sfc_phi[idx_photon] = field_sfc_phi_per_photon[idx_sfc];
            }


            // Compensating power per photon by making sure the sum of all photon's phi = total emitted phi
            double sum_field_atm_phi = kahan_sum(field_atm_phi);
            double sum_field_sfc_phi = kahan_sum(field_sfc_phi);
            double sum_arr_photons_atm_phi = kahan_sum(arr_photons_atm_phi);
            double sum_arr_photons_sfc_phi = kahan_sum(arr_photons_sfc_phi);
            for (size_t idx_photon = 0; idx_photon < Natm; idx_photon++)
            {
                arr_photons_atm_phi[idx_photon] *= (sum_field_atm_phi / sum_arr_photons_atm_phi);
            }
            for (size_t idx_photon = 0; idx_photon < Nsfc; idx_photon++)
            {
                arr_photons_sfc_phi[idx_photon] *= (sum_field_sfc_phi / sum_arr_photons_sfc_phi);
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
            // Atmosphere
            // Generating weighted choice

            AliasTable_double alias_atm(field_atm_phi);

            #pragma omp parallel
            {
                FastRNG rng_local(omp_get_thread_num() + std::chrono::high_resolution_clock::now().time_since_epoch().count());
                #pragma omp for
                for (size_t idx_photon = 0; idx_photon < Natm; idx_photon++)
                {

                    // Position
                    arr_photons_atm_pos_idx[idx_photon] = alias_atm.sample(rng_local);
                    int idx_atm = arr_photons_atm_pos_idx[idx_photon];
                                    
                    int jktot = ktot * jtot;
                    int idx_atm_z  = idx_atm / jktot;
                    int idx_atm_2D = idx_atm % jktot;
                    int idx_atm_y  = idx_atm_2D / ktot;
                    int idx_atm_x  = idx_atm_2D % ktot;

                    double random_shift_x = rng_local.uniform() * dx;
                    double random_shift_y = rng_local.uniform() * dy;
                    double random_shift_z = rng_local.uniform() * arr_dz[idx_atm_z];

                    arr_photons_atm_pos_x[idx_photon] = idx_atm_x * dx + random_shift_x;
                    arr_photons_atm_pos_y[idx_photon] = idx_atm_y * dy + random_shift_y;
                    arr_photons_atm_pos_z[idx_photon] = arr_zh[idx_atm_z] + random_shift_z;

                    // Angles and optical thickness
                    arr_photons_atm_mu[idx_photon]  = 2*rng_local.uniform() - 1;
                    arr_photons_atm_az[idx_photon]  = 2*cdouble::PI*rng_local.uniform();
                    arr_photons_atm_tau[idx_photon] = -std::log(rng_local.uniform());

                    // Countint photons per gridcell
                    #pragma omp atomic
                    field_atm_photons_per_gridcell[idx_atm] += 1;
                }
            }
            
            // Surface
            AliasTable_double alias_sfc(field_sfc_phi);
            #pragma omp parallel
            {
                FastRNG rng_local(omp_get_thread_num() + std::chrono::high_resolution_clock::now().time_since_epoch().count());
                #pragma omp for
                for (size_t idx_photon = 0; idx_photon < Nsfc; idx_photon++)
                {
                    // Position
                    arr_photons_sfc_pos_idx[idx_photon] = alias_sfc.sample(rng_local);
                    int idx_sfc = arr_photons_sfc_pos_idx[idx_photon];

                    int idx_sfc_y  = idx_sfc / ktot;
                    int idx_sfc_x  = idx_sfc % ktot;

                    double random_shift_x = rng_local.uniform() * dx;
                    double random_shift_y = rng_local.uniform() * dy;

                    arr_photons_sfc_pos_x[idx_photon] = idx_sfc_x * dx + random_shift_x;
                    arr_photons_sfc_pos_y[idx_photon] = idx_sfc_y * dy + random_shift_y;

                    // Angles and optical thickness
                    arr_photons_sfc_mu[idx_photon]  = std::sqrt(rng_local.uniform());
                    arr_photons_sfc_az[idx_photon]  = 2*cdouble::PI*rng_local.uniform();
                    arr_photons_sfc_tau[idx_photon] = -std::log(rng_local.uniform());

                    // Countint photons per gridcell
                    #pragma omp atomic
                    field_sfc_photons_per_gridcell[idx_sfc] += 1;
                }
            }
            


            // Calculating the power per photon for each cell in the field
            for (int idx_atm = 0; idx_atm < n_volumes; idx_atm++)
            {
                field_atm_phi_per_photon[idx_atm] = field_atm_phi[idx_atm] / field_atm_photons_per_gridcell[idx_atm];
            }
            for (int idx_sfc = 0; idx_sfc < n_tiles; idx_sfc++)
            {
                field_sfc_phi_per_photon[idx_sfc] = field_sfc_phi[idx_sfc] / field_sfc_photons_per_gridcell[idx_sfc];
            }


            // Determining carrying power of each photon
            for (size_t idx_photon = 0; idx_photon < Natm; idx_photon++)
            {
                int idx_atm = arr_photons_atm_pos_idx[idx_photon];
                arr_photons_atm_phi[idx_photon] = field_atm_phi_per_photon[idx_atm];
                arr_photons_atm_phi_db[idx_photon] = static_cast<double>(field_atm_phi[idx_atm]) / field_atm_photons_per_gridcell[idx_atm];
            }
            for (size_t idx_photon = 0; idx_photon < Nsfc; idx_photon++)
            {
                int idx_sfc = arr_photons_sfc_pos_idx[idx_photon];
                arr_photons_sfc_phi[idx_photon] = field_sfc_phi_per_photon[idx_sfc];
                arr_photons_sfc_phi_db[idx_photon] = static_cast<double>(field_sfc_phi[idx_sfc]) / field_sfc_photons_per_gridcell[idx_sfc];
            }


            // Compensating power per photon by making sure the sum of all photon's phi = total emitted phi
            double sum_field_atm_phi = kahan_sum(field_atm_phi);
            double sum_field_sfc_phi = kahan_sum(field_sfc_phi);
            double sum_arr_photons_atm_phi = kahan_sum(arr_photons_atm_phi);
            double sum_arr_photons_sfc_phi = kahan_sum(arr_photons_sfc_phi);
            double sum_arr_photons_atm_phi_db = kahan_sum(arr_photons_atm_phi_db);
            double sum_arr_photons_sfc_phi_db = kahan_sum(arr_photons_sfc_phi_db);
            for (size_t idx_photon = 0; idx_photon < Natm; idx_photon++)
            {
                arr_photons_atm_phi[idx_photon]    *= (sum_field_atm_phi / sum_arr_photons_atm_phi);
                arr_photons_atm_phi_db[idx_photon] *= (static_cast<double>(sum_field_atm_phi) / sum_arr_photons_atm_phi_db);
            }
            for (size_t idx_photon = 0; idx_photon < Nsfc; idx_photon++)
            {
                arr_photons_sfc_phi[idx_photon] *= (sum_field_sfc_phi / sum_arr_photons_sfc_phi);
                arr_photons_sfc_phi_db[idx_photon] *= (static_cast<double>(sum_field_sfc_phi) / sum_arr_photons_sfc_phi_db);
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
    if (INTERCELL_TECHNIQUE == "power-gradient")
    {
        if (INTRACELL_TECHNIQUE == "naive")
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
            AliasTable_double alias_atm(field_atm_phi_gradient);

            #pragma omp parallel
            {
                FastRNG rng_local(omp_get_thread_num() + std::chrono::high_resolution_clock::now().time_since_epoch().count());
                #pragma omp for
                for (size_t idx_photon = 0; idx_photon < Natm; idx_photon++)
                {

                    // Position
                    arr_photons_atm_pos_idx[idx_photon] = alias_atm.sample(rng_local);
                    int idx_atm = arr_photons_atm_pos_idx[idx_photon];
                                    
                    int jktot = ktot * jtot;
                    int idx_atm_z  = idx_atm / jktot;
                    int idx_atm_2D = idx_atm % jktot;
                    int idx_atm_y  = idx_atm_2D / ktot;
                    int idx_atm_x  = idx_atm_2D % ktot;

                    double random_shift_x = rng_local.uniform() * dx;
                    double random_shift_y = rng_local.uniform() * dy;
                    double random_shift_z = rng_local.uniform() * arr_dz[idx_atm_z];

                    arr_photons_atm_pos_x[idx_photon] = idx_atm_x * dx + random_shift_x;
                    arr_photons_atm_pos_y[idx_photon] = idx_atm_y * dy + random_shift_y;
                    arr_photons_atm_pos_z[idx_photon] = arr_zh[idx_atm_z] + random_shift_z;

                    // Angles and optical thickness
                    arr_photons_atm_mu[idx_photon]  = 2*rng_local.uniform() - 1;
                    arr_photons_atm_az[idx_photon]  = 2*cdouble::PI*rng_local.uniform();
                    arr_photons_atm_tau[idx_photon] = -std::log(rng_local.uniform());

                    // Countint photons per gridcell
                    #pragma omp atomic
                    field_atm_photons_per_gridcell[idx_atm] += 1;
                }
            }
            
            // Surface
            AliasTable_double alias_sfc(field_sfc_phi_gradient);
            #pragma omp parallel
            {
                FastRNG rng_local(omp_get_thread_num() + std::chrono::high_resolution_clock::now().time_since_epoch().count());
                #pragma omp for
                for (size_t idx_photon = 0; idx_photon < Nsfc; idx_photon++)
                {
                    // Position
                    arr_photons_sfc_pos_idx[idx_photon] = alias_sfc.sample(rng_local);
                    int idx_sfc = arr_photons_sfc_pos_idx[idx_photon];

                    int idx_sfc_y  = idx_sfc / ktot;
                    int idx_sfc_x  = idx_sfc % ktot;

                    double random_shift_x = rng_local.uniform() * dx;
                    double random_shift_y = rng_local.uniform() * dy;

                    arr_photons_sfc_pos_x[idx_photon] = idx_sfc_x * dx + random_shift_x;
                    arr_photons_sfc_pos_y[idx_photon] = idx_sfc_y * dy + random_shift_y;

                    // Angles and optical thickness
                    arr_photons_sfc_mu[idx_photon]  = std::sqrt(rng_local.uniform());
                    arr_photons_sfc_az[idx_photon]  = 2*cdouble::PI*rng_local.uniform();
                    arr_photons_sfc_tau[idx_photon] = -std::log(rng_local.uniform());

                    // Countint photons per gridcell
                    #pragma omp atomic
                    field_sfc_photons_per_gridcell[idx_sfc] += 1;
                }
            }
            


            // Calculating the power per photon for each cell in the field
            for (int idx_atm = 0; idx_atm < n_volumes; idx_atm++)
            {
                field_atm_phi_per_photon[idx_atm] = field_atm_phi[idx_atm] / field_atm_photons_per_gridcell[idx_atm];
            }
            for (int idx_sfc = 0; idx_sfc < n_tiles; idx_sfc++)
            {
                field_sfc_phi_per_photon[idx_sfc] = field_sfc_phi[idx_sfc] / field_sfc_photons_per_gridcell[idx_sfc];
            }


            // Determining carrying power of each photon
            for (size_t idx_photon = 0; idx_photon < Natm; idx_photon++)
            {
                int idx_atm = arr_photons_atm_pos_idx[idx_photon];
                arr_photons_atm_phi[idx_photon] = field_atm_phi_per_photon[idx_atm];
            }
            for (size_t idx_photon = 0; idx_photon < Nsfc; idx_photon++)
            {
                int idx_sfc = arr_photons_sfc_pos_idx[idx_photon];
                arr_photons_sfc_phi[idx_photon] = field_sfc_phi_per_photon[idx_sfc];
            }


            // Compensating power per photon by making sure the sum of all photon's phi = total emitted phi
            double sum_field_atm_phi = kahan_sum(field_atm_phi);
            double sum_field_sfc_phi = kahan_sum(field_sfc_phi);
            double sum_arr_photons_atm_phi = kahan_sum(arr_photons_atm_phi);
            double sum_arr_photons_sfc_phi = kahan_sum(arr_photons_sfc_phi);
            for (size_t idx_photon = 0; idx_photon < Natm; idx_photon++)
            {
                arr_photons_atm_phi[idx_photon] *= (sum_field_atm_phi / sum_arr_photons_atm_phi);
            }
            for (size_t idx_photon = 0; idx_photon < Nsfc; idx_photon++)
            {
                arr_photons_sfc_phi[idx_photon] *= (sum_field_sfc_phi / sum_arr_photons_sfc_phi);
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
    }


    // std::cout << "SFC: total emited power:  " << kahan_sum(field_sfc_phi) << std::endl;
    // std::cout << "SFC: sum of photon power: " << kahan_sum(arr_photons_sfc_phi) << std::endl;
    // std::cout << "ATM: total emited power:  " << kahan_sum(field_atm_phi) << std::endl;
    // std::cout << "ATM: sum of photon power: " << kahan_sum(arr_photons_atm_phi) << std::endl;

    /////////////////////// PHOTON RELEASE /////////////////////// 
    if (verbose) 
    {
        std::cout << "  MC: Photon release - atm" << std::endl;
    }

    double tot_error = 0.;

    double dummy_photon_atm_phi = (kahan_sum(field_atm_phi))/(Natm);
    double dummy_photon_sfc_phi = (kahan_sum(field_sfc_phi))/(Nsfc);
    std::vector<double> arr_z34_phi;
    arr_z34_phi.reserve(Natm);



    int n_tot = n_volumes + n_tiles + 1; // sfc + atm + TOA

    std::vector<int> M_photon_origin_counter(n_tot*n_tot);

    // Photons from the atmosphere
    photon_propagation(arr_photons_atm_pos_idx,
                       arr_photons_atm_pos_x,
                       arr_photons_atm_pos_y,
                       arr_photons_atm_pos_z,
                       arr_photons_atm_mu,
                       arr_photons_atm_az,
                       arr_photons_atm_tau,
                       arr_photons_atm_phi,
                       arr_photons_atm_phi_db,
                       field_atm_kext,
                       arr_xh,
                       arr_yh,
                       arr_zh,
                       arr_x,
                       arr_y,
                       arr_z,
                       field_atm_absorbed,
                       field_sfc_absorbed,
                       field_atm_absorbed_db,
                       field_sfc_absorbed_db,
                       TOA_absorbed,
                       TOA_absorbed_db,
                       tot_error,
                       arr_z34_phi,
                       x_max,
                       y_max,
                       z_max,
                       itot,
                       jtot,
                       ktot,
                       dx,
                       dy,
                       Natm,
                       0,
                       M_photon_origin_counter);
    
    if (verbose)
    {
        std::cout << "  MC: Photon release - sfc" << std::endl;
    }
    // Photons from the surface
    photon_propagation(arr_photons_sfc_pos_idx,
                       arr_photons_sfc_pos_x,
                       arr_photons_sfc_pos_y,
                       arr_photons_sfc_pos_z,
                       arr_photons_sfc_mu,
                       arr_photons_sfc_az,
                       arr_photons_sfc_tau,
                       arr_photons_sfc_phi,
                       arr_photons_sfc_phi_db,
                       field_atm_kext,
                       arr_xh,
                       arr_yh,
                       arr_zh,
                       arr_x,
                       arr_y,
                       arr_z,
                       field_atm_absorbed,
                       field_sfc_absorbed,
                       field_atm_absorbed_db,
                       field_sfc_absorbed_db,
                       TOA_absorbed,
                       TOA_absorbed_db,
                       tot_error,
                       arr_z34_phi,
                       x_max,
                       y_max,
                       z_max,
                       itot,
                       jtot,
                       ktot,
                       dx,
                       dy,
                       Nsfc,
                       1,
                       M_photon_origin_counter);

    


    // Storing photon counter matrix
    std::ofstream fout_Mphot("/home/gijs-hogeboom/dev/lwproj/data_output/photon_matrix/photon_counter_" + CASE + 
                                 "_" + INTERCELL_TECHNIQUE + "_Natm" + std::to_string(static_cast<int>(std::log2(Natm))) + "_Nsfc" + 
                                 std::to_string(static_cast<int>(std::log2(Nsfc))) + ".dat");
    if (!fout_Mphot.is_open())
    {
        std::cerr << "Error, cannot open fout_Mphot!" << std::endl;
    }

    std::vector<std::string> counter_axis(n_tot);
    // Write axis
    for (int iaxis = 0; iaxis < n_tot; iaxis++)
    {
        if (iaxis < n_tiles)
        {
            counter_axis[iaxis] = "sfc" + std::to_string(iaxis);
        }
        else if ((iaxis >= n_tiles) && (iaxis < (n_tiles + n_volumes)))
        {
            counter_axis[iaxis] = "atm" + std::to_string(iaxis - n_tiles);
        }
        else
        {
            counter_axis[iaxis] = "TOA";
        }
    }

    // Writing header
    fout_Mphot << "received from,";
    for (int icol = 0; icol < n_tot; icol++)
    {
        if (icol < (n_tot-1))
        {
            fout_Mphot << counter_axis[icol] << ',';
        }
        else
        {
            fout_Mphot << counter_axis[icol] << std::endl;
        }
    }

    // Writing body
    for (int irow = 0; irow < n_tot; irow++)
    {
        fout_Mphot << counter_axis[irow] << ',';

        for (int icol = 0; icol < n_tot; icol++)
        {
            int idx_M = irow*n_tot + icol;
            if (icol < (n_tot-1))
            {
                fout_Mphot << M_photon_origin_counter[idx_M] << ',';
            }
            else
            {
                fout_Mphot << M_photon_origin_counter[idx_M] << std::endl;
            }
        }
        // std::cout << "row " << irow << std::endl;
    }

    fout_Mphot.close();


    // Calculating phi table
    std::vector<double> M_PhotPhi(n_tot*n_tot);
    for (int irow = 0; irow < n_tot; irow++)
    {
        for (int icol = 0; icol < n_tot; icol++)
        {
            int idx_M = irow*n_tot + icol;
            double photon_energy = 0.;
            if (irow < n_tiles)
            {
                int idx_sfc = irow;
                if (!std::isinf(field_atm_phi_per_photon[idx_sfc]))
                {
                    photon_energy = field_sfc_phi_per_photon[idx_sfc] * M_photon_origin_counter[idx_M];
                }
            }
            else
            if ((irow >= n_tiles) && (irow < (n_tiles + n_volumes)))
            {
                int idx_atm = irow - n_tiles;
                if (!std::isinf(field_atm_phi_per_photon[idx_atm]))
                {
                    photon_energy = field_atm_phi_per_photon[idx_atm] * M_photon_origin_counter[idx_M];
                }
            }
            M_PhotPhi[idx_M] = photon_energy;
        }
    }

    // Writing phi matrix
    std::ofstream fout_MphotPhi("/home/gijs-hogeboom/dev/lwproj/data_output/photon_matrix/photon_phi_matrix_" + CASE + 
                                 "_" + INTERCELL_TECHNIQUE + "_Natm" + std::to_string(static_cast<int>(std::log2(Natm))) + "_Nsfc" + 
                                 std::to_string(static_cast<int>(std::log2(Nsfc))) + ".dat");
    if (!fout_MphotPhi.is_open())
    {
        std::cerr << "Error, cannot open fout_MphotPhi!" << std::endl;
    }

    // Writing header
    fout_MphotPhi << "received from,";
    for (int icol = 0; icol < n_tot; icol++)
    {
        if (icol < (n_tot-1))
        {
            fout_MphotPhi << counter_axis[icol] << ',';
        }
        else
        {
            fout_MphotPhi << counter_axis[icol] << std::endl;
        }
    }

    // Writing body
    for (int irow = 0; irow < n_tot; irow++)
    {
        fout_MphotPhi << counter_axis[irow] << ',';

        for (int icol = 0; icol < n_tot; icol++)
        {
            int idx_M = irow*n_tot + icol;
            if (icol < (n_tot-1))
            {
                fout_MphotPhi << M_PhotPhi[idx_M] << ',';
            }
            else
            {
                fout_MphotPhi << M_PhotPhi[idx_M] << std::endl;
            }
        }
    }

    fout_MphotPhi.close();






    double tot_absorbed_z34 = kahan_sum(arr_z34_phi);
    // std::cout << std::setprecision(16) << field_atm_absorbed[34] << '\t' << tot_absorbed_z34 << '\t' << arr_z34_phi.size() << '\t' << field_atm_absorbed[34] - tot_absorbed_z34 << std::endl;
    // std::cout << std::setprecision(16) << dummy_photon_atm_phi << ',' << dummy_photon_sfc_phi << std::endl;
    if (verbose) 
    {
        std::cout << "[][][][][] Estimation of max numerical error: " << tot_error/(dx*dy) << std::endl;
    }               
    ///////////////// CALCULATING HEATING RATES /////////////////

    if (verbose)
    {
        std::cout << "  MC: Calculating heating rates" << std::endl;
    }

    std::vector<double> field_atm_numerical_differences(n_volumes);
    for (int idx_atm = 0; idx_atm < n_volumes; idx_atm++)
    {
        field_atm_numerical_differences[idx_atm] = field_atm_absorbed_db[idx_atm] - field_atm_absorbed[idx_atm];
    }

    // std::cout << "[][][][][] Total difference: " << kahan_sum(field_atm_numerical_differences)/(dx*dy) << std::endl;

    std::vector<double> field_atm_heating_rates(n_volumes);
    std::vector<double> field_sfc_heating_rates(n_tiles);
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < itot; i++)
    {
        for (int j = 0; j < jtot; j++)
        {
            for (int k = 0; k < ktot; k++)
            {
                int idx_atm =  i * ktot * jtot + j * ktot + k;
                double atm_net_phi = field_atm_absorbed[idx_atm] - field_atm_emitted[idx_atm];
                double dz = arr_dz[i];
                field_atm_heating_rates[idx_atm] = atm_net_phi / (cdouble::RHO * cdouble::CP * dx * dy * dz) * 86400;

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    double sfc_net_phi = field_sfc_absorbed[idx_sfc] - field_sfc_emitted[idx_sfc];
                    field_sfc_heating_rates[idx_sfc] = sfc_net_phi / (cdouble::RHO * cdouble::CP * dx * dy) * 86400;
                }
            }
        }
    }




    // MC energy balance
    // sfc_source          + atm_source          + TOA_source = sfc_sink                 + atm_sink                 + TOA_sink
    // field_sfc_phi.sum() + field_atm_phi.sum() + 0          = field_sfc_absorbed.sum() + field_atm_absorbed.sum() + TOA_absorbed

    double sfc_source = kahan_sum(field_sfc_emitted) / (dx*dy);
    double sfc_sink   = kahan_sum(field_sfc_absorbed) / (dx*dy);
    double atm_source = kahan_sum(field_atm_emitted) / (dx*dy);
    double atm_sink   = kahan_sum(field_atm_absorbed) / (dx*dy);
    double TOA_source = 0.;
    double TOA_sink   = TOA_absorbed / (dx*dy);

    double netto_phi = sfc_source + atm_source + TOA_source - sfc_sink - atm_sink - TOA_sink;
    double netto_phi_percentage = netto_phi/(sfc_source + atm_source + TOA_source) * 100.;


    if (print_EB)
    {
        std::cout << "+++ MC ENERGY BALANCE ++++++++++++++" << std::endl;
        std::cout << "+ -- source ------------------------" << std::endl;
        std::cout << "+ sfc source:      " << sfc_source << std::endl;
        std::cout << "+ atm source:      " << atm_source << std::endl;
        std::cout << "+ TOA source:      " << TOA_source << std::endl;
        std::cout << "+ -- sinks -------------------------" << std::endl;
        std::cout << "+ sfc sink:        " << sfc_sink << std::endl;
        std::cout << "+ atm sink:        " << atm_sink << std::endl;
        std::cout << "+ TOA sink:        " << TOA_sink << std::endl;
        std::cout << "+ -- net ---------------------------" << std::endl;
        std::cout << "+ sfc net:         " << sfc_sink - sfc_source << std::endl;
        std::cout << "+ atm net:         " << atm_sink - atm_source << std::endl;
        std::cout << "+ TOA net:         " << TOA_sink - TOA_source << std::endl;
        std::cout << "+ -- sums --------------------------" << std::endl;
        std::cout << "+ sources:         " << sfc_source + atm_source + TOA_source << std::endl;
        std::cout << "+ sinks:           " << sfc_sink + atm_sink + TOA_sink << std::endl;
        std::cout << "+ sources - sinks: " << netto_phi << std::endl;
        std::cout << "+ as percentage:   " << netto_phi_percentage << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    }


    

    return field_atm_heating_rates;
}