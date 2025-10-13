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
#include <cstdint>
#include <queue>
#include <omp.h>
#include <chrono>

#include "util.h"




struct Xoshiro256ss {
    uint64_t s[4];

    Xoshiro256ss(uint64_t seed = 1) {
        // SplitMix64 to initialize the state
        uint64_t z = seed + 0x9e3779b97f4a7c15ULL;
        for (int i = 0; i < 4; ++i) {
            z ^= (z >> 30); z *= 0xbf58476d1ce4e5b9ULL;
            z ^= (z >> 27); z *= 0x94d049bb133111ebULL;
            z ^= (z >> 31);
            s[i] = z;
            z += 0x9e3779b97f4a7c15ULL;
        }
    }

    inline uint64_t next() {
        const uint64_t result = rotl(s[1] * 5, 7) * 9;
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    inline double next_double() {
        // Take upper 53 bits of next() and convert to double in [0,1)
        return (next() >> 11) * (1.0 / 9007199254740992.0);
    }

private:
    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
};


struct FastRNG {
    Xoshiro256ss rng;
    FastRNG(uint64_t seed) : rng(seed) {}

    inline double uniform() { return rng.next_double(); }

    inline int uniform_int(int max_exclusive) {
        return static_cast<int>(((unsigned __int128)rng.next() * (unsigned __int128)max_exclusive) >> 64);
    }
};


struct AliasTable {
    std::vector<double> prob;
    std::vector<int> alias;
    int n;

    AliasTable(const std::vector<double>& weights) {
        n = weights.size();
        prob.resize(n);
        alias.resize(n);

        std::vector<double> scaled(weights);
        double sum = std::accumulate(scaled.begin(), scaled.end(), 0.0);
        for (auto& w : scaled) w *= n / sum;

        std::queue<int> small, large;
        for (int i = 0; i < n; ++i)
            (scaled[i] < 1.0 ? small : large).push(i);

        while (!small.empty() && !large.empty()) {
            int s = small.front(); small.pop();
            int l = large.front(); large.pop();
            prob[s] = scaled[s];
            alias[s] = l;
            scaled[l] = scaled[l] + scaled[s] - 1.0;
            (scaled[l] < 1.0 ? small : large).push(l);
        }

        while (!large.empty()) { prob[large.front()] = 1.0; large.pop(); }
        while (!small.empty()) { prob[small.front()] = 1.0; small.pop(); }
    }

    inline int sample(FastRNG& rng) const {
        int i = rng.rng.next() % n;
        double r = rng.rng.next_double();
        return (r < prob[i]) ? i : alias[i];
    }
};







void photon_propagation(const std::vector<int>& arr_photons_pos_idx,
                        const std::vector<double>& arr_photons_pos_x,
                        const std::vector<double>& arr_photons_pos_y,
                        const std::vector<double>& arr_photons_pos_z,
                        const std::vector<double>& arr_photons_mu,
                        const std::vector<double>& arr_photons_az,
                        const std::vector<double>& arr_photons_tau,
                        const std::vector<double>& arr_photons_phi,
                        const std::vector<double>& field_kext,
                        const std::vector<double>& arr_xh,
                        const std::vector<double>& arr_yh,
                        const std::vector<double>& arr_zh,
                        const std::vector<double>& arr_x,
                        const std::vector<double>& arr_y,
                        const std::vector<double>& arr_z,
                        std::vector<double>& field_atm_absorbed,
                        std::vector<double>& field_sfc_absorbed,
                        double& TOA_absorbed,
                        double x_max,
                        double y_max,
                        double z_max,
                        int itot,
                        int jtot,
                        int ktot,
                        double dx,
                        double dy,
                        int N,
                        int domain_section)
{
    const double eps = 1e-5;
    const int jktot = jtot * ktot;

    // #pragma omp for schedule(dynamic)
    for (int idx_photon = 0; idx_photon < N; idx_photon++)
    {
        
        // Loading initial photon variables
        int idx_flat = arr_photons_pos_idx[idx_photon];
        double x      = arr_photons_pos_x[idx_photon];
        double y      = arr_photons_pos_y[idx_photon];
        double z      = arr_photons_pos_z[idx_photon];
        double mu     = arr_photons_mu[idx_photon];
        double az     = arr_photons_az[idx_photon];
        double tau    = arr_photons_tau[idx_photon];

        // Calculating cartesian direction vector
        double s = sqrt(1 - mu*mu);
        double dx = s*cos(az);
        double dy = s*sin(az);
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

                // #pragma omp atomic
                TOA_absorbed += arr_photons_phi[idx_photon];
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
                break;
            }

        
            // retrieving kext from flattened grid
            int idx_flat = idx_z * jktot + idx_y * ktot + idx_x;
            double current_kext = field_kext[idx_flat];



            // Scanning collision with cell boundaries
            double time_x, time_y, time_z;
            double dn = 0.;
            double f = 1.;

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
            double dist_x = f * dx;
            double dist_y = f * dy;
            double dist_z = f * dz;
            
            double ds = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
            double max_s = tau/current_kext;
            double tau_absorbed = current_kext*ds;

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
            }
        }
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
                          int Natm,
                          int Nsfc,
                          std::vector<double>& EB_MC)
{

    ///////////////////// INITIALIZING DOMAIN /////////////////////// 
    std::cout << "  MC: Initializing domain" << std::endl;
 
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
    std::vector<double> field_sfc_eps(n_tiles, 1.0f);

    std::vector<double> field_atm_absorbed(n_volumes);
    std::vector<double> field_atm_emitted(n_volumes);
    std::vector<double> field_sfc_absorbed(n_tiles);
    std::vector<double> field_sfc_emitted(n_tiles);

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
    std::vector<double> arr_photons_atm_pos_x(Natm);
    std::vector<double> arr_photons_atm_pos_y(Natm);
    std::vector<double> arr_photons_atm_pos_z(Natm);
    std::vector<double> arr_photons_atm_mu(Natm);
    std::vector<double> arr_photons_atm_az(Natm);
    std::vector<double> arr_photons_atm_tau(Natm);
    std::vector<double> arr_photons_atm_phi(Natm);
    std::vector<int>   field_atm_photons_per_gridcell(n_volumes, 0);

    std::vector<int>   arr_photons_sfc_pos_idx(Nsfc);
    std::vector<double> arr_photons_sfc_pos_x(Nsfc);
    std::vector<double> arr_photons_sfc_pos_y(Nsfc);
    std::vector<double> arr_photons_sfc_pos_z(Nsfc, 0.0f);
    std::vector<double> arr_photons_sfc_mu(Nsfc);
    std::vector<double> arr_photons_sfc_az(Nsfc);
    std::vector<double> arr_photons_sfc_tau(Nsfc);
    std::vector<double> arr_photons_sfc_phi(Nsfc);
    std::vector<int>   field_sfc_photons_per_gridcell(n_tiles, 0);

    double TOA_absorbed = 0.0f; // This variable keeps track of photons lost to TOA
    

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
                arr_photons_atm_az[idx_photon]  = 2*cf::PI*rng.uniform();
                arr_photons_atm_tau[idx_photon] = -logf(rng.uniform());

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
                arr_photons_sfc_az[idx_photon]  = 2*cf::PI*rng.uniform();
                arr_photons_sfc_tau[idx_photon] = -logf(rng.uniform());

                // Countint photons per gridcell
                field_sfc_photons_per_gridcell[idx_sfc] += 1;
                

            }

            // Calculating the power per photon for each cell in the field
            std::vector<float> field_atm_phi_per_photon(n_volumes);
            std::vector<float> field_sfc_phi_per_photon(n_tiles);
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

            AliasTable alias_atm(field_atm_phi);

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
                    arr_photons_atm_az[idx_photon]  = 2*cf::PI*rng_local.uniform();
                    arr_photons_atm_tau[idx_photon] = -logf(rng_local.uniform());

                    // Countint photons per gridcell
                    #pragma omp atomic
                    field_atm_photons_per_gridcell[idx_atm] += 1;
                }
            }
            
            // Surface
            AliasTable alias_sfc(field_sfc_phi);
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
                    arr_photons_sfc_az[idx_photon]  = 2*cf::PI*rng_local.uniform();
                    arr_photons_sfc_tau[idx_photon] = -logf(rng_local.uniform());

                    // Countint photons per gridcell
                    #pragma omp atomic
                    field_sfc_photons_per_gridcell[idx_sfc] += 1;
                }
            }
            


            // Calculating the power per photon for each cell in the field
            std::vector<float> field_atm_phi_per_photon(n_volumes);
            std::vector<float> field_sfc_phi_per_photon(n_tiles);
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

                        size_t idx_atm = i*jtot*ktot + j*ktot + k;
                        double center = field_atm_phi[idx_atm];
                        double diffs = 0.;

                        for (size_t iidx = 0; iidx < 27; iidx++)
                        {
                            int ii = i + iidx / 9 - 1;
                            int jj = j + (iidx / 3) % 3 - 1;
                            int kk = k + iidx % 3 - 1;

                            if      (jj == -1)   {jj = jtot - 1;}
                            else if (jj == jtot) {jj = 0;}
                            if      (kk == -1)   {kk = ktot - 1;}
                            else if (kk == ktot) {kk = 0;}

                            if (ii == -1)
                            {
                                size_t idx_sfc = j*ktot + k;
                                diffs += pow(center - field_sfc_phi[idx_sfc], 2);
                            }
                            else if (ii == itot)
                            {
                                diffs += pow(center,2);
                            }
                            else
                            {
                                size_t idx_atm_diffs = ii*jtot*ktot + jj*ktot + kk;
                                diffs += pow(center - field_atm_phi[idx_atm_diffs], 2);
                            }
                        }

                        field_atm_phi_gradient[idx_atm] = std::sqrt(diffs / 26.); // Taking the root mean of the differences (26 neighbors)

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
            AliasTable alias_atm(field_atm_phi_gradient);

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
                    arr_photons_atm_az[idx_photon]  = 2*cf::PI*rng_local.uniform();
                    arr_photons_atm_tau[idx_photon] = -logf(rng_local.uniform());

                    // Countint photons per gridcell
                    #pragma omp atomic
                    field_atm_photons_per_gridcell[idx_atm] += 1;
                }
            }
            
            // Surface
            AliasTable alias_sfc(field_sfc_phi_gradient);
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
                    arr_photons_sfc_az[idx_photon]  = 2*cf::PI*rng_local.uniform();
                    arr_photons_sfc_tau[idx_photon] = -logf(rng_local.uniform());

                    // Countint photons per gridcell
                    #pragma omp atomic
                    field_sfc_photons_per_gridcell[idx_sfc] += 1;
                }
            }
            


            // Calculating the power per photon for each cell in the field
            std::vector<float> field_atm_phi_per_photon(n_volumes);
            std::vector<float> field_sfc_phi_per_photon(n_tiles);
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
                       TOA_absorbed,
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
                       TOA_absorbed,
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
                field_atm_heating_rates[idx_atm] = atm_net_phi / (cf::RHO * cf::CP * dx * dy * dz) * 86400;

                if (i == 0)
                {
                    int idx_sfc = j * ktot + k;
                    double sfc_net_phi = field_sfc_absorbed[idx_sfc] - field_sfc_emitted[idx_sfc];
                    field_sfc_heating_rates[idx_sfc] = sfc_net_phi / (cf::RHO * cf::CP * dx * dy) * 86400;
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

    // Storing output
    EB_MC[0] = sfc_source;
    EB_MC[1] = sfc_sink;
    EB_MC[2] = atm_source;
    EB_MC[3] = atm_sink;
    EB_MC[4] = TOA_source;
    EB_MC[5] = TOA_sink;
    

    return field_atm_heating_rates;
}