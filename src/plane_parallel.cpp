#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <boost/math/quadrature/gauss.hpp>

#include "util.h"


using namespace boost::math::quadrature;


class LinearInterpTau {
    std::vector<float> arr_x, arr_y;
public:
    LinearInterpTau(const std::vector<float>& arr_x_in, const std::vector<float>& arr_y_in)
        : arr_x(arr_x_in), arr_y(arr_y_in) {};

    float operator()(float x) const {

        if (x >= arr_x.front()) return arr_y.front();
        if (x <= arr_x.back())  return arr_y.back();


        // Find interval safely
        size_t i = 0;
        while (i + 1 < arr_x.size() && arr_x[i+1] >= x) ++i;

        float x0 = arr_x[i];
        float x1 = arr_x[i+1];
        float y0 = arr_y[i];
        float y1 = arr_y[i+1];

        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }
};


std::vector<float> linspace(float start, float end, int N)
{
    std::vector<float> arr_result(N);

    float step_size = (end - start)/(N - 1);

    for (size_t i = 0; i < N; i++)
    {
        arr_result[i] = start + (static_cast<float>(i) * step_size);
    }

    return arr_result;
}




std::vector<float> cumulative_trapezoid(const std::vector<float>& arr_x,
                                        const std::vector<float>& arr_y) 
{
    size_t n = arr_x.size();

    std::vector<float> arr_cY(n - 1, 0.0);

    for (size_t i = 0; i < (n - 1); i++) 
    {
        float dx = arr_x[i+1] - arr_x[i];
        float Y  = arr_y[i] * dx + (arr_y[i+1] - arr_y[i]) * dx / 2.0;
        if (i == 0)
        {
            arr_cY[i] = Y;
        } 
        else
        {
            arr_cY[i] = arr_cY[i-1] + Y;
        }
    }

    return arr_cY;
}

float trapezoid(const std::vector<float>& arr_x,
                const std::vector<float>& arr_y)
{
    size_t n = arr_x.size();
    std::vector<float> values((n - 1));

    for (size_t i = 0; i < (n - 1); i++)
    {
        float dx = arr_x[i+1] - arr_x[i];
        float Y  = arr_y[i] * dx + (arr_y[i+1] - arr_y[i]) * dx / 2.0;
        values[i] = Y;
    }

    float result = std::accumulate(values.begin(), values.end(), 0.0);
    return result;
}

float trapezoid(const std::vector<float>& arr_y, float dx)
{
    size_t n = arr_y.size();
    std::vector<float> values((n - 1));

    for (size_t i = 0; i < (n - 1); i++)
    {
        float Y  = arr_y[i] * dx + (arr_y[i+1] - arr_y[i]) * dx / 2.0;
        values[i] = Y;
    }

    float result = std::accumulate(values.begin(), values.end(), 0.0);
    return result;
}


float I_upwards_emission(float tau_prime, float B_at_tau_prime, float tau, float mu)
{
    return B_at_tau_prime * std::exp( -(tau_prime - tau)/mu );
}

float I_upwards_exctinction(float tau, float mu, float I_at_sfc, float tau_at_sfc)
{
    return I_at_sfc * std::exp( -(tau_at_sfc - tau)/mu );
}


float I_downwards_emission(float tau_prime, float B_at_tau_prime, float tau, float mu)
{
    return B_at_tau_prime * std::exp( -(tau - tau_prime)/mu );
}

float I_downwards_exctinction(float tau, float mu, float I_at_TOA, float tau_at_TOA)
{
    return I_at_TOA * std::exp( -(tau - tau_at_TOA)/mu );
}


float I_complete(float term_extinction, const std::vector<float>& vec_emission_values, float dtau, float mu)
{
    float term_emission = trapezoid(vec_emission_values, dtau) / mu;
    return term_extinction + term_emission;
}


std::vector<float> run_plane_parallel(const std::vector<float>& arr_z,
                                      const std::vector<float>& arr_zh,
                                      const std::vector<float>& arr_dz,
                                      const std::vector<float>& arr_kext,
                                      const std::vector<float>& arr_Batm,
                                      const std::vector<float>& arr_Batmh,
                                      const std::string& CASE,
                                      float Bsfc,
                                      int N_mu)
{

    // Initializing domain
    std::cout << "  PP: Initializing domain" << std::endl;

    int itot  = arr_z.size();
    int itoth = arr_zh.size();
    int itotg = arr_z.size() + 2; // for ghost-cell including arrays

    // Initializing arr_zg (arr_z with 2 ghost cells)
    std::vector<float> arr_zg(itotg);


    arr_zg[0] = arr_z[0] - arr_dz[0];
    arr_zg[itotg-1] = arr_z[itot-1] + arr_dz[itot-1];
    for (int i = 1; i < (itotg-1); i++) { arr_zg[i] = arr_z[i-1]; }

    // Initializing arr_kextg (arr_kext with 2 ghost cells)
    std::vector<float> arr_kextg(itotg);
    arr_kextg[0] = arr_kext[0] - (arr_kext[1] - arr_kext[0])/(arr_z[1] - arr_z[0]);
    arr_kextg[itotg-1] = arr_kext[itot-1] + (arr_kext[itot-1] - arr_kext[itot-2])/(arr_z[itot-1] - arr_z[itot-2]);
    for (int i = 1; i < (itotg-1); i++) { arr_kextg[i] = arr_kext[i-1]; }

    if (arr_kextg[0] <= 0) { arr_kextg[0] = 0.0; }
    if (arr_kextg[itotg-1] <= 0) { arr_kextg[itotg-1] = 0.0; }

    // Initializing arr_Batmg (arr_Batm with 2 ghost cells)
    std::vector<float> arr_Batmg(itotg);
    arr_Batmg[0] = arr_Batm[0] - (arr_Batm[1] - arr_Batm[0])/(arr_z[1] - arr_z[0]);
    arr_Batmg[itotg-1] = arr_Batm[itot-1] + (arr_Batm[itot-1] - arr_Batm[itot-2])/(arr_z[itot-1] - arr_z[itot-2]);
    for (int i = 1; i < (itotg-1); i++) { arr_Batmg[i] = arr_Batm[i-1]; }

    if (arr_Batmg[0] <= 0) { arr_Batmg[0] = 0.0; }
    if (arr_Batmg[itotg-1] <= 0) { arr_Batmg[itotg-1] = 0.0; }

    // Generating tauh (cumulative trapezoid from TOA to bottom)
    std::vector<float> arr_kextg_reversed(arr_kextg.rbegin(), arr_kextg.rend());
    std::vector<float> arr_zg_reversed(arr_zg.rbegin(), arr_zg.rend());

    std::vector<float> arr_tauh = cumulative_trapezoid(arr_zg_reversed,arr_kextg_reversed);
    std::reverse(arr_tauh.begin(), arr_tauh.end());
    for (auto& v: arr_tauh) { v = -v; }

    // Creating B(tau) function
    LinearInterpTau B_tauh(arr_tauh, arr_Batmh);
    auto f_B_tauh = [&](float tau) { return B_tauh(tau); };

    float tau_at_sfc = arr_tauh[0];
    float tau_at_TOA = arr_tauh[arr_tauh.size()-1];

    float I_at_sfc = Bsfc;
    float I_at_TOA = 0.;
    

    // Angles
    std::vector<float> arr_mu = linspace(0.01, 1.0, N_mu);
    int jtot = N_mu;
    float ktot = 200.;


    // Calculating emission terms
    std::cout << "  PP: Calculating emission terms" << std::endl;

    std::vector<float> M_I_emission_uph(itoth*jtot*ktot);
    std::vector<float> M_I_emission_downh(itoth*jtot*ktot);

    std::vector<float> arr_dtau_uph(itoth);
    std::vector<float> arr_dtau_downh(itoth);
    

    for (size_t i = 0; i < itoth; i++)
    {
        // Loading tau
        float tau = arr_tauh[i];

        // Storing dtau at this level
        arr_dtau_uph[i]   = (tau_at_sfc - tau)/(ktot-1);
        arr_dtau_downh[i] = (tau - tau_at_TOA)/(ktot-1);

        for (size_t j = 0; j < jtot; j++)
        {
            // Loading mu
            float mu = arr_mu[j];
            for (float k = 0.; k < ktot; k++)
            {
                int idx = i*jtot*ktot + j*ktot + k;

                float tau_prime_up        = (k/(ktot-1)*tau        + (1 - k/(ktot-1))*tau_at_sfc);
                float B_at_tau_prime_up   = f_B_tauh(tau_prime_up);
                M_I_emission_uph[idx]     = I_upwards_emission(tau_prime_up, B_at_tau_prime_up, tau, mu);
                
                float tau_prime_down      = (k/(ktot-1)*tau_at_TOA + (1 - k/(ktot-1))*tau);
                float B_at_tau_prime_down = f_B_tauh(tau_prime_down);
                M_I_emission_downh[idx]   = I_downwards_emission(tau_prime_down, B_at_tau_prime_down, tau, mu);
            }
        }
    }

    // Calculating I(tau, mu)
    std::cout << "  PP: Calculating I(tau, mu)" << std::endl;

    std::vector<float> M_I_uph(itoth*jtot);
    std::vector<float> M_I_downh(itoth*jtot);

    for (size_t i = 0; i < itoth; i++)
    {
        float tau = arr_tauh[i];
        float dtau_up = arr_dtau_uph[i];
        float dtau_down = arr_dtau_downh[i];
        for (size_t j = 0; j < jtot; j++)
        {
            float mu = arr_mu[j];
            int idx = i*jtot + j;
            std::vector<float> I_uph_kaxis(ktot);
            std::vector<float> I_downh_kaxis(ktot);
            for (size_t k = 0; k < ktot; k++)
            {
                int idx_k = i*jtot*ktot + j*ktot + k;
                I_uph_kaxis[k]   = M_I_emission_uph[idx_k];
                I_downh_kaxis[k] = M_I_emission_downh[idx_k];
            }

            M_I_uph[idx] = I_complete(I_upwards_exctinction(tau, mu, I_at_sfc, tau_at_sfc),
                                      I_uph_kaxis, 
                                      dtau_up,
                                      mu);
            M_I_downh[idx] = I_complete(I_downwards_exctinction(tau, mu, I_at_TOA, tau_at_TOA),
                                        I_downh_kaxis, 
                                        dtau_down,
                                        mu);

        }
    }


    // Calculating F(z)
    std::cout << "  PP: Calculating F(z)" << std::endl;

    std::vector<float> arr_F_uph(itoth);
    std::vector<float> arr_F_downh(itoth);

    for (size_t i = 0; i < itoth; i++)
    {
        // Calculating and storing I(tau, mu) * mu for each angle per level
        std::vector<float> arr_Imu_uph(jtot);
        std::vector<float> arr_Imu_downh(jtot);
        for (size_t j = 0; j < jtot; j++)
        {
            float mu = arr_mu[j];

            int idx = i*jtot + j;

            arr_Imu_uph[j]   = M_I_uph[idx]*mu;
            arr_Imu_downh[j] = M_I_downh[idx]*mu;
        }
        arr_F_uph[i]   = 2*cf::PI*trapezoid(arr_mu, arr_Imu_uph);
        arr_F_downh[i] = 2*cf::PI*trapezoid(arr_mu, arr_Imu_downh);
    }

    // Calculating net flux at each cell

    std::vector<float> arr_F_net(itot);
    for (size_t i = 0; i < itoth; i++)
    {
        arr_F_net[i] = arr_F_uph[i] + arr_F_downh[i+1] - arr_F_uph[i+1] - arr_F_downh[i];
    }


    // Calculating heating rates at each cell
    std::cout << "  PP: Calculating heating rates" << std::endl;

    std::vector<float> arr_heating_rates(itot);
    for (size_t i = 0; i < itot; i++)
    {
        arr_heating_rates[i] = 1 / (cf::RHO * cf::CP * arr_dz[i]) * arr_F_net[i] * 86400;
    }

    // PP energy balance
    // sfc_source          + TOA_source                       = sfc_sink                 + TOA_sink                 + atm_netto
    // pi*I_at_sfc         + 0                                = F_downh[0]               + F_uph[-1]                + F_net.sum()

    float sfc_source = cf::PI * I_at_sfc;
    float sfc_sink   = arr_F_downh[0];
    float atm_netto  = kahan_sum(arr_F_net);
    float TOA_source = 0.;
    float TOA_sink   = arr_F_uph[itoth-1];


    float netto_phi = sfc_source + TOA_source - sfc_sink - TOA_sink - atm_netto;
    float netto_phi_percentage = netto_phi/(sfc_source + TOA_source) * 100.;


    std::cout << "+++ PP ENERGY BALANCE ++++++++++++++" << std::endl;
    std::cout << "+ sfc source:      " << sfc_source << std::endl;
    std::cout << "+ sfc sink:        " << sfc_sink << std::endl;
    std::cout << "+ atm netto:       " << atm_netto << std::endl;
    std::cout << "+ TOA source:      " << TOA_source << std::endl;
    std::cout << "+ TOA sink:        " << TOA_sink << std::endl;
    std::cout << "+                  " << std::endl;
    std::cout << "+ sources - sinks: " << netto_phi << std::endl;
    std::cout << "+ as percentage:   " << netto_phi_percentage << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;

    
    return arr_heating_rates;
}