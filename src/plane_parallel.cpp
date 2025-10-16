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
                                      float Bsfc,
                                      int N_mu,
                                      bool print_EB,
                                      bool verbose)
{

    // Initializing domain
    if (verbose)
    {
        std::cout << "  PP: Initializing domain" << std::endl;
    }

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

    // Generating tauh (cumulative sum of kext*dz from TOA to bottom)
    std::vector<float> arr_tauh(itoth);
    arr_tauh[itoth-1] = 0.;
    for (int i = (itoth-2); i >= 0; i--)
    {
        arr_tauh[i] = arr_tauh[i+1] + (arr_dz[i]*arr_kext[i]);
    }


    // Creating B(tau) function
    float tau_at_sfc = arr_tauh[0];
    float tau_at_TOA = arr_tauh[itoth-1];

    float I_at_sfc = Bsfc;
    float I_at_TOA = 0.;
    

    // Angles
    std::vector<float> arr_mu = linspace(0.01, 1.0, N_mu);
    int jtot = N_mu;
    float ktot = 200.;


    // Calculating I(tau, mu)
    if (verbose)
    {
        std::cout << "  PP: Calculating I(tau, mu)" << std::endl;
    }
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
            float extinction_term_uph   = I_upwards_exctinction(tau, mu, I_at_sfc, tau_at_sfc);
            float extinction_term_downh = I_downwards_exctinction(tau, mu, I_at_TOA, tau_at_TOA);

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

            int idx = ih*jtot + j;
            M_I_uph[idx] = extinction_term_uph + emission_term_uph;
            M_I_downh[idx] = extinction_term_downh + emission_term_downh;
        }
    }


    // Calculating F(z)
    if (verbose)
    {
        std::cout << "  PP: Calculating F(z)" << std::endl;
    }

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
        arr_F_uph[i]   = 2*cfloat::PI*trapezoid(arr_mu, arr_Imu_uph);
        arr_F_downh[i] = 2*cfloat::PI*trapezoid(arr_mu, arr_Imu_downh);
    }

    // Calculating net flux at each cell
    std::vector<float> arr_F_net(itot);
    for (size_t i = 0; i < itot; i++)
    {
        arr_F_net[i] = arr_F_uph[i] + arr_F_downh[i+1] - arr_F_uph[i+1] - arr_F_downh[i];
    }

    // LOGvec(arr_F_uph, "F up", true);
    // LOGvec(arr_F_downh, "F down", true);

    // Calculating heating rates at each cell
    if (verbose)
    {
        std::cout << "  PP: Calculating heating rates" << std::endl;
    }
    
    std::vector<float> arr_heating_rates(itot);
    for (size_t i = 0; i < itot; i++)
    {
        arr_heating_rates[i] = 1 / (cfloat::RHO * cfloat::CP * arr_dz[i]) * arr_F_net[i] * 86400;
    }
    // PP energy balance
    // sfc_source          + TOA_source                       = sfc_sink                 + TOA_sink                 + atm_netto
    // pi*I_at_sfc         + 0                                = F_downh[0]               + F_uph[-1]                + F_net.sum()

    float sfc_source = cfloat::PI * I_at_sfc;
    float sfc_sink   = arr_F_downh[0];
    float atm_netto  = kahan_sum(arr_F_net);
    float TOA_source = 0.;
    float TOA_sink   = arr_F_uph[itoth-1];


    float netto_phi = sfc_source + TOA_source - sfc_sink - TOA_sink - atm_netto;
    float netto_phi_percentage = netto_phi/(sfc_source + TOA_source) * 100.;


    if (print_EB)
    {
        std::cout << "+++ PP ENERGY BALANCE ++++++++++++++" << std::endl;
        std::cout << "+ -- source ------------------------" << std::endl;
        std::cout << "+ sfc source:      " << sfc_source << std::endl;
        std::cout << "+ atm source:      -" << std::endl;
        std::cout << "+ TOA source:      " << TOA_source << std::endl;
        std::cout << "+ -- sinks -------------------------" << std::endl;
        std::cout << "+ sfc sink:        " << sfc_sink << std::endl;
        std::cout << "+ atm sink:        -" << std::endl;
        std::cout << "+ TOA sink:        " << TOA_sink << std::endl;
        std::cout << "+ -- net ---------------------------" << std::endl;
        std::cout << "+ sfc net:         " << sfc_sink - sfc_source << std::endl;
        std::cout << "+ atm net:         " << atm_netto << std::endl;
        std::cout << "+ TOA net:         " << TOA_sink - TOA_source << std::endl;
        std::cout << "+ -- sums --------------------------" << std::endl;
        std::cout << "+ sources:          unknown" << std::endl;
        std::cout << "+ sinks:            unknown" << std::endl;
        std::cout << "+ sources - sinks: " << netto_phi << std::endl;
        std::cout << "+ as percentage:   " << netto_phi_percentage << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
    }


    return arr_heating_rates;
}