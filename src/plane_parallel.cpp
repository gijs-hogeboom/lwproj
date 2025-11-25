#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <boost/math/quadrature/gauss.hpp>

#include "util.h"


using namespace boost::math::quadrature;




std::vector<double> run_plane_parallel(const std::vector<double>& arr_z,
                                      const std::vector<double>& arr_zh,
                                      const std::vector<double>& arr_dz,
                                      const std::vector<double>& arr_kext,
                                      const std::vector<double>& arr_Batm,
                                      double Bsfc,
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


    // Generating tauh (cumulative sum of kext*dz from TOA to bottom)
    std::vector<double> arr_tauh(itoth);
    arr_tauh[itoth-1] = 0.;
    for (int i = (itoth-2); i >= 0; i--)
    {
        arr_tauh[i] = arr_tauh[i+1] + (arr_dz[i]*arr_kext[i]);
    }


    double tau_at_sfc = arr_tauh[0];
    double tau_at_TOA = arr_tauh[itoth-1];

    double I_at_sfc = Bsfc;
    double I_at_TOA = 0.;
    

    // Angles
    std::vector<double> arr_mu = linspace(0.01, 1.0, N_mu);
    int jtot = N_mu;


    // Calculating I(tau, mu) (p.32 of K. N. Liou, 2002)
    if (verbose)
    {
        std::cout << "  PP: Calculating I(tau, mu)" << std::endl;
    }
    std::vector<double> M_I_uph(itoth*jtot);
    std::vector<double> M_I_downh(itoth*jtot);

    for (int ih = 0; ih < itoth; ih++)
    {
        // Loading tau
        double tau = arr_tauh[ih];

        for (int j = 0; j < jtot; j++)
        {
            // Loading mu
            double mu = arr_mu[j];

            // Calculating extinction term
            double extinction_term_uph   = I_at_sfc * std::exp( -(tau_at_sfc - tau)/mu );
            double extinction_term_downh = I_at_TOA * std::exp( -(tau - tau_at_TOA)/mu );

            // Calculating emission term
            double emission_term_uph = 0.;
            double emission_term_downh = 0.;
            for (int k = 0; k < ih; k++)
            {
                double temp = (arr_Batm[k] * (std::exp(-(arr_tauh[k+1] - tau)/mu) - std::exp(-(arr_tauh[k] - tau)/mu))); 
                emission_term_uph += temp;
            }
            for (int k = (itot-1); k > (ih - 1); k--)
            {
                double temp = (arr_Batm[k] * (std::exp(-(tau - arr_tauh[k])/mu) - std::exp(-(tau - arr_tauh[k+1])/mu)));
                emission_term_downh += temp;
            }    

            // Combining extinction + emission term
            int idx        = ih*jtot + j;
            M_I_uph[idx]   = extinction_term_uph + emission_term_uph;
            M_I_downh[idx] = extinction_term_downh + emission_term_downh;
        }
    }


    // Calculating F(z)
    if (verbose)
    {
        std::cout << "  PP: Calculating F(z)" << std::endl;
    }

    std::vector<double> arr_F_uph(itoth);
    std::vector<double> arr_F_downh(itoth);

    for (size_t ih = 0; ih < itoth; ih++)
    {
        // Calculating and storing I(tau, mu) * mu for each angle per level
        std::vector<double> arr_Imu_uph(jtot);
        std::vector<double> arr_Imu_downh(jtot);
        for (size_t j = 0; j < jtot; j++)
        {
            double mu = arr_mu[j];

            int idx = ih*jtot + j;

            arr_Imu_uph[j]   = M_I_uph[idx]*mu;
            arr_Imu_downh[j] = M_I_downh[idx]*mu;
        }
        arr_F_uph[ih]   = 2.*cdouble::PI*trapezoid(arr_mu, arr_Imu_uph);
        arr_F_downh[ih] = 2.*cdouble::PI*trapezoid(arr_mu, arr_Imu_downh);
    }

    // Calculating net flux at each cell
    std::vector<double> arr_F_net(itot);
    for (size_t i = 0; i < itot; i++)
    {
        arr_F_net[i] = arr_F_uph[i] + arr_F_downh[i+1] - arr_F_uph[i+1] - arr_F_downh[i];
    }


    // Calculating heating rates at each cell
    if (verbose)
    {
        std::cout << "  PP: Calculating heating rates" << std::endl;
    }
    
    std::vector<double> arr_heating_rates(itot);
    for (size_t i = 0; i < itot; i++)
    {
        arr_heating_rates[i] = 1 / (cdouble::RHO * cdouble::CP * arr_dz[i]) * arr_F_net[i] * 86400.;
    }


    // PP energy balance
    // sfc_source          + TOA_source                       = sfc_sink                 + TOA_sink                 + atm_netto
    // pi*I_at_sfc         + 0                                = F_downh[0]               + F_uph[-1]                + F_net.sum()

    double sfc_source = cdouble::PI * I_at_sfc;
    double sfc_sink   = arr_F_downh[0];
    double atm_netto  = kahan_sum(arr_F_net);
    double TOA_source = 0.;
    double TOA_sink   = arr_F_uph[itoth-1];


    double netto_phi = sfc_source + TOA_source - sfc_sink - TOA_sink - atm_netto;
    double netto_phi_percentage = netto_phi/(sfc_source + TOA_source) * 100.;


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