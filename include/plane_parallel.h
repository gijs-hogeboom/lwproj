#pragma once

#include <vector>
#include <string>



std::vector<double> linspace(double start, double end, int N);

std::vector<double> cumulative_trapezoid(const std::vector<double>& arr_x,
                                        const std::vector<double>& arr_y);

double trapezoid(const std::vector<double>& arr_y, double dx);

double trapezoid(const std::vector<double>& arr_x,
                const std::vector<double>& arr_y);

std::vector<double> run_plane_parallel(const std::vector<double>& arr_z,
                                      const std::vector<double>& arr_zh,
                                      const std::vector<double>& arr_dz,
                                      const std::vector<double>& arr_kext,
                                      const std::vector<double>& arr_Batm,
                                      const std::vector<double>& arr_Batmh,
                                      double Bsfc,
                                      int N_mu,
                                      std::vector<double>& EB_PP);


double I_upwards_emission(double tau_prime, double B_at_tau_prime, double tau, double mu);
double I_upwards_exctinction(double tau, double mu, double I_at_sfc, double tau_at_sfc);
double I_downwards_emission(double tau_prime, double B_at_tau_prime, double tau, double mu);
double I_downwards_exctinction(double tau, double mu, double I_at_TOA, double tau_at_TOA);
double I_complete(double term_extinction, const std::vector<double>& vec_emission_values, double dtau, double mu);