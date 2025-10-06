#pragma once

#include <vector>
#include <string>



std::vector<float> linspace(float start, float end, int N);

float kahan_sum(const std::vector<float>& values);

std::vector<float> cumulative_trapezoid(const std::vector<float>& arr_x,
                                        const std::vector<float>& arr_y);

float trapezoid(const std::vector<float>& arr_y, float dx);

float trapezoid(const std::vector<float>& arr_x,
                const std::vector<float>& arr_y);

std::vector<float> run_plane_parallel(const std::vector<float>& arr_z,
                                      const std::vector<float>& arr_zh,
                                      const std::vector<float>& arr_dz,
                                      const std::vector<float>& arr_kext,
                                      const std::vector<float>& arr_Batm,
                                      const std::vector<float>& arr_Batmh,
                                      const std::string& CASE,
                                      float Bsfc,
                                      int N_mu);


float I_upwards_emission(float tau_prime, float B_at_tau_prime, float tau, float mu);
float I_upwards_exctinction(float tau, float mu, float I_at_sfc, float tau_at_sfc);
float I_downwards_emission(float tau_prime, float B_at_tau_prime, float tau, float mu);
float I_downwards_exctinction(float tau, float mu, float I_at_TOA, float tau_at_TOA);
float I_complete(float term_extinction, const std::vector<float>& vec_emission_values, float dtau, float mu);