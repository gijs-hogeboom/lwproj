#pragma once

#include <vector>
#include <string>

std::vector<float> run_plane_parallel(const std::vector<float>& arr_z,
                                      const std::vector<float>& arr_zh,
                                      const std::vector<float>& arr_dz,
                                      const std::vector<float>& arr_kext,
                                      const std::vector<float>& arr_Batm,
                                      const std::vector<float>& arr_Batmh,
                                      const std::string& CASE,
                                      int N_mu);
