#pragma once

#include <vector>
#include <string>

std::vector<float> run_MC(const std::vector<float>& arr_z,
                          const std::vector<float>& arr_zh,
                          const std::vector<float>& arr_dz,
                          const std::vector<float>& arr_kext,
                          const std::vector<float>& arr_Batm,
                          float Bsfc,
                          float dx,
                          float dy,
                          int ktot,
                          int jtot,
                          int itot,
                          const std::string& CASE,
                          const std::string& INTERCELL_TECHNIQUE,
                          const std::string& INTRACELL_TECHNIQUE,
                          int Natm,
                          int Nsfc);