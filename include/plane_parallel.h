#include <vector>
#include <string>


std::vector<float> run_plane_parallel(const std::vector<float>& arr_z,
                                      const std::vector<float>& arr_zh,
                                      const std::vector<float>& arr_dz,
                                      const std::vector<float>& field_atm_kext,
                                      const std::vector<float>& field_atm_B,
                                      const std::vector<float>& field_sfc_B,
                                      const std::string& CASE,
                                      const float dx,
                                      const float dy,
                                      const int jtot,
                                      const int ktot,
                                      const bool print_EB,
                                      const bool verbose,
                                      const bool OUTPUT_3D);