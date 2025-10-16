#include <vector>


std::vector<float> run_plane_parallel(const std::vector<float>& arr_z,
                                      const std::vector<float>& arr_zh,
                                      const std::vector<float>& arr_dz,
                                      const std::vector<float>& arr_kext,
                                      const std::vector<float>& arr_Batm,
                                      const std::vector<float>& arr_Batmh,
                                      float Bsfc,
                                      int N_mu,
                                      bool print_EB,
                                      bool verbose);