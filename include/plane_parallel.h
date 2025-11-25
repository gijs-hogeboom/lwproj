#include <vector>


std::vector<double> run_plane_parallel(const std::vector<double>& arr_z,
                                      const std::vector<double>& arr_zh,
                                      const std::vector<double>& arr_dz,
                                      const std::vector<double>& arr_kext,
                                      const std::vector<double>& arr_Batm,
                                      double Bsfc,
                                      int N_mu,
                                      bool print_EB,
                                      bool verbose);