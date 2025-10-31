#include <vector>
#include <string>


std::vector<double> run_MC(const std::vector<double>& arr_z,
                          const std::vector<double>& arr_zh,
                          const std::vector<double>& arr_dz,
                          const std::vector<double>& arr_kext,
                          const std::vector<double>& arr_Batm,
                          double Bsfc,
                          double dx,
                          double dy,
                          int ktot,
                          int jtot,
                          int itot,
                          const std::string& INTERCELL_TECHNIQUE,
                          const std::string& INTRACELL_TECHNIQUE,
                          const std::string& CASE,
                          int Natm,
                          int Nsfc,
                          bool print_EB,
                          bool verbose,
                          bool enable_full_counter_matrix);