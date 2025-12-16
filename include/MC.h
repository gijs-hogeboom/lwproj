#include <vector>
#include <string>


std::vector<double> run_MC(const std::vector<double>& arr_z,
                           const std::vector<double>& arr_zh,
                           const std::vector<double>& arr_dz,
                           const std::vector<double>& field_atm_kext,
                           const std::vector<double>& field_atm_B,
                           const std::vector<double>& field_sfc_B,
                           const std::vector<double>& field_atm_SSA,
                           const std::vector<double>& field_atm_ASY,
                           const double dx,
                           const double dy,
                           const int ktot,
                           const int jtot,
                           const int itot,
                           const std::string& INTERCELL_TECHNIQUE,
                           const std::string& CASE,
                           const int Nphot,
                           const bool print_EB,
                           const bool verbose,
                           const bool enable_full_counter_matrix,
                           const bool Pesc_mode,
                           const bool OUTPUT_3D,
                           const bool enable_scattering);