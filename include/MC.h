#include <vector>
#include <string>


std::vector<float> run_MC(const std::vector<float>& arr_z,
                           const std::vector<float>& arr_zh,
                           const std::vector<float>& arr_dz,
                           const std::vector<float>& field_atm_kext,
                           const std::vector<float>& field_atm_B,
                           const std::vector<float>& field_sfc_B,
                           const std::vector<float>& field_atm_SSA,
                           const std::vector<float>& field_atm_ASY,
                           const float dx,
                           const float dy,
                           const int ktot,
                           const int jtot,
                           const int itot,
                           const std::string& INTERCELL_TECHNIQUE,
                           const std::string& CASE,
                           const long int Nphot,
                           const bool print_EB,
                           const bool verbose,
                           const bool enable_full_counter_matrix,
                           const bool Pesc_mode,
                           const bool OUTPUT_3D,
                           const bool enable_scattering);