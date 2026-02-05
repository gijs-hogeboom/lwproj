#include <vector>
#include <string>

#include "precision.h"


std::vector<Real> run_MC(const std::vector<Real>& arr_z,
                           const std::vector<Real>& arr_zh,
                           const std::vector<Real>& arr_dz,
                           const std::vector<Real>& field_atm_kext,
                           const std::vector<Real>& field_atm_B,
                           const std::vector<Real>& field_sfc_B,
                           const std::vector<Real>& field_atm_SSA,
                           const std::vector<Real>& field_atm_ASY,
                           const Real dx,
                           const Real dy,
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