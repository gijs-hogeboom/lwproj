#include <vector>
#include <string>

#include "precision.h"


std::vector<Real> run_plane_parallel(const std::vector<Real>& arr_z,
                                     const std::vector<Real>& arr_zh,
                                     const std::vector<Real>& arr_dz,
                                     const std::vector<Real>& field_atm_kext,
                                     const std::vector<Real>& field_atm_B,
                                     const std::vector<Real>& field_sfc_B,
                                     const std::string& CASE,
                                     const Real dx,
                                     const Real dy,
                                     const int jtot,
                                     const int ktot,
                                     const bool print_EB,
                                     const bool verbose,
                                     const bool OUTPUT_3D);