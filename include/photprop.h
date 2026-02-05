#include <vector>
#include <string>

#include "util.h"
#include "precision.h"


void photon_propagation(const AliasTable& aliastable,
                        const std::vector<Real>& field_kext,
                        const std::vector<Real>& field_sfc_eps,
                        const std::vector<Real>& field_SSA,
                        const std::vector<Real>& field_ASY,
                        const std::vector<Real>& arr_xh,
                        const std::vector<Real>& arr_yh,
                        const std::vector<Real>& arr_zh,
                        const std::vector<Real>& arr_x,
                        const std::vector<Real>& arr_y,
                        const std::vector<Real>& arr_z,
                        const std::vector<Real>& arr_dz,
                        const std::vector<Real>& field_phi,
                        std::vector<Real>& field_atm_net_phi,
                        std::vector<Real>& field_sfc_net_phi,
                        std::vector<Real>& field_TOA_net_phi,
                        const long int N,
                        const int domain_section,
                        const std::string& INTERCELL_TECHNIQUE,
                        const bool Pesc_mode,
                        const bool enable_scattering);