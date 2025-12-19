#include <vector>
#include <string>

#include "util.h"


void photon_propagation(const AliasTable_float& aliastable,
                        const std::vector<float>& field_kext,
                        const std::vector<float>& field_sfc_eps,
                        const std::vector<float>& field_SSA,
                        const std::vector<float>& field_ASY,
                        const std::vector<float>& arr_xh,
                        const std::vector<float>& arr_yh,
                        const std::vector<float>& arr_zh,
                        const std::vector<float>& arr_x,
                        const std::vector<float>& arr_y,
                        const std::vector<float>& arr_z,
                        const std::vector<float>& arr_dz,
                        const std::vector<float>& field_phi,
                        std::vector<float>& field_atm_net_phi,
                        std::vector<float>& field_sfc_net_phi,
                        std::vector<float>& field_TOA_net_phi,
                        const long int N,
                        const int domain_section,
                        const std::string& INTERCELL_TECHNIQUE,
                        const bool Pesc_mode,
                        const bool enable_scattering);