#include <vector>
#include <string>

#include "util.h"


void photon_propagation_incl_scattering(const AliasTable_double& aliastable,
                        FastRNG& rng,
                        const std::vector<double>& field_kext,
                        const std::vector<double>& field_sfc_eps,
                        const std::vector<double>& field_SSA,
                        const std::vector<double>& field_ASY,
                        const std::vector<double>& arr_xh,
                        const std::vector<double>& arr_yh,
                        const std::vector<double>& arr_zh,
                        const std::vector<double>& arr_x,
                        const std::vector<double>& arr_y,
                        const std::vector<double>& arr_z,
                        const std::vector<double>& arr_dz,
                        std::vector<double>& field_phi,
                        std::vector<double>& field_atm_net_phi,
                        std::vector<double>& field_sfc_net_phi,
                        std::vector<double>& field_TOA_net_phi,
                        const int N,
                        const int domain_section,
                        const std::string& INTERCELL_TECHNIQUE,
                        const bool Pesc_mode);


void photon_propagation(const AliasTable_double& aliastable,
                        FastRNG& rng,
                        const std::vector<double>& field_kext,
                        const std::vector<double>& field_sfc_eps,
                        const std::vector<double>& arr_xh,
                        const std::vector<double>& arr_yh,
                        const std::vector<double>& arr_zh,
                        const std::vector<double>& arr_x,
                        const std::vector<double>& arr_y,
                        const std::vector<double>& arr_z,
                        const std::vector<double>& arr_dz,
                        std::vector<double>& field_phi,
                        std::vector<double>& field_atm_net_phi,
                        std::vector<double>& field_sfc_net_phi,
                        std::vector<double>& field_TOA_net_phi,
                        const int N,
                        const int domain_section,
                        const std::string& INTERCELL_TECHNIQUE,
                        const bool Pesc_mode);