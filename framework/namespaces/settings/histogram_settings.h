//
// Created by Alon Sportes on 07/04/2025.
//

#ifndef HISTOGRAM_FUNCTIONS_H
#define HISTOGRAM_FUNCTIONS_H

#include <iostream>

#include "../../classes/VariableHistograms/VariableHistograms.h"

namespace histogram_settings {

std::vector<YVarConfig> yvars = {
    {"P_n", "Neutron Momentum", "P_{n} [GeV/c]", 0.0, 1.5, 50},
    {"theta_n", "#theta_{n}", "#theta_{n} [#circ]", 0.0, 180.0, 50},
    {"phi_n", "#phi_{n}", "#phi_{n} [#circ]", -180.0, 180.0, 50},

    {"P_miss", "Missing Momentum", "P_{miss} [GeV/c]", 0.0, 1.5, 50},
    {"theta_miss", "#theta_{miss}", "#theta_{miss} [#circ]", 0.0, 180.0, 50},
    {"phi_miss", "#phi_{miss}", "#phi_{miss} [#circ]", -180.0, 180.0, 50},

    {"E_miss", "E_{miss}", "E_{miss} = E_{beam} + m_{D} - |#font[62]{P}_{e}| - E_{p} [GeV]", 0.5, 1.5, 50},
    {"M_miss", "M_{miss}", "M_{miss} = #left[#left(E_{beam} + m_{D} - |#font[62]{P}_{e}| - E_{p}#right)^{2} - #font[62]{P}_{miss}^{2}#right]^{1/2}  #left[GeV/c^{2}#right]", 0.65, 1.25, 50},

    {"dpp", "(|#font[62]{P}_{miss}|-|#font[62]{P}_{n}|)/P_{miss}", "(|#font[62]{P}_{miss}|-|#font[62]{P}_{n}|)/P_{miss}", -3.0, 1.0, 50},
    {"theta_n_miss", "#theta_{n,miss}", "#theta_{n,miss} [#circ]", 0.0, 180.0, 50},

    {"beta_n", "#beta_{n}", "#beta_{n}", -0.1, 1.1, 50},
    {"path_n", "Path length", "Path length [cm]", 0.0, 100.0, 50},
    {"ToF_n", "Neutron ToF", "t_{ToF}^{n} [ns]", 0.0, 20.0, 50},

    {"E_p", "E_{p}", "E_{p} = #left[m^{2} + #font[62]{P}_{p}^{2}#right]^{1/2} [GeV]", 0.5, 1.5, 50},
    {"theta_p", "E_{p}", "E_{p} = #left[m^{2} + #font[62]{P}_{p}^{2}#right]^{1/2} [GeV]", 0.5, 1.5, 50},
    {"E_p", "E_{p}", "E_{p} = #left[m^{2} + #font[62]{P}_{p}^{2}#right]^{1/2} [GeV]", 0.5, 1.5, 50},

    {"nSector", "Neutron Sector Number", "Sector Number", 0.5, 24.5, 24},

    {"Edep_CND1", "E^{CND,1}_{dep}", "E^{CND,1}_{dep} [MeV]", 0.0, 100.0, 50},
    {"Edep_CND2", "E^{CND,2}_{dep}", "E^{CND,2}_{dep} [MeV]", 0.0, 100.0, 50},
    {"Edep_CND3", "E^{CND,3}_{dep}", "E^{CND,3}_{dep} [MeV]", 0.0, 100.0, 50},
    {"Edep_CND", "E^{CND}_{dep}", "E^{CND}_{dep} [MeV]", 0.0, 100.0, 50},

    {"Size_CND1", "Size(CND1)", "Size(CND1)", -0.5, 4.5, 5},
    {"Size_CND2", "Size(CND2)", "Size(CND2)", -0.5, 4.5, 5},
    {"Size_CND3", "Size(CND3)", "Size(CND3)", -0.5, 4.5, 5},

    {"LayerMult_CND1", "LayerMult(CND1)", "LayerMult(CND1)", -0.5, 4.5, 5},
    {"LayerMult_CND2", "LayerMult(CND2)", "LayerMult(CND2)", -0.5, 4.5, 5},
    {"LayerMult_CND3", "LayerMult(CND3)", "LayerMult(CND3)", -0.5, 4.5, 5},
    {"LayerMult_CND", "LayerMult(CND)", "LayerMult(CND)", -0.5, 4.5, 5}

};

// VariableHistograms vh(
//     "Edep_CND",
//     "Total Neutron Energy Deposition in the CND",
//     "E^{CND}_{dep} [MeV]",
//     "epCDn",
//     yvars
// );

// double Edep = 42.0;
// std::map<std::string, double> values = {
//     {"P_n", 0.8},
//     {"theta_n", 24},
//     {"phi_n", -135},
//     {"dpp", -0.2},
//     {"beta_n", 0.85}
// };

// vh.FillHistograms(Edep, values);

};  // namespace histogram_settings

#endif  // HISTOGRAM_FUNCTIONS_H
