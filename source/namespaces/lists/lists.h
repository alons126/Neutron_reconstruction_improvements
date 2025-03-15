//
// Created by Alon Sportes on 13/03/2025.
//

#ifndef LISTS_H
#define LISTS_H

#include <map>

namespace lists {
std::map<std::string, std::map<std::string, std::string>> VariableNames = {
    //
    {"P_e", {{"VarName", "P_e"}, {"VarLabel", "|#font[62]{P}_{e}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"theta_e", {{"VarName", "theta_e"}, {"VarLabel", "#theta_{e}"}, {"VarDim", "#left[#circ#right]"}}},
    {"phi_e", {{"VarName", "phi_e"}, {"VarLabel", "#phi_{e}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"P_p", {{"VarName", "P_p"}, {"VarLabel", "|#font[62]{P}_{p}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"theta_p", {{"VarName", "theta_p"}, {"VarLabel", "#theta_{p}"}, {"VarDim", "#left[#circ#right]"}}},
    {"phi_p", {{"VarName", "phi_p"}, {"VarLabel", "#phi_{p}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"P_n", {{"VarName", "P_n"}, {"VarLabel", "|#font[62]{P}_{n}|"}, {"VarDim", "#left[GeV/c##right]"}}},
    {"theta_n", {{"VarName", "theta_n"}, {"VarLabel", "#theta_{n}"}, {"VarDim", "#left[#circ#right]"}}},
    {"phi_n", {{"VarName", "phi_n"}, {"VarLabel", "#phi_{n}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"P_pFD", {{"VarName", "P_pFD"}, {"VarLabel", "|#font[62]{P}_{pFD}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"theta_pFD", {{"VarName", "theta_pFD"}, {"VarLabel", "#theta_{pFD}"}, {"VarDim", "#left[#circ#right]"}}},
    {"phi_pFD", {{"VarName", "phi_pFD"}, {"VarLabel", "#phi_{pFD}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"P_nFD", {{"VarName", "P_nFD"}, {"VarLabel", "|#font[62]{P}_{nFD}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"theta_nFD", {{"VarName", "theta_nFD"}, {"VarLabel", "#theta_{nFD}"}, {"VarDim", "#left[#circ#right]"}}},
    {"phi_nFD", {{"VarName", "phi_nFD"}, {"VarLabel", "#phi_{nFD}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"P_pCD", {{"VarName", "P_pCD"}, {"VarLabel", "|#font[62]{P}_{pCD}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"theta_pCD", {{"VarName", "theta_pCD"}, {"VarLabel", "#theta_{pCD}"}, {"VarDim", "#left[#circ#right]"}}},
    {"phi_pCD", {{"VarName", "phi_pCD"}, {"VarLabel", "#phi_{pCD}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"P_nCD", {{"VarName", "P_nCD"}, {"VarLabel", "|#font[62]{P}_{nCD}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"theta_nCD", {{"VarName", "theta_nCD"}, {"VarLabel", "#theta_{nCD}"}, {"VarDim", "#left[#circ#right]"}}},
    {"phi_nCD", {{"VarName", "phi_nCD"}, {"VarLabel", "#phi_{nCD}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"momentum", {{"VarName", "momentum"}, {"VarLabel", "momentum"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"leading_nuc_momentum", {{"VarName", "leading_nuc_momentum"}, {"VarLabel", "|#font[62]{P}_{nucL}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"recoil_nuc_momentum", {{"VarName", "recoil_nuc_momentum"}, {"VarLabel", "|#font[62]{P}_{nucR}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    //
    {"total_3momentum", {{"VarName", "total_3momentum"}, {"VarLabel", "Total 3-momentum"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"theta_tot", {{"VarName", "theta_tot"}, {"VarLabel", "#theta_{tot}"}, {"VarDim", "#left[#circ#right]"}}},
    {"phi_tot", {{"VarName", "phi_tot"}, {"VarLabel", "#phi_{tot}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"relative_3momentum", {{"VarName", "relative_3momentum"}, {"VarLabel", "Relative 3-momentum"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"theta_rel", {{"VarName", "theta_rel"}, {"VarLabel", "#theta_{rel}"}, {"VarDim", "#left[#circ#right]"}}},
    {"phi_rel", {{"VarName", "phi_rel"}, {"VarLabel", "#phi_{rel}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"total_4momentum", {{"VarName", "total_4momentum"}, {"VarLabel", "Total 4-momentum"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"relative_4momentum", {{"VarName", "relative_4momentum"}, {"VarLabel", "Relative 4-momentum"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"P_tot_minus_q", {{"VarName", "P_tot_minus_q"}, {"VarLabel", "|#font[62]{P}_{tot}-#font[62]{q}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"Opening_ang_P_nucL_minus_q_nucR",
     {{"VarName", "Opening_ang_P_nucL_minus_q_nucR"}, {"VarLabel", "#theta_{#font[62]{P}_{nucL}-#font[62]{q},#font[62]{P}_{nucR}}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"P_tot_minus_q", {{"VarName", "P_tot_minus_q"}, {"VarLabel", "#vec{P}_{tot}-#vec{q}"}, {"VarDim", "#left[#circ#right]"}}},
    //
    {"q", {{"VarName", "q"}, {"VarLabel", "|#font[62]{q}|"}, {"VarDim", "#left[GeV/c#right]"}}},
    {"omega", {{"VarName", "omega"}, {"VarLabel", "#omega"}, {"VarDim", "#left[GeV#right]"}}},
    {"Q2", {{"VarName", "Q2"}, {"VarLabel", "Q^{2}"}, {"VarDim", "#left[GeV^{2}/c^{2}#right]"}}},
    {"xB", {{"VarName", "xB"}, {"VarLabel", "x_{B}"}, {"VarDim", ""}}},
    {"W", {{"VarName", "W"}, {"VarLabel", "W"}, {"VarDim", "#left[GeV/c^{2}#right]"}}},
    //
    {"tof", {{"VarName", "ToF"}, {"VarLabel", "t_{ToF}"}, {"VarDim", "#left[ns#right]"}}},
    //
    {"E_e", {{"VarName", "E_e"}, {"VarLabel", "E_{e}"}, {"VarDim", "#left[GeV#right]"}}},
    {"Ecal", {{"VarName", "Ecal"}, {"VarLabel", "E_{cal}"}, {"VarDim", "#left[GeV#right]"}}},
    {"Ecal_ext_1N", {{"VarName", "Ecal"}, {"VarLabel", "E_{cal} = E_{e} + T_{nuc}"}, {"VarDim", "#left[GeV#right]"}}},
    {"Ecal_ext_2N", {{"VarName", "Ecal"}, {"VarLabel", "E_{cal} = E_{e} + T_{nuc,1} + T_{nuc,2}"}, {"VarDim", "#left[GeV#right]"}}},
    //
    {"deltaP_T_tot", {{"VarName", "deltaP_T_tot"}, {"VarLabel", "#deltaP_{T,tot}"}, {"VarDim", ""}}},
    {"deltaP_T_L", {{"VarName", "deltaP_T_L"}, {"VarLabel", "#deltaP_{T,L}"}, {"VarDim", ""}}},
    {"deltaAlpha_T_tot", {{"VarName", "deltaAlpha_T_tot"}, {"VarLabel", "#delta#alpha_{T,tot}"}, {"VarDim", ""}}},
};

std::map<std::string, std::map<std::string, double>> HistogramLimits;

}  // namespace lists

#endif  // LISTS_H
