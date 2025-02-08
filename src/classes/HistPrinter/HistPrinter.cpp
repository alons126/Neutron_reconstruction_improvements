//
// Created by Alon Sportes on 22/01/2025.
//

#include "HistPrinter.h"

// PrintPage function --------------------------------------------------------------------------------------------------------------------------------------------------------

void HistPrinter::PrintPage(vector<TH1 *> &HistoList, const std::string &PageTitle, TCanvas *myText, char fileName[100], TLatex titles, TLatex text, const std::string &Constraint1,
                            const std::string &Constraint2) {
    vector<TString> summary_table_Line;

    if (PageTitle == "Manual Veto Plots") {
        if (Constraint2 == "") {
            titles.DrawLatex(0.05, 0.9, "Manual Veto Plots");
        } else {
            titles.DrawLatex(0.05, 0.9, ("Manual Veto Plots - " + Constraint2).c_str());
        }

        text.DrawLatex(0.05, 0.7, "#diamond  #font[12]{(e,e'p)} Cuts:");

        if (Constraint1 == "") {
            text.DrawLatex(0.10, 0.6, ("#bullet  " + to_string_with_precision(Num_of_e_cut, 0) + " electron").c_str());
            text.DrawLatex(0.10, 0.5, ("#bullet  " + to_string_with_precision(Num_of_p_cut, 0) + " proton in CD or FD").c_str());
            text.DrawLatex(0.10, 0.4, "#bullet  Any number of neutrons in CND");
            text.DrawLatex(0.10, 0.3, "#bullet  Only particles with pdg=2112,11,2212,0,22 in event");
        } else if (Constraint1 == "CD") {
            text.DrawLatex(0.10, 0.6, ("#bullet  " + to_string_with_precision(Num_of_e_cut, 0) + " electron").c_str());
            text.DrawLatex(0.10, 0.5, ("#bullet  " + to_string_with_precision(Num_of_p_cut, 0) + " proton in CD").c_str());
            text.DrawLatex(0.10, 0.4, "#bullet  Any number of neutrons in CND");
            text.DrawLatex(0.10, 0.3, "#bullet  Only particles with pdg=2112,11,2212,0,22 in event");
        } else if (Constraint1 == "FD") {
            text.DrawLatex(0.10, 0.6, ("#bullet  " + to_string_with_precision(Num_of_e_cut, 0) + " electron").c_str());
            text.DrawLatex(0.10, 0.5, ("#bullet  " + to_string_with_precision(Num_of_p_cut, 0) + " proton in FD").c_str());
            text.DrawLatex(0.10, 0.4, "#bullet  Any number of neutrons in CND");
            text.DrawLatex(0.10, 0.3, "#bullet  Only particles with pdg=2112,11,2212,0,22 in event");
        }

        myText->Print(fileName, "pdf");
        myText->Clear();
    } else if (PageTitle == "PID Plots") {
        titles.DrawLatex(0.05, 0.9, "PID Plots");

        if (Constraint1 == "") {
            text.DrawLatex(0.05, 0.8, "#diamond  CD protons:");
            text.DrawLatex(0.10, 0.75, ("#bullet  #font[12]{#lbarV_{z}^{p} - V_{z}^{e}#lbar #leq " + to_string_with_precision(dVz_pCD_cut, 1) + "} cm").c_str());
            text.DrawLatex(0.10, 0.7, ("#bullet  #font[12]{" + to_string_with_precision(P_pCD_lcut, 1) + " #leq P_{p} #leq " + to_string_with_precision(P_pCD_ucut, 1) + "} GeV/c").c_str());
            text.DrawLatex(0.10, 0.65, ("#bullet  #font[12]{" + to_string_with_precision(pCD_chi2_lcut, 1) + " #leq #chi^{2} #leq " + to_string_with_precision(pCD_chi2_ucut, 1) + "}").c_str());

            text.DrawLatex(0.05, 0.55, "#diamond  FD protons:");
            text.DrawLatex(0.10, 0.5, ("#bullet  #font[12]{#lbarV_{z}^{p} - V_{z}^{e}#lbar #leq " + to_string_with_precision(dVz_pFD_cut, 1) + "} cm").c_str());
            text.DrawLatex(0.10, 0.45, ("#bullet  #font[12]{" + to_string_with_precision(P_pFD_lcut, 1) + " #leq P_{p} #leq " + to_string_with_precision(P_pFD_ucut, 1) + "} GeV/c").c_str());
            text.DrawLatex(0.10, 0.4, ("#bullet  #font[12]{" + to_string_with_precision(pFD_chi2_lcut, 1) + " #leq #chi^{2} #leq " + to_string_with_precision(pFD_chi2_ucut, 1) + "}").c_str());
        } else if (Constraint1 == "CD") {
            text.DrawLatex(0.05, 0.8, "#diamond  CD protons:");
            text.DrawLatex(0.10, 0.7, ("#bullet  #font[12]{#lbarV_{z}^{p} - V_{z}^{e}#lbar #leq " + to_string_with_precision(dVz_pCD_cut, 1) + "} cm").c_str());
            text.DrawLatex(0.10, 0.6, ("#bullet  #font[12]{" + to_string_with_precision(P_pCD_lcut, 1) + " #leq P_{p} #leq " + to_string_with_precision(P_pCD_ucut, 1) + "} GeV/c").c_str());
            text.DrawLatex(0.10, 0.5, ("#bullet  #font[12]{" + to_string_with_precision(pCD_chi2_lcut, 1) + " #leq #chi^{2} #leq " + to_string_with_precision(pCD_chi2_ucut, 1) + "}").c_str());
        } else if (Constraint1 == "FD") {
            text.DrawLatex(0.05, 0.8, "#diamond  FD protons:");
            text.DrawLatex(0.10, 0.7, ("#bullet  #font[12]{#lbarV_{z}^{p} - V_{z}^{e}#lbar #leq " + to_string_with_precision(dVz_pFD_cut, 1) + "} cm").c_str());
            text.DrawLatex(0.10, 0.6, ("#bullet  #font[12]{" + to_string_with_precision(P_pFD_lcut, 1) + " #leq P_{p} #leq " + to_string_with_precision(P_pFD_ucut, 1) + "} GeV/c").c_str());
            text.DrawLatex(0.10, 0.5, ("#bullet  #font[12]{" + to_string_with_precision(pFD_chi2_lcut, 1) + " #leq #chi^{2} #leq " + to_string_with_precision(pFD_chi2_ucut, 1) + "}").c_str());
        }

        myText->Print(fileName, "pdf");
        myText->Clear();
    } else if (PageTitle == "Plots with basic cuts") {
        titles.DrawLatex(0.05, 0.9, "Plots with basic cuts");
        text.DrawLatex(0.05, 0.8, "#diamond  Missing variables cuts:");
        text.DrawLatex(0.10, 0.75, ("#bullet  #font[12]{" + to_string_with_precision(P_miss_lcut, 2) + " #leq P_{miss} #leq " + to_string_with_precision(P_miss_ucut, 2) + "} GeV/c").c_str());
        text.DrawLatex(0.10, 0.7,
                       ("#bullet  #font[12]{" + to_string_with_precision(Theta_miss_lcut, 0) + "#circ #leq #theta_{miss} #leq " + to_string_with_precision(Theta_miss_ucut, 0) + "#circ}").c_str());
        text.DrawLatex(0.10, 0.65, ("#bullet  #font[12]{" + to_string_with_precision(M_miss_lcut, 2) + " #leq M_{miss} #leq " + to_string_with_precision(M_miss_ucut, 2) + "} GeV/c^{2}").c_str());

        text.DrawLatex(0.05, 0.55, "#diamond  Neutron PID cuts:");
        text.DrawLatex(0.10, 0.5, ("#bullet  #font[12]{" + to_string_with_precision(Beta_n_lcut, 2) + " #leq #beta_{n} #leq " + to_string_with_precision(Beta_n_ucut, 2) + "}").c_str());
        text.DrawLatex(0.10, 0.45, ("#bullet  #font[12]{" + to_string_with_precision(Theta_n_lcut, 0) + "#circ #leq #theta_{n} #leq " + to_string_with_precision(Theta_n_ucut, 0) + "#circ}").c_str());
        text.DrawLatex(0.10, 0.4, ("#bullet  Status = " + to_string_with_precision(Status_n_cut, 0) + " (no double-hits)").c_str());
        text.DrawLatex(0.10, 0.35, "#bullet  CTOF veto (neutron cluster does not have a CTOF hit)");

        myText->Print(fileName, "pdf");
        myText->Clear();
    } else if (PageTitle == "Definition of neutrons in veto steps") {
        titles.DrawLatex(0.05, 0.9, "Definition of neutrons in veto steps");
        text.DrawLatex(0.05, 0.8, "#diamond  Good neutrons definition:");
        text.DrawLatex(0.10, 0.7, ("#bullet  #font[12]{#theta_{n,miss} #leq " + to_string_with_precision(GN_theta_n_miss_ucut, 0) + "#circ}, and").c_str());
        text.DrawLatex(0.10, 0.6,
                       ("#bullet  #font[12]{" + to_string_with_precision(GN_dpp_lcut, 1) + " #leq #left(#lbar#vec{P}_{miss}#lbar - #lbar#vec{P}_{n}#lbar#right)/P_{miss} #leq " +
                        to_string_with_precision(GN_dpp_ucut, 1) + "}")
                           .c_str());
        text.DrawLatex(0.05, 0.5, "#diamond  Bad neutrons definition:");
        text.DrawLatex(0.10, 0.4, ("#bullet  #font[12]{#theta_{n,miss} #geq " + to_string_with_precision(BN_theta_n_miss_lcut, 0) + "#circ}, or").c_str());
        text.DrawLatex(0.10, 0.3, ("#bullet  #font[12]{#left(#lbar#vec{P}_{miss}#lbar - #lbar#vec{P}_{n}#lbar#right)/P_{miss} #leq " + to_string_with_precision(BN_dpp_ucut, 1) + "}").c_str());

        myText->Print(fileName, "pdf");
        myText->Clear();
    } else if (PageTitle == "Step0 Plots") {
        text.DrawLatex(0.05, 0.8, "#diamond  Step0 cuts (included in Step1):");
        text.DrawLatex(0.10, 0.7, ("#bullet  #font[12]{#lbar#beta_{n} - L/(t_{ToF,n}c)#lbar #leq " + to_string_with_precision(dBeta_n_cut) + "}").c_str());
        text.DrawLatex(0.10, 0.6, ("#bullet  #font[12]{" + to_string_with_precision(Vz_n_lcut) + " #leq V_{hit,z} #leq " + to_string_with_precision(Vz_n_ucut) + "} cm").c_str());
        text.DrawLatex(0.10, 0.5, ("#bullet  #font[12]{" + to_string_with_precision(ToF_n_lcut) + " #leq t_{ToF,n} #leq " + to_string_with_precision(ToF_n_ucut) + "} ns").c_str());

        myText->Print(fileName, "pdf");
        myText->Clear();
    } else if (PageTitle == "Step1 Plots") {
        text.DrawLatex(0.05, 0.8, "#diamond  Step0 cuts (included in Step1):");
        text.DrawLatex(0.10, 0.7, ("#bullet  #font[12]{#lbar#beta_{n} - L/(t_{ToF,n}c)#lbar #leq " + to_string_with_precision(dBeta_n_cut) + "}").c_str());
        text.DrawLatex(0.10, 0.6, ("#bullet  #font[12]{" + to_string_with_precision(Vz_n_lcut) + " #leq V_{hit,z} #leq " + to_string_with_precision(Vz_n_ucut) + "} cm").c_str());
        text.DrawLatex(0.10, 0.5, ("#bullet  #font[12]{" + to_string_with_precision(ToF_n_lcut) + " #leq t_{ToF,n} #leq " + to_string_with_precision(ToF_n_ucut) + "} ns").c_str());

        text.DrawLatex(0.05, 0.4, "#diamond  Step1 cuts:");
        text.DrawLatex(0.10, 0.3, ("#bullet  #font[12]{" + to_string_with_precision(Edep_CND_lcut, 0) + " #leq E_{dep}^{CND} #leq (#gamma_{n} - 1) m_{n}} MeV").c_str());

        myText->Print(fileName, "pdf");
        myText->Clear();
    } else if (PageTitle == "Step2 Plots") {
        text.DrawLatex(0.05, 0.8, "#diamond  Step0 cuts (included in Step2):");
        text.DrawLatex(0.10, 0.75, ("#bullet  #font[12]{#lbar#beta_{n} - L/(t_{ToF,n}c)#lbar #leq " + to_string_with_precision(dBeta_n_cut) + "}").c_str());
        text.DrawLatex(0.10, 0.7, ("#bullet  #font[12]{" + to_string_with_precision(Vz_n_lcut) + " #leq V_{hit,z} #leq " + to_string_with_precision(Vz_n_ucut) + "} cm").c_str());
        text.DrawLatex(0.10, 0.65, ("#bullet  #font[12]{" + to_string_with_precision(ToF_n_lcut) + " #leq t_{ToF,n} #leq " + to_string_with_precision(ToF_n_ucut) + "} ns").c_str());

        text.DrawLatex(0.05, 0.55, "#diamond  Step1 cuts (included in Step2):");
        text.DrawLatex(0.10, 0.5, ("#bullet  #font[12]{" + to_string_with_precision(Edep_CND_lcut, 0) + " #leq E_{dep}^{CND} #leq (#gamma_{n} - 1) m_{n}} MeV").c_str());

        text.DrawLatex(0.05, 0.4, "#diamond  Step2 cuts:");
        text.DrawLatex(0.10, 0.3, ("#bullet  Cluster size (= width) is " + to_string_with_precision(Cluster_size_cut, 0) + " hit").c_str());
        text.DrawLatex(0.10, 0.25, "#bullet  Layer multiplicity:");
        text.DrawLatex(0.15, 0.2, ("#Box  Hit in CND1 #rightarrow layer multiplicity = " + to_string_with_precision(CND1_LayerMult_cut, 0)).c_str());
        text.DrawLatex(0.15, 0.15, ("#Box  Hit in CND2 or CND3 #rightarrow layer multiplicity = up to " + to_string_with_precision(CND2andCND3_LayerMult_ucut, 0)).c_str());
        text.DrawLatex(0.10, 0.35, "#bullet  No nearby hits associated with the charged particle track");

        myText->Print(fileName, "pdf");
        myText->Clear();
    }
}

// SummaryTablePlotter function ----------------------------------------------------------------------------------------------------------------------------------------------

void HistPrinter::GenerateSummaryTable(TCanvas *myTable, vector<TH1 *> HistoList, string Constraint1, string Constraint2) {
    if (First_table_epCDn_generation && (Constraint1 == "" || Constraint1 == "CD")) {
        myTable->SetTopMargin(0.15);

        /* Before Step Cuts */
        // Before Step Cuts overall:
        vector<TString> summary_table_bfSteps_epCDn_1stLine = {"#splitline{Before}{step cuts}", to_string_with_precision(Num_of_goodN_bfSteps_epCDn, 0),
                                                               to_string_with_precision(Num_of_badN_bfSteps_epCDn, 0), "--", "--"};
        summary_table_bfSteps_epCDn.push_back(summary_table_bfSteps_epCDn_1stLine);

        for (int i = 0; i < summary_table_bfSteps_epCDn.size(); i++) { table_epCDn.push_back(summary_table_bfSteps_epCDn.at(i)); }

        /* Step0 */
        if (Apply_Step0_Cuts) {
            // Step0 overall:
            double Num_of_goodN_Step0_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step0_epCDn");
            double Num_of_badN_Step0_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step0_epCDn");

            vector<TString> summary_table_Step0_epCDn_1stLine = {"Step0", to_string_with_precision(Num_of_goodN_Step0_epCDn, 0), to_string_with_precision(Num_of_badN_Step0_epCDn, 0),
                                                                 to_string_with_precision(Num_of_goodN_Step0_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                                 to_string_with_precision(Num_of_goodN_Step0_epCDn / (Num_of_goodN_Step0_epCDn + Num_of_badN_Step0_epCDn), PurAccuracy)};
            summary_table_Step0_epCDn.push_back(summary_table_Step0_epCDn_1stLine);

            // dBeta_n test (Step0):
            double Num_of_goodN_Step0_dBeta_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_dBeta_n_Step0_goodN_BCTest_epCDn");
            // double Num_of_badN_Step0_dBeta_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_dBeta_n_Step0_badN_BCTest_epCDn");
            double Num_of_goodN_Step0_dBeta_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_dBeta_n_Step0_goodN_ACTest_epCDn");
            double Num_of_badN_Step0_dBeta_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_dBeta_n_Step0_badN_ACTest_epCDn");

            summary_table_Step0_epCDn.push_back(
                {"#splitline{#Delta#beta_{n} cuts}{(Step0)}", to_string_with_precision(Num_of_goodN_Step0_dBeta_n_ACTest_epCDn, 0), to_string_with_precision(Num_of_badN_Step0_dBeta_n_ACTest_epCDn, 0),
                 //  to_string_with_precision(Num_of_goodN_Step0_dBeta_n_ACTest_epCDn / Num_of_goodN_Step0_dBeta_n_BCTest_epCDn, EffAccuracy),
                 to_string_with_precision(Num_of_goodN_Step0_dBeta_n_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                 to_string_with_precision(Num_of_goodN_Step0_dBeta_n_ACTest_epCDn / (Num_of_goodN_Step0_dBeta_n_ACTest_epCDn + Num_of_badN_Step0_dBeta_n_ACTest_epCDn), PurAccuracy)});

            // Vz_n test (Step0):
            double Num_of_goodN_Step0_Vz_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Vz_n_Step0_goodN_BCTest_epCDn");
            // double Num_of_badN_Step0_Vz_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Vz_n_Step0_badN_BCTest_epCDn");
            double Num_of_goodN_Step0_Vz_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Vz_n_Step0_goodN_ACTest_epCDn");
            double Num_of_badN_Step0_Vz_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Vz_n_Step0_badN_ACTest_epCDn");

            summary_table_Step0_epCDn.push_back(
                {"#splitline{V_{hit,z} cuts}{(Step0)}", to_string_with_precision(Num_of_goodN_Step0_Vz_n_ACTest_epCDn, 0), to_string_with_precision(Num_of_badN_Step0_Vz_n_ACTest_epCDn, 0),
                 //  to_string_with_precision(Num_of_goodN_Step0_Vz_n_ACTest_epCDn / Num_of_goodN_Step0_Vz_n_BCTest_epCDn, EffAccuracy),
                 to_string_with_precision(Num_of_goodN_Step0_Vz_n_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                 to_string_with_precision(Num_of_goodN_Step0_Vz_n_ACTest_epCDn / (Num_of_goodN_Step0_Vz_n_ACTest_epCDn + Num_of_badN_Step0_Vz_n_ACTest_epCDn), PurAccuracy)});

            // ToF_n test (Step0):
            double Num_of_goodN_Step0_ToF_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_ToF_n_Step0_goodN_BCTest_epCDn");
            // double Num_of_badN_Step0_ToF_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_ToF_n_Step0_badN_BCTest_epCDn");
            double Num_of_goodN_Step0_ToF_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_ToF_n_Step0_goodN_ACTest_epCDn");
            double Num_of_badN_Step0_ToF_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_ToF_n_Step0_badN_ACTest_epCDn");

            summary_table_Step0_epCDn.push_back(
                {"#splitline{t_{ToF,n} cuts}{(Step0)}", to_string_with_precision(Num_of_goodN_Step0_ToF_n_ACTest_epCDn, 0), to_string_with_precision(Num_of_badN_Step0_ToF_n_ACTest_epCDn, 0),
                 //  to_string_with_precision(Num_of_goodN_Step0_ToF_n_ACTest_epCDn / Num_of_goodN_Step0_ToF_n_BCTest_epCDn, EffAccuracy),
                 to_string_with_precision(Num_of_goodN_Step0_ToF_n_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                 to_string_with_precision(Num_of_goodN_Step0_ToF_n_ACTest_epCDn / (Num_of_goodN_Step0_ToF_n_ACTest_epCDn + Num_of_badN_Step0_ToF_n_ACTest_epCDn), PurAccuracy)});

            for (int i = 0; i < summary_table_Step0_epCDn.size(); i++) { table_epCDn.push_back(summary_table_Step0_epCDn.at(i)); }
        }

        /* Step1 */
        // Step1 overall:
        double Num_of_goodN_Step1_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step1_epCDn");
        double Num_of_badN_Step1_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step1_epCDn");

        vector<TString> summary_table_Step1_epCDn_1stLine = {"Step1", to_string_with_precision(Num_of_goodN_Step1_epCDn, 0), to_string_with_precision(Num_of_badN_Step1_epCDn, 0),
                                                             // to_string_with_precision(Num_of_goodN_Step1_epCDn / Num_of_goodN_Step0_epCDn, EffAccuracy),
                                                             to_string_with_precision(Num_of_goodN_Step1_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                             to_string_with_precision(Num_of_goodN_Step1_epCDn / (Num_of_goodN_Step1_epCDn + Num_of_badN_Step1_epCDn), PurAccuracy)};
        summary_table_Step1_epCDn.push_back(summary_table_Step1_epCDn_1stLine);

        // Edep_CND test (Step1):
        double Num_of_goodN_Step1_Edep_CND_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Edep_CND_Step1_goodN_BCTest_epCDn");
        // double Num_of_badN_Step1_Edep_CND_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Edep_CND_Step1_badN_BCTest_epCDn");
        double Num_of_goodN_Step1_Edep_CND_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Edep_CND_Step1_goodN_ACTest_epCDn");
        double Num_of_badN_Step1_Edep_CND_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Edep_CND_Step1_badN_ACTest_epCDn");

        summary_table_Step1_epCDn.push_back(
            {"#splitline{E^{CND}_{dep} cuts}{(Step1)}", to_string_with_precision(Num_of_goodN_Step1_Edep_CND_ACTest_epCDn, 0), to_string_with_precision(Num_of_badN_Step1_Edep_CND_ACTest_epCDn, 0),
             //  to_string_with_precision(Num_of_goodN_Step1_Edep_CND_ACTest_epCDn / Num_of_goodN_Step1_Edep_CND_BCTest_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step1_Edep_CND_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step1_Edep_CND_ACTest_epCDn / (Num_of_goodN_Step1_Edep_CND_ACTest_epCDn + Num_of_badN_Step1_Edep_CND_ACTest_epCDn), PurAccuracy)});

        for (int i = 0; i < summary_table_Step1_epCDn.size(); i++) { table_epCDn.push_back(summary_table_Step1_epCDn.at(i)); }

        /* Step2 */
        // Step2 overall:
        double Num_of_goodN_Step2_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step2_epCDn");
        double Num_of_badN_Step2_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step2_epCDn");

        vector<TString> summary_table_Step2_epCDn_1stLine = {"Step2", to_string_with_precision(Num_of_goodN_Step2_epCDn, 0), to_string_with_precision(Num_of_badN_Step2_epCDn, 0),
                                                             // to_string_with_precision(Num_of_goodN_Step2_epCDn / Num_of_goodN_Step1_epCDn, EffAccuracy),
                                                             to_string_with_precision(Num_of_goodN_Step2_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                             to_string_with_precision(Num_of_goodN_Step2_epCDn / (Num_of_goodN_Step2_epCDn + Num_of_badN_Step2_epCDn), PurAccuracy)};
        // summary_table_Step2_epCDn.push_back(summary_table_Step2_epCDn_1stLine);
        summary_table2_Step2_epCDn.push_back(summary_table_Step2_epCDn_1stLine);
        summary_table3_Step2_epCDn.push_back(summary_table_Step2_epCDn_1stLine);

        // Size_CND1 test (Step2):
        double Num_of_goodN_Step2_Size_CND1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND1_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_Size_CND1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND1_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_Size_CND1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND1_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_Size_CND1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND1_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{Size(CND1) cuts}{(Step2)}", to_string_with_precision(Num_of_goodN_Step2_Size_CND1_ACTest_epCDn, 0), to_string_with_precision(Num_of_badN_Step2_Size_CND1_ACTest_epCDn, 0),
             //  to_string_with_precision(Num_of_goodN_Step2_Size_CND1_ACTest_epCDn / Num_of_goodN_Step2_Size_CND1_BCTest_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step2_Size_CND1_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step2_Size_CND1_ACTest_epCDn / (Num_of_goodN_Step2_Size_CND1_ACTest_epCDn + Num_of_badN_Step2_Size_CND1_ACTest_epCDn), PurAccuracy)});

        // Size_CND2 test (Step2):
        double Num_of_goodN_Step2_Size_CND2_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND2_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_Size_CND2_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND2_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_Size_CND2_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND2_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_Size_CND2_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND2_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{Size(CND2) cuts}{(Step2)}", to_string_with_precision(Num_of_goodN_Step2_Size_CND2_ACTest_epCDn, 0), to_string_with_precision(Num_of_badN_Step2_Size_CND2_ACTest_epCDn, 0),
             //  to_string_with_precision(Num_of_goodN_Step2_Size_CND2_ACTest_epCDn / Num_of_goodN_Step2_Size_CND2_BCTest_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step2_Size_CND2_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step2_Size_CND2_ACTest_epCDn / (Num_of_goodN_Step2_Size_CND2_ACTest_epCDn + Num_of_badN_Step2_Size_CND2_ACTest_epCDn), PurAccuracy)});

        // Size_CND3 test (Step2):
        double Num_of_goodN_Step2_Size_CND3_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND3_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_Size_CND3_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND3_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_Size_CND3_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND3_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_Size_CND3_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND3_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{Size(CND3) cuts}{(Step2)}", to_string_with_precision(Num_of_goodN_Step2_Size_CND3_ACTest_epCDn, 0), to_string_with_precision(Num_of_badN_Step2_Size_CND3_ACTest_epCDn, 0),
             //  to_string_with_precision(Num_of_goodN_Step2_Size_CND3_ACTest_epCDn / Num_of_goodN_Step2_Size_CND3_BCTest_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step2_Size_CND3_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step2_Size_CND3_ACTest_epCDn / (Num_of_goodN_Step2_Size_CND3_ACTest_epCDn + Num_of_badN_Step2_Size_CND3_ACTest_epCDn), PurAccuracy)});

        // LayerMult_CND1 test (Step2):
        double Num_of_goodN_Step2_LayerMult_CND1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND1_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_LayerMult_CND1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND1_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND1_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_LayerMult_CND1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND1_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{LayerMult(CND1)}{cuts (Step2)}", to_string_with_precision(Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn, 0),
             to_string_with_precision(Num_of_badN_Step2_LayerMult_CND1_ACTest_epCDn, 0),
             //  to_string_with_precision(Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn / Num_of_goodN_Step2_LayerMult_CND1_BCTest_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn / (Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn + Num_of_badN_Step2_LayerMult_CND1_ACTest_epCDn), PurAccuracy)});

        // LayerMult_CND2andCND3 test (Step2):
        double Num_of_goodN_Step2_LayerMult_CND2andCND3_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND2andCND3_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_LayerMult_CND2andCND3_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND2andCND3_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND2andCND3_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_LayerMult_CND2andCND3_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND2andCND3_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{LayerMult(CND1+2)}{cuts (Step2)}", to_string_with_precision(Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn, 0),
             to_string_with_precision(Num_of_badN_Step2_LayerMult_CND2andCND3_ACTest_epCDn, 0),
             //  to_string_with_precision(Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn / Num_of_goodN_Step2_LayerMult_CND2andCND3_BCTest_epCDn, EffAccuracy),
             to_string_with_precision(Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             to_string_with_precision(
                 Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn / (Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn + Num_of_badN_Step2_LayerMult_CND2andCND3_ACTest_epCDn), PurAccuracy)});

        for (int i = 0; i < summary_table2_Step2_epCDn.size(); i++) { table2_epCDn.push_back(summary_table2_Step2_epCDn.at(i)); }

        for (int k = 0; k < 7; k++) {
            double Num_of_goodN_of1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of1_pos_Step2_layer_" + to_string_with_precision(k - 3, 0) + "_goodN_BCTest_epCDn");
            double Num_of_goodN_of1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of1_pos_Step2_layer_" + to_string_with_precision(k - 3, 0) + "_goodN_ACTest_epCDn");
            double Num_of_badN_of1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of1_pos_Step2_layer_" + to_string_with_precision(k - 3, 0) + "_badN_ACTest_epCDn");

            summary_table3_Step2_epCDn.push_back({"#splitline{|#DeltaS_{n,+}|>1 & #DeltaL_{n,+} = " + to_string_with_precision(k - 3, 0) + " cuts}{(Step2)}",
                                                  to_string_with_precision(Num_of_goodN_of1_ACTest_epCDn, 0), to_string_with_precision(Num_of_badN_of1_ACTest_epCDn, 0),
                                                  to_string_with_precision(Num_of_goodN_of1_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                  to_string_with_precision(Num_of_goodN_of1_ACTest_epCDn / (Num_of_goodN_of1_ACTest_epCDn + Num_of_badN_of1_ACTest_epCDn), PurAccuracy)});

            double Num_of_goodN_of2_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of2_pos_Step2_layer_" + to_string_with_precision(k - 3, 0) + "_goodN_BCTest_epCDn");
            double Num_of_goodN_of2_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of2_pos_Step2_layer_" + to_string_with_precision(k - 3, 0) + "_goodN_ACTest_epCDn");
            double Num_of_badN_of2_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of2_pos_Step2_layer_" + to_string_with_precision(k - 3, 0) + "_badN_ACTest_epCDn");

            summary_table3_Step2_epCDn.push_back({"#splitline{|#DeltaS_{n,+}|>2 & #DeltaL_{n,+} = " + to_string_with_precision(k - 3, 0) + " cuts}{(Step2)}",
                                                  to_string_with_precision(Num_of_goodN_of2_ACTest_epCDn, 0), to_string_with_precision(Num_of_badN_of2_ACTest_epCDn, 0),
                                                  to_string_with_precision(Num_of_goodN_of2_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                  to_string_with_precision(Num_of_goodN_of2_ACTest_epCDn / (Num_of_goodN_of2_ACTest_epCDn + Num_of_badN_of2_ACTest_epCDn), PurAccuracy)});
        }

        for (int i = 0; i < summary_table3_Step2_epCDn.size(); i++) { table3_epCDn.push_back(summary_table3_Step2_epCDn.at(i)); }

        // Prevent the generation of multiple lines:
        First_table_epCDn_generation = false;
    }
}

// SummaryTablePlotter function ----------------------------------------------------------------------------------------------------------------------------------------------

void HistPrinter::SummaryTablePlotter(int n_col, int n_row, TCanvas *myCanvas, TCanvas *myText, TCanvas *myTable, vector<TH1 *> HistoList, TLatex titles, TLatex text, char fileName[100],
                                      string PDFFile, string Constraint1, string Constraint2, bool LogScale2D) {
    if (Constraint1 == "" || Constraint1 == "CD") {
        myTable->SetTopMargin(0.15);

        titles.SetTextSize(0.05);
        // titles.SetTextSize(0.065);

        GenerateSummaryTable(myTable, HistoList, Constraint1, Constraint2);

        // Draw a frame without axis numbers and ticks
        TH2F *frame1_epCDn = new TH2F("frame1_epCDn", "", summary_table_title.size(), 0, summary_table_title.size(), table_epCDn.size(), 0, table_epCDn.size());
        frame1_epCDn->SetStats(0);                   // Disable statistics box
        frame1_epCDn->GetXaxis()->SetLabelSize(0);   // Remove x-axis labels
        frame1_epCDn->GetXaxis()->SetTickLength(0);  // Remove x-axis ticks
        frame1_epCDn->GetYaxis()->SetLabelSize(0);   // Remove y-axis labels
        frame1_epCDn->GetYaxis()->SetTickLength(0);  // Remove y-axis ticks
        frame1_epCDn->Draw();

        titles.DrawLatexNDC(0.05, 0.9, "Step by step statistics - CD proton (table 1 of 3)");

        TPave *row1 = new TPave(0, (table_epCDn.size() - 1), summary_table_title.size(), table_epCDn.size(), 0, "br");  // Row 1 (top row)
        row1->SetFillColor(kGray + 1);
        // row1->SetFillColor(kAzure - 9);
        row1->SetFillStyle(1001);
        row1->SetLineColor(0);
        row1->Draw();

        // Loop over rows and columns to position text
        for (int i = 0; i < table_epCDn.size(); i++) {
            std::string cellTemp = table_epCDn.at(i).at(0).Data();

            if ((!findSubstring(cellTemp, "}{(Step") && !findSubstring(cellTemp, "}{cuts")) && (i > 0)) {
                TPave *row = new TPave(0, (table_epCDn.size() - i - 1), summary_table_title.size(), (table_epCDn.size() - i), 0, "br");
                row->SetFillColor(kAzure - 9);
                row->SetFillStyle(1001);
                row->SetLineColor(0);
                row->Draw();
            }
        }

        // Create an instance of TLatex
        TLatex latex_epCDn;

        // Set text alignment and font size
        latex_epCDn.SetTextAlign(22);  // Centered
        // latex_epCDn.SetTextSize(0.01);
        latex_epCDn.SetTextSize(0.02);

        // Loop over rows and columns to position text
        for (int i = 0; i < table_epCDn.size(); i++) {
            // table_epCDn.size() rows
            for (int j = 0; j < summary_table_title.size(); j++) {
                // summary_table_title.size() columns
                latex_epCDn.DrawLatex(j + 0.5, table_epCDn.size() - i - 0.5, table_epCDn[i][j]);  // Adjust positioning
            }
        }

        // Add gridlines for clarity (optional)
        for (int i = 0; i <= table_epCDn.size(); i++) {
            // Horizontal lines
            TLine *line_epCDn = new TLine(0, i, summary_table_title.size(), i);
            line_epCDn->SetLineStyle(2);
            line_epCDn->Draw();
        }
        for (int j = 0; j <= summary_table_title.size(); j++) {
            // Vertical lines
            TLine *line_epCDn = new TLine(j, 0, j, table_epCDn.size());
            line_epCDn->SetLineStyle(2);
            line_epCDn->Draw();
        }

        myTable->Print(fileName, "pdf");
        myTable->Clear();

        TH2F *frame2_epCDn = new TH2F("frame2_epCDn", "", summary_table_title.size(), 0, summary_table_title.size(), table2_epCDn.size(), 0, table2_epCDn.size());
        frame2_epCDn->SetStats(0);                   // Disable statistics box
        frame2_epCDn->GetXaxis()->SetLabelSize(0);   // Remove x-axis labels
        frame2_epCDn->GetXaxis()->SetTickLength(0);  // Remove x-axis ticks
        frame2_epCDn->GetYaxis()->SetLabelSize(0);   // Remove y-axis labels
        frame2_epCDn->GetYaxis()->SetTickLength(0);  // Remove y-axis ticks
        frame2_epCDn->Draw();

        titles.DrawLatexNDC(0.05, 0.9, "Step by step statistics - CD proton (table 2 of 3)");

        TPave *row2 = new TPave(0, (table2_epCDn.size() - 1), summary_table_title.size(), table2_epCDn.size(), 0, "br");  // Row 1 (top row)
        row2->SetFillColor(kGray + 1);
        // row2->SetFillColor(kAzure - 9);
        row2->SetFillStyle(1001);
        row2->SetLineColor(0);
        row2->Draw();

        // Loop over rows and columns to position text
        for (int i = 0; i < table2_epCDn.size(); i++) {
            std::string cellTemp = table2_epCDn.at(i).at(0).Data();

            if ((!findSubstring(cellTemp, "}{(Step") && !findSubstring(cellTemp, "}{cuts")) && (i > 0)) {
                TPave *row = new TPave(0, (table2_epCDn.size() - i - 1), summary_table_title.size(), (table2_epCDn.size() - i), 0, "br");
                row->SetFillColor(kAzure - 9);
                row->SetFillStyle(1001);
                row->SetLineColor(0);
                row->Draw();
            }
        }

        // Create an instance of TLatex
        // TLatex latex_epCDn;

        // Set text alignment and font size
        latex_epCDn.SetTextAlign(22);  // Centered
        // latex_epCDn.SetTextSize(0.01);
        latex_epCDn.SetTextSize(0.02);

        // Loop over rows and columns to position text
        for (int i = 0; i < table2_epCDn.size(); i++) {
            // table2_epCDn.size() rows
            for (int j = 0; j < summary_table_title.size(); j++) {
                // summary_table_title.size() columns
                latex_epCDn.DrawLatex(j + 0.5, table2_epCDn.size() - i - 0.5, table2_epCDn[i][j]);  // Adjust positioning
            }
        }

        // Add gridlines for clarity (optional)
        for (int i = 0; i <= table2_epCDn.size(); i++) {
            // Horizontal lines
            TLine *line_epCDn = new TLine(0, i, summary_table_title.size(), i);
            line_epCDn->SetLineStyle(2);
            line_epCDn->Draw();
        }
        for (int j = 0; j <= summary_table_title.size(); j++) {
            // Vertical lines
            TLine *line_epCDn = new TLine(j, 0, j, table2_epCDn.size());
            line_epCDn->SetLineStyle(2);
            line_epCDn->Draw();
        }

        myTable->Print(fileName, "pdf");
        myTable->Clear();

        TH2F *frame3_epCDn = new TH2F("frame3_epCDn", "", summary_table_title.size(), 0, summary_table_title.size(), table3_epCDn.size(), 0, table3_epCDn.size());
        frame3_epCDn->SetStats(0);                   // Disable statistics box
        frame3_epCDn->GetXaxis()->SetLabelSize(0);   // Remove x-axis labels
        frame3_epCDn->GetXaxis()->SetTickLength(0);  // Remove x-axis ticks
        frame3_epCDn->GetYaxis()->SetLabelSize(0);   // Remove y-axis labels
        frame3_epCDn->GetYaxis()->SetTickLength(0);  // Remove y-axis ticks
        frame3_epCDn->Draw();

        titles.DrawLatexNDC(0.05, 0.9, "Step by step statistics - CD proton (table 3 of 3)");

        TPave *row3 = new TPave(0, (table3_epCDn.size() - 1), summary_table_title.size(), table3_epCDn.size(), 0, "br");  // Row 1 (top row)
        row3->SetFillColor(kGray + 1);
        // row3->SetFillColor(kAzure - 9);
        row3->SetFillStyle(1001);
        row3->SetLineColor(0);
        row3->Draw();

        // Loop over rows and columns to position text
        for (int i = 0; i < table3_epCDn.size(); i++) {
            std::string cellTemp = table3_epCDn.at(i).at(0).Data();

            if ((!findSubstring(cellTemp, "}{(Step") && !findSubstring(cellTemp, "}{cuts")) && (i > 0)) {
                TPave *row = new TPave(0, (table3_epCDn.size() - i - 1), summary_table_title.size(), (table3_epCDn.size() - i), 0, "br");
                row->SetFillColor(kAzure - 9);
                row->SetFillStyle(1001);
                row->SetLineColor(0);
                row->Draw();
            }
        }

        // Create an instance of TLatex
        // TLatex latex_epCDn;

        // Set text alignment and font size
        latex_epCDn.SetTextAlign(22);  // Centered
        latex_epCDn.SetTextSize(0.015);
        // latex_epCDn.SetTextSize(0.02);

        // Loop over rows and columns to position text
        for (int i = 0; i < table3_epCDn.size(); i++) {
            // table3_epCDn.size() rows
            for (int j = 0; j < summary_table_title.size(); j++) {
                // summary_table_title.size() columns
                latex_epCDn.DrawLatex(j + 0.5, table3_epCDn.size() - i - 0.5, table3_epCDn[i][j]);  // Adjust positioning
            }
        }

        // Add gridlines for clarity (optional)
        for (int i = 0; i <= table3_epCDn.size(); i++) {
            // Horizontal lines
            TLine *line_epCDn = new TLine(0, i, summary_table_title.size(), i);
            line_epCDn->SetLineStyle(2);
            line_epCDn->Draw();
        }
        for (int j = 0; j <= summary_table_title.size(); j++) {
            // Vertical lines
            TLine *line_epCDn = new TLine(j, 0, j, table3_epCDn.size());
            line_epCDn->SetLineStyle(2);
            line_epCDn->Draw();
        }

        myTable->Print(fileName, "pdf");
        myTable->Clear();
        titles.SetTextSize(0.065);
    }

    if (Constraint1 == "" || Constraint1 == "FD") {
        myTable->SetTopMargin(0.15);

        double Num_of_goodN_bfSteps_epFDn = GetHistogramEntries(HistoList, "dpp_goodN_epFDn");
        double Num_of_badN_bfSteps_epFDn = GetHistogramEntries(HistoList, "dpp_badN_epFDn");
        string Num_of_allN_bfSteps_epFDn_str = to_string_with_precision(Num_of_goodN_bfSteps_epFDn + Num_of_badN_bfSteps_epFDn, 0);
        string Num_of_goodN_bfSteps_epFDn_str = to_string_with_precision(Num_of_goodN_bfSteps_epFDn, 0);
        string Num_of_badN_bfSteps_epFDn_str = to_string_with_precision(Num_of_badN_bfSteps_epFDn, 0);
        string Single_eff_bfSteps_epFDn_str = "--";
        string Single_purity_bfSteps_epFDn_str = "--";
        const char *Num_of_allN_bfSteps_epFDn_char = Num_of_allN_bfSteps_epFDn_str.c_str();
        const char *Num_of_goodN_bfSteps_epFDn_char = Num_of_goodN_bfSteps_epFDn_str.c_str();
        const char *Num_of_badN_bfSteps_epFDn_char = Num_of_badN_bfSteps_epFDn_str.c_str();
        const char *Single_eff_bfSteps_epFDn_char = Single_eff_bfSteps_epFDn_str.c_str();
        const char *Single_purity_bfSteps_epFDn_char = Single_purity_bfSteps_epFDn_str.c_str();

        double Num_of_goodN_Step0_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step0_epFDn");
        double Num_of_badN_Step0_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step0_epFDn");
        string Num_of_allN_Step0_epFDn_str = to_string_with_precision(Num_of_goodN_Step0_epFDn + Num_of_badN_Step0_epFDn, 0);
        string Num_of_goodN_Step0_epFDn_str = to_string_with_precision(Num_of_goodN_Step0_epFDn, 0);
        string Num_of_badN_Step0_epFDn_str = to_string_with_precision(Num_of_badN_Step0_epFDn, 0);
        string Single_eff_Step0_epFDn_str = to_string_with_precision(Num_of_goodN_Step0_epFDn / Num_of_goodN_bfSteps_epFDn);
        string Single_purity_Step0_epFDn_str = to_string_with_precision(Num_of_goodN_Step0_epFDn / (Num_of_goodN_Step0_epFDn + Num_of_badN_Step0_epFDn));
        const char *Num_of_allN_Step0_epFDn_char = Num_of_allN_Step0_epFDn_str.c_str();
        const char *Num_of_goodN_Step0_epFDn_char = Num_of_goodN_Step0_epFDn_str.c_str();
        const char *Num_of_badN_Step0_epFDn_char = Num_of_badN_Step0_epFDn_str.c_str();
        const char *Single_eff_Step0_epFDn_char = Single_eff_Step0_epFDn_str.c_str();
        const char *Single_purity_Step0_epFDn_char = Single_purity_Step0_epFDn_str.c_str();

        double Num_of_goodN_Step1_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step1_epFDn");
        double Num_of_badN_Step1_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step1_epFDn");
        string Num_of_allN_Step1_epFDn_str = to_string_with_precision(Num_of_goodN_Step1_epFDn + Num_of_badN_Step1_epFDn, 0);
        string Num_of_goodN_Step1_epFDn_str = to_string_with_precision(Num_of_goodN_Step1_epFDn, 0);
        string Num_of_badN_Step1_epFDn_str = to_string_with_precision(Num_of_badN_Step1_epFDn, 0);
        string Single_eff_Step1_epFDn_str = to_string_with_precision(Num_of_goodN_Step1_epFDn / Num_of_goodN_Step0_epFDn);
        string Single_purity_Step1_epFDn_str = to_string_with_precision(Num_of_goodN_Step1_epFDn / (Num_of_goodN_Step1_epFDn + Num_of_badN_Step1_epFDn));
        const char *Num_of_allN_Step1_epFDn_char = Num_of_allN_Step1_epFDn_str.c_str();
        const char *Num_of_goodN_Step1_epFDn_char = Num_of_goodN_Step1_epFDn_str.c_str();
        const char *Num_of_badN_Step1_epFDn_char = Num_of_badN_Step1_epFDn_str.c_str();
        const char *Single_eff_Step1_epFDn_char = Single_eff_Step1_epFDn_str.c_str();
        const char *Single_purity_Step1_epFDn_char = Single_purity_Step1_epFDn_str.c_str();

        double Num_of_goodN_Step2_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step2_epFDn");
        double Num_of_badN_Step2_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step2_epFDn");
        string Num_of_allN_Step2_epFDn_str = to_string_with_precision(Num_of_goodN_Step2_epFDn + Num_of_badN_Step2_epFDn, 0);
        string Num_of_goodN_Step2_epFDn_str = to_string_with_precision(Num_of_goodN_Step2_epFDn, 0);
        string Num_of_badN_Step2_epFDn_str = to_string_with_precision(Num_of_badN_Step2_epFDn, 0);
        string Single_eff_Step2_epFDn_str = to_string_with_precision(Num_of_goodN_Step2_epFDn / Num_of_goodN_Step1_epFDn);
        string Single_purity_Step2_epFDn_str = to_string_with_precision(Num_of_goodN_Step2_epFDn / (Num_of_goodN_Step2_epFDn + Num_of_badN_Step2_epFDn));
        const char *Num_of_allN_Step2_epFDn_char = Num_of_allN_Step2_epFDn_str.c_str();
        const char *Num_of_goodN_Step2_epFDn_char = Num_of_goodN_Step2_epFDn_str.c_str();
        const char *Num_of_badN_Step2_epFDn_char = Num_of_badN_Step2_epFDn_str.c_str();
        const char *Single_eff_Step2_epFDn_char = Single_eff_Step2_epFDn_str.c_str();
        const char *Single_purity_Step2_epFDn_char = Single_purity_Step2_epFDn_str.c_str();

        double Num_of_goodN_Step3_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step3_epFDn");
        double Num_of_badN_Step3_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step3_epFDn");
        string Num_of_allN_Step3_epFDn_str = to_string_with_precision(Num_of_goodN_Step3_epFDn + Num_of_badN_Step3_epFDn, 0);
        string Num_of_goodN_Step3_epFDn_str = to_string_with_precision(Num_of_goodN_Step3_epFDn, 0);
        string Num_of_badN_Step3_epFDn_str = to_string_with_precision(Num_of_badN_Step3_epFDn, 0);
        string Single_eff_Step3_epFDn_str = to_string_with_precision(Num_of_goodN_Step3_epFDn / Num_of_goodN_Step2_epFDn);
        string Single_purity_Step3_epFDn_str = to_string_with_precision(Num_of_goodN_Step3_epFDn / (Num_of_goodN_Step3_epFDn + Num_of_badN_Step3_epFDn));
        const char *Num_of_allN_Step3_epFDn_char = Num_of_allN_Step3_epFDn_str.c_str();
        const char *Num_of_goodN_Step3_epFDn_char = Num_of_goodN_Step3_epFDn_str.c_str();
        const char *Num_of_badN_Step3_epFDn_char = Num_of_badN_Step3_epFDn_str.c_str();
        const char *Single_eff_Step3_epFDn_char = Single_eff_Step3_epFDn_str.c_str();
        const char *Single_purity_Step3_epFDn_char = Single_purity_Step3_epFDn_str.c_str();

        double Num_of_goodN_Step4_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step4_epFDn");
        double Num_of_badN_Step4_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step4_epFDn");
        string Num_of_allN_Step4_epFDn_str = to_string_with_precision(Num_of_goodN_Step4_epFDn + Num_of_badN_Step4_epFDn, 0);
        string Num_of_goodN_Step4_epFDn_str = to_string_with_precision(Num_of_goodN_Step4_epFDn, 0);
        string Num_of_badN_Step4_epFDn_str = to_string_with_precision(Num_of_badN_Step4_epFDn, 0);
        string Single_eff_Step4_epFDn_str = to_string_with_precision(Num_of_goodN_Step4_epFDn / Num_of_goodN_Step3_epFDn);
        string Single_purity_Step4_epFDn_str = to_string_with_precision(Num_of_goodN_Step4_epFDn / (Num_of_goodN_Step4_epFDn + Num_of_badN_Step4_epFDn));
        const char *Num_of_allN_Step4_epFDn_char = Num_of_allN_Step4_epFDn_str.c_str();
        const char *Num_of_goodN_Step4_epFDn_char = Num_of_goodN_Step4_epFDn_str.c_str();
        const char *Num_of_badN_Step4_epFDn_char = Num_of_badN_Step4_epFDn_str.c_str();
        const char *Single_eff_Step4_epFDn_char = Single_eff_Step4_epFDn_str.c_str();
        const char *Single_purity_Step4_epFDn_char = Single_purity_Step4_epFDn_str.c_str();

        double Num_of_goodN_Step5_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step5_epFDn");
        double Num_of_badN_Step5_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step5_epFDn");
        string Num_of_allN_Step5_epFDn_str = to_string_with_precision(Num_of_goodN_Step5_epFDn + Num_of_badN_Step5_epFDn, 0);
        string Num_of_goodN_Step5_epFDn_str = to_string_with_precision(Num_of_goodN_Step5_epFDn, 0);
        string Num_of_badN_Step5_epFDn_str = to_string_with_precision(Num_of_badN_Step5_epFDn, 0);
        string Num_of_goodN_Step5_loss_epFDn_str = to_string_with_precision(Num_of_goodN_Step5_epFDn / Num_of_goodN_Step4_epFDn);
        string Num_of_badN_Step5_loss_epFDn_str = to_string_with_precision(Num_of_goodN_Step5_epFDn / (Num_of_goodN_Step5_epFDn + Num_of_badN_Step5_epFDn));
        const char *Num_of_allN_Step5_epFDn_char = Num_of_allN_Step5_epFDn_str.c_str();
        const char *Num_of_goodN_Step5_epFDn_char = Num_of_goodN_Step5_epFDn_str.c_str();
        const char *Num_of_badN_Step5_epFDn_char = Num_of_badN_Step5_epFDn_str.c_str();
        const char *Num_of_goodN_Step5_loss_epFDn_char = Num_of_goodN_Step5_loss_epFDn_str.c_str();
        const char *Num_of_badN_Step5_loss_epFDn_char = Num_of_badN_Step5_loss_epFDn_str.c_str();

        // Draw a frame without axis numbers and ticks
        TH2F *frame_epFDn = new TH2F("frame_epFDn", "", 5, 0, 5, 8, 0, 8);
        frame_epFDn->SetStats(0);                   // Disable statistics box
        frame_epFDn->GetXaxis()->SetLabelSize(0);   // Remove x-axis labels
        frame_epFDn->GetXaxis()->SetTickLength(0);  // Remove x-axis ticks
        frame_epFDn->GetYaxis()->SetLabelSize(0);   // Remove y-axis labels
        frame_epFDn->GetYaxis()->SetTickLength(0);  // Remove y-axis ticks
        frame_epFDn->Draw();

        titles.DrawLatexNDC(0.05, 0.9, "Step by step statistics - FD proton");

        // Create an instance of TLatex
        TLatex latex_epFDn;

        // Set text alignment and font size
        latex_epFDn.SetTextAlign(22);  // Centered
        latex_epFDn.SetTextSize(0.02);

        // Define table content
        // const char *table_epFDn[8][5] = {
        table_epFDn = {{"", "#(goodN)", "#(badN)", "#splitline{Signal}{Efficiency}", "#splitline{Signal}{Purity}"},
                       {"#splitline{Before}{Step Cuts}", Num_of_goodN_bfSteps_epFDn_char, Num_of_badN_bfSteps_epFDn_char, Single_eff_bfSteps_epFDn_char, Single_purity_bfSteps_epFDn_char},
                       {"Step 0", Num_of_goodN_Step0_epFDn_char, Num_of_badN_Step0_epFDn_char, Single_eff_Step0_epFDn_char, Single_purity_Step0_epFDn_char},
                       {"Step 1", Num_of_goodN_Step1_epFDn_char, Num_of_badN_Step1_epFDn_char, Single_eff_Step1_epFDn_char, Single_purity_Step1_epFDn_char},
                       {"Step 2", Num_of_goodN_Step2_epFDn_char, Num_of_badN_Step2_epFDn_char, Single_eff_Step2_epFDn_char, Single_purity_Step2_epFDn_char},
                       {"Step 3", Num_of_goodN_Step3_epFDn_char, Num_of_badN_Step3_epFDn_char, Single_eff_Step3_epFDn_char, Single_purity_Step3_epFDn_char},
                       {"Step 4", Num_of_goodN_Step4_epFDn_char, Num_of_badN_Step4_epFDn_char, Single_eff_Step4_epFDn_char, Single_purity_Step4_epFDn_char},
                       {"Step 5", Num_of_goodN_Step5_epFDn_char, Num_of_badN_Step5_epFDn_char, Num_of_goodN_Step5_loss_epFDn_char, Num_of_badN_Step5_loss_epFDn_char}};

        // Loop over rows and columns to position text
        for (int i = 0; i < 8; i++) {
            // 8 rows
            for (int j = 0; j < 5; j++) {
                // 5 columns
                latex_epFDn.DrawLatex(j + 0.5, 8 - i - 0.5, table_epFDn[i][j]);  // Adjust positioning
            }
        }

        // Add gridlines for clarity (optional)
        for (int i = 0; i <= 8; i++) {
            // Horizontal lines
            TLine *line_epFDn = new TLine(0, i, 5, i);
            line_epFDn->SetLineStyle(2);
            line_epFDn->Draw();
        }
        for (int j = 0; j <= 5; j++) {
            // Vertical lines
            TLine *line_epFDn = new TLine(j, 0, j, 8);
            line_epFDn->SetLineStyle(2);
            line_epFDn->Draw();
        }

        myTable->Print(fileName, "pdf");
        myTable->Clear();
    }
}

// GetHistogramEntries function ----------------------------------------------------------------------------------------------------------------------------------------------

double HistPrinter::GetHistogramEntries(const std::vector<TH1 *> &HistoList, const std::string &histName) {
    for (const auto &hist : HistoList) {
        if (hist && hist->GetName() == histName) {
            // Check if the histogram exists and the name matches
            return hist->GetEntries();  // Return the number of entries
        }
    }
    return -1;  // Return -1 if no match is found
}

// extractStep function ------------------------------------------------------------------------------------------------------------------------------------------------------

std::string HistPrinter::extractStep(const std::string &input) {
    std::regex stepRegex(R"(Step\d+)");  // Regex to match "Step" followed by digits
    std::smatch match;

    if (std::regex_search(input, match, stepRegex)) {
        return match.str();  // Return the matched substring
    }
    return "";  // Return an empty string if no match is found
}

// SkippingCondition function ------------------------------------------------------------------------------------------------------------------------------------------------

bool HistPrinter::SkippingCondition(string HistoName, int canvas_ind) {
    bool PrintOut = true;

    // TODO: fix this in the all plots file!
    if (HistoName == "Chi2pid_p_APID_epCD" || HistoName == "Chi2pid_p_APID_epFD"                                                                       // Last PID plot
        || HistoName == "nSector_VS_ToF_n_epCDn" || HistoName == "nSector_VS_ToF_n_epFDn"                                                              // Last miss cuts plot
        || HistoName == "n_multiplicity_badN_epCDn_Step0" || HistoName == "n_multiplicity_badN_epFDn_Step0"                                            // Last Step0 multiplicity plot
        || HistoName == "beta_n_badN_Step0_epCDn" || HistoName == "beta_n_badN_Step0_epFDn"                                                            // Last Step0 plot
        || HistoName == "n_multiplicity_badN_epCDn_Step1" || HistoName == "n_multiplicity_badN_epFDn_Step1"                                            // Last Step1 multiplicity plot
        || HistoName == "diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epCDn" || HistoName == "diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epFDn"                  // Last Step1 plot
        || HistoName == "n_multiplicity_badN_epCDn_Step2" || HistoName == "n_multiplicity_badN_epFDn_Step2"                                            // Last Step2 multiplicity plot
        || HistoName == "Edep_CND_badN_withNearbyPos_Step2prep_epCDn" || HistoName == "Edep_CND_badN_withNearbyPos_Step2prep_epFDn"                    // Last Step2 multiplicity plot
        || HistoName == "sdiff_pos_VS_dToF_rel_n_badN_Step2prep_layer_3_epCDn" || HistoName == "sdiff_pos_VS_dToF_rel_n_badN_Step2prep_layer_3_epFDn"  // Last Step2
        || (HistoName == "Edep_CND_badN_withNearbyPos_Step2_epCDn" || HistoName == "Edep_CND_badN_withNearbyPos_Step2_epFDn")                          // Last before ldiff = -3 plots
    ) {
        if (PrintOut) { cout << "\n\nHistoName = '" << HistoName << "'; Skipped!\n\n"; }

        return true;
    }

    return false;
}

// replaceSubstring function -------------------------------------------------------------------------------------------------------------------------------------------------

// Function to replace one substring with another
std::string HistPrinter::replaceSubstring(const std::string &input, const std::string &toReplace, const std::string &replaceWith) {
    size_t pos = input.find(toReplace);

    if (pos == std::string::npos) {
        // If 'toReplace' is not found, return the original string
        return input;
    }
    return input.substr(0, pos) + replaceWith + input.substr(pos + toReplace.length());
}

// SectionPlotter function ---------------------------------------------------------------------------------------------------------------------------------------------------

void HistPrinter::SectionPlotter(int n_col, int n_row, TCanvas *myCanvas, TCanvas *myText, TCanvas *myTable, vector<TH1 *> HistoList, string PDFFile, string Constraint1, string Constraint2,
                                 bool LogScale2D) {
    Num_of_goodN_bfSteps_epCDn = GetHistogramEntries(HistoList, "dpp_goodN_epCDn");
    Num_of_badN_bfSteps_epCDn = GetHistogramEntries(HistoList, "dpp_badN_epCDn");

    TLatex titles, text;
    titles.SetTextSize(0.065);
    text.SetTextSize(0.04);

    string pdfFile0;

    if (!LogScale2D) {
        if (Constraint1 == "" && Constraint2 == "") {
            pdfFile0 = PDFFile;
        } else if (Constraint1 != "" && Constraint2 == "") {
            string pdfFile1 = ConfigOutPutName(PDFFile, Constraint1);
            pdfFile0 = pdfFile1;
        } else if (Constraint1 == "" && Constraint2 != "") {
            string pdfFile1 = ConfigOutPutName(PDFFile, Constraint2);
            pdfFile0 = pdfFile1;
        } else if (Constraint1 != "" && Constraint2 != "") {
            string pdfFile2 = ConfigOutPutName(PDFFile, Constraint1);
            string pdfFile1 = ConfigOutPutName(pdfFile2, Constraint2);
            pdfFile0 = pdfFile1;
        }
    } else {
        if (Constraint1 == "" && Constraint2 == "") {
            string pdfFile1 = ConfigOutPutName(PDFFile, "LogScale2D");
            pdfFile0 = pdfFile1;
        } else if (Constraint1 != "" && Constraint2 == "") {
            string pdfFile2 = ConfigOutPutName(PDFFile, Constraint1);
            string pdfFile1 = ConfigOutPutName(pdfFile2, "LogScale2D");
            pdfFile0 = pdfFile1;
        } else if (Constraint1 == "" && Constraint2 != "") {
            string pdfFile2 = ConfigOutPutName(PDFFile, Constraint2);
            string pdfFile1 = ConfigOutPutName(pdfFile2, "LogScale2D");
            pdfFile0 = pdfFile1;
        } else if (Constraint1 != "" && Constraint2 != "") {
            string pdfFile3 = ConfigOutPutName(PDFFile, Constraint1);
            string pdfFile2 = ConfigOutPutName(pdfFile3, Constraint2);
            string pdfFile1 = ConfigOutPutName(pdfFile2, "LogScale2D");
            pdfFile0 = pdfFile1;
        }
    }

    const char *pdfFile = pdfFile0.c_str();

    // My plots root file
    TList *plots = new TList();
    string listName = replaceSubstring(pdfFile0, ".pdf", ".root");
    const char *TListName = listName.c_str();

    char fileName[100];
    sprintf(fileName, "%s[", pdfFile);
    myText->SaveAs(fileName);
    sprintf(fileName, "%s", pdfFile);

    myText->cd();

    PrintPage(HistoList, "Manual Veto Plots", myText, fileName, titles, text, Constraint1, Constraint2);

    myCanvas->cd();
    myCanvas->Divide(n_col, n_row);

    double x_1 = 0.2, y_1 = 0.3, x_2 = 0.86, y_2 = 0.7;
    double diplayTextSize = 0.1;

    int canvas_ind = 1;

    bool FilledConstraint1Bookmark = false;

    bool FirstPIDPlot = true;

    bool FirstOnlyMissCutsPlot = true;

    map<string, bool> FirstStepPlot;
    FirstStepPlot["Step0"] = true, FirstStepPlot["Step1"] = true, FirstStepPlot["Step2"] = true, FirstStepPlot["Step3"] = true, FirstStepPlot["Step4"] = true, FirstStepPlot["Step5"] = true;

    for (int i = 0; i < HistoList.size(); i++) {
        string TempHistName = HistoList[i]->GetName();

        bool GoodHistogram;

        if (Constraint1 == "" && Constraint2 == "") {
            GoodHistogram = true;
        } else if (Constraint1 != "" && Constraint2 == "") {
            GoodHistogram = findSubstring(TempHistName, Constraint1);
        } else if (Constraint1 == "" && Constraint2 != "") {
            GoodHistogram = findSubstring(TempHistName, Constraint2);
        } else {
            GoodHistogram = (findSubstring(TempHistName, Constraint1) && findSubstring(TempHistName, Constraint2));
        }

        if (GoodHistogram) {
            if (findSubstring(TempHistName, "BPID")) {
                if (FirstPIDPlot) {
                    myText->cd();

                    PrintPage(HistoList, "PID Plots", myText, fileName, titles, text, Constraint1, Constraint2);

                    FirstPIDPlot = false;
                }
            } else if (findSubstring(TempHistName, "BmissC")) {
                if (FirstOnlyMissCutsPlot) {
                    myText->cd();

                    PrintPage(HistoList, "Plots with basic cuts", myText, fileName, titles, text, Constraint1, Constraint2);
                    PrintPage(HistoList, "Definition of neutrons in veto steps", myText, fileName, titles, text, Constraint1, Constraint2);

                    myTable->cd();
                    SummaryTablePlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, titles, text, fileName, PDFFile, Constraint1, Constraint2, LogScale2D);
                    myText->cd();

                    FirstOnlyMissCutsPlot = false;
                }
            } else if (findSubstring(TempHistName, "Step")) {
                string Step = extractStep(TempHistName);

                if (FirstStepPlot[Step] == true) {
                    myText->cd();

                    titles.DrawLatex(0.05, 0.9, (Step + " Plots").c_str());

                    if (Step == "Step0") {
                        PrintPage(HistoList, "Step0 Plots", myText, fileName, titles, text, Constraint1, Constraint2);
                        PrintPage(HistoList, "Definition of neutrons in veto steps", myText, fileName, titles, text, Constraint1, Constraint2);
                    } else if (Step == "Step1") {
                        PrintPage(HistoList, "Step1 Plots", myText, fileName, titles, text, Constraint1, Constraint2);
                        PrintPage(HistoList, "Definition of neutrons in veto steps", myText, fileName, titles, text, Constraint1, Constraint2);
                    } else if (Step == "Step2") {
                        PrintPage(HistoList, "Step2 Plots", myText, fileName, titles, text, Constraint1, Constraint2);
                        PrintPage(HistoList, "Definition of neutrons in veto steps", myText, fileName, titles, text, Constraint1, Constraint2);
                    }

                    //                    myText->Print(fileName, "pdf");
                    //                    myText->Clear();

                    myTable->cd();
                    SummaryTablePlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, titles, text, fileName, PDFFile, Constraint1, Constraint2, LogScale2D);
                    myText->cd();

                    FirstStepPlot[Step] = false;
                }
            }

            myCanvas->cd(canvas_ind);
            myCanvas->cd(canvas_ind)->SetBottomMargin(0.14), myCanvas->cd(canvas_ind)->SetLeftMargin(0.16), myCanvas->cd(canvas_ind)->SetRightMargin(0.16),
                myCanvas->cd(canvas_ind)->SetTopMargin(0.12);

            HistoList[i]->GetYaxis()->SetTitleOffset(1.5);
            HistoList[i]->GetXaxis()->SetTitleOffset(1.1);

            // gPad->SetRightMargin(0.23);

            gPad->SetGrid();
            gPad->SetFrameLineWidth(1);  // Reset frame line width to 1

            gStyle->SetOptStat("ourmen");

            if (HistoList[i]->InheritsFrom("TH1D")) {
                HistoList[i]->SetMinimum(0);
                HistoList[i]->SetLineWidth(1);
                HistoList[i]->SetLineColor(kRed);
            }

            if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0) {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize * 0.6), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);

                if (HistoList[i]->InheritsFrom("TH1D")) {
                    HistoList[i]->Draw(), displayText->Draw("same");
                    plots->Add(HistoList[i]);
                } else if (HistoList[i]->InheritsFrom("TH2D")) {
                    if (LogScale2D) {
                        TH2D *HistoList_i_LogScale = dynamic_cast<TH2D *>(HistoList[i]);

                        if (HistoList_i_LogScale) {
                            // HistoList_i_LogScale->SetLogz(1);
                            gPad->SetLogz(1);

                            if (findSubstring(TempHistName, "Size_CND1_VS_Size_CND2") || findSubstring(TempHistName, "Size_CND1_VS_Size_CND3") ||
                                findSubstring(TempHistName, "Size_CND2_VS_Size_CND3") || findSubstring(TempHistName, "LayerMult_CND1_VS_LayerMult_CND2") ||
                                findSubstring(TempHistName, "LayerMult_CND1_VS_LayerMult_CND3") || findSubstring(TempHistName, "LayerMult_CND2_VS_LayerMult_CND3")) {
                                HistoList_i_LogScale->Draw("text colz"), displayText->Draw("same");
                                HistoList[i]->SetMarkerSize(3.0);  // Increase marker size, which scales the text
                                HistoList[i]->SetMarkerColor(kMagenta);
                                plots->Add(HistoList_i_LogScale);
                            } else {
                                HistoList_i_LogScale->Draw("colz"), displayText->Draw("same");
                                plots->Add(HistoList_i_LogScale);
                            }
                        }
                    } else {
                        if (findSubstring(TempHistName, "Size_CND1_VS_Size_CND2") || findSubstring(TempHistName, "Size_CND1_VS_Size_CND3") || findSubstring(TempHistName, "Size_CND2_VS_Size_CND3") ||
                            findSubstring(TempHistName, "LayerMult_CND1_VS_LayerMult_CND2") || findSubstring(TempHistName, "LayerMult_CND1_VS_LayerMult_CND3") ||
                            findSubstring(TempHistName, "LayerMult_CND2_VS_LayerMult_CND3")) {
                            HistoList[i]->Draw("text colz"), displayText->Draw("same");
                            HistoList[i]->SetMarkerSize(3.0);  // Increase marker size, which scales the text
                            HistoList[i]->SetMarkerColor(kMagenta);
                            plots->Add(HistoList[i]);
                        } else {
                            HistoList[i]->Draw("colz"), displayText->Draw("same");
                            plots->Add(HistoList[i]);
                        }
                    }
                } else {
                    cout << "ERROR! could not determine histogram class! Exiting...";
                    exit(0);
                }
            } else {
                if (HistoList[i]->InheritsFrom("TH1D")) {
                    HistoList[i]->Draw();
                    plots->Add(HistoList[i]);
                } else if (HistoList[i]->InheritsFrom("TH2D")) {
                    if (findSubstring(TempHistName, "Size_CND1_VS_Size_CND2") || findSubstring(TempHistName, "Size_CND1_VS_Size_CND3") || findSubstring(TempHistName, "Size_CND2_VS_Size_CND3") ||
                        findSubstring(TempHistName, "LayerMult_CND1_VS_LayerMult_CND2") || findSubstring(TempHistName, "LayerMult_CND1_VS_LayerMult_CND3") ||
                        findSubstring(TempHistName, "LayerMult_CND2_VS_LayerMult_CND3")) {
                        HistoList[i]->Draw("text colz");
                        HistoList[i]->SetMarkerSize(3.0);  // Increase marker size, which scales the text
                        HistoList[i]->SetMarkerColor(kMagenta);

                        gPad->Update();
                        TPaletteAxis *palette = (TPaletteAxis *)HistoList[i]->GetListOfFunctions()->FindObject("palette");
                        palette->SetY2NDC(0.55);
                        gPad->Modified();
                        gPad->Update();

                        plots->Add(HistoList[i]);
                    } else {
                        HistoList[i]->Draw("colz");

                        gPad->Update();
                        TPaletteAxis *palette = (TPaletteAxis *)HistoList[i]->GetListOfFunctions()->FindObject("palette");
                        palette->SetY2NDC(0.55);
                        gPad->Modified();
                        gPad->Update();

                        plots->Add(HistoList[i]);
                    }
                } else {
                    cout << "ERROR! could not determine histogram class! Exiting...";
                    exit(0);
                }
            }

            // Save the canvas to a PDF page after filling 12 pads or processing the last histogram
            if (canvas_ind == n_col * n_row || SkippingCondition(TempHistName, canvas_ind)) {
                myCanvas->Print(fileName);       // Save the current page
                myCanvas->Clear();               // Clear the canvas for the next page
                myCanvas->Divide(n_col, n_row);  // Reset the grid layout

                canvas_ind = 0;
            }

            ++canvas_ind;
        }
    }

    sprintf(fileName, "%s]", pdfFile);
    myCanvas->Print(fileName, "pdf");

    myCanvas->Clear();
    myText->Clear();

    // Saving histogram TList
    TFile *plots_fout = new TFile(TListName, "recreate");
    plots_fout->cd();
    plots->Write();
    plots_fout->Write();
    plots_fout->Close();

    delete plots_fout;
}

// PlotHistograms function ---------------------------------------------------------------------------------------------------------------------------------------------------

void HistPrinter::PlotHistograms(const vector<TH1 *> HistoList, const string &PDFFile, bool LogScale2D) {
    /////////////////////////////////////////////////////
    // Now create the output PDFs
    /////////////////////////////////////////////////////

    int n_col = 2, n_row = 2;

    // int pixelx = 1980, pixely = 1530;
    // int pixelx = 1980 * n_col, pixely = 1530 * 4;
    // int pixelx = 1980 * n_col * 1.5 * 2, pixely = 1530 * 4 * 1.5 * 2;
    // int pixelx = 1000 * n_col * 1.5 * 2, pixely = 750 * 3 * 1.5 * 2;
    // int pixelx = 1000 * n_col * 5, pixely = 750 * 3 * 4;
    int pixelx = 1000 * n_col * 5, pixely = 750 * n_row * 5;

    TCanvas *myCanvas = new TCanvas("myPage", "myPage", pixelx, pixely);
    TCanvas *myText = new TCanvas("myText", "myText", pixelx, pixely);
    TCanvas *myTable = new TCanvas("myTable", "myTable", pixelx, pixely);

    /* Saving all plots */
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile);

    /* Saving only CD proton plots */
    SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD");
    if (Apply_Step0_Cuts) { SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step0"); }
    SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step1");
    SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step2");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step3");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step4");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step5");

    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "", true);
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step0", true);
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step1", true);
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step2", true);
    // // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step3", true);
    // // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step4", true);
    // // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "CD", "Step5", true);

    // /* Saving only FD proton plots */
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step0");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step1");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step2");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step3");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step4");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step5");
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "", true);
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step0", true);
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step1", true);
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step2", true);
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step3", true);
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step4", true);
    // SectionPlotter(n_col, n_row, myCanvas, myText, myTable, HistoList, PDFFile, "FD", "Step5", true);

    delete myCanvas;
    delete myText;
    delete myTable;
}
