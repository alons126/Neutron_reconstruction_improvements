// SummaryTablePlotter function ----------------------------------------------------------------------------------------------------------------------------------------------

void HistPrinter::GenerateSummaryTable(TCanvas *myTable, vector<TH1 *> HistoList, string Constraint1, string Constraint2) {
    if (First_table_epCDn_generation && (Constraint1 == "" || Constraint1 == "CD")) {
        myTable->SetTopMargin(0.15);

        /* Before Step Cuts */
        // Before Step Cuts overall:
        vector<TString> summary_table_bfSteps_epCDn_1stLine = {"#splitline{Before}{step cuts}", basic_tools::ToStringWithPrecision(Num_of_goodN_bfSteps_epCDn, 0),
                                                               basic_tools::ToStringWithPrecision(Num_of_badN_bfSteps_epCDn, 0), "--", "--"};
        summary_table_bfSteps_epCDn.push_back(summary_table_bfSteps_epCDn_1stLine);

        for (int i = 0; i < summary_table_bfSteps_epCDn.size(); i++) { table_epCDn.push_back(summary_table_bfSteps_epCDn.at(i)); }

        /* Step0 */
        if (Apply_Step0_Cuts) {
            // Step0 overall:
            double Num_of_goodN_Step0_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step0_epCDn");
            double Num_of_badN_Step0_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step0_epCDn");

            vector<TString> summary_table_Step0_epCDn_1stLine = {"Step0", basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step0_epCDn, 0),
                                                                 basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                                 basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_epCDn / (Num_of_goodN_Step0_epCDn + Num_of_badN_Step0_epCDn), PurAccuracy)};
            summary_table_Step0_epCDn.push_back(summary_table_Step0_epCDn_1stLine);

            // dBeta_n test (Step0):
            double Num_of_goodN_Step0_dBeta_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_dBeta_n_Step0_goodN_BCTest_epCDn");
            // double Num_of_badN_Step0_dBeta_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_dBeta_n_Step0_badN_BCTest_epCDn");
            double Num_of_goodN_Step0_dBeta_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_dBeta_n_Step0_goodN_ACTest_epCDn");
            double Num_of_badN_Step0_dBeta_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_dBeta_n_Step0_badN_ACTest_epCDn");

            summary_table_Step0_epCDn.push_back(
                {"#splitline{#Delta#beta_{n} cuts}{(Step0)}", basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_dBeta_n_ACTest_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step0_dBeta_n_ACTest_epCDn, 0),
                 //  basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_dBeta_n_ACTest_epCDn / Num_of_goodN_Step0_dBeta_n_BCTest_epCDn, EffAccuracy),
                 basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_dBeta_n_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                 basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_dBeta_n_ACTest_epCDn / (Num_of_goodN_Step0_dBeta_n_ACTest_epCDn + Num_of_badN_Step0_dBeta_n_ACTest_epCDn), PurAccuracy)});

            // Vz_n test (Step0):
            double Num_of_goodN_Step0_Vz_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Vz_n_Step0_goodN_BCTest_epCDn");
            // double Num_of_badN_Step0_Vz_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Vz_n_Step0_badN_BCTest_epCDn");
            double Num_of_goodN_Step0_Vz_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Vz_n_Step0_goodN_ACTest_epCDn");
            double Num_of_badN_Step0_Vz_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Vz_n_Step0_badN_ACTest_epCDn");

            summary_table_Step0_epCDn.push_back(
                {"#splitline{V_{hit,z} cuts}{(Step0)}", basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_Vz_n_ACTest_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step0_Vz_n_ACTest_epCDn, 0),
                 //  basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_Vz_n_ACTest_epCDn / Num_of_goodN_Step0_Vz_n_BCTest_epCDn, EffAccuracy),
                 basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_Vz_n_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                 basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_Vz_n_ACTest_epCDn / (Num_of_goodN_Step0_Vz_n_ACTest_epCDn + Num_of_badN_Step0_Vz_n_ACTest_epCDn), PurAccuracy)});

            // ToF_n test (Step0):
            double Num_of_goodN_Step0_ToF_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_ToF_n_Step0_goodN_BCTest_epCDn");
            // double Num_of_badN_Step0_ToF_n_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_ToF_n_Step0_badN_BCTest_epCDn");
            double Num_of_goodN_Step0_ToF_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_ToF_n_Step0_goodN_ACTest_epCDn");
            double Num_of_badN_Step0_ToF_n_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_ToF_n_Step0_badN_ACTest_epCDn");

            summary_table_Step0_epCDn.push_back(
                {"#splitline{t_{ToF,n} cuts}{(Step0)}", basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_ToF_n_ACTest_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step0_ToF_n_ACTest_epCDn, 0),
                 //  basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_ToF_n_ACTest_epCDn / Num_of_goodN_Step0_ToF_n_BCTest_epCDn, EffAccuracy),
                 basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_ToF_n_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                 basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_ToF_n_ACTest_epCDn / (Num_of_goodN_Step0_ToF_n_ACTest_epCDn + Num_of_badN_Step0_ToF_n_ACTest_epCDn), PurAccuracy)});

            for (int i = 0; i < summary_table_Step0_epCDn.size(); i++) { table_epCDn.push_back(summary_table_Step0_epCDn.at(i)); }
        }

        /* Step1 */
        // Step1 overall:
        double Num_of_goodN_Step1_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step1_epCDn");
        double Num_of_badN_Step1_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step1_epCDn");

        vector<TString> summary_table_Step1_epCDn_1stLine = {"Step1", basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step1_epCDn, 0),
                                                             // basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_epCDn / Num_of_goodN_Step0_epCDn, EffAccuracy),
                                                             basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                             basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_epCDn / (Num_of_goodN_Step1_epCDn + Num_of_badN_Step1_epCDn), PurAccuracy)};
        summary_table_Step1_epCDn.push_back(summary_table_Step1_epCDn_1stLine);

        // Edep_CND test (Step1):
        double Num_of_goodN_Step1_Edep_CND_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Edep_CND_Step1_goodN_BCTest_epCDn");
        // double Num_of_badN_Step1_Edep_CND_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Edep_CND_Step1_badN_BCTest_epCDn");
        double Num_of_goodN_Step1_Edep_CND_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Edep_CND_Step1_goodN_ACTest_epCDn");
        double Num_of_badN_Step1_Edep_CND_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Edep_CND_Step1_badN_ACTest_epCDn");

        summary_table_Step1_epCDn.push_back(
            {"#splitline{E^{CND}_{dep} cuts}{(Step1)}", basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_Edep_CND_ACTest_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step1_Edep_CND_ACTest_epCDn, 0),
             //  basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_Edep_CND_ACTest_epCDn / Num_of_goodN_Step1_Edep_CND_BCTest_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_Edep_CND_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_Edep_CND_ACTest_epCDn / (Num_of_goodN_Step1_Edep_CND_ACTest_epCDn + Num_of_badN_Step1_Edep_CND_ACTest_epCDn), PurAccuracy)});

        for (int i = 0; i < summary_table_Step1_epCDn.size(); i++) { table_epCDn.push_back(summary_table_Step1_epCDn.at(i)); }

        /* Step2 */
        // Step2 overall:
        double Num_of_goodN_Step2_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step2_epCDn");
        double Num_of_badN_Step2_epCDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step2_epCDn");

        vector<TString> summary_table_Step2_epCDn_1stLine = {"Step2", basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step2_epCDn, 0),
                                                             // basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_epCDn / Num_of_goodN_Step1_epCDn, EffAccuracy),
                                                             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_epCDn / (Num_of_goodN_Step2_epCDn + Num_of_badN_Step2_epCDn), PurAccuracy)};
        // summary_table_Step2_epCDn.push_back(summary_table_Step2_epCDn_1stLine);
        summary_table2_Step2_epCDn.push_back(summary_table_Step2_epCDn_1stLine);
        summary_table3_Step2_epCDn.push_back(summary_table_Step2_epCDn_1stLine);

        // Size_CND1 test (Step2):
        double Num_of_goodN_Step2_Size_CND1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND1_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_Size_CND1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND1_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_Size_CND1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND1_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_Size_CND1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND1_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{Size(CND1) cuts}{(Step2)}", basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND1_ACTest_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step2_Size_CND1_ACTest_epCDn, 0),
             //  basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND1_ACTest_epCDn / Num_of_goodN_Step2_Size_CND1_BCTest_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND1_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND1_ACTest_epCDn / (Num_of_goodN_Step2_Size_CND1_ACTest_epCDn + Num_of_badN_Step2_Size_CND1_ACTest_epCDn), PurAccuracy)});

        // Size_CND2 test (Step2):
        double Num_of_goodN_Step2_Size_CND2_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND2_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_Size_CND2_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND2_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_Size_CND2_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND2_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_Size_CND2_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND2_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{Size(CND2) cuts}{(Step2)}", basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND2_ACTest_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step2_Size_CND2_ACTest_epCDn, 0),
             //  basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND2_ACTest_epCDn / Num_of_goodN_Step2_Size_CND2_BCTest_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND2_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND2_ACTest_epCDn / (Num_of_goodN_Step2_Size_CND2_ACTest_epCDn + Num_of_badN_Step2_Size_CND2_ACTest_epCDn), PurAccuracy)});

        // Size_CND3 test (Step2):
        double Num_of_goodN_Step2_Size_CND3_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND3_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_Size_CND3_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND3_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_Size_CND3_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND3_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_Size_CND3_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_Size_CND3_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{Size(CND3) cuts}{(Step2)}", basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND3_ACTest_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_Step2_Size_CND3_ACTest_epCDn, 0),
             //  basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND3_ACTest_epCDn / Num_of_goodN_Step2_Size_CND3_BCTest_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND3_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_Size_CND3_ACTest_epCDn / (Num_of_goodN_Step2_Size_CND3_ACTest_epCDn + Num_of_badN_Step2_Size_CND3_ACTest_epCDn), PurAccuracy)});

        // LayerMult_CND1 test (Step2):
        double Num_of_goodN_Step2_LayerMult_CND1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND1_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_LayerMult_CND1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND1_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND1_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_LayerMult_CND1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND1_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{LayerMult(CND1)}{cuts (Step2)}", basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn, 0),
             basic_tools::ToStringWithPrecision(Num_of_badN_Step2_LayerMult_CND1_ACTest_epCDn, 0),
             //  basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn / Num_of_goodN_Step2_LayerMult_CND1_BCTest_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn / (Num_of_goodN_Step2_LayerMult_CND1_ACTest_epCDn + Num_of_badN_Step2_LayerMult_CND1_ACTest_epCDn), PurAccuracy)});

        // LayerMult_CND2andCND3 test (Step2):
        double Num_of_goodN_Step2_LayerMult_CND2andCND3_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND2andCND3_Step2_goodN_BCTest_epCDn");
        // double Num_of_badN_Step2_LayerMult_CND2andCND3_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND2andCND3_Step2_badN_BCTest_epCDn");
        double Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND2andCND3_Step2_goodN_ACTest_epCDn");
        double Num_of_badN_Step2_LayerMult_CND2andCND3_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_LayerMult_CND2andCND3_Step2_badN_ACTest_epCDn");

        summary_table2_Step2_epCDn.push_back(
            {"#splitline{LayerMult(CND1+2)}{cuts (Step2)}", basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn, 0),
             basic_tools::ToStringWithPrecision(Num_of_badN_Step2_LayerMult_CND2andCND3_ACTest_epCDn, 0),
             //  basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn / Num_of_goodN_Step2_LayerMult_CND2andCND3_BCTest_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
             basic_tools::ToStringWithPrecision(
                 Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn / (Num_of_goodN_Step2_LayerMult_CND2andCND3_ACTest_epCDn + Num_of_badN_Step2_LayerMult_CND2andCND3_ACTest_epCDn), PurAccuracy)});

        for (int i = 0; i < summary_table2_Step2_epCDn.size(); i++) { table2_epCDn.push_back(summary_table2_Step2_epCDn.at(i)); }

        for (int k = 0; k < 7; k++) {
            double Num_of_goodN_of1_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of1_pos_Step2_layer_" + basic_tools::ToStringWithPrecision(k - 3, 0) + "_goodN_BCTest_epCDn");
            double Num_of_goodN_of1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of1_pos_Step2_layer_" + basic_tools::ToStringWithPrecision(k - 3, 0) + "_goodN_ACTest_epCDn");
            double Num_of_badN_of1_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of1_pos_Step2_layer_" + basic_tools::ToStringWithPrecision(k - 3, 0) + "_badN_ACTest_epCDn");

            summary_table3_Step2_epCDn.push_back({"#splitline{|#DeltaS_{n,+}|>1 & #DeltaL_{n,+} = " + basic_tools::ToStringWithPrecision(k - 3, 0) + " cuts}{(Step2)}",
                                                  basic_tools::ToStringWithPrecision(Num_of_goodN_of1_ACTest_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_of1_ACTest_epCDn, 0),
                                                  basic_tools::ToStringWithPrecision(Num_of_goodN_of1_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                  basic_tools::ToStringWithPrecision(Num_of_goodN_of1_ACTest_epCDn / (Num_of_goodN_of1_ACTest_epCDn + Num_of_badN_of1_ACTest_epCDn), PurAccuracy)});

            double Num_of_goodN_of2_BCTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of2_pos_Step2_layer_" + basic_tools::ToStringWithPrecision(k - 3, 0) + "_goodN_BCTest_epCDn");
            double Num_of_goodN_of2_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of2_pos_Step2_layer_" + basic_tools::ToStringWithPrecision(k - 3, 0) + "_goodN_ACTest_epCDn");
            double Num_of_badN_of2_ACTest_epCDn = GetHistogramEntries(HistoList, "Test_sdiff_of2_pos_Step2_layer_" + basic_tools::ToStringWithPrecision(k - 3, 0) + "_badN_ACTest_epCDn");

            summary_table3_Step2_epCDn.push_back({"#splitline{|#DeltaS_{n,+}|>2 & #DeltaL_{n,+} = " + basic_tools::ToStringWithPrecision(k - 3, 0) + " cuts}{(Step2)}",
                                                  basic_tools::ToStringWithPrecision(Num_of_goodN_of2_ACTest_epCDn, 0), basic_tools::ToStringWithPrecision(Num_of_badN_of2_ACTest_epCDn, 0),
                                                  basic_tools::ToStringWithPrecision(Num_of_goodN_of2_ACTest_epCDn / Num_of_goodN_bfSteps_epCDn, EffAccuracy),
                                                  basic_tools::ToStringWithPrecision(Num_of_goodN_of2_ACTest_epCDn / (Num_of_goodN_of2_ACTest_epCDn + Num_of_badN_of2_ACTest_epCDn), PurAccuracy)});
        }

        for (int i = 0; i < summary_table3_Step2_epCDn.size(); i++) { table3_epCDn.push_back(summary_table3_Step2_epCDn.at(i)); }

        // Prevent the generation of multiple lines:
        First_table_epCDn_generation = false;
    }
}

