// SummaryTablePlotter function ----------------------------------------------------------------------------------------------------------------------------------------------

void HistPrinter::SummaryTablePlotter(int n_col, int n_row, TCanvas *myCanvas, TCanvas *myText, TCanvas *myTable, vector<TH1 *> HistoList, TLatex titles, TLatex text, char fileName[100],
                                      string PDFFile, string Constraint1, string Constraint2, bool LogScale2D) {
    if (Constraint1 == "" || Constraint1 == "CD") {
        myTable.SetTopMargin(0.15);

        titles.SetTextSize(0.05);
        // titles.SetTextSize(0.065);

        GenerateSummaryTable(myTable, HistoList, Constraint1, Constraint2);

        // Draw a frame without axis numbers and ticks
        TH2F frame1_epCDn = TH2F("frame1_epCDn", "", summary_table_title.size(), 0, summary_table_title.size(), table_epCDn.size(), 0, table_epCDn.size());
        frame1_epCDn.SetDirectory(0);
        frame1_epCDn.SetStats(0);                   // Disable statistics box
        frame1_epCDn.GetXaxis().SetLabelSize(0);   // Remove x-axis labels
        frame1_epCDn.GetXaxis().SetTickLength(0);  // Remove x-axis ticks
        frame1_epCDn.GetYaxis().SetLabelSize(0);   // Remove y-axis labels
        frame1_epCDn.GetYaxis().SetTickLength(0);  // Remove y-axis ticks
        frame1_epCDn.Draw();

        titles.DrawLatexNDC(0.05, 0.9, "Step by step statistics - CD proton (table 1 of 3)");

        TPave row1 = TPave(0, (table_epCDn.size() - 1), summary_table_title.size(), table_epCDn.size(), 0, "br");  // Row 1 (top row)
        row1.SetFillColor(kGray + 1);
        // row1.SetFillColor(kAzure - 9);
        row1.SetFillStyle(1001);
        row1.SetLineColor(0);
        row1.Draw();

        // Loop over rows and columns to position text
        for (int i = 0; i < table_epCDn.size(); i++) {
            std::string cellTemp = table_epCDn.at(i).at(0).Data();

            if ((!basic_tools::FindSubstring(cellTemp, "}{(Step") && !basic_tools::FindSubstring(cellTemp, "}{cuts")) && (i > 0)) {
                TPave row = TPave(0, (table_epCDn.size() - i - 1), summary_table_title.size(), (table_epCDn.size() - i), 0, "br");
                row.SetFillColor(kAzure - 9);
                row.SetFillStyle(1001);
                row.SetLineColor(0);
                row.Draw();
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
            TLine line_epCDn = TLine(0, i, summary_table_title.size(), i);
            line_epCDn.SetLineStyle(2);
            line_epCDn.Draw();
        }
        for (int j = 0; j <= summary_table_title.size(); j++) {
            // Vertical lines
            TLine line_epCDn = TLine(j, 0, j, table_epCDn.size());
            line_epCDn.SetLineStyle(2);
            line_epCDn.Draw();
        }

        myTable.Print(fileName, "pdf");
        myTable.Clear();

        TH2F frame2_epCDn = TH2F("frame2_epCDn", "", summary_table_title.size(), 0, summary_table_title.size(), table2_epCDn.size(), 0, table2_epCDn.size());
        frame2_epCDn.SetDirectory(0);
        frame2_epCDn.SetStats(0);                   // Disable statistics box
        frame2_epCDn.GetXaxis().SetLabelSize(0);   // Remove x-axis labels
        frame2_epCDn.GetXaxis().SetTickLength(0);  // Remove x-axis ticks
        frame2_epCDn.GetYaxis().SetLabelSize(0);   // Remove y-axis labels
        frame2_epCDn.GetYaxis().SetTickLength(0);  // Remove y-axis ticks
        frame2_epCDn.Draw();

        titles.DrawLatexNDC(0.05, 0.9, "Step by step statistics - CD proton (table 2 of 3)");

        TPave row2 = TPave(0, (table2_epCDn.size() - 1), summary_table_title.size(), table2_epCDn.size(), 0, "br");  // Row 1 (top row)
        row2.SetFillColor(kGray + 1);
        // row2.SetFillColor(kAzure - 9);
        row2.SetFillStyle(1001);
        row2.SetLineColor(0);
        row2.Draw();

        // Loop over rows and columns to position text
        for (int i = 0; i < table2_epCDn.size(); i++) {
            std::string cellTemp = table2_epCDn.at(i).at(0).Data();

            if ((!basic_tools::FindSubstring(cellTemp, "}{(Step") && !basic_tools::FindSubstring(cellTemp, "}{cuts")) && (i > 0)) {
                TPave row = TPave(0, (table2_epCDn.size() - i - 1), summary_table_title.size(), (table2_epCDn.size() - i), 0, "br");
                row.SetFillColor(kAzure - 9);
                row.SetFillStyle(1001);
                row.SetLineColor(0);
                row.Draw();
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
            TLine line_epCDn = TLine(0, i, summary_table_title.size(), i);
            line_epCDn.SetLineStyle(2);
            line_epCDn.Draw();
        }
        for (int j = 0; j <= summary_table_title.size(); j++) {
            // Vertical lines
            TLine line_epCDn = TLine(j, 0, j, table2_epCDn.size());
            line_epCDn.SetLineStyle(2);
            line_epCDn.Draw();
        }

        myTable.Print(fileName, "pdf");
        myTable.Clear();

        TH2F frame3_epCDn = TH2F("frame3_epCDn", "", summary_table_title.size(), 0, summary_table_title.size(), table3_epCDn.size(), 0, table3_epCDn.size());
        frame3_epCDn.SetDirectory(0);
        frame3_epCDn.SetStats(0);                   // Disable statistics box
        frame3_epCDn.GetXaxis().SetLabelSize(0);   // Remove x-axis labels
        frame3_epCDn.GetXaxis().SetTickLength(0);  // Remove x-axis ticks
        frame3_epCDn.GetYaxis().SetLabelSize(0);   // Remove y-axis labels
        frame3_epCDn.GetYaxis().SetTickLength(0);  // Remove y-axis ticks
        frame3_epCDn.Draw();

        titles.DrawLatexNDC(0.05, 0.9, "Step by step statistics - CD proton (table 3 of 3)");

        TPave row3 = TPave(0, (table3_epCDn.size() - 1), summary_table_title.size(), table3_epCDn.size(), 0, "br");  // Row 1 (top row)
        row3.SetFillColor(kGray + 1);
        // row3.SetFillColor(kAzure - 9);
        row3.SetFillStyle(1001);
        row3.SetLineColor(0);
        row3.Draw();

        // Loop over rows and columns to position text
        for (int i = 0; i < table3_epCDn.size(); i++) {
            std::string cellTemp = table3_epCDn.at(i).at(0).Data();

            if ((!basic_tools::FindSubstring(cellTemp, "}{(Step") && !basic_tools::FindSubstring(cellTemp, "}{cuts")) && (i > 0)) {
                TPave row = TPave(0, (table3_epCDn.size() - i - 1), summary_table_title.size(), (table3_epCDn.size() - i), 0, "br");
                row.SetFillColor(kAzure - 9);
                row.SetFillStyle(1001);
                row.SetLineColor(0);
                row.Draw();
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
            TLine line_epCDn =  TLine(0, i, summary_table_title.size(), i);
            line_epCDn.SetLineStyle(2);
            line_epCDn.Draw();
        }
        for (int j = 0; j <= summary_table_title.size(); j++) {
            // Vertical lines
            TLine line_epCDn = TLine(j, 0, j, table3_epCDn.size());
            line_epCDn.SetLineStyle(2);
            line_epCDn.Draw();
        }

        myTable.Print(fileName, "pdf");
        myTable.Clear();
        titles.SetTextSize(0.065);
    }

    if (Constraint1 == "" || Constraint1 == "FD") {
        myTable.SetTopMargin(0.15);

        double Num_of_goodN_bfSteps_epFDn = GetHistogramEntries(HistoList, "dpp_goodN_epFDn");
        double Num_of_badN_bfSteps_epFDn = GetHistogramEntries(HistoList, "dpp_badN_epFDn");
        string Num_of_allN_bfSteps_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_bfSteps_epFDn + Num_of_badN_bfSteps_epFDn, 0);
        string Num_of_goodN_bfSteps_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_bfSteps_epFDn, 0);
        string Num_of_badN_bfSteps_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_badN_bfSteps_epFDn, 0);
        string Single_eff_bfSteps_epFDn_str = "--";
        string Single_purity_bfSteps_epFDn_str = "--";
        const char *Num_of_allN_bfSteps_epFDn_char = Num_of_allN_bfSteps_epFDn_str.c_str();
        const char *Num_of_goodN_bfSteps_epFDn_char = Num_of_goodN_bfSteps_epFDn_str.c_str();
        const char *Num_of_badN_bfSteps_epFDn_char = Num_of_badN_bfSteps_epFDn_str.c_str();
        const char *Single_eff_bfSteps_epFDn_char = Single_eff_bfSteps_epFDn_str.c_str();
        const char *Single_purity_bfSteps_epFDn_char = Single_purity_bfSteps_epFDn_str.c_str();

        double Num_of_goodN_Step0_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step0_epFDn");
        double Num_of_badN_Step0_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step0_epFDn");
        string Num_of_allN_Step0_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_epFDn + Num_of_badN_Step0_epFDn, 0);
        string Num_of_goodN_Step0_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_epFDn, 0);
        string Num_of_badN_Step0_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_badN_Step0_epFDn, 0);
        string Single_eff_Step0_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_epFDn / Num_of_goodN_bfSteps_epFDn);
        string Single_purity_Step0_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step0_epFDn / (Num_of_goodN_Step0_epFDn + Num_of_badN_Step0_epFDn));
        const char *Num_of_allN_Step0_epFDn_char = Num_of_allN_Step0_epFDn_str.c_str();
        const char *Num_of_goodN_Step0_epFDn_char = Num_of_goodN_Step0_epFDn_str.c_str();
        const char *Num_of_badN_Step0_epFDn_char = Num_of_badN_Step0_epFDn_str.c_str();
        const char *Single_eff_Step0_epFDn_char = Single_eff_Step0_epFDn_str.c_str();
        const char *Single_purity_Step0_epFDn_char = Single_purity_Step0_epFDn_str.c_str();

        double Num_of_goodN_Step1_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step1_epFDn");
        double Num_of_badN_Step1_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step1_epFDn");
        string Num_of_allN_Step1_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_epFDn + Num_of_badN_Step1_epFDn, 0);
        string Num_of_goodN_Step1_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_epFDn, 0);
        string Num_of_badN_Step1_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_badN_Step1_epFDn, 0);
        string Single_eff_Step1_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_epFDn / Num_of_goodN_Step0_epFDn);
        string Single_purity_Step1_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step1_epFDn / (Num_of_goodN_Step1_epFDn + Num_of_badN_Step1_epFDn));
        const char *Num_of_allN_Step1_epFDn_char = Num_of_allN_Step1_epFDn_str.c_str();
        const char *Num_of_goodN_Step1_epFDn_char = Num_of_goodN_Step1_epFDn_str.c_str();
        const char *Num_of_badN_Step1_epFDn_char = Num_of_badN_Step1_epFDn_str.c_str();
        const char *Single_eff_Step1_epFDn_char = Single_eff_Step1_epFDn_str.c_str();
        const char *Single_purity_Step1_epFDn_char = Single_purity_Step1_epFDn_str.c_str();

        double Num_of_goodN_Step2_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step2_epFDn");
        double Num_of_badN_Step2_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step2_epFDn");
        string Num_of_allN_Step2_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_epFDn + Num_of_badN_Step2_epFDn, 0);
        string Num_of_goodN_Step2_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_epFDn, 0);
        string Num_of_badN_Step2_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_badN_Step2_epFDn, 0);
        string Single_eff_Step2_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_epFDn / Num_of_goodN_Step1_epFDn);
        string Single_purity_Step2_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step2_epFDn / (Num_of_goodN_Step2_epFDn + Num_of_badN_Step2_epFDn));
        const char *Num_of_allN_Step2_epFDn_char = Num_of_allN_Step2_epFDn_str.c_str();
        const char *Num_of_goodN_Step2_epFDn_char = Num_of_goodN_Step2_epFDn_str.c_str();
        const char *Num_of_badN_Step2_epFDn_char = Num_of_badN_Step2_epFDn_str.c_str();
        const char *Single_eff_Step2_epFDn_char = Single_eff_Step2_epFDn_str.c_str();
        const char *Single_purity_Step2_epFDn_char = Single_purity_Step2_epFDn_str.c_str();

        double Num_of_goodN_Step3_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step3_epFDn");
        double Num_of_badN_Step3_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step3_epFDn");
        string Num_of_allN_Step3_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step3_epFDn + Num_of_badN_Step3_epFDn, 0);
        string Num_of_goodN_Step3_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step3_epFDn, 0);
        string Num_of_badN_Step3_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_badN_Step3_epFDn, 0);
        string Single_eff_Step3_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step3_epFDn / Num_of_goodN_Step2_epFDn);
        string Single_purity_Step3_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step3_epFDn / (Num_of_goodN_Step3_epFDn + Num_of_badN_Step3_epFDn));
        const char *Num_of_allN_Step3_epFDn_char = Num_of_allN_Step3_epFDn_str.c_str();
        const char *Num_of_goodN_Step3_epFDn_char = Num_of_goodN_Step3_epFDn_str.c_str();
        const char *Num_of_badN_Step3_epFDn_char = Num_of_badN_Step3_epFDn_str.c_str();
        const char *Single_eff_Step3_epFDn_char = Single_eff_Step3_epFDn_str.c_str();
        const char *Single_purity_Step3_epFDn_char = Single_purity_Step3_epFDn_str.c_str();

        double Num_of_goodN_Step4_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step4_epFDn");
        double Num_of_badN_Step4_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step4_epFDn");
        string Num_of_allN_Step4_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step4_epFDn + Num_of_badN_Step4_epFDn, 0);
        string Num_of_goodN_Step4_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step4_epFDn, 0);
        string Num_of_badN_Step4_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_badN_Step4_epFDn, 0);
        string Single_eff_Step4_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step4_epFDn / Num_of_goodN_Step3_epFDn);
        string Single_purity_Step4_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step4_epFDn / (Num_of_goodN_Step4_epFDn + Num_of_badN_Step4_epFDn));
        const char *Num_of_allN_Step4_epFDn_char = Num_of_allN_Step4_epFDn_str.c_str();
        const char *Num_of_goodN_Step4_epFDn_char = Num_of_goodN_Step4_epFDn_str.c_str();
        const char *Num_of_badN_Step4_epFDn_char = Num_of_badN_Step4_epFDn_str.c_str();
        const char *Single_eff_Step4_epFDn_char = Single_eff_Step4_epFDn_str.c_str();
        const char *Single_purity_Step4_epFDn_char = Single_purity_Step4_epFDn_str.c_str();

        double Num_of_goodN_Step5_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_goodN_Step5_epFDn");
        double Num_of_badN_Step5_epFDn = GetHistogramEntries(HistoList, "beta_n_VS_Edep_CND_badN_Step5_epFDn");
        string Num_of_allN_Step5_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step5_epFDn + Num_of_badN_Step5_epFDn, 0);
        string Num_of_goodN_Step5_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step5_epFDn, 0);
        string Num_of_badN_Step5_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_badN_Step5_epFDn, 0);
        string Num_of_goodN_Step5_loss_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step5_epFDn / Num_of_goodN_Step4_epFDn);
        string Num_of_badN_Step5_loss_epFDn_str = basic_tools::ToStringWithPrecision(Num_of_goodN_Step5_epFDn / (Num_of_goodN_Step5_epFDn + Num_of_badN_Step5_epFDn));
        const char *Num_of_allN_Step5_epFDn_char = Num_of_allN_Step5_epFDn_str.c_str();
        const char *Num_of_goodN_Step5_epFDn_char = Num_of_goodN_Step5_epFDn_str.c_str();
        const char *Num_of_badN_Step5_epFDn_char = Num_of_badN_Step5_epFDn_str.c_str();
        const char *Num_of_goodN_Step5_loss_epFDn_char = Num_of_goodN_Step5_loss_epFDn_str.c_str();
        const char *Num_of_badN_Step5_loss_epFDn_char = Num_of_badN_Step5_loss_epFDn_str.c_str();

        // Draw a frame without axis numbers and ticks
        TH2F frame_epFDn = TH2F("frame_epFDn", "", 5, 0, 5, 8, 0, 8);
        frame_epFDn.SetDirectory(0);
        frame_epFDn.SetStats(0);                   // Disable statistics box
        frame_epFDn.GetXaxis().SetLabelSize(0);   // Remove x-axis labels
        frame_epFDn.GetXaxis().SetTickLength(0);  // Remove x-axis ticks
        frame_epFDn.GetYaxis().SetLabelSize(0);   // Remove y-axis labels
        frame_epFDn.GetYaxis().SetTickLength(0);  // Remove y-axis ticks
        frame_epFDn.Draw();

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
            TLine line_epFDn = TLine(0, i, 5, i);
            line_epFDn.SetLineStyle(2);
            line_epFDn.Draw();
        }
        for (int j = 0; j <= 5; j++) {
            // Vertical lines
            TLine line_epFDn = TLine(j, 0, j, 8);
            line_epFDn.SetLineStyle(2);
            line_epFDn.Draw();
        }

        myTable.Print(fileName, "pdf");
        myTable.Clear();
    }
}

