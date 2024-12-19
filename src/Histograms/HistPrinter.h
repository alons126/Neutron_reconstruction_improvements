#ifndef HISTPRINTER_H
#define HISTPRINTER_H

#include <cstdlib>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "../functions/GeneralFunctions.h"

using namespace std;

// extractStep function -------------------------------------------------------------------------------------------------------------------------------------------------------

std::string extractStep(const std::string &input)
{
    std::regex stepRegex(R"(Step\d+)"); // Regex to match "Step" followed by digits
    std::smatch match;

    if (std::regex_search(input, match, stepRegex))
    {
        return match.str(); // Return the matched substring
    }
    return ""; // Return an empty string if no match is found
}

// SkippingCondition function -------------------------------------------------------------------------------------------------------------------------------------------------

bool SkippingCondition(string HistoName)
{
    // TODO: fix this in the all plots file!
    if (HistoName == "Chi2pid_p_APID_epCD" || HistoName == "Chi2pid_p_APID_epFD"                                                      // Last PID plot
        || HistoName == "nSector_VS_ToF_epCDn" || HistoName == "nSector_VS_ToF_epFDn"                                                 // Last miss cuts plot
        || HistoName == "beta_n_badN_Step0_epCDn" || HistoName == "beta_n_badN_Step0_epFDn"                                           // Last Step0 plot
        || HistoName == "diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epCDn" || HistoName == "diff_ToFc_z_VS_Edep_yesNear_badN_Step1_epFDn" // Last Step1 plot
    )
    {
        return true;
    }

    return false;
}

// SectionPlotter function ----------------------------------------------------------------------------------------------------------------------------------------------------

void SectionPlotter(int n_col, int n_row, TCanvas *myCanvas, TCanvas *myText, vector<TH1 *> HistoList, string PDFFile, string Constraint1 = "", string Constraint2 = "")
{
    TLatex titles, text;
    titles.SetTextSize(0.065);
    text.SetTextSize(0.04);

    string pdfFile0;

    if (Constraint1 == "" && Constraint2 == "")
    {
        pdfFile0 = PDFFile;
    }
    else if (Constraint1 != "" && Constraint2 == "")
    {
        string pdfFile1 = ConfigOutPutName(PDFFile, Constraint1);
        pdfFile0 = pdfFile1;
    }
    else if (Constraint1 == "" && Constraint2 != "")
    {
        string pdfFile1 = ConfigOutPutName(PDFFile, Constraint2);
        pdfFile0 = pdfFile1;
    }
    else if (Constraint1 != "" && Constraint2 != "")
    {
        string pdfFile2 = ConfigOutPutName(PDFFile, Constraint1);
        string pdfFile1 = ConfigOutPutName(pdfFile2, Constraint2);
        pdfFile0 = pdfFile1;
    }

    const char *pdfFile = pdfFile0.c_str();

    char fileName[100];
    sprintf(fileName, "%s[", pdfFile);
    myText->SaveAs(fileName);
    sprintf(fileName, "%s", pdfFile);

    myText->cd();

    if (Constraint2 == "")
    {
        titles.DrawLatex(0.05, 0.9, "Manual Veto Plots");
    }
    else
    {
        titles.DrawLatex(0.05, 0.9, ("Manual Veto Plots - " + Constraint2).c_str());
    }

    text.DrawLatex(0.1, 0.7, "(e,e'p) Cuts:");

    if (Constraint1 == "")
    {
        text.DrawLatex(0.15, 0.6, "1 electron");
        text.DrawLatex(0.15, 0.5, "1 proton in CD or FD");
        text.DrawLatex(0.15, 0.4, "Any number of neutrons in CND");
        text.DrawLatex(0.15, 0.3, "Only particles with pdg=2112,11,2212,0,22 in event");
    }
    else if (Constraint1 == "CD")
    {
        text.DrawLatex(0.15, 0.6, "1 electron");
        text.DrawLatex(0.15, 0.5, "1 proton in CD");
        text.DrawLatex(0.15, 0.4, "Any number of neutrons in CND");
        text.DrawLatex(0.15, 0.3, "Only particles with pdg=2112,11,2212,0,22 in event");
    }
    else if (Constraint1 == "FD")
    {
        text.DrawLatex(0.15, 0.6, "1 electron");
        text.DrawLatex(0.15, 0.5, "1 proton in FD");
        text.DrawLatex(0.15, 0.4, "Any number of neutrons in CND");
        text.DrawLatex(0.15, 0.3, "Only particles with pdg=2112,11,2212,0,22 in event");
    }

    myText->Print(fileName, "pdf");
    myText->Clear();

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

    for (int i = 0; i < HistoList.size(); i++)
    {
        string TempHistName = HistoList[i]->GetName();

        bool GoodHistogram;

        if (Constraint1 == "" && Constraint2 == "")
        {
            GoodHistogram = true;
        }
        else if (Constraint1 != "" && Constraint2 == "")
        {
            GoodHistogram = findSubstring(TempHistName, Constraint1);
        }
        else if (Constraint1 == "" && Constraint2 != "")
        {
            GoodHistogram = findSubstring(TempHistName, Constraint2);
        }
        else
        {
            GoodHistogram = (findSubstring(TempHistName, Constraint1) && findSubstring(TempHistName, Constraint2));
        }

        if (GoodHistogram)
        {
            if (findSubstring(TempHistName, "BPID"))
            {
                if (FirstPIDPlot)
                {
                    myText->cd();

                    titles.DrawLatex(0.05, 0.9, "PID Plots");

                    if (Constraint1 == "")
                    {
                        text.DrawLatex(0.1, 0.8, "CD protons:");
                        text.DrawLatex(0.2, 0.7, "#fons[25]{#lbarV_{z}^{p} - V_{z}^{e}#lbar #leq 4} cm");
                        text.DrawLatex(0.2, 0.6, "0.3 #leq P_{p} #leq 1.5 GeV/c");
                        text.DrawLatex(0.2, 0.5, "#lbar#Delta#beta_{p}#lbar #leq 0.05");

                        text.DrawLatex(0.1, 0.4, "FD protons:");
                        text.DrawLatex(0.2, 0.3, "#lbarV_{z}^{p} - V_{z}^{e}#lbar #leq 5 cm");
                        text.DrawLatex(0.2, 0.2, "0.4 #leq P_{p} #leq 3.0 GeV/c");
                        text.DrawLatex(0.2, 0.1, "#lbar#Delta#beta_{p}#lbar #leq 0.03");
                    }
                    else if (Constraint1 == "CD")
                    {
                        text.DrawLatex(0.1, 0.7, "CD protons:");
                        text.DrawLatex(0.2, 0.6, "#lbarV_{z}^{p} - V_{z}^{e}#lbar #leq 4 cm");
                        text.DrawLatex(0.2, 0.5, "0.3 #leq P_{p} #leq 1.5 GeV/c");
                        text.DrawLatex(0.2, 0.4, "#lbar#Delta#beta_{p}#lbar #leq 0.05");
                    }
                    else if (Constraint1 == "FD")
                    {
                        text.DrawLatex(0.1, 0.7, "FD protons:");
                        text.DrawLatex(0.2, 0.6, "#lbarV_{z}^{p} - V_{z}^{e}#lbar #leq 5 cm");
                        text.DrawLatex(0.2, 0.5, "0.4 #leq P_{p} #leq 3.0 GeV/c");
                        text.DrawLatex(0.2, 0.4, "#lbar#Delta#beta_{p}#lbar #leq 0.03");
                    }

                    myText->Print(fileName, "pdf");
                    myText->Clear();

                    // titles.DrawLatex(0.05, 0.9, "Neutron cuts and definitions");
                    // text.DrawLatex(0.1, 0.8, "Neutron PID cuts:");
                    // text.DrawLatex(0.2, 0.7, "0.15 #leq #beta_{n} #leq 0.8");
                    // text.DrawLatex(0.2, 0.6, "#theta_{n} #leq 160#circ");
                    // text.DrawLatex(0.2, 0.5, "Status = 0 (no double-hits)");

                    // text.DrawLatex(0.1, 0.3, "Good neutron definition:");
                    // text.DrawLatex(0.2, 0.2, "#theta_{n,miss} #leq 25#circ");
                    // text.DrawLatex(0.2, 0.1, "#lbar#left(#lbar#vec{P}_{miss}#lbar - #lbar#vec{P}_{n}#lbar#right)/P_{miss}#lbar #leq 0.3");

                    // myText->Print(fileName, "pdf");
                    // myText->Clear();

                    FirstPIDPlot = false;
                }
            }
            else if (findSubstring(TempHistName, "BmissC"))
            {
                if (FirstOnlyMissCutsPlot)
                {
                    myText->cd();

                    titles.DrawLatex(0.05, 0.9, "Basic cuts & definitions");
                    text.DrawLatex(0.1, 0.8, "Missing variables cuts:");
                    text.DrawLatex(0.15, 0.7, "0.2 #leq P_{miss} #leq 1.5 GeV/c; 40#circ #leq #theta_{miss} #leq 135#circ; 0.7 #leq M_{miss} #leq 1.2 GeV/c^{2}");
                    // text.DrawLatex(0.2, 0.6, "40#circ #leq #theta_{miss} #leq 135#circ");
                    // text.DrawLatex(0.2, 0.5, "0.7 #leq M_{miss} #leq 1.2 GeV/c^{2}");

                    text.DrawLatex(0.1, 0.6, "Neutron PID cuts:");
                    text.DrawLatex(0.15, 0.5, "0.15 #leq #beta_{n} #leq 0.8; #theta_{n} #leq 160#circ; Status = 0 (no double-hits)");
                    // text.DrawLatex(0.2, 0.3, "#theta_{n} #leq 160#circ");
                    // text.DrawLatex(0.2, 0.1, "Status = 0 (no double-hits)");

                    text.DrawLatex(0.1, 0.4, "Good neutrons definition:");
                    text.DrawLatex(0.15, 0.3, "#theta_{n,miss} #leq 25#circ");
                    text.DrawLatex(0.15, 0.2, "#lbar#left(#lbar#vec{P}_{miss}#lbar - #lbar#vec{P}_{n}#lbar#right)/P_{miss}#lbar #leq 0.3");
                    text.DrawLatex(0.1, 0.1, "Bad neutrons definition: not good neutrons (TEMP!)");

                    // titles.DrawLatex(0.05, 0.9, "Before and After P_{miss}, #theta_{miss}, and M_{miss} Cuts Plots");
                    // text.DrawLatex(0.1, 0.7, "Used cuts:");
                    // text.DrawLatex(0.2, 0.6, "0.2 #leq P_{miss} #leq 1.5 GeV/c");
                    // text.DrawLatex(0.2, 0.5, "40#circ #leq #theta_{miss} #leq 135#circ");
                    // text.DrawLatex(0.2, 0.4, "0.7 #leq M_{miss} #leq 1.2 GeV/c^{2}");

                    myText->Print(fileName, "pdf");
                    myText->Clear();

                    FirstOnlyMissCutsPlot = false;
                }
            }
            else if (findSubstring(TempHistName, "Step"))
            {
                string Step = extractStep(TempHistName);

                if (FirstStepPlot[Step] == true)
                {
                    myText->cd();

                    titles.DrawLatex(0.05, 0.9, (Step + " Cuts").c_str());

                    if (Step == "Step0")
                    {
                        text.DrawLatex(0.15, 0.7, "#lbar#beta_{n} - L/(t_{ToF,n} * c)#lbar #leq 0.01");
                        text.DrawLatex(0.15, 0.6, "-40 #leq V_{hit,z} #leq 45 cm");
                        text.DrawLatex(0.15, 0.5, "0 #leq t_{ToF,n} #leq 20 ns");
                    }
                    else if (Step == "Step1")
                    {
                        text.DrawLatex(0.1, 0.4, "Step0 cuts:");
                        text.DrawLatex(0.15, 0.7, "#lbar#beta_{n} - L/(t_{ToF,n} * c)#lbar #leq 0.01");
                        text.DrawLatex(0.15, 0.6, "-40 #leq V_{hit,z} #leq 45 cm");
                        text.DrawLatex(0.15, 0.5, "0 #leq t_{ToF,n} #leq 20 ns");

                        text.DrawLatex(0.1, 0.4, "Step1 cuts:");
                        text.DrawLatex(0.15, 0.3, "5 #leq E_{dep}^{CND} #leq (#gamma_{n} - 1) * m_{n} MeV");
                    }
                    else if (Step == "Step2")
                    {
                        text.DrawLatex(0.15, 0.8, "Step0 cuts: #lbar#beta_{n} - L/(t_{ToF,n} * c)#lbar #leq 0.01; -40 #leq V_{hit,z} #leq 45 cm; 0 #leq t_{ToF,n} #leq 20 ns");
                        text.DrawLatex(0.15, 0.7, "Step1 cuts: 5 #leq E_{dep}^{CND} #leq (#gamma_{n} - 1) * m_{n})");

                        text.DrawLatex(0.15, 0.6, "No nearby hits associated with the charged particle track");
                        text.DrawLatex(0.15, 0.5, "Cluster width is 1 hit");
                        text.DrawLatex(0.15, 0.4, "Layer multiplicity:");
                        text.DrawLatex(0.2, 0.3, "Hit in CND1 #rightarrow layer multiplicity = 1");
                        text.DrawLatex(0.2, 0.2, "Hit in CND2 or CND3 #rightarrow layer multiplicity = 1 or 2");
                    }

                    myText->Print(fileName, "pdf");
                    myText->Clear();

                    FirstStepPlot[Step] = false;
                }
            }

            myCanvas->cd(canvas_ind);
            myCanvas->cd(canvas_ind)->SetBottomMargin(0.14), myCanvas->cd(canvas_ind)->SetLeftMargin(0.16), myCanvas->cd(canvas_ind)->SetRightMargin(0.16), myCanvas->cd(canvas_ind)->SetTopMargin(0.12);
            gPad->SetGrid();
            gPad->SetFrameLineWidth(1); // Reset frame line width to 1

            gStyle->SetOptStat("ourmen");
            // gStyle->SetOptStat(000111111);

            if (HistoList[i]->InheritsFrom("TH1D"))
            {
                HistoList[i]->SetMinimum(0);
                HistoList[i]->SetLineWidth(2);
                HistoList[i]->SetLineColor(kRed);
                // HistoList[i]->SetLineColor(kBlue);
            }

            if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize * 0.6), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);

                if (HistoList[i]->InheritsFrom("TH1D"))
                {
                    HistoList[i]->Draw(), displayText->Draw("same");
                }
                else if (HistoList[i]->InheritsFrom("TH2D"))
                {
                    HistoList[i]->Draw("COLZ"), displayText->Draw("same");
                }
                else
                {
                    cout << "ERROR! could not determine histogram class! Exiting...";
                    exit(0);
                }
            }
            else
            {
                if (HistoList[i]->InheritsFrom("TH1D"))
                // if (HistoList[i]->GetClassName() == "TH1D")
                {
                    HistoList[i]->Draw();
                }
                else if (HistoList[i]->InheritsFrom("TH2D"))
                // else if (HistoList[i]->GetClassName() == "TH2D")
                {
                    HistoList[i]->Draw("COLZ");
                }
                else
                {
                    cout << "ERROR! could not determine histogram class! Exiting...";
                    exit(0);
                }
            }

            // Save the canvas to a PDF page after filling 12 pads or processing the last histogram
            if (canvas_ind == n_col * n_row || SkippingCondition(TempHistName))
            {
                myCanvas->Print(fileName);      // Save the current page
                myCanvas->Clear();              // Clear the canvas for the next page
                myCanvas->Divide(n_col, n_row); // Reset the grid layout

                canvas_ind = 0;
            }

            ++canvas_ind;
        }
    }

    sprintf(fileName, "%s]", pdfFile);
    myCanvas->Print(fileName, "pdf");

    myCanvas->Clear();
    myText->Clear();
}

// HistPrinter function -------------------------------------------------------------------------------------------------------------------------------------------------------

void HistPrinter(vector<TH1 *> HistoList, string PDFFile)
{
#pragma region /* Andrew's wrap up - start */

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

    /* Saving all plots - start */
    SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile);

    /* Saving only CD proton plots - start */
    SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "CD");
    SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "CD", "Step0");
    SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "CD", "Step1");
    SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "CD", "Step2");
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "CD", "Step3");
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "CD", "Step4");
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "CD", "Step5");

    // /* Saving only FD proton plots */
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "FD");
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "FD", "Step0");
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "FD", "Step1");
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "FD", "Step2");
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "FD", "Step3");
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "FD", "Step4");
    // SectionPlotter(n_col, n_row, myCanvas, myText, HistoList, PDFFile, "FD", "Step5");

#pragma endregion /* Andrew's wrap up - end */
}

#endif // HISTPRINTER_H
