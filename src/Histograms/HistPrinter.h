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
    if (findSubstring(HistoName, "Chi2pid_p_APID_ep") // Last PID plot
                                                      // || findSubstring(HistoName, "Z_badN_Step0_ep") // Last Step0 plot
                                                      // || findSubstring(HistoName, "Z_badN_Step1_ep") // Last Step1 plot
    )
    {
        return true;
    }

    return false;
}

// SectionPlotter function ----------------------------------------------------------------------------------------------------------------------------------------------------

void SectionPlotter(TCanvas *myCanvas, TCanvas *myText, vector<TH1 *> HistoList, string PDFFile, string Constraint1 = "", string Constraint2 = "")
{
    TLatex text;
    text.SetTextSize(0.05);

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

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text.DrawLatex(0.1, 0.9, "(e,e'p) Cuts:");
    text.DrawLatex(0.1, 0.8, "(e,e') Cuts");
    text.DrawLatex(0.1, 0.7, "Neutrons in CND");

    myText->Print(fileName, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->Divide(4, 3);

    double x_1 = 0.2, y_1 = 0.3, x_2 = 0.86, y_2 = 0.7;
    double diplayTextSize = 0.1;

    int canvas_ind = 1;

    bool FilledConstraint1Bookmark = false;

    bool FirstPIDPlot = true;

    bool FirstOnlyMissCutsPlot = true;

    map<string, bool> FirstStepPlot;
    FirstStepPlot["Step0"] = true;
    FirstStepPlot["Step1"] = true;
    FirstStepPlot["Step2"] = true;
    FirstStepPlot["Step3"] = true;
    FirstStepPlot["Step4"] = true;
    FirstStepPlot["Step5"] = true;

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

                    text.DrawLatex(0.1, 0.9, "PID Plots");
                    text.DrawLatex(0.1, 0.8, "");
                    text.DrawLatex(0.1, 0.7, "(e,e'p) Cuts:");
                    text.DrawLatex(0.1, 0.6, "(e,e') Cuts");
                    text.DrawLatex(0.1, 0.5, "Neutrons in CND");

                    myText->Print(fileName, "pdf");
                    myText->Clear();

                    FirstPIDPlot = false;
                }
            }
            else if (findSubstring(TempHistName, "BmissC"))
            {
                if (FirstOnlyMissCutsPlot)
                {
                    myText->cd();

                    text.DrawLatex(0.1, 0.9, "Before and after P_{miss}, #theta_{miss}, and M_{miss} Cuts Plots");
                    text.DrawLatex(0.1, 0.8, "");
                    text.DrawLatex(0.1, 0.7, "(e,e'p) Cuts:");
                    text.DrawLatex(0.1, 0.6, "(e,e') Cuts");
                    text.DrawLatex(0.1, 0.5, "Neutrons in CND");

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

                    text.DrawLatex(0.1, 0.9, (Step + " Plots").c_str());
                    text.DrawLatex(0.1, 0.8, "");
                    text.DrawLatex(0.1, 0.7, "(e,e'p) Cuts:");
                    text.DrawLatex(0.1, 0.6, "(e,e') Cuts");
                    text.DrawLatex(0.1, 0.5, "Neutrons in CND");

                    myText->Print(fileName, "pdf");
                    myText->Clear();

                    FirstStepPlot[Step] = false;
                }
            }

            myCanvas->cd(canvas_ind);
            myCanvas->cd(canvas_ind)->SetBottomMargin(0.14), myCanvas->cd(canvas_ind)->SetLeftMargin(0.16), myCanvas->cd(canvas_ind)->SetRightMargin(0.16), myCanvas->cd(canvas_ind)->SetTopMargin(0.12);
            gPad->SetGrid();

            HistoList[i]->SetLineWidth(1);
            HistoList[i]->SetLineColor(kBlue);

            if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize * 0.6), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);

                if (HistoList[i]->InheritsFrom("TH1D"))
                // if (HistoList[i]->GetClassName() == "TH1D")
                {
                    HistoList[i]->Draw(), displayText->Draw("same");
                }
                else if (HistoList[i]->InheritsFrom("TH2D"))
                // else if (HistoList[i]->GetClassName() == "TH2D")
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
            if (canvas_ind == 12 || SkippingCondition(TempHistName))
            {
                // if (findSubstring(TempHistName, Constraint1) && !FilledConstraint1Bookmark)
                // {
                //     string bookmarkOption = "pdfBookmark=" + Constraint1 + "_Proton_Plots";
                //     myCanvas->Print(fileName, bookmarkOption.c_str());
                //     FilledConstraint1Bookmark = true;
                // }

                myCanvas->Print(fileName); // Save the current page
                myCanvas->Clear();         // Clear the canvas for the next page
                myCanvas->Divide(4, 3);    // Reset the grid layout

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

    // int pixelx = 1980, pixely = 1530;
    // int pixelx = 1980 * 4, pixely = 1530 * 4;
    // int pixelx = 1980 * 4 * 1.5 * 2, pixely = 1530 * 4 * 1.5 * 2;
    // int pixelx = 1000 * 4 * 1.5 * 2, pixely = 750 * 3 * 1.5 * 2;
    int pixelx = 1000 * 4 * 4, pixely = 750 * 3 * 4;

    TCanvas *myCanvas = new TCanvas("myPage", "myPage", pixelx, pixely);
    TCanvas *myText = new TCanvas("myText", "myText", pixelx, pixely);

    /* Saving all plots - start */
    SectionPlotter(myCanvas, myText, HistoList, PDFFile);

    /* Saving only CD proton plots - start */
    SectionPlotter(myCanvas, myText, HistoList, PDFFile, "CD");
    SectionPlotter(myCanvas, myText, HistoList, PDFFile, "CD", "Step0");
    SectionPlotter(myCanvas, myText, HistoList, PDFFile, "CD", "Step1");
    SectionPlotter(myCanvas, myText, HistoList, PDFFile, "CD", "Step2");
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "CD", "Step3");
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "CD", "Step4");
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "CD", "Step5");

    // /* Saving only FD proton plots */
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "FD");
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "FD", "Step0");
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "FD", "Step1");
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "FD", "Step2");
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "FD", "Step3");
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "FD", "Step4");
    // SectionPlotter(myCanvas, myText, HistoList, PDFFile, "FD", "Step5");

#pragma endregion /* Andrew's wrap up - end */
}

#endif // HISTPRINTER_H
