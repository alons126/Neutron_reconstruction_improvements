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

void HistPrinter(vector<TH1 *> hist_list_1_A,vector<TH2 *> hist_list_2_A,string PDFFile)
{
#pragma region /* Andrew's wrap up - start */

    /////////////////////////////////////////////////////
    // Now create the output PDFs
    /////////////////////////////////////////////////////

    int pixelx = 1980, pixely = 1530;
    // int pixelx = 1000 * 4 * 1.5 * 2, pixely = 750 * 3 * 1.5 * 2;

    TCanvas *myCanvas = new TCanvas("myPage", "myPage", pixelx, pixely);
    TCanvas *myText = new TCanvas("myText", "myText", pixelx, pixely);

#pragma region /* Saving all plots - start */

    TLatex text;
    text.SetTextSize(0.05);

    const char *pdfFile = PDFFile.c_str();

    char fileName[100];
    sprintf(fileName, "%s[", pdfFile);
    myText->SaveAs(fileName);
    sprintf(fileName, "%s", pdfFile);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text.DrawLatex(0.2, 0.7, "Neutrons in CND");

    myText->Print(fileName, "pdf");
    myText->Clear();

    myCanvas->cd();
    // myCanvas->SetGrid();
    myCanvas->Divide(4, 3);
    // myCanvas->SetGrid(), myCanvas->cd()->SetBottomMargin(0.14), myCanvas->cd()->SetLeftMargin(0.16), myCanvas->cd()->SetRightMargin(0.16), myCanvas->cd()->SetTopMargin(0.12);

    double x_1 = 0.2, y_1 = 0.3, x_2 = 0.86, y_2 = 0.7;
    double diplayTextSize = 0.1;

    // int canvas_ind = 1;

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        int canvas_ind = (i % 12) + 1; // Determine the pad number (1 to 12)

        myCanvas->cd(canvas_ind);
        gPad->SetGrid();

        hist_list_1_A[i]->SetLineWidth(2);
        hist_list_1_A[i]->SetLineColor(kBlue);

        if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
        {
            TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
            displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
            hist_list_1_A[i]->Draw(), displayText->Draw("same");
        }
        else
        {
            hist_list_1_A[i]->Draw();
        }

        // Save the canvas to a PDF page after filling 12 pads or processing the last histogram
        if (canvas_ind == 12 || i == hist_list_1_A.size() - 1)
        {
            myCanvas->Print(fileName); // Save the current page
            if (i != hist_list_1_A.size() - 1)
            {
                myCanvas->Clear();      // Clear the canvas for the next page
                myCanvas->Divide(4, 3); // Reset the grid layout
            }
        }
    }

    // for (int i = 0; i < hist_list_1_A.size(); i++)
    // {
    //     // myCanvas->cd(canvas_ind);
    //     // gPad->SetGrid();

    //     hist_list_1_A[i]->SetLineWidth(2);
    //     hist_list_1_A[i]->SetLineColor(kBlue);

    //     if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
    //     {
    //         TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //         displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //         hist_list_1_A[i]->Draw(), displayText->Draw("same");
    //     }
    //     else
    //     {
    //         hist_list_1_A[i]->Draw();
    //     }

    //     myCanvas->Print(fileName, "pdf");
    //     myCanvas->Clear();

    //     // ++canvas_ind;

    //     // if (i > 12 && 12 % i == 0)
    //     // {
    //     //     myCanvas->Print(fileName, "pdf");
    //     //     myCanvas->Clear();
    //     //     myCanvas->Divide(4, 3);

    //     //     canvas_ind = 1;
    //     // }
    // }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        myCanvas->cd(1);

        if (hist_list_2_A[i]->GetEntries() == 0 || hist_list_2_A[i]->Integral() == 0)
        {
            TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
            displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
            hist_list_2_A[i]->Draw("colz"), displayText->Draw("same");
        }
        else
        {
            hist_list_2_A[i]->Draw("colz");
        }

        myCanvas->Print(fileName, "pdf");
        myCanvas->Clear();
    }

    sprintf(fileName, "%s]", pdfFile);
    myCanvas->Print(fileName, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving all plots - end */

#pragma region /* Saving only CD proton plots - start */

    TLatex text_CD;
    text_CD.SetTextSize(0.05);

    string pdfFile_CD_0 = ConfigOutPutName(PDFFile, "pCD_only").c_str();
    const char *pdfFile_CD = pdfFile_CD_0.c_str();

    char fileName_CD[100];
    sprintf(fileName_CD, "%s[", pdfFile_CD);
    myText->SaveAs(fileName_CD);
    sprintf(fileName_CD, "%s", pdfFile_CD);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text_CD.DrawLatex(0.2, 0.9, "(e,e'pCD) Cuts:");
    text_CD.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text_CD.DrawLatex(0.2, 0.7, "Neutrons in CND");

    myText->Print(fileName_CD, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->SetGrid();
    // myCanvas->Divide(4, 3);
    // myCanvas->SetGrid(), myCanvas->cd()->SetBottomMargin(0.14), myCanvas->cd()->SetLeftMargin(0.16), myCanvas->cd()->SetRightMargin(0.16), myCanvas->cd()->SetTopMargin(0.12);

    // double x_1 = 0.2, y_1 = 0.3, x_2 = 0.86, y_2 = 0.7;
    // double diplayTextSize = 0.1;

    // int canvas_ind = 1;

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        string TempHistName = hist_list_1_A[i]->GetName();

        if (findSubstring(TempHistName, "CD"))
        {
            // myCanvas->cd(canvas_ind);
            // gPad->SetGrid();

            hist_list_1_A[i]->SetLineWidth(2);
            hist_list_1_A[i]->SetLineColor(kBlue);

            if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                hist_list_1_A[i]->Draw(), displayText->Draw("same");
            }
            else
            {
                hist_list_1_A[i]->Draw();
            }

            myCanvas->Print(fileName_CD, "pdf");
            myCanvas->Clear();

            // ++canvas_ind;

            // if (i > 12 && 12 % i == 0)
            // {
            //     myCanvas->Print(fileName_CD, "pdf");
            //     myCanvas->Clear();
            //     myCanvas->Divide(4, 3);

            //     canvas_ind = 1;
            // }
        }
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        string TempHistName = hist_list_2_A[i]->GetName();

        if (findSubstring(TempHistName, "CD"))
        {
            myCanvas->cd(1);

            if (hist_list_2_A[i]->GetEntries() == 0 || hist_list_2_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                hist_list_2_A[i]->Draw("colz"), displayText->Draw("same");
            }
            else
            {
                hist_list_2_A[i]->Draw("colz");
            }

            myCanvas->Print(fileName_CD, "pdf");
            myCanvas->Clear();
        }
    }

    sprintf(fileName_CD, "%s]", pdfFile_CD);
    myCanvas->Print(fileName_CD, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving only CD proton plots - end */

#pragma region /* Saving Step0 plots - start */

    TLatex text_Step0;
    text_Step0.SetTextSize(0.05);

    string pdfFile_Step0_0 = ConfigOutPutName(PDFFile, "Step0").c_str();
    const char *pdfFile_Step0 = pdfFile_Step0_0.c_str();

    cout << "\nPDFFile = " << PDFFile << "\n";
    cout << "\npdfFile_Step0_0 = " << pdfFile_Step0_0 << "\n";
    cout << "\npdfFile_Step0 = " << pdfFile_Step0 << "\n";

    char fileName_Step0[100];
    sprintf(fileName_Step0, "%s[", pdfFile_Step0);
    myText->SaveAs(fileName_Step0);
    sprintf(fileName_Step0, "%s", pdfFile_Step0);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text_Step0.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text_Step0.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text_Step0.DrawLatex(0.2, 0.7, "Neutrons in CND - step 0");

    myText->Print(fileName_Step0, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->SetGrid();

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        string TempHistName = hist_list_1_A[i]->GetName();

        if (findSubstring(TempHistName, "Step0"))
        {
            myCanvas->cd(1);
            hist_list_1_A[i]->SetLineWidth(2);
            hist_list_1_A[i]->SetLineColor(kBlue);

            if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                displayText->SetBorderSize(6);
                hist_list_1_A[i]->Draw(), displayText->Draw("same");
            }
            else
            {
                hist_list_1_A[i]->Draw();
            }

            myCanvas->Print(fileName_Step0, "pdf");
            myCanvas->Clear();
        }
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        string TempHistName = hist_list_2_A[i]->GetName();

        if (findSubstring(TempHistName, "Step0"))
        {
            myCanvas->cd(1);

            if (hist_list_2_A[i]->GetEntries() == 0 || hist_list_2_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                hist_list_2_A[i]->Draw("colz"), displayText->Draw("same");
            }
            else
            {
                hist_list_2_A[i]->Draw("colz");
            }

            myCanvas->Print(fileName_Step0, "pdf");
            myCanvas->Clear();
        }
    }

    sprintf(fileName_Step0, "%s]", pdfFile_Step0);
    myCanvas->Print(fileName_Step0, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving Step0 plots - end */

#pragma region /* Saving Step1 plots - start */

    TLatex text_Step1;
    text_Step1.SetTextSize(0.05);

    string pdfFile_Step1_0 = ConfigOutPutName(PDFFile, "Step1").c_str();
    const char *pdfFile_Step1 = pdfFile_Step1_0.c_str();

    char fileName_Step1[100];
    sprintf(fileName_Step1, "%s[", pdfFile_Step1);
    myText->SaveAs(fileName_Step1);
    sprintf(fileName_Step1, "%s", pdfFile_Step1);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text_Step1.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text_Step1.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text_Step1.DrawLatex(0.2, 0.7, "Neutrons in CND - step 1");

    myText->Print(fileName_Step1, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->SetGrid();

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        string TempHistName = hist_list_1_A[i]->GetName();

        if (findSubstring(TempHistName, "Step1"))
        {
            myCanvas->cd(1);
            hist_list_1_A[i]->SetLineWidth(2);
            hist_list_1_A[i]->SetLineColor(kBlue);

            if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                displayText->SetBorderSize(6);
                hist_list_1_A[i]->Draw(), displayText->Draw("same");
            }
            else
            {
                hist_list_1_A[i]->Draw();
            }

            myCanvas->Print(fileName_Step1, "pdf");
            myCanvas->Clear();
        }
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        string TempHistName = hist_list_2_A[i]->GetName();

        if (findSubstring(TempHistName, "Step1"))
        {
            myCanvas->cd(1);

            if (hist_list_2_A[i]->GetEntries() == 0 || hist_list_2_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                hist_list_2_A[i]->Draw("colz"), displayText->Draw("same");
            }
            else
            {
                hist_list_2_A[i]->Draw("colz");
            }

            myCanvas->Print(fileName_Step1, "pdf");
            myCanvas->Clear();
        }
    }

    sprintf(fileName_Step1, "%s]", pdfFile_Step1);
    myCanvas->Print(fileName_Step1, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving Step1 plots - end */

#pragma region /* Saving Step2 plots - start */

    TLatex text_Step2;
    text_Step2.SetTextSize(0.05);

    string pdfFile_Step2_0 = ConfigOutPutName(PDFFile, "Step2").c_str();
    const char *pdfFile_Step2 = pdfFile_Step2_0.c_str();

    char fileName_Step2[100];
    sprintf(fileName_Step2, "%s[", pdfFile_Step2);
    myText->SaveAs(fileName_Step2);
    sprintf(fileName_Step2, "%s", pdfFile_Step2);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text_Step2.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text_Step2.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text_Step2.DrawLatex(0.2, 0.7, "Neutrons in CND - step 2");

    myText->Print(fileName_Step2, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->SetGrid();

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        string TempHistName = hist_list_1_A[i]->GetName();

        if (findSubstring(TempHistName, "Step2"))
        {
            myCanvas->cd(1);
            hist_list_1_A[i]->SetLineWidth(2);
            hist_list_1_A[i]->SetLineColor(kBlue);

            if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                displayText->SetBorderSize(6);
                hist_list_1_A[i]->Draw(), displayText->Draw("same");
            }
            else
            {
                hist_list_1_A[i]->Draw();
            }

            myCanvas->Print(fileName_Step2, "pdf");
            myCanvas->Clear();
        }
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        string TempHistName = hist_list_2_A[i]->GetName();

        if (findSubstring(TempHistName, "Step2"))
        {
            myCanvas->cd(1);

            if (hist_list_2_A[i]->GetEntries() == 0 || hist_list_2_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                hist_list_2_A[i]->Draw("colz"), displayText->Draw("same");
            }
            else
            {
                hist_list_2_A[i]->Draw("colz");
            }

            myCanvas->Print(fileName_Step2, "pdf");
            myCanvas->Clear();
        }
    }

    sprintf(fileName_Step2, "%s]", pdfFile_Step2);
    myCanvas->Print(fileName_Step2, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving Step2 plots - end */

#pragma region /* Saving Step3 plots - start */

    TLatex text_Step3;
    text_Step3.SetTextSize(0.05);

    string pdfFile_Step3_0 = ConfigOutPutName(PDFFile, "Step3").c_str();
    const char *pdfFile_Step3 = pdfFile_Step3_0.c_str();

    char fileName_Step3[100];
    sprintf(fileName_Step3, "%s[", pdfFile_Step3);
    myText->SaveAs(fileName_Step3);
    sprintf(fileName_Step3, "%s", pdfFile_Step3);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text_Step3.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text_Step3.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text_Step3.DrawLatex(0.2, 0.7, "Neutrons in CND - step 3");

    myText->Print(fileName_Step3, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->SetGrid();

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        string TempHistName = hist_list_1_A[i]->GetName();

        if (findSubstring(TempHistName, "Step3"))
        {
            myCanvas->cd(1);
            hist_list_1_A[i]->SetLineWidth(2);
            hist_list_1_A[i]->SetLineColor(kBlue);

            if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                displayText->SetBorderSize(6);
                hist_list_1_A[i]->Draw(), displayText->Draw("same");
            }
            else
            {
                hist_list_1_A[i]->Draw();
            }

            myCanvas->Print(fileName_Step3, "pdf");
            myCanvas->Clear();
        }
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        string TempHistName = hist_list_2_A[i]->GetName();

        if (findSubstring(TempHistName, "Step3"))
        {
            myCanvas->cd(1);

            if (hist_list_2_A[i]->GetEntries() == 0 || hist_list_2_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                hist_list_2_A[i]->Draw("colz"), displayText->Draw("same");
            }
            else
            {
                hist_list_2_A[i]->Draw("colz");
            }

            myCanvas->Print(fileName_Step3, "pdf");
            myCanvas->Clear();
        }
    }

    sprintf(fileName_Step3, "%s]", pdfFile_Step3);
    myCanvas->Print(fileName_Step3, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving Step3 plots - end */

#pragma region /* Saving Step4 plots - start */

    TLatex text_Step4;
    text_Step4.SetTextSize(0.05);

    string pdfFile_Step4_0 = ConfigOutPutName(PDFFile, "Step4").c_str();
    const char *pdfFile_Step4 = pdfFile_Step4_0.c_str();

    char fileName_Step4[100];
    sprintf(fileName_Step4, "%s[", pdfFile_Step4);
    myText->SaveAs(fileName_Step4);
    sprintf(fileName_Step4, "%s", pdfFile_Step4);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text_Step4.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text_Step4.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text_Step4.DrawLatex(0.2, 0.7, "Neutrons in CND - step 4");

    myText->Print(fileName_Step4, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->SetGrid();

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        string TempHistName = hist_list_1_A[i]->GetName();

        if (findSubstring(TempHistName, "Step4"))
        {
            myCanvas->cd(1);
            hist_list_1_A[i]->SetLineWidth(2);
            hist_list_1_A[i]->SetLineColor(kBlue);

            if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                displayText->SetBorderSize(6);
                hist_list_1_A[i]->Draw(), displayText->Draw("same");
            }
            else
            {
                hist_list_1_A[i]->Draw();
            }

            myCanvas->Print(fileName_Step4, "pdf");
            myCanvas->Clear();
        }
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        string TempHistName = hist_list_2_A[i]->GetName();

        if (findSubstring(TempHistName, "Step4"))
        {
            myCanvas->cd(1);
            if (hist_list_2_A[i]->GetEntries() == 0 || hist_list_2_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                hist_list_2_A[i]->Draw("colz"), displayText->Draw("same");
            }
            else
            {
                hist_list_2_A[i]->Draw("colz");
            }

            myCanvas->Print(fileName_Step4, "pdf");
            myCanvas->Clear();
        }
    }

    sprintf(fileName_Step4, "%s]", pdfFile_Step4);
    myCanvas->Print(fileName_Step4, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving Step4 plots - end */

#pragma region /* Saving Step5 plots - start */

    TLatex text_Step5;
    text_Step5.SetTextSize(0.05);

    string pdfFile_Step5_0 = ConfigOutPutName(PDFFile, "Step5").c_str();
    const char *pdfFile_Step5 = pdfFile_Step5_0.c_str();

    char fileName_Step5[100];
    sprintf(fileName_Step5, "%s[", pdfFile_Step5);
    myText->SaveAs(fileName_Step5);
    sprintf(fileName_Step5, "%s", pdfFile_Step5);

    /////////////////////////////////////
    // CND Neutron Information
    /////////////////////////////////////

    myText->cd();

    text_Step5.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    text_Step5.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    text_Step5.DrawLatex(0.2, 0.7, "Neutrons in CND - step 5");

    myText->Print(fileName_Step5, "pdf");
    myText->Clear();

    myCanvas->cd();
    myCanvas->SetGrid();

    for (int i = 0; i < hist_list_1_A.size(); i++)
    {
        string TempHistName = hist_list_1_A[i]->GetName();

        if (findSubstring(TempHistName, "Step5"))
        {
            myCanvas->cd(1);
            hist_list_1_A[i]->SetLineWidth(2);
            hist_list_1_A[i]->SetLineColor(kBlue);

            if (hist_list_1_A[i]->GetEntries() == 0 || hist_list_1_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                displayText->SetBorderSize(6);
                hist_list_1_A[i]->Draw(), displayText->Draw("same");
            }
            else
            {
                hist_list_1_A[i]->Draw();
            }

            myCanvas->Print(fileName_Step5, "pdf");
            myCanvas->Clear();
        }
    }

    for (int i = 0; i < hist_list_2_A.size(); i++)
    {
        string TempHistName = hist_list_2_A[i]->GetName();

        if (findSubstring(TempHistName, "Step5"))
        {
            myCanvas->cd(1);

            if (hist_list_2_A[i]->GetEntries() == 0 || hist_list_2_A[i]->Integral() == 0)
            {
                TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
                displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
                hist_list_2_A[i]->Draw("colz"), displayText->Draw("same");
            }
            else
            {
                hist_list_2_A[i]->Draw("colz");
            }

            myCanvas->Print(fileName_Step5, "pdf");
            myCanvas->Clear();
        }
    }

    sprintf(fileName_Step5, "%s]", pdfFile_Step5);
    myCanvas->Print(fileName_Step5, "pdf");

    myCanvas->Clear();
    myText->Clear();

#pragma endregion /* Saving Step5 plots - end */

#pragma endregion /* Andrew's wrap up - end */
}

#endif // HISTPRINTER_H
