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

#include "clas12reader.h"
#include "HipoChain.h"

#include "src/constants.h"
// #include "src/Histograms/ManualVetoHistograms.h"
#include "src/Histograms/HistPrinter.h"
#include "src/functions/GeneralFunctions.h"
#include "src/functions/NeutronFunctions.h"
#include "src/functions/neutron-veto/veto_functions.cpp"
#include "src/functions/Andrews_functions/Andrews_functions.cpp"
#include "src/functions/HipoChain_config.cpp"
#include "src/classes/clas12ana/clas12ana.cpp"

using namespace std;
using namespace clas12;

#pragma region /* Erin main function - start */

int D_getfeatures_Phase6(                                                                             //
    const string OutDir, string output_pdf_Erin,                                                      // My arguments
    double Ebeam, bool keep_good, string output_root_Erin, string output_txt_Erin, string input_hipo, // Erin's arguments
    string PDFFile, int isMC = 0                                                                      // Andrew's arguments
)
// int main(int argc, char **argv)
{
    auto Code_start_time = std::chrono::system_clock::now(); // Start counting running time

    // ======================================================================================================================================================================
    // Printouts
    // ======================================================================================================================================================================

#pragma region /* Printouts - start */

    cout << "\033[33m\n\033[0m";
    // cout << "\033[33mHIPO_FILES:\033[0m\t\t" << gSystem->Getenv("HIPO_FILES") << "\n";
    cout << "\033[33minput_hipo:\033[0m\t\t" << input_hipo << "\n";
    cout << "\033[33m\n\033[0m";
    cout << "\033[33mOUTDIR:\033[0m\t\t\t" << gSystem->Getenv("OUTDIR") << "\n";
    cout << "\033[33mOutDir:\033[0m\t\t\t" << OutDir << "\n";
    cout << "\033[33moutput_pdf_Erin:\033[0m\t" << output_pdf_Erin << "\n";
    cout << "\033[33moutput_root_Erin:\033[0m\t" << output_root_Erin << "\n";
    cout << "\033[33moutput_txt_Erin:\033[0m\t" << output_txt_Erin << "\n";
    cout << "\033[33mPDFFile:\033[0m\t\t" << PDFFile << "\n\n";

#pragma endregion /* Printouts - end */

    // ======================================================================================================================================================================
    // Initial setup
    // ======================================================================================================================================================================

#pragma region /* Initial setup - start */

    // Delete old output folder
    cout << "\033[33m\nClearing\033[0m '" << OutDir << "'\n";
    system(("rm -r " + OutDir).c_str());
    cout << "\n";

    // Remake old output folder
    cout << "\033[33m\nRemaking\033[0m '" << OutDir << "'\n";
    system(("mkdir -p " + OutDir).c_str());
    cout << "\n\n";

    // Erin's output file names
    TFile *f = new TFile(output_root_Erin.c_str(), "RECREATE");
    TTree *ntree = new TTree("T", "NeutronTree");
    std::ofstream outtxt(output_txt_Erin);

    // Input hipo file
    clas12root::HipoChain chain;
    HipoChain_config(chain, input_hipo);

    auto config_c12 = chain.GetC12Reader();
    chain.SetReaderTags({0});
    const std::unique_ptr<clas12::clas12reader> &c12 = chain.C12ref();
    chain.db()->turnOffQADB();

    int numevent = 0;

    // Set up root tree for TMVA
    Int_t nhits;
    double px, py, pz, momentum;
    Int_t sec[100] = {-1};
    Int_t lay[100] = {-1};
    int event;
    double energy, cnd_energy, ctof_energy, angle_diff;
    int layermult, size, cnd_hits, ctof_hits;
    bool is_CTOF, is_CND1, is_CND2, is_CND3;

    int counter = 0;
    // cout << endl;

    // set up instance of clas12ana
    clas12ana *clasAna = new clas12ana();

    clasAna->readEcalSFPar("src/cuts/paramsSF_LD2_x2.dat"); // TODO: check if applied
    clasAna->readEcalPPar("src/cuts/paramsPI_LD2_x2.dat");  // TODO: check if applied

    clasAna->setProtonPidCuts(true);

#pragma endregion /* Initial setup - end */

    // ======================================================================================================================================================================
    // Andrew's histograms
    // ======================================================================================================================================================================


#pragma region /* Chain loop - start */

    int counter_A = 0; /* From Andrew */

    int counter_epXn = 0;
    int counter_pass_step0_cuts = 0, counter_pass_step1_cuts = 0, counter_pass_step2_cuts = 0, counter_pass_step3_cuts = 0, counter_pass_step4_cuts = 0, counter_pass_step5_cuts = 0;
    int counter_n_multiplicity_allN_epCDn = 0, counter_n_multiplicity_goodN_epCDn = 0, counter_n_multiplicity_badN_epCDn = 0;
    int counter_n_multiplicity_allN_epCDn_Step0 = 0, counter_n_multiplicity_goodN_epCDn_Step0 = 0, counter_n_multiplicity_badN_epCDn_Step0 = 0;
    int counter_n_multiplicity_allN_epCDn_Step1 = 0, counter_n_multiplicity_goodN_epCDn_Step1 = 0, counter_n_multiplicity_badN_epCDn_Step1 = 0;
    int counter_n_multiplicity_allN_epCDn_Step2 = 0, counter_n_multiplicity_goodN_epCDn_Step2 = 0, counter_n_multiplicity_badN_epCDn_Step2 = 0;
    int counter_n_multiplicity_allN_epCDn_Step3 = 0, counter_n_multiplicity_goodN_epCDn_Step3 = 0, counter_n_multiplicity_badN_epCDn_Step3 = 0;
    int counter_n_multiplicity_allN_epCDn_Step4 = 0, counter_n_multiplicity_goodN_epCDn_Step4 = 0, counter_n_multiplicity_badN_epCDn_Step4 = 0;
    int counter_n_multiplicity_allN_epCDn_Step5 = 0, counter_n_multiplicity_goodN_epCDn_Step5 = 0, counter_n_multiplicity_badN_epCDn_Step5 = 0;
    int counter_n_multiplicity_allN_epFDn = 0, counter_n_multiplicity_goodN_epFDn = 0, counter_n_multiplicity_badN_epFDn = 0;
    int counter_n_multiplicity_allN_epFDn_Step0 = 0, counter_n_multiplicity_goodN_epFDn_Step0 = 0, counter_n_multiplicity_badN_epFDn_Step0 = 0;
    int counter_n_multiplicity_allN_epFDn_Step1 = 0, counter_n_multiplicity_goodN_epFDn_Step1 = 0, counter_n_multiplicity_badN_epFDn_Step1 = 0;
    int counter_n_multiplicity_allN_epFDn_Step2 = 0, counter_n_multiplicity_goodN_epFDn_Step2 = 0, counter_n_multiplicity_badN_epFDn_Step2 = 0;
    int counter_n_multiplicity_allN_epFDn_Step3 = 0, counter_n_multiplicity_goodN_epFDn_Step3 = 0, counter_n_multiplicity_badN_epFDn_Step3 = 0;
    int counter_n_multiplicity_allN_epFDn_Step4 = 0, counter_n_multiplicity_goodN_epFDn_Step4 = 0, counter_n_multiplicity_badN_epFDn_Step4 = 0;
    int counter_n_multiplicity_allN_epFDn_Step5 = 0, counter_n_multiplicity_goodN_epFDn_Step5 = 0, counter_n_multiplicity_badN_epFDn_Step5 = 0;

    while (chain.Next())
    {
        // Display completed (from Andrew)
        counter_A++;

        if ((counter_A % 1000000) == 0)
        {
            cerr << "\n\n";
            cerr << "\033[33m" << counter_A / 1000000 << " million completed\033[0m\n\n";
        }

    } // closes event loop

#pragma endregion /* Chain loop - end */

    // ======================================================================================================================================================================
    // Andrew's wrap up
    // ======================================================================================================================================================================


#pragma region /* Andrew's wrap up - start */

    //     /////////////////////////////////////////////////////
    //     // Now create the output PDFs
    //     /////////////////////////////////////////////////////

    //     int pixelx = 1980, pixely = 1530;
    //     // int pixelx = 1000 * 4 * 1.5 * 2, pixely = 750 * 3 * 1.5 * 2;

    //     TCanvas *myCanvas = new TCanvas("myPage", "myPage", pixelx, pixely);
    //     TCanvas *myText = new TCanvas("myText", "myText", pixelx, pixely);

    // #pragma region /* Saving all plots - start */

    //     TLatex text;
    //     text.SetTextSize(0.05);

    //     const char *pdfFile = PDFFile.c_str();

    //     char fileName[100];
    //     sprintf(fileName, "%s[", pdfFile);
    //     myText->SaveAs(fileName);
    //     sprintf(fileName, "%s", pdfFile);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text.DrawLatex(0.2, 0.7, "Neutrons in CND");

    //     myText->Print(fileName, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     // myCanvas->SetGrid();
    //     myCanvas->Divide(4, 3);
    //     // myCanvas->SetGrid(), myCanvas->cd()->SetBottomMargin(0.14), myCanvas->cd()->SetLeftMargin(0.16), myCanvas->cd()->SetRightMargin(0.16), myCanvas->cd()->SetTopMargin(0.12);

    //     double x_1 = 0.2, y_1 = 0.3, x_2 = 0.86, y_2 = 0.7;
    //     double diplayTextSize = 0.1;

    //     // int canvas_ind = 1;

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         int canvas_ind = (i % 12) + 1; // Determine the pad number (1 to 12)

    //         myCanvas->cd(canvas_ind);
    //         gPad->SetGrid();

    //         HistoList[i]->SetLineWidth(2);
    //         HistoList[i]->SetLineColor(kBlue);

    //         if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //         {
    //             TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //             displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //             HistoList[i]->Draw(), displayText->Draw("same");
    //         }
    //         else
    //         {
    //             HistoList[i]->Draw();
    //         }

    //         // Save the canvas to a PDF page after filling 12 pads or processing the last histogram
    //         if (canvas_ind == 12 || i == HistoList.size() - 1)
    //         {
    //             myCanvas->Print(fileName); // Save the current page
    //             if (i != HistoList.size() - 1)
    //             {
    //                 myCanvas->Clear();      // Clear the canvas for the next page
    //                 myCanvas->Divide(4, 3); // Reset the grid layout
    //             }
    //         }
    //     }

    //     // for (int i = 0; i < HistoList.size(); i++)
    //     // {
    //     //     // myCanvas->cd(canvas_ind);
    //     //     // gPad->SetGrid();

    //     //     HistoList[i]->SetLineWidth(2);
    //     //     HistoList[i]->SetLineColor(kBlue);

    //     //     if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //     //     {
    //     //         TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //     //         displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //     //         HistoList[i]->Draw(), displayText->Draw("same");
    //     //     }
    //     //     else
    //     //     {
    //     //         HistoList[i]->Draw();
    //     //     }

    //     //     myCanvas->Print(fileName, "pdf");
    //     //     myCanvas->Clear();

    //     //     // ++canvas_ind;

    //     //     // if (i > 12 && 12 % i == 0)
    //     //     // {
    //     //     //     myCanvas->Print(fileName, "pdf");
    //     //     //     myCanvas->Clear();
    //     //     //     myCanvas->Divide(4, 3);

    //     //     //     canvas_ind = 1;
    //     //     // }
    //     // }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         myCanvas->cd(1);

    //         if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //         {
    //             TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //             displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //             HistoList[i]->Draw("colz"), displayText->Draw("same");
    //         }
    //         else
    //         {
    //             HistoList[i]->Draw("colz");
    //         }

    //         myCanvas->Print(fileName, "pdf");
    //         myCanvas->Clear();
    //     }

    //     sprintf(fileName, "%s]", pdfFile);
    //     myCanvas->Print(fileName, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving all plots - end */

    // #pragma region /* Saving only CD proton plots - start */

    //     TLatex text_CD;
    //     text_CD.SetTextSize(0.05);

    //     string pdfFile_CD_0 = ConfigOutPutName(PDFFile, "pCD_only").c_str();
    //     const char *pdfFile_CD = pdfFile_CD_0.c_str();

    //     char fileName_CD[100];
    //     sprintf(fileName_CD, "%s[", pdfFile_CD);
    //     myText->SaveAs(fileName_CD);
    //     sprintf(fileName_CD, "%s", pdfFile_CD);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_CD.DrawLatex(0.2, 0.9, "(e,e'pCD) Cuts:");
    //     text_CD.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_CD.DrawLatex(0.2, 0.7, "Neutrons in CND");

    //     myText->Print(fileName_CD, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();
    //     // myCanvas->Divide(4, 3);
    //     // myCanvas->SetGrid(), myCanvas->cd()->SetBottomMargin(0.14), myCanvas->cd()->SetLeftMargin(0.16), myCanvas->cd()->SetRightMargin(0.16), myCanvas->cd()->SetTopMargin(0.12);

    //     // double x_1 = 0.2, y_1 = 0.3, x_2 = 0.86, y_2 = 0.7;
    //     // double diplayTextSize = 0.1;

    //     // int canvas_ind = 1;

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "CD"))
    //         {
    //             // myCanvas->cd(canvas_ind);
    //             // gPad->SetGrid();

    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_CD, "pdf");
    //             myCanvas->Clear();

    //             // ++canvas_ind;

    //             // if (i > 12 && 12 % i == 0)
    //             // {
    //             //     myCanvas->Print(fileName_CD, "pdf");
    //             //     myCanvas->Clear();
    //             //     myCanvas->Divide(4, 3);

    //             //     canvas_ind = 1;
    //             // }
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "CD"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_CD, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_CD, "%s]", pdfFile_CD);
    //     myCanvas->Print(fileName_CD, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving only CD proton plots - end */

    // #pragma region /* Saving Step0 plots - start */

    //     TLatex text_Step0;
    //     text_Step0.SetTextSize(0.05);

    //     string pdfFile_Step0_0 = ConfigOutPutName(PDFFile, "Step0").c_str();
    //     const char *pdfFile_Step0 = pdfFile_Step0_0.c_str();

    //     cout << "\nPDFFile = " << PDFFile << "\n";
    //     cout << "\npdfFile_Step0_0 = " << pdfFile_Step0_0 << "\n";
    //     cout << "\npdfFile_Step0 = " << pdfFile_Step0 << "\n";

    //     char fileName_Step0[100];
    //     sprintf(fileName_Step0, "%s[", pdfFile_Step0);
    //     myText->SaveAs(fileName_Step0);
    //     sprintf(fileName_Step0, "%s", pdfFile_Step0);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step0.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step0.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step0.DrawLatex(0.2, 0.7, "Neutrons in CND - step 0");

    //     myText->Print(fileName_Step0, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step0"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step0, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step0"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step0, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step0, "%s]", pdfFile_Step0);
    //     myCanvas->Print(fileName_Step0, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step0 plots - end */

    // #pragma region /* Saving Step1 plots - start */

    //     TLatex text_Step1;
    //     text_Step1.SetTextSize(0.05);

    //     string pdfFile_Step1_0 = ConfigOutPutName(PDFFile, "Step1").c_str();
    //     const char *pdfFile_Step1 = pdfFile_Step1_0.c_str();

    //     char fileName_Step1[100];
    //     sprintf(fileName_Step1, "%s[", pdfFile_Step1);
    //     myText->SaveAs(fileName_Step1);
    //     sprintf(fileName_Step1, "%s", pdfFile_Step1);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step1.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step1.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step1.DrawLatex(0.2, 0.7, "Neutrons in CND - step 1");

    //     myText->Print(fileName_Step1, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step1"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step1, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step1"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step1, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step1, "%s]", pdfFile_Step1);
    //     myCanvas->Print(fileName_Step1, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step1 plots - end */

    // #pragma region /* Saving Step2 plots - start */

    //     TLatex text_Step2;
    //     text_Step2.SetTextSize(0.05);

    //     string pdfFile_Step2_0 = ConfigOutPutName(PDFFile, "Step2").c_str();
    //     const char *pdfFile_Step2 = pdfFile_Step2_0.c_str();

    //     char fileName_Step2[100];
    //     sprintf(fileName_Step2, "%s[", pdfFile_Step2);
    //     myText->SaveAs(fileName_Step2);
    //     sprintf(fileName_Step2, "%s", pdfFile_Step2);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step2.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step2.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step2.DrawLatex(0.2, 0.7, "Neutrons in CND - step 2");

    //     myText->Print(fileName_Step2, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step2"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step2, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step2"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step2, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step2, "%s]", pdfFile_Step2);
    //     myCanvas->Print(fileName_Step2, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step2 plots - end */

    // #pragma region /* Saving Step3 plots - start */

    //     TLatex text_Step3;
    //     text_Step3.SetTextSize(0.05);

    //     string pdfFile_Step3_0 = ConfigOutPutName(PDFFile, "Step3").c_str();
    //     const char *pdfFile_Step3 = pdfFile_Step3_0.c_str();

    //     char fileName_Step3[100];
    //     sprintf(fileName_Step3, "%s[", pdfFile_Step3);
    //     myText->SaveAs(fileName_Step3);
    //     sprintf(fileName_Step3, "%s", pdfFile_Step3);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step3.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step3.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step3.DrawLatex(0.2, 0.7, "Neutrons in CND - step 3");

    //     myText->Print(fileName_Step3, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step3"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step3, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step3"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step3, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step3, "%s]", pdfFile_Step3);
    //     myCanvas->Print(fileName_Step3, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step3 plots - end */

    // #pragma region /* Saving Step4 plots - start */

    //     TLatex text_Step4;
    //     text_Step4.SetTextSize(0.05);

    //     string pdfFile_Step4_0 = ConfigOutPutName(PDFFile, "Step4").c_str();
    //     const char *pdfFile_Step4 = pdfFile_Step4_0.c_str();

    //     char fileName_Step4[100];
    //     sprintf(fileName_Step4, "%s[", pdfFile_Step4);
    //     myText->SaveAs(fileName_Step4);
    //     sprintf(fileName_Step4, "%s", pdfFile_Step4);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step4.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step4.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step4.DrawLatex(0.2, 0.7, "Neutrons in CND - step 4");

    //     myText->Print(fileName_Step4, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step4"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step4, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step4"))
    //         {
    //             myCanvas->cd(1);
    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step4, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step4, "%s]", pdfFile_Step4);
    //     myCanvas->Print(fileName_Step4, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step4 plots - end */

    // #pragma region /* Saving Step5 plots - start */

    //     TLatex text_Step5;
    //     text_Step5.SetTextSize(0.05);

    //     string pdfFile_Step5_0 = ConfigOutPutName(PDFFile, "Step5").c_str();
    //     const char *pdfFile_Step5 = pdfFile_Step5_0.c_str();

    //     char fileName_Step5[100];
    //     sprintf(fileName_Step5, "%s[", pdfFile_Step5);
    //     myText->SaveAs(fileName_Step5);
    //     sprintf(fileName_Step5, "%s", pdfFile_Step5);

    //     /////////////////////////////////////
    //     // CND Neutron Information
    //     /////////////////////////////////////

    //     myText->cd();

    //     text_Step5.DrawLatex(0.2, 0.9, "(e,e'p) Cuts:");
    //     text_Step5.DrawLatex(0.2, 0.8, "(e,e') Cuts");
    //     text_Step5.DrawLatex(0.2, 0.7, "Neutrons in CND - step 5");

    //     myText->Print(fileName_Step5, "pdf");
    //     myText->Clear();

    //     myCanvas->cd();
    //     myCanvas->SetGrid();

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step5"))
    //         {
    //             myCanvas->cd(1);
    //             HistoList[i]->SetLineWidth(2);
    //             HistoList[i]->SetLineColor(kBlue);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 displayText->SetBorderSize(6);
    //                 HistoList[i]->Draw(), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw();
    //             }

    //             myCanvas->Print(fileName_Step5, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     for (int i = 0; i < HistoList.size(); i++)
    //     {
    //         string TempHistName = HistoList[i]->GetName();

    //         if (findSubstring(TempHistName, "Step5"))
    //         {
    //             myCanvas->cd(1);

    //             if (HistoList[i]->GetEntries() == 0 || HistoList[i]->Integral() == 0)
    //             {
    //                 TPaveText *displayText = new TPaveText(x_1, y_1, x_2, y_2, "NDC");
    //                 displayText->SetTextSize(diplayTextSize), displayText->SetFillColor(0), displayText->AddText("Empty histogram"), displayText->SetTextAlign(22);
    //                 HistoList[i]->Draw("colz"), displayText->Draw("same");
    //             }
    //             else
    //             {
    //                 HistoList[i]->Draw("colz");
    //             }

    //             myCanvas->Print(fileName_Step5, "pdf");
    //             myCanvas->Clear();
    //         }
    //     }

    //     sprintf(fileName_Step5, "%s]", pdfFile_Step5);
    //     myCanvas->Print(fileName_Step5, "pdf");

    //     myCanvas->Clear();
    //     myText->Clear();

    // #pragma endregion /* Saving Step5 plots - end */

#pragma endregion /* Andrew's wrap up - end */

    // ======================================================================================================================================================================
    // Save log file
    // ======================================================================================================================================================================

    // ======================================================================================================================================================================
    // Printouts
    // ======================================================================================================================================================================

#pragma region /* Printouts - start */

    cout << "\033[33m\n\033[0m";
    // cout << "\033[33mHIPO_FILES:\033[0m\t\t" << gSystem->Getenv("HIPO_FILES") << "\n";
    cout << "\033[33minput_hipo:\033[0m\t\t" << input_hipo << "\n";
    cout << "\033[33m\n\033[0m";
    cout << "\033[33mOUTDIR:\033[0m\t\t\t" << gSystem->Getenv("OUTDIR") << "\n";
    cout << "\033[33mOutDir:\033[0m\t\t\t" << OutDir << "\n";
    cout << "\033[33moutput_pdf_Erin:\033[0m\t" << output_pdf_Erin << "\n";
    cout << "\033[33moutput_root_Erin:\033[0m\t" << output_root_Erin << "\n";
    cout << "\033[33moutput_txt_Erin:\033[0m\t" << output_txt_Erin << "\n";
    cout << "\033[33mPDFFile:\033[0m\t\t" << PDFFile << "\n\n\n";

    /* Timing output */
    auto Code_end_time = std::chrono::system_clock::now();
    auto Elapsed_time_seconds = std::chrono::duration_cast<std::chrono::seconds>(Code_end_time - Code_start_time);
    double Elapsed_time_minutes = Elapsed_time_seconds.count() / 60;

    if (Elapsed_time_seconds.count() < 60)
    {
        std::cout << "\033[33mRunning time:\033[0m\t\t" << Elapsed_time_seconds.count() << " seconds\n\n";
    }
    else
    {
        std::cout << "\033[33mRunning time:\033[0m\t\t" << to_string_with_precision(Elapsed_time_minutes, 3) << " minutes\n\n";
    }

#pragma endregion /* Printouts - end */

    return 0;

} // closes main function

#pragma endregion /* Erin main function - end */