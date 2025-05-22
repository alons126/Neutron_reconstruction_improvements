#include <iostream>

#include "ManualNeutronVeto.cpp"

int main() {
    CodeDirectories codeDirectories;  // Get the directories

    std::string OUTDIR = codeDirectories.plots_path + "/MnCDV_" + "Output_data_P11";

    // Only run 015050 of D2 @ 6GeV data (1 file) ---------------------------------------
    std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    double Ebeam = 5.98636;
    bool keep_good = true;
    std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015050/rec_clas_015050.evio.00210-00214.hipo";
    // ----------------------------------------------------------------------------------

    // // Only run 015449 of D2 @ 6GeV data (57 files) ------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015449/*.hipo";
    // // ----------------------------------------------------------------------------------

    // // Only run 015050 of D2 @ 6GeV data (127 files) ------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015050/*.hipo";
    // // ----------------------------------------------------------------------------------

    // // Only run 015443 of D2 @ 6GeV data (370 files) ------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo";
    // // ----------------------------------------------------------------------------------

    // // Only run 015443 of D2 @ 6GeV data (1292 files) ------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/333";
    // // ----------------------------------------------------------------------------------

    // // All D2 @ 6GeV data ---------------------------------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*";
    // // ----------------------------------------------------------------------------------

    ManualNeutronVeto();
    
    return 0;
}
