#include <iostream>

#include "ManualNeutronVeto.cpp"

int main() {
    CodeDirectories codeDirectories;  // Get the directories

    std::string OUTDIR_prefix = "Output_data_P11_2";
    
    std::string OUTDIR = codeDirectories.plots_path + "/" + codeDirectories.plots_path_prefix + OUTDIR_prefix;

    // Only run 015050 of D2 @ 6GeV data (1 file) ---------------------------------------
    std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    double Ebeam = 5.98636;
    bool keep_good = true;
    std::string Erin_plots_root = OUTDIR + "/Erin_plots.root";
    std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015050/rec_clas_015050.evio.00210-00214.hipo";
    std::string ManualVeto_plots_pdf = OUTDIR + "/ManualVeto_plots.pdf";
    // // ----------------------------------------------------------------------------------

    // // Only run 015449 of D2 @ 6GeV data (57 files) ------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_root = OUTDIR + "/Erin_plots.root";
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015449/*.hipo";
    // std::string ManualVeto_plots_pdf = OUTDIR + "/ManualVeto_plots.pdf";
    // // ----------------------------------------------------------------------------------

    // // Only run 015050 of D2 @ 6GeV data (127 files) ------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_root = OUTDIR + "/Erin_plots.root";
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015050/*.hipo";
    // std::string ManualVeto_plots_pdf = OUTDIR + "/ManualVeto_plots.pdf";
    // // ----------------------------------------------------------------------------------

    // // Only run 015443 of D2 @ 6GeV data (370 files) ------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_root = OUTDIR + "/Erin_plots.root";
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo";
    // std::string ManualVeto_plots_pdf = OUTDIR + "/ManualVeto_plots.pdf";
    // // ----------------------------------------------------------------------------------

    // // Only run 015443 of D2 @ 6GeV data (1292 files) ------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_root = OUTDIR + "/Erin_plots.root";
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/333";
    // std::string ManualVeto_plots_pdf = OUTDIR + "/ManualVeto_plots.pdf";
    // // ----------------------------------------------------------------------------------

    // // All D2 @ 6GeV data ---------------------------------------------------------------
    // std::string Erin_plots_pdf = OUTDIR + "/Erin_plots.pdf";
    // double Ebeam = 5.98636;
    // bool keep_good = true;
    // std::string Erin_plots_root = OUTDIR + "/Erin_plots.root";
    // std::string Erin_plots_txt = OUTDIR + "/Erin_plots.txt";
    // std::string Data_dir = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*";
    // std::string ManualVeto_plots_pdf = OUTDIR + "/ManualVeto_plots.pdf";
    // // ----------------------------------------------------------------------------------

    ManualNeutronVeto(OUTDIR, Erin_plots_pdf, Ebeam, keep_good, Erin_plots_root, Erin_plots_txt, Data_dir, ManualVeto_plots_pdf);

    return 0;
}
