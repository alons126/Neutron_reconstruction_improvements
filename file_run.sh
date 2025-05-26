#!/bin/tcsh

# Set the output directory (/w/hallb-scshelf2102/clas12/asportes/Neutron_reconstruction_improvements/)
setenv OUTDIR Output_data_P11_4

 # Only run 015050 of D2 @ 6GeV data (1 file) ---------------------------------------
 clas12root -l -q 'ManualNeutronVeto.cpp("${OUTDIR}", "Erin_plots.pdf", 5.98636, true, "Erin_plots.root", "Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015050/rec_clas_015050.evio.00210-00214.hipo", "ManualVeto_plots.pdf")'
 # ----------------------------------------------------------------------------------

# # Only run 015449 of D2 @ 6GeV data (57 files) ------------------------------------
# clas12root -l -q 'ManualNeutronVeto.cpp("${OUTDIR}", "Erin_plots.pdf", 5.98636, true, "Erin_plots.root", "Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015449/*.hipo", "ManualVeto_plots.pdf")'
# # ----------------------------------------------------------------------------------

# # Only run 015050 of D2 @ 6GeV data (131 files) ------------------------------------
# clas12root -l -q 'ManualNeutronVeto.cpp("${OUTDIR}", "Erin_plots.pdf", 5.98636, true, "Erin_plots.root", "Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015050/*.hipo", "ManualVeto_plots.pdf")'
# # ----------------------------------------------------------------------------------

# # Only run 015443 of D2 @ 6GeV data (370 files) ------------------------------------
# clas12root -l -q 'ManualNeutronVeto.cpp("${OUTDIR}", "Erin_plots.pdf", 5.98636, true, "Erin_plots.root", "Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo", "ManualVeto_plots.pdf")'
# # ----------------------------------------------------------------------------------

# # Only run 015443 of D2 @ 6GeV data (1292 files) ------------------------------------
# clas12root -l -q 'ManualNeutronVeto.cpp("${OUTDIR}", "Erin_plots.pdf", 5.98636, true, "Erin_plots.root", "Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/333", "ManualVeto_plots.pdf")'
# # ----------------------------------------------------------------------------------

# # All D2 @ 6GeV data ---------------------------------------------------------------
# clas12root -l -q 'ManualNeutronVeto.cpp("${OUTDIR}", "Erin_plots.pdf", 5.98636, true, "Erin_plots.root", "Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "ManualVeto_plots.pdf")'
# # ----------------------------------------------------------------------------------
