#!/bin/tcsh

# Set the output directory
setenv OUTDIR Output_full_data_3_t

# All D2 @ 6GeV data:
clas12root -l -q 'D_getfeatures_Phase5.cpp("${OUTDIR}", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "Output_full_data_3_t/Andrew_plots.pdf")'
# clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_3_t", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "Output_full_data_3_t/Andrew_plots.pdf")'

## Only run 015045 of D2 @ 6GeV data:
# clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_3_t", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/*.hipo", "Output_full_data_3_t/Andrew_plots.pdf")'

# # Only run 015443 of D2 @ 6GeV data:
# clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_3_t", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo", "Output_full_data_3_t/Andrew_plots.pdf")'

# # Only one file of run 015045 of D2 @ 6GeV data:
# clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_3_t", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo", "Output_full_data_3_t/Andrew_plots.pdf")'




###!/bin/bash

# # OUTDIR="Output_full_data_3_t"
# export OUTDIR="Output_full_data_3_t"

# # All D2 @ 6GeV data:
# clas12root -l -q 'D_getfeatures_Phase5.cpp("${OUTDIR}", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "Output_full_data_3_t/Andrew_plots.pdf")'
# # clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_3_t", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "Output_full_data_3_t/Andrew_plots.pdf")'

# ## Only run 015045 of D2 @ 6GeV data:
# # clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_3_t", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/*.hipo", "Output_full_data_3_t/Andrew_plots.pdf")'

# # # Only run 015443 of D2 @ 6GeV data:
# # clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_3_t", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo", "Output_full_data_3_t/Andrew_plots.pdf")'

# # # Only one file of run 015045 of D2 @ 6GeV data:
# # clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_3_t", "Output_full_data_3_t/Erin_plots.pdf", 5.98636, true, "Output_full_data_3_t/Erin_plots.root", "Output_full_data_3_t/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo", "Output_full_data_3_t/Andrew_plots.pdf")'
