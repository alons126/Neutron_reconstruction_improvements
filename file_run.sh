#!/bin/bash

# OUTDIR="Output_full_data_2_ErinOnly"

# All D2 @ 6GeV data:
clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_2_ErinOnly", "Output_full_data_2_ErinOnly/Erin_plots.pdf", 5.98636, true, "Output_full_data_2_ErinOnly/Erin_plots.root", "Output_full_data_2_ErinOnly/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "Output_full_data_2_ErinOnly/Andrew_plots.pdf")'

## Only run 015045 of D2 @ 6GeV data:
# clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_2_ErinOnly", "Output_full_data_2_ErinOnly/Erin_plots.pdf", 5.98636, true, "Output_full_data_2_ErinOnly/Erin_plots.root", "Output_full_data_2_ErinOnly/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/*.hipo", "Output_full_data_2_ErinOnly/Andrew_plots.pdf")'

# # Only run 015443 of D2 @ 6GeV data:
# clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_2_ErinOnly", "Output_full_data_2_ErinOnly/Erin_plots.pdf", 5.98636, true, "Output_full_data_2_ErinOnly/Erin_plots.root", "Output_full_data_2_ErinOnly/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo", "Output_full_data_2_ErinOnly/Andrew_plots.pdf")'

# # Only one file of run 015045 of D2 @ 6GeV data:
# clas12root -l -q 'D_getfeatures_Phase5.cpp("Output_full_data_2_ErinOnly", "Output_full_data_2_ErinOnly/Erin_plots.pdf", 5.98636, true, "Output_full_data_2_ErinOnly/Erin_plots.root", "Output_full_data_2_ErinOnly/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo", "Output_full_data_2_ErinOnly/Andrew_plots.pdf")'
