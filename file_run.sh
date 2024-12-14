#!/bin/tcsh

# Set the output directory
unset OUTDIR
setenv OUTDIR Output_data_11_Andrew_131

# # All D2 @ 6GeV data ---------------------------------------------------------------
# clas12root -l -q 'D_getfeatures_Phase6.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "${OUTDIR}/Andrew_plots.pdf")'
# # ----------------------------------------------------------------------------------

# # Only run 015449 of D2 @ 6GeV data (57 files) ------------------------------------
# clas12root -l -q 'D_getfeatures_Phase6.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015449/*.hipo", "${OUTDIR}/Andrew_plots.pdf")'
# # ----------------------------------------------------------------------------------

# Only run 015045 of D2 @ 6GeV data (131 files) ------------------------------------
clas12root -l -q 'D_getfeatures_Phase6.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/*.hipo", "${OUTDIR}/Andrew_plots.pdf")'
# ----------------------------------------------------------------------------------

# # Only run 015443 of D2 @ 6GeV data (370 files) ------------------------------------
# clas12root -l -q 'D_getfeatures_Phase6.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo", "${OUTDIR}/Andrew_plots.pdf")'
# # ----------------------------------------------------------------------------------

# # Only one file of run 015045 of D2 @ 6GeV data ------------------------------------
# clas12root -l -q 'D_getfeatures_Phase6.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo", "${OUTDIR}/Andrew_plots.pdf")'
# # ----------------------------------------------------------------------------------