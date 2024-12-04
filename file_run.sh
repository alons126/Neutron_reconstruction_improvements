#!/bin/tcsh

# Set the output directory
unset OUTDIR
setenv OUTDIR Output_full_data_3_mmiss_test_OnlydbetaCut

# # All D2 @ 6GeV data ---------------------------------------------------------------
# clas12root -l -q 'D_getfeatures_Phase5.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "${OUTDIR}/Andrew_plots.pdf")'
# # ----------------------------------------------------------------------------------

# Only run 015045 of D2 @ 6GeV data ------------------------------------------------
clas12root -l -q 'D_getfeatures_Phase5.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/*.hipo", "${OUTDIR}/Andrew_plots.pdf")'
# ----------------------------------------------------------------------------------

# # Only run 015443 of D2 @ 6GeV data ------------------------------------------------
# clas12root -l -q 'D_getfeatures_Phase5.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo", "${OUTDIR}/Andrew_plots.pdf")'
# # ----------------------------------------------------------------------------------

# # Only one file of run 015045 of D2 @ 6GeV data ------------------------------------
# clas12root -l -q 'D_getfeatures_Phase5.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo", "${OUTDIR}/Andrew_plots.pdf")'
# # ----------------------------------------------------------------------------------