#!/bin/tcsh

# Set the output directory
unset OUTDIR
setenv OUTDIR Output_full_data_3_test2

# Set the hipo files to run over
unset HIPO_FILES
# setenv HIPO_FILES /cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/* # All D2 @ 6GeV data
setenv HIPO_FILES "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/*.hipo" # Only run 015045 of D2 @ 6GeV data
# setenv HIPO_FILES /cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo # Only run 015443 of D2 @ 6GeV data
# setenv HIPO_FILES /cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo # Only one file of run 015045 of D2 @ 6GeV data

clas12root -l -q 'D_getfeatures_Phase5.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", ${HIPO_FILES}, "${OUTDIR}/Andrew_plots.pdf")'