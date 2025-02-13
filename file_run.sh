#!/bin/tcsh

# Set the output directory (/w/hallb-scshelf2102/clas12/asportes/Neutron_reconstruction_improvements/)
setenv OUTDIR Output_data_P9_run9_131_newBetaAndToF_NoCTOF
# setenv OUTDIR Output_data_P9_run9_131_newBetaAndToF
# setenv OUTDIR Output_data_P9_run9_57_newBetaAndToF_NoCTOF
# setenv OUTDIR Output_data_P9_run9_57_newBetaAndToF
# setenv OUTDIR Output_data_P9_run9_full_REDO_new_LargeMmissCuts_2
# setenv OUTDIR Output_data_P9_run9_full_REDO_new_2

#  # Only run 015045 of D2 @ 6GeV data (1 file) ---------------------------------------
#  clas12root -l -q 'ManualVeto_Phase9.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo", "${OUTDIR}/ManualVeto_plots.pdf")'
#  # ----------------------------------------------------------------------------------

# # Only run 015449 of D2 @ 6GeV data (57 files) ------------------------------------
# clas12root -l -q 'ManualVeto_Phase9.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015449/*.hipo", "${OUTDIR}/ManualVeto_plots.pdf")'
# # ----------------------------------------------------------------------------------

# Only run 015045 of D2 @ 6GeV data (131 files) ------------------------------------
clas12root -l -q 'ManualVeto_Phase9.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/*.hipo", "${OUTDIR}/ManualVeto_plots.pdf")'
# ----------------------------------------------------------------------------------

# # Only run 015443 of D2 @ 6GeV data (370 files) ------------------------------------
# clas12root -l -q 'ManualVeto_Phase9.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015443/*.hipo", "${OUTDIR}/ManualVeto_plots.pdf")'
# # ----------------------------------------------------------------------------------

# # Only run 015443 of D2 @ 6GeV data (1292 files) ------------------------------------
# clas12root -l -q 'ManualVeto_Phase9.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/333", "${OUTDIR}/ManualVeto_plots.pdf")'
# # ----------------------------------------------------------------------------------

# # All D2 @ 6GeV data ---------------------------------------------------------------
# clas12root -l -q 'ManualVeto_Phase9.cpp("${OUTDIR}", "${OUTDIR}/Erin_plots.pdf", 5.98636, true, "${OUTDIR}/Erin_plots.root", "${OUTDIR}/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "${OUTDIR}/ManualVeto_plots.pdf")'
# # ----------------------------------------------------------------------------------


# setenv BaseName Output_data
# setenv Run run7
# setenv Phase P9
# setenv Ending P9