#!/bin/bash

OUTDIR=Output_trtrtrtr

# # All D2 @ 6GeV data:
# clas12root -l -q 'D_getfeatures_Phase4.cpp($OUTDIR, 5.98636, true, "$OUTDIR/Erin_plots.root", "$OUTDIR/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*", "$OUTDIR/Andrew_plots.pdf")'

## Only run 015045 of D2 @ 6GeV data:
# clas12root -l -q 'D_getfeatures_Phase4.cpp($OUTDIR, 5.98636, true, "$OUTDIR/Erin_plots.root", "$OUTDIR/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/*.hipo", "$OUTDIR/Andrew_plots.pdf")'

# Only one file of run 015045 of D2 @ 6GeV data:
clas12root -l -q 'D_getfeatures_Phase4.cpp('"$OUTDIR"', 5.98636, true, "'"$OUTDIR"'/Erin_plots.root", "'"$OUTDIR"'/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo", "'"$OUTDIR"'/Andrew_plots.pdf")'