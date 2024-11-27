#!/bin/bash
echo ""
git pull
echo ""
echo "- Re-pulling repository -----------------------------------------------"
echo ""
git reset --hard
git pull
echo ""
echo "- Lunching Erin's D_getfeatures code ----------------------------------"
echo ""
# clas12root -l -q 'D_getfeatures_Phase4.cpp(5.98636, true, "Output/Erin_plots.root", "Output/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/*.hipo", "Output/Andrew_plots.pdf")'
clas12root -l -q 'D_getfeatures_Phase4.cpp(5.98636, true, "Output/Erin_plots.root", "Output/Erin_plots.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo", "Output/Andrew_plots.pdf")'
echo ""
echo "- Operation finished --------------------------------------------------"