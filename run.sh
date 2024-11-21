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
clas12root -l -q 'D_getfeatures_Phase1.cpp(5.98636, true, "test.root", "test.txt", "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/015045/rec_clas_015045.evio.00650-00654.hipo")'
echo ""
echo "- Operation finished --------------------------------------------------"