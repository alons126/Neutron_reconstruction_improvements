#!/bin/bash

echo ""
git pull
echo ""
echo "- Re-pulling repository -----------------------------------------------"
echo ""
git reset --hard
git pull
echo ""
echo "- Lunching manual neutron veto code -----------------------------------"
echo ""
source file_run.sh
echo ""
echo "- Operation finished --------------------------------------------------"