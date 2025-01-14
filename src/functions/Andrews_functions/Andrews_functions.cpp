#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "clas12reader.h"
#include "HipoChain.h"
// #include "eventcut.h"
// #include "functions.h"

using namespace std;
using namespace clas12;

/*
Possible ldiff:
n=1,2,3; c=0,1,2,3; ldiff=n-c
ldiff == -3: (n,c)=(0,3)x
ldiff == -2: (n,c)=(0,2)x, (1,3)
ldiff == -1: (n,c)=(0,1)x, (1,2), (2,3)
ldiff ==  0: (n,c)=(1,1), (2,2), (3,3)
ldiff ==  1: (n,c)=(1,0), (2,1), (3,2)
ldiff ==  2: (n,c)=(2,0), (3,1)
ldiff ==  3: (n,c)=(3,0)
*/

// ==========================================================================================================================================================================
// Andrew's functions
// ==========================================================================================================================================================================

// printProgress function ---------------------------------------------------------------------------------------------------------------------------------------------------

void printProgress(double percentage);

// isPosNear function -------------------------------------------------------------------------------------------------------------------------------------------------------

bool isPosNear(int sdiff, int ldiff)
{
    if ((ldiff == -2) && (sdiff >= -1) && (sdiff <= 0))
    {
        return true;
    }

    if ((ldiff == -1) && (sdiff >= -1) && (sdiff <= 2))
    {
        return true;
    }

    if ((ldiff == 0) && (sdiff >= -1) && (sdiff <= 2))
    {
        return true;
    }

    if ((ldiff == 1) && (sdiff >= -1) && (sdiff <= 2))
    {
        return true;
    }

    if ((ldiff == 2) && (sdiff >= -1) && (sdiff <= 2))
    {
        return true;
    }

    if ((ldiff == 3) && (sdiff >= -1) && (sdiff <= 2))
    {
        return true;
    }

    return false;
}

// isPosNear_PhiCut function ------------------------------------------------------------------------------------------------------------------------------------------------

bool isPosNear_PhiCut(int sdiff, int ldiff, double Phi_n)
{

    if (ldiff == -2)
    {
        bool Phi_Range = ((30. <= Phi_n) && (Phi_n <= 60.));

        if (Phi_Range && ((sdiff >= 1) && (sdiff <= 2)))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    if (ldiff == -1)
    {
        bool Phi_Range = ((30. <= Phi_n) && (Phi_n <= 60.));

        if (Phi_Range && ((sdiff >= 1) && (sdiff <= 2)))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    if (ldiff == 0)
    {
        bool Phi_Range = ((40. <= Phi_n) && (Phi_n <= 60.));

        if (Phi_Range && (sdiff == 2))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    if (ldiff == 1)
    {
        if (((30. <= Phi_n) && (Phi_n <= 50.)) && ((sdiff >= 1) && (sdiff <= 2)))
        {
            return true;
        }
        else if (((120. <= Phi_n) && (Phi_n <= 140.)) && ((sdiff >= -1) && (sdiff <= 0)))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    if (ldiff == 2)
    {
        if (((30. <= Phi_n) && (Phi_n <= 50.)) && ((sdiff >= 1) && (sdiff <= 2)))
        {
            return true;
        }
        else if (((120. <= Phi_n) && (Phi_n <= 140.)) && ((sdiff >= -1) && (sdiff <= 0)))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    if (ldiff == 3)
    {
        if (((30. <= Phi_n) && (Phi_n <= 50.)) && ((sdiff >= 1) && (sdiff <= 2)))
        {
            return true;
        }
        else if (((120. <= Phi_n) && (Phi_n <= 140.)) && ((sdiff >= -1) && (sdiff <= 0)))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    return false;
}

// isPosNear_dToF function --------------------------------------------------------------------------------------------------------------------------------------------------

bool isPosNear_dToF(int sdiff, int ldiff, double dToF)
{

    if (ldiff == -2)
    {
        return false;
    }

    if (ldiff == -1)
    {
        if (((dToF >= 0) && (dToF <= 2)) && (abs(sdiff) <= 2))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    if (ldiff == 0)
    {
        if (((dToF >= 0) && (dToF <= 2)) && (abs(sdiff) <= 2))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    if (ldiff == 1)
    {
        if (((dToF >= 0) && (dToF <= 1.25)) && (abs(sdiff) <= 1))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    if (ldiff == 2)
    {
        if (((dToF >= 0) && (dToF <= 1.25)) && (abs(sdiff) <= 1))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    if (ldiff == 3)
    {
        if (((dToF >= 0) && (dToF <= 1.25)) && (abs(sdiff) <= 1))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    return false;
}

// isNear function ----------------------------------------------------------------------------------------------------------------------------------------------------------

bool isNear(int sdiff, int ldiff)
{
    /*
    //if((ldiff== 2) && (sdiff==-2)){return true;}
    //if((ldiff== 2) && (sdiff==-1)){return true;}
    if((ldiff== 2) && (sdiff== 0)){return true;}
    if((ldiff== 2) && (sdiff== 1)){return true;}
    if((ldiff== 2) && (sdiff== 2)){return true;}

    //if((ldiff== 1) && (sdiff== 1)){return true;}
    if((ldiff== 1) && (sdiff== 2)){return true;}

    //if((ldiff== 0) && (sdiff== 1)){return true;}
    if((ldiff== 0) && (sdiff== 2)){return true;}

    if((ldiff==-1) && (sdiff== -1)){return true;}
    */

    if ((ldiff == -2) && (sdiff == -2))
    {
        return true;
    }
    if ((ldiff == -2) && (sdiff == -1))
    {
        return true;
    }
    if ((ldiff == -2) && (sdiff == 0))
    {
        return true;
    }
    if ((ldiff == -2) && (sdiff == 1))
    {
        return true;
    }
    if ((ldiff == -2) && (sdiff == 2))
    {
        return true;
    }

    if ((ldiff == -1) && (sdiff == -2))
    {
        return true;
    }
    if ((ldiff == -1) && (sdiff == -1))
    {
        return true;
    }
    // if((ldiff==-1) && (sdiff== 0)){return true;}
    if ((ldiff == -1) && (sdiff == 1))
    {
        return true;
    }
    if ((ldiff == -1) && (sdiff == 2))
    {
        return true;
    }

    if ((ldiff == 0) && (sdiff == -2))
    {
        return true;
    }
    // if((ldiff== 0) && (sdiff==-1)){return true;}
    // if((ldiff== 0) && (sdiff== 0)){return true;}
    // if((ldiff== 0) && (sdiff== 1)){return true;}
    if ((ldiff == 0) && (sdiff == 2))
    {
        return true;
    }

    if ((ldiff == 1) && (sdiff == -2))
    {
        return true;
    }
    if ((ldiff == 1) && (sdiff == -1))
    {
        return true;
    }
    // if((ldiff== 1) && (sdiff== 0)){return true;}
    if ((ldiff == 1) && (sdiff == 1))
    {
        return true;
    }
    if ((ldiff == 1) && (sdiff == 2))
    {
        return true;
    }

    if ((ldiff == 2) && (sdiff == -2))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == -1))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == 0))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == 1))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == 2))
    {
        return true;
    }

    /*
    //if((ldiff==-1) && (sdiff==-2)){return true;}
    if((ldiff==-1) && (sdiff==-1)){return true;}
    if((ldiff==-1) && (sdiff== 0)){return true;}
    if((ldiff==-1) && (sdiff== 1)){return true;}
    //if((ldiff==-1) && (sdiff== 2)){return true;}

    //if((ldiff== 0) && (sdiff==-4)){return true;}
    //if((ldiff== 0) && (sdiff==-3)){return true;}
    //if((ldiff== 0) && (sdiff==-2)){return true;}
    if((ldiff== 0) && (sdiff==-1)){return true;}
    if((ldiff== 0) && (sdiff== 1)){return true;}
    //if((ldiff== 0) && (sdiff== 2)){return true;}
    //if((ldiff== 0) && (sdiff== 3)){return true;}

    //if((ldiff== 1) && (sdiff==-2)){return true;}
    if((ldiff== 1) && (sdiff==-1)){return true;}
    if((ldiff== 1) && (sdiff== 0)){return true;}
    if((ldiff== 1) && (sdiff== 1)){return true;}
    //if((ldiff== 1) && (sdiff== 2)){return true;}
    */
    return false;
}

// isNearCTOF function ------------------------------------------------------------------------------------------------------------------------------------------------------

bool isNearCTOF(int sdiff, int ldiff)
{
    if ((ldiff == 1) && (sdiff == -3))
    {
        return true;
    }
    if ((ldiff == 1) && (sdiff == -2))
    {
        return true;
    }
    if ((ldiff == 1) && (sdiff == -1))
    {
        return true;
    }
    if ((ldiff == 1) && (sdiff == 1))
    {
        return true;
    }
    if ((ldiff == 1) && (sdiff == 2))
    {
        return true;
    }
    if ((ldiff == 1) && (sdiff == 3))
    {
        return true;
    }

    if ((ldiff == 2) && (sdiff == -3))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == -2))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == -1))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == 0))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == 1))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == 2))
    {
        return true;
    }
    if ((ldiff == 2) && (sdiff == 3))
    {
        return true;
    }

    if ((ldiff == 3) && (sdiff == -3))
    {
        return true;
    }
    if ((ldiff == 3) && (sdiff == -2))
    {
        return true;
    }
    if ((ldiff == 3) && (sdiff == -1))
    {
        return true;
    }
    if ((ldiff == 3) && (sdiff == 0))
    {
        return true;
    }
    if ((ldiff == 3) && (sdiff == 1))
    {
        return true;
    }
    if ((ldiff == 3) && (sdiff == 2))
    {
        return true;
    }
    if ((ldiff == 3) && (sdiff == 3))
    {
        return true;
    }

    return false;
}

// printProgress function ---------------------------------------------------------------------------------------------------------------------------------------------------

void printProgress(double percentage)
{
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
