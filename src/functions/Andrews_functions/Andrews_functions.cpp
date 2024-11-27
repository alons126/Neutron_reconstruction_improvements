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

// ==========================================================================================================================================================================
// Andrew's functions
// ==========================================================================================================================================================================

void printProgress(double percentage);

bool isPosNear(int sdiff, int ldiff)
{
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

    // TODO: ask Andrew for motivations for these sdiff cuts
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

void printProgress(double percentage)
{
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
