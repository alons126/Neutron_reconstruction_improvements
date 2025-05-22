//
// Created by Alon Sportes on 15/03/2025.
//

#ifndef NEUTRON_VETO_FUNCTIONS_H
#define NEUTRON_VETO_FUNCTIONS_H

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <TLorentzVector.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Include libraries:
#include "../basic_tools.h"
#include "../constants.h"

// Include CLAS12 libraries:
#include "../../../includes/clas12_include.h"

namespace neutron_veto_functions {

// SetNeutronCounters function ------------------------------------------------------------------------------------------------------------------------------------------------

void SetNeutronCounters(const bool isGN, int &counter_n_multiplicity_allN, int &counter_n_multiplicity_goodN, int &counter_n_multiplicity_badN) {
    ++counter_n_multiplicity_allN;
    isGN ? ++counter_n_multiplicity_goodN : ++counter_n_multiplicity_badN;
}

// SetNeutronCounters function ------------------------------------------------------------------------------------------------------------------------------------------------

void SetNeutronCounters(const bool pInCD, const bool pInFD, const bool isGN, int &counter_n_multiplicity_allN_epCD, int &counter_n_multiplicity_goodN_epCD, int &counter_n_multiplicity_badN_epCD,
                        int &counter_n_multiplicity_allN_epFD, int &counter_n_multiplicity_goodN_epFD, int &counter_n_multiplicity_badN_epFD) {
    if (pInCD || pInFD) {
        SetNeutronCounters(isGN, pInCD ? counter_n_multiplicity_allN_epCD : counter_n_multiplicity_allN_epFD, pInCD ? counter_n_multiplicity_goodN_epCD : counter_n_multiplicity_goodN_epFD,
                           pInCD ? counter_n_multiplicity_badN_epCD : counter_n_multiplicity_badN_epFD);
    }
}

// GetVzHitLocation function ------------------------------------------------------------------------------------------------------------------------------------------------

TVector3 GetVzHitLocation(region_part_ptr Electron) {
    TVector3 V_nvtx_3v;  // Neutron's vertex location -> set as the electron vertex

    V_nvtx_3v.SetXYZ(Electron->par()->getVx(), Electron->par()->getVy(), Electron->par()->getVz());

    return V_nvtx_3v;
}

// GetVzHitInCND function ------------------------------------------------------------------------------------------------------------------------------------------------

TVector3 GetVzHitInCND(region_part_ptr NeutronCD) {
    TVector3 V_hit_3v;  // Neutron's hit location in CND

    // Andrew's response checks:
    bool nCD_C1_hit = (NeutronCD->sci(CND1)->getDetector() == 3);
    bool nCD_C2_hit = (NeutronCD->sci(CND2)->getDetector() == 3);
    bool nCD_C3_hit = (NeutronCD->sci(CND3)->getDetector() == 3);
    auto nCD_detlayer = nCD_C1_hit ? CND1 : nCD_C2_hit ? CND2 : CND3;  // CND layer with hit

    V_hit_3v.SetXYZ(NeutronCD->sci(nCD_detlayer)->getX(), NeutronCD->sci(nCD_detlayer)->getY(), NeutronCD->sci(nCD_detlayer)->getZ());

    return V_hit_3v;
}

// GetnCDPath function ------------------------------------------------------------------------------------------------------------------------------------------------

TVector3 GetnCDPath(region_part_ptr NeutronCD, region_part_ptr Electron) {
    TVector3 V_nvtx_3v = GetVzHitLocation(Electron);  // Neutron's vertex location -> set as the electron vertex
    TVector3 V_hit_3v = GetVzHitInCND(NeutronCD);     // Neutron's hit location in CND
    TVector3 V_path_3v = V_hit_3v - V_nvtx_3v;        // Direct calculation of neutron's path (in vector form)

    return V_path_3v;
}

// GetnCDToF function ------------------------------------------------------------------------------------------------------------------------------------------------

double GetnCDToF(region_part_ptr NeutronCD, double starttime) {
    // Andrew's response checks:
    bool nCD_C1_hit = (NeutronCD->sci(CND1)->getDetector() == 3);
    bool nCD_C2_hit = (NeutronCD->sci(CND2)->getDetector() == 3);
    bool nCD_C3_hit = (NeutronCD->sci(CND3)->getDetector() == 3);
    auto nCD_detlayer = nCD_C1_hit ? CND1 : nCD_C2_hit ? CND2 : CND3;  // CND layer with hit

    return NeutronCD->sci(nCD_detlayer)->getTime() - starttime;
}

// GetnCDBeta function ------------------------------------------------------------------------------------------------------------------------------------------------

double GetnCDBeta(region_part_ptr NeutronCD, region_part_ptr Electron, double starttime) {
    double Beta_n;

    TVector3 V_path_3v = GetnCDPath(NeutronCD, Electron);  // Direct calculation of neutron's path (in vector form)
    double ToF_n = GetnCDToF(NeutronCD, starttime);

    Beta_n = V_path_3v.Mag() / (ToF_n * constants::c);

    return Beta_n;
}

//////

void HipoChain_config(HipoChain &chain, const string &AnalyseFilePath) {
    if (AnalyseFilePath == "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/*") {
        const bool PrintOut = true;

        string D2_6GeV_Data_Path = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/";

        /* Data in cache/clas12/rg-m/production/pass1/6gev/D/dst/recon */
        vector<string> Runs = {"015045", "015052", "015058", "015066", "015077", "015094", "015100", "015106", "015442", "015449", "015456", "015046", "015053", "015059", "015067", "015078",
                               "015095", "015101", "015435", "015443", "015450", "015047", "015054", "015060", "015072", "015079", "015096", "015102", "015436", "015444", "015451", "015049",
                               "015055", "015061", "015073", "015081", "015097", "015103", "015437", "015445", "015452", "015050", "015056", "015062", "015074", "015082", "015098", "015104",
                               "015439", "015447", "015454", "015051", "015057", "015065", "015075", "015093", "015099", "015105", "015441", "015448", "015455"};

        for (int i = 0; i < Runs.size(); i++) {
            string TempAnalyseFile = D2_6GeV_Data_Path + Runs.at(i) + "/*.hipo";
            chain.Add(TempAnalyseFile.c_str());

            if (PrintOut) { cout << TempAnalyseFile << " directory added to HipoChain!\n"; }
        }

        if (PrintOut) { cout << "\n"; }
    } else if (AnalyseFilePath == "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/333") {
        const bool PrintOut = true;

        string D2_6GeV_Data_Path = "/cache/clas12/rg-m/production/pass1/6gev/D/dst/recon/";

        /* Data in runs with at least 300 HIPO files */
        vector<string> Runs = {"015443", "015437", "015448", "015455"};

        for (int i = 0; i < Runs.size(); i++) {
            string TempAnalyseFile = D2_6GeV_Data_Path + Runs.at(i) + "/*.hipo";
            chain.Add(TempAnalyseFile.c_str());

            if (PrintOut) { cout << TempAnalyseFile << " directory added to HipoChain!\n"; }
        }

        if (PrintOut) { cout << "\n"; }
    } else {
        chain.Add(AnalyseFilePath);
    }
}

//////

// ConfigOutPutName function --------------------------------------------------------------------------------------------------------------------------------------------------

std::string ConfigOutPutName(const std::string &original, const std::string &toInsert) {
    size_t pos = original.find(".pdf");
    if (pos != std::string::npos) { return original.substr(0, pos) + "_" + toInsert + original.substr(pos); }
    // If ".pdf" is not found, return the original string
    return original;
}

//////

// printProgress function ---------------------------------------------------------------------------------------------------------------------------------------------------

void printProgress(double percentage);

// isPosNear function -------------------------------------------------------------------------------------------------------------------------------------------------------

bool isPosNear(int sdiff, int ldiff) {
    if ((ldiff == -2) && (sdiff >= -1) && (sdiff <= 0)) { return true; }

    if ((ldiff == -1) && (sdiff >= -1) && (sdiff <= 2)) { return true; }

    if ((ldiff == 0) && (sdiff >= -1) && (sdiff <= 2)) { return true; }

    if ((ldiff == 1) && (sdiff >= -1) && (sdiff <= 2)) { return true; }

    if ((ldiff == 2) && (sdiff >= -1) && (sdiff <= 2)) { return true; }

    if ((ldiff == 3) && (sdiff >= -1) && (sdiff <= 2)) { return true; }

    return false;
}

// isPosNear_PhiCut function ------------------------------------------------------------------------------------------------------------------------------------------------

bool isPosNear_PhiCut(int sdiff, int ldiff, double Phi_n) {
    if (ldiff == -2) {
        bool Phi_Range = ((30. <= Phi_n) && (Phi_n <= 60.));

        if (Phi_Range && ((sdiff >= 1) && (sdiff <= 2))) {
            return true;
        } else {
            return false;
        }
    }

    if (ldiff == -1) {
        bool Phi_Range = ((30. <= Phi_n) && (Phi_n <= 60.));

        if (Phi_Range && ((sdiff >= 1) && (sdiff <= 2))) {
            return true;
        } else {
            return false;
        }
    }

    if (ldiff == 0) {
        bool Phi_Range = ((40. <= Phi_n) && (Phi_n <= 60.));

        if (Phi_Range && (sdiff == 2)) {
            return true;
        } else {
            return false;
        }
    }

    if (ldiff == 1) {
        if (((30. <= Phi_n) && (Phi_n <= 50.)) && ((sdiff >= 1) && (sdiff <= 2))) {
            return true;
        } else if (((120. <= Phi_n) && (Phi_n <= 140.)) && ((sdiff >= -1) && (sdiff <= 0))) {
            return true;
        } else {
            return false;
        }
    }

    if (ldiff == 2) {
        if (((30. <= Phi_n) && (Phi_n <= 50.)) && ((sdiff >= 1) && (sdiff <= 2))) {
            return true;
        } else if (((120. <= Phi_n) && (Phi_n <= 140.)) && ((sdiff >= -1) && (sdiff <= 0))) {
            return true;
        } else {
            return false;
        }
    }

    if (ldiff == 3) {
        if (((30. <= Phi_n) && (Phi_n <= 50.)) && ((sdiff >= 1) && (sdiff <= 2))) {
            return true;
        } else if (((120. <= Phi_n) && (Phi_n <= 140.)) && ((sdiff >= -1) && (sdiff <= 0))) {
            return true;
        } else {
            return false;
        }
    }

    return false;
}

// isPosNear_dToF function --------------------------------------------------------------------------------------------------------------------------------------------------

bool isPosNear_dToF(int sdiff, int ldiff, double dToF) {
    if (ldiff == -2) { return false; }

    if (ldiff == -1) {
        if (((dToF >= 0) && (dToF <= 2)) && (abs(sdiff) <= 2)) {
            return true;
        } else {
            return false;
        }
    }

    if (ldiff == 0) {
        if (((dToF >= 0) && (dToF <= 2)) && (abs(sdiff) <= 2)) {
            return true;
        } else {
            return false;
        }
    }

    if (ldiff == 1) {
        if (((dToF >= 0) && (dToF <= 1.25)) && (abs(sdiff) <= 1)) {
            return true;
        } else {
            return false;
        }
    }

    if (ldiff == 2) {
        if (((dToF >= 0) && (dToF <= 1.25)) && (abs(sdiff) <= 1)) {
            return true;
        } else {
            return false;
        }
    }

    if (ldiff == 3) {
        if (((dToF >= 0) && (dToF <= 1.25)) && (abs(sdiff) <= 1)) {
            return true;
        } else {
            return false;
        }
    }

    return false;
}

// isNear function ----------------------------------------------------------------------------------------------------------------------------------------------------------

bool isNear(int sdiff, int ldiff) {
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

    if ((ldiff == -2) && (sdiff == -2)) { return true; }
    if ((ldiff == -2) && (sdiff == -1)) { return true; }
    if ((ldiff == -2) && (sdiff == 0)) { return true; }
    if ((ldiff == -2) && (sdiff == 1)) { return true; }
    if ((ldiff == -2) && (sdiff == 2)) { return true; }

    if ((ldiff == -1) && (sdiff == -2)) { return true; }
    if ((ldiff == -1) && (sdiff == -1)) { return true; }
    // if((ldiff==-1) && (sdiff== 0)){return true;}
    if ((ldiff == -1) && (sdiff == 1)) { return true; }
    if ((ldiff == -1) && (sdiff == 2)) { return true; }

    if ((ldiff == 0) && (sdiff == -2)) { return true; }
    // if((ldiff== 0) && (sdiff==-1)){return true;}
    // if((ldiff== 0) && (sdiff== 0)){return true;}
    // if((ldiff== 0) && (sdiff== 1)){return true;}
    if ((ldiff == 0) && (sdiff == 2)) { return true; }

    if ((ldiff == 1) && (sdiff == -2)) { return true; }
    if ((ldiff == 1) && (sdiff == -1)) { return true; }
    // if((ldiff== 1) && (sdiff== 0)){return true;}
    if ((ldiff == 1) && (sdiff == 1)) { return true; }
    if ((ldiff == 1) && (sdiff == 2)) { return true; }

    if ((ldiff == 2) && (sdiff == -2)) { return true; }
    if ((ldiff == 2) && (sdiff == -1)) { return true; }
    if ((ldiff == 2) && (sdiff == 0)) { return true; }
    if ((ldiff == 2) && (sdiff == 1)) { return true; }
    if ((ldiff == 2) && (sdiff == 2)) { return true; }

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

bool isNearCTOF(int sdiff, int ldiff) {
    if ((ldiff == 1) && (sdiff == -3)) { return true; }
    if ((ldiff == 1) && (sdiff == -2)) { return true; }
    if ((ldiff == 1) && (sdiff == -1)) { return true; }
    if ((ldiff == 1) && (sdiff == 1)) { return true; }
    if ((ldiff == 1) && (sdiff == 2)) { return true; }
    if ((ldiff == 1) && (sdiff == 3)) { return true; }

    if ((ldiff == 2) && (sdiff == -3)) { return true; }
    if ((ldiff == 2) && (sdiff == -2)) { return true; }
    if ((ldiff == 2) && (sdiff == -1)) { return true; }
    if ((ldiff == 2) && (sdiff == 0)) { return true; }
    if ((ldiff == 2) && (sdiff == 1)) { return true; }
    if ((ldiff == 2) && (sdiff == 2)) { return true; }
    if ((ldiff == 2) && (sdiff == 3)) { return true; }

    if ((ldiff == 3) && (sdiff == -3)) { return true; }
    if ((ldiff == 3) && (sdiff == -2)) { return true; }
    if ((ldiff == 3) && (sdiff == -1)) { return true; }
    if ((ldiff == 3) && (sdiff == 0)) { return true; }
    if ((ldiff == 3) && (sdiff == 1)) { return true; }
    if ((ldiff == 3) && (sdiff == 2)) { return true; }
    if ((ldiff == 3) && (sdiff == 3)) { return true; }

    return false;
}

// printProgress function ---------------------------------------------------------------------------------------------------------------------------------------------------

void printProgress(double percentage) {
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

};  // namespace neutron_veto_functions

#endif  // NEUTRON_VETO_FUNCTIONS_H
