#ifndef NEUTRONFUNCTIONS_H
#define NEUTRONFUNCTIONS_H

#include <cstdlib>
#include <iostream>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TTree.h"
//
#include "../constants.h"
//
#include "HipoChain.h"
#include "clas12reader.h"

using namespace std;
using namespace clas12;

// SetNeutronCounters function ------------------------------------------------------------------------------------------------------------------------------------------------

void SetNeutronCounters(const bool isGN, int &counter_n_multiplicity_allN, int &counter_n_multiplicity_goodN, int &counter_n_multiplicity_badN) {
    ++counter_n_multiplicity_allN;

    if (isGN) {
        ++counter_n_multiplicity_goodN;
    } else {
        ++counter_n_multiplicity_goodN;
    }
}

// SetNeutronCounters function ------------------------------------------------------------------------------------------------------------------------------------------------

void SetNeutronCounters(const bool pInCD, const bool pInFD, const bool isGN, int &counter_n_multiplicity_allN_epCD, int &counter_n_multiplicity_goodN_epCD, int &counter_n_multiplicity_badN_epCD,
                        int &counter_n_multiplicity_allN_epFD, int &counter_n_multiplicity_goodN_epFD, int &counter_n_multiplicity_badN_epFD) {
    if (pInCD) {
        SetNeutronCounters(isGN, counter_n_multiplicity_allN_epCD, counter_n_multiplicity_goodN_epCD, counter_n_multiplicity_badN_epCD);
    } else if (pInFD) {
        SetNeutronCounters(isGN, counter_n_multiplicity_allN_epFD, counter_n_multiplicity_goodN_epFD, counter_n_multiplicity_badN_epFD);
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
    double ToF_n;

    // Andrew's response checks:
    bool nCD_C1_hit = (NeutronCD->sci(CND1)->getDetector() == 3);
    bool nCD_C2_hit = (NeutronCD->sci(CND2)->getDetector() == 3);
    bool nCD_C3_hit = (NeutronCD->sci(CND3)->getDetector() == 3);
    auto nCD_detlayer = nCD_C1_hit ? CND1 : nCD_C2_hit ? CND2 : CND3;  // CND layer with hit

    ToF_n = NeutronCD->sci(nCD_detlayer)->getTime() - starttime;

    return ToF_n;
}

// GetnCDBeta function ------------------------------------------------------------------------------------------------------------------------------------------------

double GetnCDBeta(region_part_ptr NeutronCD, region_part_ptr Electron, double starttime) {
    double Beta_n;

    TVector3 V_path_3v = GetnCDPath(NeutronCD, Electron);  // Direct calculation of neutron's path (in vector form)
    double ToF_n = GetnCDToF(NeutronCD, starttime);

    Beta_n = V_path_3v.Mag() / (ToF_n * c);

    return Beta_n;
}

#endif  // NEUTRONFUNCTIONS_H
