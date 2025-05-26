//
// Created by Alon Sportes on 27/03/2025.
//

#ifndef NEUTRONVETOVARIABLES_H
#define NEUTRONVETOVARIABLES_H

#include <TVector3.h>

#include <cmath>
#include <iostream>
#include <string>

// Include libraries:
// #include "../../namespaces/general_utilities/basic_tools.h"
#include "../../namespaces/general_utilities/constants.h"

// Include CLAS12 libraries:
#include "../../includes/clas12_include.h"

struct NeutronVetoVariables {
    // Andrew's response checks:
    bool CT, C1, C2, C3;

    // Neutron position and momentum
    TVector3 v_nvtx_3v;
    TVector3 v_hit_3v;
    TVector3 v_path_3v;
    double ToF_n;
    double theta_n;
    double beta_n;
    double gamma_n;
    double P_n_mag;

    int detINTlayer;
    int detlayer;

    int Size_CND1, Size_CND2, Size_CND3, Size_CND;
    int LayerMult_CND1, LayerMult_CND2, LayerMult_CND3, LayerMult_CND;
    double Edep_CND1, Edep_CND2, Edep_CND3, Edep_CND;
    double Edep_single;
    double Edep_CTOF;

    // TVector3 v_nvtx_3v(nvtx_x, nvtx_y, nvtx_z);  // Neutron's vertex location
    // TVector3 v_hit_3v;  // Neutron's hit location in CND
    // v_hit_3v.SetXYZ(AllParticles[itr1]->sci(detlayer)->getX(), AllParticles[itr1]->sci(detlayer)->getY(), AllParticles[itr1]->sci(detlayer)->getZ());
    // TVector3 v_path_3v = v_hit_3v - v_nvtx_3v;  // Direct calculation of neutron's path (in vector form)

    TVector3 P_n_3v;

    // Why "v_path_3v.Mag() / 100"? unit conversion.
    // TODO: check if this unit conversion is needed!
    double path_n;
    double theta_n_miss;
    double dpp;
    int nSector;

    NeutronVetoVariables(std::vector<region_part_ptr>& AllParticles, std::vector<region_part_ptr>& Electrons, int itr1, double starttime, const TVector3& P_miss_3v) {
        // Andrew's response checks:
        CT = (AllParticles[itr1]->sci(clas12::CTOF)->getDetector() == 4);
        C1 = (AllParticles[itr1]->sci(clas12::CND1)->getDetector() == 3);
        C2 = (AllParticles[itr1]->sci(clas12::CND2)->getDetector() == 3);
        C3 = (AllParticles[itr1]->sci(clas12::CND3)->getDetector() == 3);

        // Explicit calculation of the neutron's momentum (to bypass cases where P_n is E_dep)
        v_nvtx_3v = neutron_veto_functions::GetVzHitLocation(Electrons[0]);                // Neutron's vertex location -> set as the electron vertex
        v_hit_3v = neutron_veto_functions::GetVzHitInCND(AllParticles[itr1]);              // Neutron's hit location in CND
        v_path_3v = neutron_veto_functions::GetnCDPath(AllParticles[itr1], Electrons[0]);  // Direct calculation of neutron's path (in vector form)
        ToF_n = neutron_veto_functions::GetnCDToF(AllParticles[itr1], starttime);
        theta_n = AllParticles[itr1]->getTheta() * 180 / M_PI;
        beta_n = neutron_veto_functions::GetnCDBeta(AllParticles[itr1], Electrons[0], starttime);
        // double beta_n = AllParticles[itr1]->par()->getBeta();
        gamma_n = 1 / sqrt(1 - (beta_n * beta_n));
        P_n_mag = gamma_n * beta_n * constants::m_n;
        // double ToF_n = AllParticles[itr1]->getTime() - starttime;

        detINTlayer = C1 ? 1 : C2 ? 2 : 3;
        detlayer = C1 ? clas12::CND1 : C2 ? clas12::CND2 : clas12::CND3;  // CND layer with hit

        Size_CND1 = AllParticles[itr1]->sci(clas12::CND1)->getSize();
        Size_CND2 = AllParticles[itr1]->sci(clas12::CND2)->getSize();
        Size_CND3 = AllParticles[itr1]->sci(clas12::CND3)->getSize();
        Size_CND = Size_CND1 + Size_CND2 + Size_CND3;

        LayerMult_CND1 = AllParticles[itr1]->sci(clas12::CND1)->getLayermulti();
        LayerMult_CND2 = AllParticles[itr1]->sci(clas12::CND2)->getLayermulti();
        LayerMult_CND3 = AllParticles[itr1]->sci(clas12::CND3)->getLayermulti();
        LayerMult_CND = LayerMult_CND1 + LayerMult_CND2 + LayerMult_CND3;

        Edep_CND1 = AllParticles[itr1]->sci(clas12::CND1)->getEnergy();
        Edep_CND2 = AllParticles[itr1]->sci(clas12::CND2)->getEnergy();
        Edep_CND3 = AllParticles[itr1]->sci(clas12::CND3)->getEnergy();
        Edep_CND = Edep_CND1 + Edep_CND2 + Edep_CND3;

        Edep_single = AllParticles[itr1]->sci(detlayer)->getEnergy();
        Edep_CTOF = AllParticles[itr1]->sci(clas12::CTOF)->getEnergy();

        // TVector3 v_nvtx_3v(nvtx_x, nvtx_y, nvtx_z);  // Neutron's vertex location
        // TVector3 v_hit_3v;  // Neutron's hit location in CND
        // v_hit_3v.SetXYZ(AllParticles[itr1]->sci(detlayer)->getX(), AllParticles[itr1]->sci(detlayer)->getY(), AllParticles[itr1]->sci(detlayer)->getZ());
        // TVector3 v_path_3v = v_hit_3v - v_nvtx_3v;  // Direct calculation of neutron's path (in vector form)

        P_n_3v.SetMagThetaPhi(P_n_mag, v_path_3v.Theta(), v_path_3v.Phi());

        // Why "v_path_3v.Mag() / 100"? unit conversion.
        // TODO: check if this unit conversion is needed!
        path_n = v_path_3v.Mag() / 100;
        // double path_n = v_path_3v.Mag();
        theta_n_miss = P_n_3v.Angle(P_miss_3v) * 180 / M_PI;       // Opening angle between calculated neutron's momentum and predicted neutron momentum (= missing momentum)
        dpp = (P_miss_3v.Mag() - P_n_3v.Mag()) / P_miss_3v.Mag();  // Momentum resolution
        nSector = AllParticles[itr1]->sci(detlayer)->getSector();  // Number of CND sector with a neutron hit in the layer detlayer
    }
};

#endif  // NEUTRONVETOVARIABLES_H
