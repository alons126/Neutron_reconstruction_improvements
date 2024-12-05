
#ifndef NEUTRONFUNCTIONS_H
#define NEUTRONFUNCTIONS_H

// ConfigOutPutName function --------------------------------------------------------------------------------------------------------------------------------------------------

void SetNeutronCounters(const bool isGN, int &counter_n_multiplicity_allN, int &counter_n_multiplicity_goodN, int &counter_n_multiplicity_badN)
{
    ++counter_n_multiplicity_allN;

    if (isGN)
    {
        ++counter_n_multiplicity_goodN;
    }
    else
    {
        ++counter_n_multiplicity_goodN;
    }
}

void SetNeutronCounters(const bool pInCD, const bool pInFD,
                        const bool isGN,
                        int &counter_n_multiplicity_allN_epCD, int &counter_n_multiplicity_goodN_epCD, int &counter_n_multiplicity_badN_epCD,
                        int &counter_n_multiplicity_allN_epFD, int &counter_n_multiplicity_goodN_epFD, int &counter_n_multiplicity_badN_epFD)
{
    if (pInCD)
    {
        SetNeutronCounters(isGN, counter_n_multiplicity_allN_epCD, counter_n_multiplicity_goodN_epCD, counter_n_multiplicity_badN_epCD);
    }
    else if (pInFD)
    {
        SetNeutronCounters(isGN, counter_n_multiplicity_allN_epFD, counter_n_multiplicity_goodN_epFD, counter_n_multiplicity_badN_epFD);
    }
}

#endif // NEUTRONFUNCTIONS_H
