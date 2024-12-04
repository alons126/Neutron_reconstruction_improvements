
#ifndef NEUTRONFUNCTIONS_H
#define NEUTRONFUNCTIONS_H

// ConfigOutPutName function --------------------------------------------------------------------------------------------------------------------------------------------------

void SetNeutronCounters(const bool isGN, int &counter_nsize_allN, int &counter_nsize_goodN, int &counter_nsize_badN)
{
    ++counter_nsize_allN;

    if (isGN)
    {
        ++counter_nsize_goodN;
    }
    else
    {
        ++counter_nsize_goodN;
    }
}

#endif // NEUTRONFUNCTIONS_H
