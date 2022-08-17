//
// ********************************************************************
// *                                                                  *
// ********************************************************************
//
//
/// \file HT1SteAct.hh
/// \brief Definition of the HT1SteAct class

#ifndef HT1SteAct_h
#define HT1SteAct_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include "HT1EveAct.hh"
#include "HT1DetCon.hh"

class G4LogicalVolume;

/// Stepping action class
/// 

class HT1SteAct : public G4UserSteppingAction
{
    public:
        HT1SteAct(HT1EveAct* EA, HT1PatMap* PM, HT1SphDetMap* SDM);
        virtual ~HT1SteAct();

        // method from the base class
        virtual void UserSteppingAction(const G4Step*);

    private:
        HT1EveAct* m_EA;
        HT1PatMap* m_PM;
        HT1SphDetMap* m_SDM;

        int m_currTrackID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
