//
// ********************************************************************
// *                                                                  *
// ********************************************************************
//
//
/// \file HT1RunAct.hh
/// \brief Definition of the HT1RunAct class

#ifndef HT1RunAct_h
#define HT1RunAct_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

#include "HT1Analysis.hh"
#include "HT1ConMan.hh"
#include "HT1PriGenAct.hh"
#include "HT1DetCon.hh"

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class HT1RunAct : public G4UserRunAction
{
  public:
    HT1RunAct(HT1ConMan* CM);
    virtual ~HT1RunAct();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

  private:
      //config object
      HT1ConMan* m_CM;
      
      // limits for 1D histogram
      G4double m_BeamKE;
      
      // number of patient parts
      G4int      m_PatiDV;
      
      // number of bins
      G4int      m_bn = 1000;
      
    // number of azimutal and polar divisions of spherical detector 
    G4int    m_PolaDV;
    G4int    m_AzziDV;
};

#endif

