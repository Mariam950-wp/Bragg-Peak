/// \file HT1PriGenAct.hh
/// \brief Definition of the HT1PriGenAct class

#ifndef HT1PriGenAct_h
#define HT1PriGenAct_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

#include "HT1ConMan.hh"

class G4ParticleGun;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued 
/// in front of the phantom across 80% of the (X,Y) phantom size.

class HT1PriGenAct : public G4VUserPrimaryGeneratorAction
{
    public:
        HT1PriGenAct(HT1ConMan* CM);    
        virtual ~HT1PriGenAct();

        // method from the base class
        virtual void GeneratePrimaries(G4Event*);         
  
        // method to access particle gun
        const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
  
    private:
        G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
    
        // config file
        HT1ConMan* m_CM;
        
        // beam
        G4double      m_BeamKE;
        G4double      m_BeamOR;
        G4String      m_SParName;
        G4ThreeVector m_ParPos;
        G4double      m_Z0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
