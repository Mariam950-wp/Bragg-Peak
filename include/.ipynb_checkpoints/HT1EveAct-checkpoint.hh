/// \file HT1EveAct.hh
/// \brief Definition of the HT1EveAct class

#include "HT1RunAct.hh"

#ifndef HT1EveAct_h
#define HT1EveAct_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "HT1RunAct.hh"

/// Event action class
///

class HT1EveAct : public G4UserEventAction
{
  public:
  	HT1EveAct(HT1RunAct* RA, HT1ConMan* CM);
	virtual ~HT1EveAct();

	virtual void BeginOfEventAction(const G4Event* event);
	virtual void EndOfEventAction(const G4Event* event);
        
	// energy deposition to degrader, collimator nad patient
	void ProtonEDepToPati(G4double Edep, G4int ID){m_ProEdepFromID[ID] += Edep;};
	void PhotonEDepToPati(G4double Edep, G4int ID){m_PhoEdepFromID[ID] += Edep;};
	void EletonEDepToPati(G4double Edep, G4int ID){m_EleEdepFromID[ID] += Edep;}; 
        
    // photon energy deposition to spherical detectors
    void PhotonEDepToSPHD(G4double Edep, G4int ID){m_IDtoPhoEdepSD[ID] += Edep;}; 

	private:
	HT1RunAct* m_RA;
        
	//config object
	HT1ConMan* m_CM;
        
	// energy deposition variables
	std::map<G4int, G4double> m_ProEdepFromID;
	std::map<G4int, G4double> m_PhoEdepFromID;
	std::map<G4int, G4double> m_EleEdepFromID;

    // energy deposition maps for spherical detectors 
    std::map<G4int, G4double> m_IDtoPhoEdepSD;

	// number of patient parts
	G4int      m_PatiDV;
        
    // number of angular divisions 
    G4int    m_AzziDV;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
