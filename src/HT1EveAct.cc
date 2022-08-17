/// \file HT1EveAct.cc
/// \brief Implementation of the HT1EveAct class

#include "HT1EveAct.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1EveAct::HT1EveAct(HT1RunAct* RA, HT1ConMan* CM)
: G4UserEventAction(),m_RA(RA), m_CM(CM)
{
    m_PatiDV = m_CM->GetPatiDV();
    m_AzziDV = m_CM->GetAzziDV();
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1EveAct::~HT1EveAct()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HT1EveAct::BeginOfEventAction(const G4Event*)
{
	// initial enegry deposition is zero
	for (G4int i = 0; i<m_PatiDV; i++){
		m_ProEdepFromID[i] = 0.0;
		m_PhoEdepFromID[i] = 0.0;
		m_EleEdepFromID[i] = 0.0;
	}
    
    for (G4int j = 0; j<m_AzziDV; j++)
        m_IDtoPhoEdepSD[j] = 0.0;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HT1EveAct::EndOfEventAction(const G4Event*)
{   
  G4AnalysisManager* m_AM = G4AnalysisManager::Instance();
			
	// 2D bragg peak
	for (G4int i = 0; i<m_PatiDV; i++)
		if (m_ProEdepFromID[i] > 0.0)
			m_AM->FillH2(0, i, m_ProEdepFromID[i]);

	// 2D gamaa distribution
	for (G4int i = 0; i<m_PatiDV; i++)
		if (m_PhoEdepFromID[i] > 0.0)
			m_AM->FillH2(1, i, m_PhoEdepFromID[i]);

	// 2D electron distribution
	for (G4int i = 0; i<m_PatiDV; i++)
		if (m_EleEdepFromID[i] > 0.0)
			m_AM->FillH2(2, i, m_EleEdepFromID[i]);
    
    // angular Gamma distribution
    for (G4int j = 0; j<m_AzziDV; j++)
        if (m_IDtoPhoEdepSD[j] > 0.0)
            m_AM->FillH2(3, j, m_IDtoPhoEdepSD[j]);
    
}



