/// \file HT1RunAct.cc
/// \brief Implementation of the HT1RunAct class

#include <ctime>

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "HT1RunAct.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1RunAct::HT1RunAct(HT1ConMan* CM)
: G4UserRunAction(), m_CM(CM)
{ 
    // Get analysis manager
    auto m_AM = G4AnalysisManager::Instance();
    
    // limits for 1D histogram
    m_BeamKE = m_CM -> GetBeamKE();
    // number of patient parts 
    m_PatiDV = m_CM->GetPatiDV();
    
    // number of anguar divizions 
    m_AzziDV = m_CM -> GetAzziDV();


	// 2D histogramm for bragg peak
	m_AM->CreateH2("Proton_Edep_CV_VolID", "Proton_Energy_deposition_VC_VolID", 
		m_PatiDV, 0.0, m_PatiDV, m_bn, 0., m_BeamKE/4);

	// 2D histogramm for bragg peak
	m_AM->CreateH2("Photon_Edep_CV_VolID", "Photon_Energy_deposition_VC_VolID", 
		m_PatiDV, 0.0, m_PatiDV, m_bn, 0., m_BeamKE/10000);

	// 2D histogramm for bragg peak
	m_AM->CreateH2("Electron_Edep_CV_VolID", "Electron_Energy_deposition_VC_VolID", 
                   m_PatiDV, 0.0, m_PatiDV, m_bn, 0., m_BeamKE/100);
    
    // angular distribution of prompt gamma energy deposition
    m_AM->CreateH2("Angular_dist_of_Prompt_Gamma_EDep_on_SD", "Angular_dist_of_Prompt_Gamma_EDep_on_SD", 
                   m_AzziDV, 0.0, m_AzziDV, 1000, 0.0, m_BeamKE/10*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1RunAct::~HT1RunAct()
{
    delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HT1RunAct::BeginOfRunAction(const G4Run*)
{
    auto m_AM = G4AnalysisManager::Instance();
    
    //time_t rawtime;
    //struct tm * timeinfo;
    //char buffer [80];
    //time (&rawtime);
    //timeinfo = localtime (&rawtime);
    //strftime (buffer,80,"%Y-%m-%d-%H:%M",timeinfo);
    // Open an output file
    //m_AM->OpenFile("/home/abuladze/simulations/HadronTheraphy1/Output/" + std::string(buffer));

    m_AM->OpenFile("/home/abuladze/simulations/HadronTheraphy1/Output/HT1_Output");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HT1RunAct::EndOfRunAction(const G4Run* run)
{ 
    auto m_AM = G4AnalysisManager::Instance();
    // Save histograms
    m_AM->Write();
    // close file
    m_AM->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

