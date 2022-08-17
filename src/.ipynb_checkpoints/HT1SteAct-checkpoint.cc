//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file HT1SteAct.cc
/// \brief Implementation of the HT1SteAct class

#include "HT1SteAct.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1SteAct::HT1SteAct(HT1EveAct* EA, HT1PatMap* PM, HT1SphDetMap* SDM)
: G4UserSteppingAction(),m_EA(EA), m_PM(PM), m_SDM(SDM)
{
	m_currTrackID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1SteAct::~HT1SteAct()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HT1SteAct::UserSteppingAction(const G4Step* step)
{
    // I wrote examples of some information which can be extracted from a step.
    // Uncomment whatever you want to use.

    // Who am I? Where am I? What am I undergoing?
    // // Track ID
    //G4int trackID = step -> GetTrack() -> GetTrackID();
    // // Particle name
    G4String parName = step -> GetTrack() -> GetDefinition() -> GetParticleName();
    // // Particle ID
    //G4int parID = step -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    // // Particle charge
    //G4int parCharge = step -> GetTrack() -> GetDefinition() -> GetPDGCharge();
    // // Process name
    //const G4VProcess* creProc = step -> GetTrack() -> GetCreatorProcess();
    //	G4String procName = step -> GetPostStepPoint() -> GetProcessDefinedStep() -> GetProcessName();
    // // Physical volume
    G4VPhysicalVolume* postPV = step -> GetPostStepPoint() -> GetPhysicalVolume();
    G4String namePrePV = step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName();
    G4String namePostPV = "";
    if(postPV)namePostPV = postPV -> GetName();
    //	if ( postPV != 0 ) namePostPV = postPV -> GetName();
    //	else namePostPV = "outside";
    // // Status
    //	G4StepStatus preStat = step -> GetPreStepPoint() -> GetStepStatus();
    //	G4StepStatus postStat = step -> GetPostStepPoint() -> GetStepStatus();
    // // Position
    //G4ThreeVector prePos = step -> GetPreStepPoint() -> GetPosition();
    //G4ThreeVector postPos = step -> GetPostStepPoint() -> GetPosition();
    // // Momentum
    //	G4ThreeVector preMom = step -> GetPreStepPoint() -> GetMomentum();
    //	G4ThreeVector postMom = step -> GetPostStepPoint() -> GetMomentum();
    // // Kinetic energy
    G4double preKinEgy = step -> GetPreStepPoint() -> GetKineticEnergy();
    G4double postKinEgy = step -> GetPostStepPoint() -> GetKineticEnergy();
	// // Energy deposition
    G4double eDep = step -> GetTotalEnergyDeposit();
    
    
    // proton, gamma, eleqtron energy deposited to patient parts 
	if ( namePrePV.contains("Patient")){
		if(parName == "proton")
			m_EA -> ProtonEDepToPati(eDep,m_PM -> GetDetIDFromDetName(namePrePV));
		else if (parName == "gamma")
			m_EA -> PhotonEDepToPati(eDep,m_PM -> GetDetIDFromDetName(namePrePV));
		else if (parName == "e-" || parName == "e+")
			m_EA -> EletonEDepToPati(eDep,m_PM -> GetDetIDFromDetName(namePrePV));
    }else if ( namePrePV.contains("SphereDetector_") && parName == "gamma"){
        m_EA -> PhotonEDepToSPHD(eDep,m_SDM ->GetDetIDFromDetName(namePrePV));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

