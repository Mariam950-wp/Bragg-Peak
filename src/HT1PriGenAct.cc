/// \file HT1PriGenAct.cc
/// \brief Implementation of the HT1PriGenAct class

#include "HT1PriGenAct.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1PriGenAct::HT1PriGenAct(HT1ConMan* CM)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),m_CM(CM)
{
    // beam parameters 
    m_BeamKE   = m_CM -> GetBeamKE();
    m_BeamOR   = m_CM -> GetBeamOR();
    m_SParName = m_CM -> GetParName();
    
    // particle innitial Z position  
    m_Z0 = 0.0;
    
    fParticleGun  = new G4ParticleGun();

    // default particle kinematic
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(m_SParName);
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(m_BeamKE);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1PriGenAct::~HT1PriGenAct()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HT1PriGenAct::GeneratePrimaries(G4Event* anEvent)
{
    // choose radius and angle
    G4double r = abs(G4RandGauss::shoot(0,m_BeamOR/3));
    G4double theta = 2*M_PI*(G4UniformRand()-0.5);

    // set values to position vector
    m_ParPos.setX(r*sin(theta));
    m_ParPos.setY(r*cos(theta));
    m_ParPos.setZ(m_Z0);

    //std::cout << "Event ID: " << anEvent->GetEventID() << std:: endl;
    // set position to particlegun
    fParticleGun->SetParticlePosition(m_ParPos);

    // generate primaries 
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

