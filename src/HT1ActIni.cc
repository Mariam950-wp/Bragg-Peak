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
/// \file HT1ActIni.cc
/// \brief Implementation of the HT1ActIni class

#include "HT1ActIni.hh"

#include "HT1RunAct.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1ActIni::HT1ActIni(HT1ConMan* CM)
 : G4VUserActionInitialization(), m_CM(CM)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1ActIni::~HT1ActIni()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HT1ActIni::BuildForMaster() const
{
  HT1RunAct* RA = new HT1RunAct(m_CM);
  SetUserAction(RA);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HT1ActIni::Build() const
{
    G4cout << "built ***********************" << G4endl;
  SetUserAction(new HT1PriGenAct(m_CM));

  HT1RunAct* RA = new HT1RunAct(m_CM);
  SetUserAction(RA);
  
  HT1EveAct* EA = new HT1EveAct(RA, m_CM);
  SetUserAction(EA);
  
  SetUserAction(new HT1SteAct(EA, m_CM->GetPatiMap(), m_CM->GetSpDeMap()));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
