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
/// \file HT1ActIni.hh
/// \brief Definition of the HT1ActIni class

#ifndef HT1ActIni_h
#define HT1ActIni_h 1

#include "G4VUserActionInitialization.hh"

#include "HT1PriGenAct.hh"
#include "HT1EveAct.hh"
#include "HT1SteAct.hh"

/// Action initialization class.

class HT1ActIni : public G4VUserActionInitialization
{
    public:
        HT1ActIni(HT1ConMan* CM);
        virtual ~HT1ActIni();

        virtual void BuildForMaster() const;
        virtual void Build() const;
    
    private:
        HT1ConMan* m_CM;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
