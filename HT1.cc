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
/// \file HT1.cc
/// \brief Main program of the HT1 example

#include "HT1DetCon.hh"
#include "HT1ConMan.hh"
#include "HT1ActIni.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "QBBC.hh"
#include "QGSP_BERT_HP.hh"
//#include "QGSC_CHIPS.hh"
//#include "CHIPS.hh"
//#include "QGSC_BERT.hh"
//#include "QGSP.hh"
#include "QGSP_BERT.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

	int option; 				// options 
	int flag_b = 0;			// bach mode flag 
	char * macro;				// macro file for bash mode
	int flag_d = 0;			// division mode (it showld be folowed by an int argumets to devide patient material)
	int m_PatiDV = 0;
	int flag_ZL = 0;			// division mode (it showld be folowed by an int argumets to devide patient material)
	G4double m_PatiZL = 0.0; // patient Z length

	// switch cases
	while ((option = getopt(argc, argv, "b:d:l:")) != -1) {
		switch (option) {
			case 'b': 
				flag_b = 1;
				macro = optarg;
				break;
			case 'd':
				flag_d = 1;
				m_PatiDV = atoi(optarg);
				break;
			case 'l':
				flag_ZL = 1;
				m_PatiZL = std::stod(optarg);
				break;
			default:
				fprintf(stderr, "Usage: %s \n",  argv[0]);
				exit(EXIT_FAILURE);
        }
    }

  // Detect interactive mode and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( ! flag_b )ui = new G4UIExecutive(argc, argv);

    const G4String config = "/home/abuladze/simulations/HadronTheraphy1/HT1/config.cfg";

    // read config file 
    HT1ConMan* CM = new HT1ConMan(config);

		// change  value in config file object
		if(flag_d){
			CM->SetPatiDV(m_PatiDV); 
			CM->SetPatiMap(m_PatiDV);   
		}
		if(flag_ZL)CM->SetPatiZL(m_PatiZL);

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new HT1DetCon(CM));

  // Physics list
  G4VModularPhysicsList* physicsList = new QGSP_BERT;
	//G4VModularPhysicsList* physicsList = new HT1PhyList();
	std::cout << "physics list" << std::endl;
	physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
    
  // User action initialization
  runManager->SetUserInitialization(new HT1ActIni(CM));
  
    // initialize G4 kernel
    runManager->Initialize();
    
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = macro;
    UImanager->ApplyCommand(command+fileName);
  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
