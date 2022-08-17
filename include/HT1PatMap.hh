#ifndef HT1PATMAP_h
#define HT1PATMAP_h 1

////////////////////////////////////////////////////////////////////////////////
//   HT1PatMap.cc for Hadron Theraphy                                         //
//                                                                            //
//   Definitions of HT1PatMap class's function definition                     //
//                                                                            //
//              - 20. Ock. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include <map>
#include "G4String.hh"

class HT1PatMap
{
  public:
	// Constructor and destructor
	HT1PatMap(G4int ND);
	~HT1PatMap();

	// Initializaion
	void Init();

	// Get
	G4String GetDetNameFromDetID(G4int id){return m_NameFromID[id];};
	G4int GetDetIDFromDetName(G4String name){return m_IDFromName[name];};

  private:
	// naming for patient parts 
	std::map<G4int, G4String> m_NameFromID;
	std::map<G4String, G4int> m_IDFromName;
    
	// number of detectors 
	G4int m_ND;
};

#endif
