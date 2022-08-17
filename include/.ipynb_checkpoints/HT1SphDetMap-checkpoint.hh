#ifndef HT1SPHDETMAP_h
#define HT1SPHDETMAP_h 1

////////////////////////////////////////////////////////////////////////////////
//   HT1SphDetMap.cc for Hadron Theraphy                                      //
//                                                                            //
//   Definitions of HT1SphDetMap class's function definition                  //
//                                                                            //
//              - 20. Ock. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include <map>
#include "G4String.hh"

class HT1SphDetMap
{
  public:
	// Constructor and destructor
	HT1SphDetMap(G4int AzziDV);
	~HT1SphDetMap();

	// Initializaion
	void Init();
    
    // Get
	G4String GetDetNameFromDetID(G4int id)  {return m_DetNameFromDetID[id];};
	G4int GetDetIDFromDetName(G4String name){return m_DetIDFromDetName[name];};
    
    // get map parameters
    G4int GetNumAziDV(){return m_AzziDV;};

  private:
    
	// naming for detector parts 
	std::map<G4int,G4String> m_DetNameFromDetID;
	std::map<G4String,G4int> m_DetIDFromDetName;
    
    G4int    m_AzziDV;
};

#endif
