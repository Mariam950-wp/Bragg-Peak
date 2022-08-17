////////////////////////////////////////////////////////////////////////////////
//   HT1SphDetMap.cc for Hadron Theraphy                                      //
//                                                                            //
//   Definitions of HT1SphDetMap class's member functions                     //
//                                                                            //
//              - 20. Ock. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include "HT1SphDetMap.hh"

//////////////////////////////////////////////////
//   Constructors and destructor                //
//////////////////////////////////////////////////
HT1SphDetMap::HT1SphDetMap(G4int AzziDV):m_AzziDV(AzziDV)
{
	Init();
}

HT1SphDetMap::~HT1SphDetMap()
{
}


//////////////////////////////////////////////////
//   Initialization                             //
//////////////////////////////////////////////////
void HT1SphDetMap::Init()
{
    	// Mapping from ID to name
    for (G4int j = 0; j<m_AzziDV; j++)
        m_DetNameFromDetID[j]   = "SphereDetector_" + std::to_string(j); 

	// Mapping from name to ID
    for (G4int j = 0; j<m_AzziDV; j++)
        m_DetIDFromDetName[m_DetNameFromDetID[j]] = j;
}