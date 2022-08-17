////////////////////////////////////////////////////////////////////////////////
//   HT1PatMap.cc for Hadron Theraphy                                         //
//                                                                            //
//   Definitions of HT1PatMap class's member functions                        //
//                                                                            //
//              - 20. Ock. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include "HT1PatMap.hh"

//////////////////////////////////////////////////
//   Constructors and destructor                //
//////////////////////////////////////////////////
HT1PatMap::HT1PatMap(G4int ND):m_ND(ND)
{
	m_NameFromID.clear();
	m_IDFromName.clear();
    
	Init();
}

HT1PatMap::~HT1PatMap()
{
}


//////////////////////////////////////////////////
//   Initialization                             //
//////////////////////////////////////////////////
void HT1PatMap::Init()
{

	// Mapping from ID to name
    for (G4int i = 0; i<m_ND; i++)
        m_NameFromID[i]   = "Patient_" + std::to_string(i); 

	// Mapping from name to ID
	for ( int i = 0; i < m_ND; i++ ) 
		m_IDFromName[m_NameFromID[i]] = i;
}