#ifndef HT1PHYLIST_h
#define HT1PHYLIST_h 1

////////////////////////////////////////////////////////////////////////////////
//   HT1PhyList.hh for hadron theraphy                                        //
//                                                                            //
//              - 20. Sep. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include "G4VModularPhysicsList.hh"

class HT1PhyList: public G4VModularPhysicsList
{
  public:
		HT1PhyList();
		virtual ~HT1PhyList();
   
		virtual void SetCuts();
    
};

#endif
