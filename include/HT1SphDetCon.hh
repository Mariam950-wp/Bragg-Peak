#ifndef HT1SPHDETCON_h
#define HT1SPHDETCON_h 1

////////////////////////////////////////////////////////////////////////////////
//   HT1SphDetCon.hh for hadron theraphy                                         //
//                                                                            //
//              - 20. Sep. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include "G4NistManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "HT1SphDetMap.hh"

/// Detector construction class to define materials and geometry.

class HT1DetCon;

class HT1SphDetCon : public G4VUserDetectorConstruction
{
  public:
    HT1SphDetCon(G4LogicalVolume* logicWorld, G4double DetIR, HT1SphDetMap* SDM);
    virtual ~HT1SphDetCon();

    virtual G4VPhysicalVolume* Construct();

    protected:
        G4NistManager*   nist;
    
    private:
        G4LogicalVolume*   m_logicWorld;
        G4double           m_DetIR;
        HT1SphDetMap*      m_SDM;
    
        G4Element*         m_ElLu;
        G4Element*         m_ElY;
        G4Element*         m_ElSi;
        G4Element*         m_ElO;
        G4Element*         m_ElCe;
    
        G4int                           m_AzziDV;
        G4Material*                     m_DetMat;
        
        G4Sphere*          m_SphDetector;
        G4LogicalVolume*   m_logicDetector;
        G4VPhysicalVolume* m_physDetector;
        
};


#endif

