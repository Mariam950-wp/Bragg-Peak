#ifndef HT1DetCon_h
#define HT1DetCon_h 1

////////////////////////////////////////////////////////////////////////////////
//   HT1DetCon.hh for hadron theraphy                                         //
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
#include "G4MagneticField.hh"

#include "HT1ConMan.hh"
#include "HT1SphDetCon.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class HT1DetCon : public G4VUserDetectorConstruction
{
  public:
    HT1DetCon(HT1ConMan* CM);
    virtual ~HT1DetCon();

    virtual G4VPhysicalVolume* Construct();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  protected:
    G4LogicalVolume*  fScoringVolume;
    G4NistManager*    nist;
    
    private:
    
        HT1ConMan* m_CM;
        
        // spherical detector
        HT1SphDetCon* m_SpDet;
        
        //reaction area
        G4double    m_InterXL;
        G4double    m_InterYL;
        G4Material* m_WorldMat;
    
        // patient
        G4double                        m_ColltoPat;
        G4double                        m_PatiZL;
        G4int                           m_PatiDV;
        G4Material*                     m_PatiMat;
        std::vector< G4Box* >           m_Patient;
        std::vector< G4LogicalVolume* > m_logicPatient;
        std::vector< G4ThreeVector >    m_PatientPos;
        
        HT1PatMap*                      m_PM;
        
        // world dimessions 
        G4double    m_WorldX;
        G4double    m_WorldY;
        G4double    m_WorldZ;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

