////////////////////////////////////////////////////////////////////////////////
//   HT1SphDetCon.hh for hadron theraphy                                      //
//                                                                            //
//              - 20. Sep. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include "G4VisAttributes.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4RunManager.hh"
#include "G4PVPlacement.hh"

#include "HT1SphDetCon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1SphDetCon::HT1SphDetCon(G4LogicalVolume* logicWorld, G4double DetIR, HT1SphDetMap* SDM):
    m_logicWorld(logicWorld), m_DetIR(DetIR), m_SDM(SDM)
{ 
    nist = G4NistManager::Instance();
        
    m_AzziDV = m_SDM->GetNumAziDV();
        
    Construct();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1SphDetCon::~HT1SphDetCon()
{ 
    delete m_ElLu;
	delete m_ElY;
    delete m_ElSi;
    delete m_ElO;
    delete m_ElCe;
    
    delete m_DetMat;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* HT1SphDetCon::Construct()
{    
    const G4double labTemp = 300.0 * kelvin;
    
    m_ElLu = new G4Element( "Lutetium", "Lu", 71, 176.944  *g/mole);
    m_ElY  = new G4Element(  "Yttrium",  "Y", 39,  88.90585*g/mole);
    m_ElSi = new G4Element(  "Silicon", "Si", 14,  28.0855 *g/mole);
    m_ElO  = new G4Element(   "Oxygen",  "O",  8,  15.9994 *g/mole);
    m_ElCe = new G4Element(   "Cerium", "Ce", 58, 140.116  *g/mole);
    
    m_DetMat = new G4Material("LYSO", 7.1*g/cm3, 5, kStateSolid, labTemp);
	m_DetMat -> AddElement(m_ElLu, 71.43/100.0);
	m_DetMat -> AddElement(m_ElY,   4.03/100.0);
	m_DetMat -> AddElement(m_ElSi,  6.37/100.0);
	m_DetMat -> AddElement(m_ElO,  18.14/100.0);
	m_DetMat -> AddElement(m_ElCe,  0.02/100.0);
    
    G4bool checkOverlaps = true; 
  
    G4double AU = M_PI/m_AzziDV;
    
    for(G4int j = 0; j < m_AzziDV; j++){
        m_SphDetector = 
            new G4Sphere(m_SDM->GetDetNameFromDetID(j),(m_DetIR-300.0)*mm,m_DetIR*mm,0.0,2*M_PI,AU*j,AU);
        
        m_logicDetector = new G4LogicalVolume(m_SphDetector,m_DetMat,"Logic_" + m_SDM->GetDetNameFromDetID(j));
        m_logicDetector -> SetVisAttributes(new G4VisAttributes(G4Colour::Brown()));
        
        new G4PVPlacement(0,G4ThreeVector(0,0,0),m_logicDetector,m_SDM->GetDetNameFromDetID(j),m_logicWorld,false,0,checkOverlaps);
    }
    return m_physDetector;
}
