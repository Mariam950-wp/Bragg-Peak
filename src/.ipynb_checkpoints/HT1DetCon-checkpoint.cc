////////////////////////////////////////////////////////////////////////////////
//   HT1DetCon.hh for hadron theraphy                                         //
//                                                                            //
//              - 20. Sep. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include "G4VisAttributes.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"

#include "HT1DetCon.hh"

#include "G4RunManager.hh"
#include "G4PVPlacement.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1DetCon::HT1DetCon(HT1ConMan* CM)
: G4VUserDetectorConstruction(),
  fScoringVolume(0), m_CM(CM)
{ 
    nist = G4NistManager::Instance();
    
    //reaction area
    m_InterXL = m_CM->GetInteXL();
    m_InterYL = m_CM->GetInteYL();
    m_WorldMat= nist->FindOrBuildMaterial(m_CM->GetLabMat());
      
    // number of detectors and detector map
    m_PatiDV    = m_CM->GetPatiDV();
    m_PM        = m_CM->GetPatiMap();
      
    // patient
    m_ColltoPat = m_CM->GetCollToPat();
		//m_ColltoPat = 0.0*mm;
    m_PatiZL    = m_CM->GetPatiZL();
      
    //m_PatiMat   = nist->FindOrBuildMaterial(m_CM->GetPatiMat());
    G4String name, symbol; 
    G4int ncomponents;
    G4Isotope* isoC12 = new G4Isotope("isoC12", 6,12,12.0*g/mole);
    G4Isotope* isoC13 = new G4Isotope("isoC13", 6,13,13.0*g/mole);
      
    G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C", 2);
    elC->AddIsotope( isoC12, 98.93*perCent);
    elC->AddIsotope( isoC13, 1.07*perCent);
      
	m_PatiMat = new G4Material(name="GRAPHITE_174", 1.74*g/cm3, ncomponents=1);
    m_PatiMat->AddElement(elC, 100.*perCent);
      
      
    // world dimessions 
    //m_WorldY = m_InterYL;
    //m_WorldX = m_ColltoPat + m_PatiZL + m_InterXL;
    //m_WorldZ = m_ColltoPat + m_PatiZL + m_InterXL + m_InterYL;
      
    m_WorldX = m_ColltoPat + m_PatiZL + m_InterXL + m_InterYL;
    m_WorldY = m_ColltoPat + m_PatiZL + m_InterXL + m_InterYL;
    m_WorldZ = m_ColltoPat + m_PatiZL + m_InterXL + m_InterYL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HT1DetCon::~HT1DetCon()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* HT1DetCon::Construct()
{     
    G4bool checkOverlaps = true;

    // World
    G4Box* solidWorld =    
        new G4Box("World",m_WorldX, m_WorldY, m_WorldZ); 
      
    G4LogicalVolume* logicWorld =                         
        new G4LogicalVolume(solidWorld,m_WorldMat,"World");    
                                   
    G4VPhysicalVolume* physWorld = 
        new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World",0,false,0,checkOverlaps);  
    
    // solid volumes 
    for (G4int i = 0; i<m_PatiDV; i++)
        m_Patient.push_back(new G4Box (m_PM->GetDetNameFromDetID(i),m_InterXL/2,m_InterYL/2,m_PatiZL/(2*m_PatiDV)));
    
    // logicals 
    for (G4int i = 0; i<m_PatiDV; i++)
        m_logicPatient.push_back(new G4LogicalVolume(m_Patient[i],m_PatiMat,"logic"+m_PM->GetDetNameFromDetID(i))); 
    
    // visattributes
    for (G4int i = 0; i<m_PatiDV; i++)
	m_logicPatient[i] -> SetVisAttributes(new G4VisAttributes(G4Colour::Cyan()));
        
    // position vectors
    for (G4int i = 0; i<m_PatiDV; i++)
        m_PatientPos.push_back(G4ThreeVector(0.0, 0.0, m_ColltoPat + (i+1.0/2)*m_PatiZL/m_PatiDV));
    
    // placement 
    for (G4int i = 0; i<m_PatiDV; i++)
        new G4PVPlacement(0,m_PatientPos[i],m_logicPatient[i],m_PM->GetDetNameFromDetID(i),logicWorld,false,0,checkOverlaps);
    
    //m_SpDet = new HT1SphDetCon(logicWorld, m_WorldZ, m_CM->GetSpDeMap());
    
    fScoringVolume = logicWorld;
    
    //always return the physical World
    return physWorld;
}

