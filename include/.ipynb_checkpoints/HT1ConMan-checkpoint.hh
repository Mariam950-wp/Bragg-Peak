#ifndef HT1CONMAN_h
#define HT1CONMAN_h 1

////////////////////////////////////////////////////////////////////////////////
//   HT1ConMan.hh for hadron theraphy                                         //
//                                                                            //
//              - 20. Sep. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <vector>
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include "G4VModularPhysicsList.hh"

#include "HT1PatMap.hh"
#include "HT1SphDetMap.hh"
#include "HT1PhyList.hh"

class HT1ConMan
{
    public:
	// ructors and Destructor
	HT1ConMan();
	HT1ConMan(const G4String fileName);
	~HT1ConMan();

	// Initialize
	void Init(const G4String fileName = "");

    // Load configuration file
	bool Load(const G4String fileName);

    // print config file
    void PrintConfiguration();
    
	//Computing
	void SetUseMTD( G4bool   useMTD){m_UseMTD = useMTD;};
	G4bool   GetUseMTD(){return m_UseMTD;};
    
	void SetNofTRD( G4int    nofTRD){m_NofTRD = nofTRD;};
	G4int    GetNofTRD(){return m_NofTRD;};

	// physics list
	void SetHTPhys(G4VModularPhysicsList* PhyList){m_PhyList = PhyList;};
	G4VModularPhysicsList*	GetHTPhys(){return m_PhyList;};

    //reaction area
    void SetInteXL(G4double ractXL){m_InterXL = ractXL * mm;};
	G4double GetInteXL(){return m_InterXL;};
    
	void SetInteYL(G4double ractYL){m_InterYL = ractYL * mm;};
	G4double GetInteYL(){return m_InterYL;};
    
    void SetLabMat(G4String matname){m_SLabMat = matname;};
	G4String GetLabMat(){return m_SLabMat;};   
    
	//Beam
	void SetBeamKE(G4double beamKE){m_BeamKE = beamKE * MeV;};
	G4double GetBeamKE(){return m_BeamKE;};
    
    void SetBeamOR(G4double beamOR){m_BeamOR = beamOR * mm;};
	G4double GetBeamOR(){return m_BeamOR;};
    
	void SetParName(G4String parname){m_SParName = parname;};
	G4String GetParName(){return m_SParName;};
    
    // patient
    void SetCollToPat(G4double Colltopat){m_ColltoPat = Colltopat * mm;};
	G4double GetCollToPat(){return m_ColltoPat;};
    
    void SetPatiZL(G4double patiZL){m_PatiZL = patiZL * mm;};
	G4double GetPatiZL(){return m_PatiZL;};
    
    void SetPatiDV(G4int patiDV){m_PatiDV = patiDV;};
	G4int GetPatiDV(){return m_PatiDV;};
    
    void SetPatiMat(G4String matname){m_SPatiMat = matname;};
	G4String GetPatiMat(){return m_SPatiMat;};

    // spherical detector
    void SetSpDeIR(G4double SpdeIR){m_SpDeIR = SpdeIR;};
	G4double GetSpDeIR(){return m_SpDeIR;};
    
    void SetPolaDV(G4int polaDV){m_PolaDV = polaDV;};
	G4int GetPolaDV(){return m_PolaDV;};
    
    void SetAzziDV(G4int azziDV){m_AzziDV = azziDV;};
	G4int GetAzziDV(){return m_AzziDV;};
    
    // Maps
	void SetPatiMap(G4int PatiDV){m_PM = new HT1PatMap(PatiDV);};
    HT1PatMap* GetPatiMap(){return m_PM;};
    
    void SetPsDeMap(G4int azziDV){m_SDM = new HT1SphDetMap(azziDV);};
    HT1SphDetMap* GetSpDeMap(){return m_SDM;};
    
	// Does the line start with specific word?
	bool StartsWith(const std::string& line, const std::string& token);

  private:
	void SetDefault();

	// Computing resource
    G4bool   m_UseMTD;
    G4int    m_NofTRD;

	// physics list	
	G4VModularPhysicsList* m_PhyList;
    
    //reaction area
    G4double m_InterXL;
    G4double m_InterYL;
    G4String m_SLabMat;

    // Beam 
    G4double m_BeamKE;
    G4double m_BeamOR;
	G4String m_SParName;
    
    // patient 
    G4double m_ColltoPat;
    G4double m_PatiZL;
    G4int    m_PatiDV;
    G4String m_SPatiMat;
    
    // spherical detector
    G4double m_SpDeIR;
    G4int    m_PolaDV = 0;
    G4int    m_AzziDV = 0;
    
     // detector maps
    HT1PatMap*    m_PM;
    HT1SphDetMap* m_SDM;
};

#endif
