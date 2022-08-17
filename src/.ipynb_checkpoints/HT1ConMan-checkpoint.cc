////////////////////////////////////////////////////////////////////////////////
//   HT1ConMan.cc for hadron tharaphy                                         //
//                                                                            //
//              - 20. Sep. 2020. Mariam Abuladze (mariam.abuladze@kiu.edu.ge) //
////////////////////////////////////////////////////////////////////////////////

#include "QGSP_BIC.hh"

#include "HT1ConMan.hh"

//////////////////////////////////////////////////
//   Constructors and destructor                //
//////////////////////////////////////////////////
HT1ConMan::HT1ConMan()
{
	Init();
}


HT1ConMan::HT1ConMan(const G4String fileName)
{
	Init(fileName);
}


HT1ConMan::~HT1ConMan()
{
}


//////////////////////////////////////////////////
//   Initialization function                    //
//////////////////////////////////////////////////
void HT1ConMan::Init(const G4String fileName)
{
	if ( !fileName ) SetDefault();
	else if ( !Load(fileName) ) SetDefault();
    
    // initialize detector maps
	SetPatiMap(m_PatiDV);
    SetPsDeMap(m_AzziDV);
    
    PrintConfiguration();
}

void HT1ConMan::SetDefault()
{
	// set default values to config features
	// Computing resource
	m_UseMTD		= true;
	m_NofTRD		= 7;

	// physics list
	m_PhyList = new QGSP_BIC;

	//reaction area
	m_InterXL		= 150.5 *  mm;
	m_InterYL		= 150.5 *  mm;
	m_SLabMat		= "G4_AIR";

	// Beam 
	m_BeamKE		= 230.0 * MeV;
	m_BeamOR		= 3.0		*  mm;
	m_SParName	= "proton";

	// patient 
	m_ColltoPat	= 200.0 *  mm;
	m_PatiZL		= 500.0 *  mm;
	m_PatiDV		= 100;
	m_SPatiMat	= "G4_WATER";
    
    // spherical detector 
    m_SpDeIR = 5000.0 * mm;
    m_PolaDV = 20;
    m_AzziDV = 10;
    
    // detector map
	m_PM  = new HT1PatMap(100);
    m_SDM = new HT1SphDetMap(10);
}


//////////////////////////////////////////////////
//   Load configuration file                    //
//////////////////////////////////////////////////
bool HT1ConMan::Load(const G4String fileName)
{
	// Open file
	std::ifstream file(fileName, std::ifstream::in);
	if ( !file.is_open() ) return false;

	// Read line by line
	std::string line;
	char hfile[1000];
	while ( std::getline(file, line) )
	{
		if ( StartsWith(line, "#") ) continue;

		// Read configurations: Computing
		if (StartsWith(line, "USEMULTITHREADS")){
			sscanf(line.data(), "USEMULTITHREADS %s", hfile);
			if      ( hfile == "ON"  ) SetUseMTD(true);
			else if ( hfile == "OFF" ) SetUseMTD(false);
            
		}else if (StartsWith(line, "NTHREADS")){
			sscanf(line.data(), "NTHREADS %s", hfile);
			SetNofTRD(std::stoi(hfile));
		
        // physics list 
		}else if (StartsWith(line, "USEHTPHYSICS")){
			sscanf(line.data(), "USEHTPHYSICS %s", hfile);
			if      ( hfile == "ON"  ) SetHTPhys(new HT1PhyList());
			else if ( hfile == "OFF" ) SetHTPhys(new QGSP_BIC);
		
        // lab parameters 
		}else if (StartsWith(line, "INTERAREXLEN")){
			sscanf(line.data(), "INTERAREXLEN %s", hfile);
            SetInteXL(std::stod(hfile));
		
		}else if (StartsWith(line, "INTERAREYLEN")){
			sscanf(line.data(), "INTERAREYLEN %s", hfile);
			SetInteYL(std::stod(hfile));
		
		}else if (StartsWith(line, "LABMATERIAL")){
			sscanf(line.data(), "LABMATERIAL %s", hfile);
			SetLabMat(hfile);
		
        // beam 
		}else if (StartsWith(line, "BEAMKINETICENERGY")){
			sscanf(line.data(), "BEAMKINETICENERGY %s", hfile);
			SetBeamKE(std::stod(hfile));
		
		}else if (StartsWith(line, "BEAMRADIUS")){
			sscanf(line.data(), "BEAMRADIUS %s", hfile);
			SetBeamOR(std::stod(hfile));
		
		}else if (StartsWith(line, "PARTICLENAME")){
			sscanf(line.data(), "PARTICLENAME %s", hfile);
			SetParName(hfile);
		
        // patient
		}else if (StartsWith(line, "CENTRETOPATIENT")){
			sscanf(line.data(), "CENTRETOPATIENT %s", hfile);
			SetCollToPat(std::stod(hfile));
		
		}else if (StartsWith(line, "PATIENTZLEN")){
			sscanf(line.data(), "PATIENTZLEN %s", hfile);
			SetPatiZL(std::stod(hfile));
		
		}else if (StartsWith(line, "DEVIDEPATGEOM")){
			sscanf(line.data(), "DEVIDEPATGEOM %s", hfile);
			SetPatiDV(std::stoi(hfile));
		
		}else if (StartsWith(line, "PATIENTMAT")){
			sscanf(line.data(), "PATIENTMAT %s", hfile);
			SetPatiMat(hfile);
		
        // spherical detector 
        }else if (StartsWith(line, "SPHERICALDETECTORRADIUS")){
			sscanf(line.data(), "SPHERICALDETECTORRADIUS %s", hfile);
			SetSpDeIR(std::stod(hfile));
        
        }else if (StartsWith(line, "POLARDIVISIONS")){
			sscanf(line.data(), "POLARDIVISIONS %s", hfile);
			SetPolaDV(std::stoi(hfile));
            
        }else if (StartsWith(line, "AZZIMUTDIVISIONS")){
			sscanf(line.data(), "AZZIMUTDIVISIONS %s", hfile);
			SetAzziDV(std::stoi(hfile));
            
        }
	}

	return true;
}

//////////////////////////////////////////////////
//   Print Configuration ...                    //
//////////////////////////////////////////////////
void HT1ConMan::PrintConfiguration()
{
    // Computing
    std::cout << "*********************Configuration*********************"<< std::endl;
    std::cout << "HT1::Load() => USEMULTITHREADS "         << m_UseMTD    << std::endl;
    std::cout << "HT1::Load() => NTHREADS "                << m_NofTRD    << std::endl;
    std::cout << "HT1::Load() => INTERAREXLEN "            << m_InterXL   << std::endl;
    std::cout << "HT1::Load() => INTERAREYLEN "            << m_InterYL   << std::endl;
    std::cout << "HT1::Load() => BEAMKINETICENERGY "       << m_BeamKE    << std::endl;
    std::cout << "HT1::Load() => BEAMRADIUS "              << m_BeamOR    << std::endl;
    std::cout << "HT1::Load() => PARTICLENAME "            << m_SParName  << std::endl;
    std::cout << "HT1::Load() => CENTRETOPATIENT "         << m_ColltoPat << std::endl;
    std::cout << "HT1::Load() => PATIENTZLEN "             << m_PatiZL    << std::endl;
    std::cout << "HT1::Load() => DEVIDEPATGEOM "           << m_PatiDV    << std::endl;
    std::cout << "HT1::Load() => PATIENTMAT "              << m_SPatiMat  << std::endl;
    std::cout << "HT1::Load() => SPHERICALDETECTORRADIUS " << m_SpDeIR    << std::endl;
    std::cout << "HT1::Load() => POLARDIVISIONS "          << m_PolaDV    << std::endl;
    std::cout << "HT1::Load() => AZZIMUTDIVISIONS "        << m_AzziDV    << std::endl;
    std::cout << "******************EndOfConfiguration*******************"<< std::endl;
}

//////////////////////////////////////////////////
//   Test whether a line starts with ...        //
//////////////////////////////////////////////////
bool HT1ConMan::StartsWith(const std::string& text, const std::string& token)
{
	if ( text.length() < token.length() ) return false;
	return ( text.compare(0, token.length(), token) == 0 );
}
