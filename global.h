//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  global.h               *// 
//*  Last Edit: 01/23/18    *//  
//***************************//

#ifndef GLOBAL_H
#define GLOBAL_H

//COLORS
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
#define CLEAR "\033[2J"  // clear screen escape code 

//Status indicators
#define EventLimit 4294967295  //Event limit to shut off the unpacker (2^32 -1)
#define ProgressInterval 1000000   //Progress bar incriments

//File Inputs
#define TMatrixFile "./Config/TMatrix.txt" 
#define	DanceMapFile "./Config/DanceMap.txt" 
#define DMatrixFile "./Config/DetectorMatrix.txt"

//File Input/Output
#define STAGE0_ROOT "./stage0_root"
#define STAGE1_ROOT "./stage1_root"
#define STAGE0_BIN "./stage0_bin"
#define STAGE1_BIN "./stage1_bin"

//PI Gates
#define GAMMAGATE "Gamma.dat"
#define ALPHAGATE "Alpha.dat"

//Physics Stuff
#define	DANCE_FlightPath 20.2407  // was 20.28437 meters
#define U235_FlightPath 22.8299  // was 22.8 //meters
#define He3_FlightPath 22.74 //meters
#define Li6_FlightPath 22.6121 // was 22.607 meters

#define DANCE_Delay 282.26  // ns (Was 523.0 in FARE, and 205.0 in last version)
#define U235_Delay 578.2201 //ns (Was 322.4 in FARE and 553.2 in previous version)
#define He3_Delay 359.3 //ns (Was -207 in previous version)
#define Li6_Delay 500.156 //ns (Was 376 in FARE, and 502.9 in previous version)

//Histogramming
#define	NeutronE_From 0.002 //Neutron energy from [eV]
#define	NeutronE_To 5e6 //Neutron energy to [eV]
#define	NeutronE_BinsPerDecade 500.0 // Number of neutron energy bins per decade

#define	GammaE_From 0.0 //Gamma energy [MeV] - low limit
#define	GammaE_To 20.0 //Gamma energy [MeV] - upper limit
#define	GammaE_NoOfBins 400.0 //Number of bins


#endif
