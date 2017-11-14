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
/* VERBOSITY */

//#define Unpacker_Verbose
//#define CheckTheDeque
//#define Eventbuilder_Verbose



//Coincidence window for the eventbuilder
//#define CoincidenceWindow 10 //in ns
#define CoincidenceWindow 500 //in ns

//Status indicators
#define EventLimit 4294967295  //Event limit to shut off the unpacker (2^32 -1)
//#define EventLimit 50000  //Event limit to shut off the unpacker 
#define ProgressInterval 1000000   //Progress bar incriments


#define FITTIMEDEV 0  //this turns on finding time deviations

#define Blocking_Time 0.0 //this is the minimum time (in ns) between crystal hits that will be valid in the unpacker

#define READ_BINARY 0  //This is basically the stage 0 or 1 flag.  If we read binary then 
#define WRITE_BINARY 1 //this will spit out the binary file of timesorted DEVT_BANK_wWF objects for easy resorting in staged analysis

#define HAVE_Threshold 0 //1 is yes 0 is no
#define Energy_Threshold 0.15 //MeV

//File Input
#define TMatrixFile "./Config/TMatrix.txt"  // File that includes ID's of crystals for which the time deviations histograms will be created
#define	DanceMapFile "./Config/DanceMap.txt"  // File that includes ID's of crystals for which the time deviations histograms will be created

//File Output
#define STAGE0_ROOT "./stage0_root"
#define STAGE1_ROOT "./stage1_root"
#define STAGE0_BIN "./stage0_bin"
#define STAGE1_BIN "./stage1_bin"

//PI Gates
#define GAMMAGATE "Gamma.dat"
#define ALPHAGATE "Alpha.dat"

//Physics Stuff
#define	DANCE_FlightPath 20.28437  // meters
#define DANCE_Delay 523.0  // ns

//Histogramming
#define	NeutronE_From 0.02 //Neutron energy from [eV]
#define	NeutronE_To 1e6 //Neutron energy to [eV]
#define	NeutronE_BinsPerDecade 500.0 // Number of neutron energy bins per decade

#define	GammaE_From 0.0 //Gamma energy [MeV] - low limit
#define	GammaE_To 20.0 //Gamma energy [MeV] - upper limit
#define	GammaE_NoOfBins 400.0 //Number of bins


#endif
