#ifndef GLOBAL_H
#define GLOBAL_H

/* THINGS THAT CONTROLL THE UNPACKER */

/* VERBOSITY */

//#define Unpacker_Verbose
//#define CheckTheDeque
//#define Eventbuilder_Verbose


/*size of the DEVT array in the unpacker*/
#define MaxDEVTArrSize 1000000 

/*size of the block buffer. The unpacker waits this many
  events before looking at the deque and time sorting.  This 
  number has huge efficiency implications so dont change it if
  you dont understand it */
#define BlockBufferSize 350000


/*Time depth (in seconds) of the buffer.  Seems like most boards
  report data at least every 30 seconds.  While you cant know if you 
  will fail the program will tell you if you did. */
#define BufferDepth 30 

//Coincidence window for the eventbuilder
#define CoincidenceWindow 500 //in ns

//Status indicators
#define EventLimit 4294967295  //Event limit to shut off the unpacker (2^32 -1)
//#define EventLimit 50000  //Event limit to shut off the unpacker 
#define ProgressInterval 1000000   //Progress bar incriments


#define TMatrixFile "./Config/TMatrix.txt"  // File that includes ID's of crystals for which the time deviations histograms will be created
#define	DanceMapFile "./Config/DanceMap.txt"  // File that includes ID's of crystals for which the time deviations histograms will be created


#define FITTIMEDEV 0  //this turns on finding time deviations

#define Blocking_Time 0.0 //this is the minimum time (in ns) between crystal hits that will be valid in the unpacker

#define READ_BINARY 0  //This is basically the stage 0 or 1 flag.  If we read binary then 
#define WRITE_BINARY 1 //this will spit out the binary file of timesorted DEVT_BANK_wWF objects for easy resorting in staged analysis

#define HAVE_Threshold 0 //1 is yes 0 is no
#define Energy_Threshold 0.15 //MeV


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
