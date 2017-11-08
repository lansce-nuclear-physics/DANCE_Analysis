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
#define CoincidenceWindow 1000 //in ns

//Status indicators
#define EventLimit 4294967295  //Event limit to shut off the unpacker (2^32 -1)
#define ProgressInterval 1000000   //Progress bar incriments


#define TMatrixFile "./Config/TMatrix.txt"  // File that includes ID's of crystals for which the time deviations histograms will be created
#define	DanceMapFile "./Config/DanceMap.txt"  // File that includes ID's of crystals for which the time deviations histograms will be created


#define FITTIMEDEV 0  //this turns on finding time deviations

#define Blocking_Time 2500 //this is the minimum time (in ns) between crystal hits that will be valid in the unpacker


#define HAVE_Threshold 1 //1 is yes 0 is no
#define Energy_Threshold 0.15 //MeV


//PI Gates
#define GAMMAGATE "Gamma.dat"
#define ALPHAGATE "Alpha.dat"


#endif
