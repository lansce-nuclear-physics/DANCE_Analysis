#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "global.h"

// C/C++ includes 
#include <stdint.h>  //uint16_t, uint32_t, uint64_t

/*size of the CEVT_BANK P Array*/
#define MaxHitsPerT0 500000  

//MIDAS Event Structures

/// Event header
struct EventHeader_t {
  uint16_t fEventId;      ///< event id
  uint16_t fTriggerMask;  ///< event trigger mask
  uint32_t fSerialNumber; ///< event serial number
  uint32_t fTimeStamp;    ///< event timestamp in seconds
  uint32_t fDataSize;     ///< event size in bytes
};

/// Bank header
struct BankHeader_t {
  uint32_t fDataSize; 
  uint32_t fFlags; 
};

/// 16-bit data bank
struct Bank_t {
  char fName[4];      ///< bank name
  uint16_t fType;     ///< type of data (see midas.h TID_xxx)
  uint16_t fDataSize;
};

/// 32-bit data bank
struct Bank32_t {
  char fName[4];      ///< bank name
  uint32_t fType;     ///< type of data (see midas.h TID_xxx)
  uint32_t fDataSize;
};

/// Scalers
struct Sclr_t {
  uint32_t mmtm;
  char fName[4];      ///< bank name
  uint32_t fType;     ///< type of data (see midas.h TID_xxx)
  uint32_t fDataSize;  
  uint32_t scaler[24];
};

//Digitizer Header
struct V1730_Header_t {
  uint32_t dataword_1;
  uint32_t dataword_2;
  uint32_t dataword_3;
  uint32_t dataword_4;
};

//Channel Aggregate Header
struct V1730_ChAgg_Header_t {
  uint32_t dataword_1;
  uint32_t dataword_2;
};

//Channel Aggregate Footer
struct V1730_ChAgg_Footer_t {
  uint32_t dataword_1;
  uint32_t dataword_2;  
};
 
typedef struct {
  uint64_t  position;
  uint32_t     extras;
  uint16_t      width;
  uint16_t      detector_id;
  uint16_t      integral[2];
} CEVT_BANK;

typedef struct{
	uint32_t N;
	CEVT_BANK P[MaxHitsPerT0];
        short wavelets[MaxHitsPerT0][1024]; // this ties to the wavelets with P
        short imported_peaks[256][16384]; // this is actually supported channels / supported length of PXXX bank
        short something[10000]; // this is actually supported channels / supported length of PXXX bank
} test_struct_cevt;


typedef struct {
    uint64_t timestamp;  // timestamp
    uint16_t Ns;         // number of samples in waveform
    uint16_t sgate;      // short gate
    uint16_t lgate;      // long gate
    uint16_t baseline;   // baseline
    uint8_t board;       // board number
    uint8_t channel;     // channel number
    //uint16_t wf[300];
    uint64_t N;
    uint64_t EVTS;
    
    
    double TOF;
    double Eslow;
    double Efast;
    int ID;
    int Valid;
} DEVT_BANK_wWF;

#endif
