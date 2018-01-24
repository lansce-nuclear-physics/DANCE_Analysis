//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  structures.h           *// 
//*  Last Edit: 01/23/18    *//  
//***************************//

#ifndef STRUCTURES_H
#define STRUCTURES_H

//File includes
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

//Digitizer Header (NOT USED YET)
struct V1730_Header_t {
  uint32_t dataword_1;
  uint32_t dataword_2;
  uint32_t dataword_3;
  uint32_t dataword_4;
};

//Channel Aggregate Header (NOT USED YET)
struct V1730_ChAgg_Header_t {
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
  uint64_t timestamp;        // timestamp
  uint16_t Ns;               // number of samples in waveform
  uint16_t Ifast;            // short integral
  uint16_t Islow;            // long integral
  uint8_t board;             // board number
  uint8_t channel;           // channel number
  double TOF;                // Time of Flight in ns
  double Eslow;              // Calibrated slow integral
  double Efast;              // Calibrated fast integral
  uint16_t ID;               // Extracted ID using DANCE Map
  uint8_t Valid;             // Valid Flag
} DEVT_BANK;


typedef struct {
  uint16_t Ifast;            // short integral
  uint16_t Islow;            // long integral
  double TOF;                // Time-Of-Flight in ns
  uint8_t ID;                // ID from DANCE Map 0 to 161 are dance //241 to 244 are beam monitors //200 is t0
} DEVT_STAGE1;


// DANCE event
typedef struct{
  double En;                 //Neutron energy from TOF
  double En_corr;            //Neutron energy from TOF corrected for moderator function
  double tof[162];           //Neutron Time-Of-Flight from each crystal
  double tof_corr[162];      //Neutron Time-Of-Flight from each crystal corrected for moderator function
  uint16_t Crystal_mult;     //DANCE crystal multiplicity 
  uint16_t Cluster_mult;     //DANCE cluster multiplicity 
  uint16_t Crystal_ID[162];  //DANCE crystal ID
  uint16_t Cluster_ID[162];  //DANCE cluster ID
  uint16_t Ifast[162];       //Uncalibrated fast integral
  uint16_t Islow[162];       //Uncalibrated slow integral
  double Ecluster[162];      //Calibrated cluster energy (From sum of crystal ESlow in the cluster)
  double Ecrystal[162];      //Calibrated crystal energy from ESlow
  double ESum;               //Total DANCE Event energy from sum of crystal ESlow in the event 
  uint16_t Valid;            //Valid event flag
} DANCE_Event;


// U235 Beam Monitor event
typedef struct{
  double En;                 //Neutron energy from TOF
  double En_corr;            //Neutron energy from TOF corrected for moderator function
  double tof;                //Neutron Time-Of-Flight 
  double tof_corr;           //Neutron Time-Of-Flight corrected for moderator function
  uint16_t Ifast;            //Uncalibrated fast integral  
  uint16_t Islow;            //Uncalibrated slow integral 
  uint16_t Valid;            //Valid event flag
} U235_Event;

// He3 Beam Monitor event
typedef struct{
  double En;                 //Neutron energy from TOF
  double En_corr;            //Neutron energy from TOF corrected for moderator function
  double tof;                //Neutron Time-Of-Flight 
  double tof_corr;           //Neutron Time-Of-Flight corrected for moderator function
  uint16_t Ifast;            //Uncalibrated fast integral  
  uint16_t Islow;            //Uncalibrated slow integral 
  uint16_t Valid;            //Valid event flag
} He3_Event;

// Li6 Beam Monitor event
typedef struct{
  double En;                 //Neutron energy from TOF
  double En_corr;            //Neutron energy from TOF corrected for moderator function
  double tof;                //Neutron Time-Of-Flight 
  double tof_corr;           //Neutron Time-Of-Flight corrected for moderator function
  uint16_t Ifast;            //Uncalibrated fast integral  
  uint16_t Islow;            //Uncalibrated slow integral 
  uint16_t Valid;            //Valid event flag
} Li6_Event;


#endif
