//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  structures.h           *// 
//*  Last Edit: 09/04/19    *//  
//***************************//

#ifndef STRUCTURES_H
#define STRUCTURES_H

//File includes
#include "global.h"
#include "message.h"

// C/C++ includes 
#include <stdint.h>  //uint16_t, uint32_t, uint64_t
#include <string>

/*size of the CEVT_BANK P Array*/
#define MaxHitsPerT0 500000  

//Number of Scalers to read
#define N_SCLR 35

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

struct Sclr_Totals_t {
  uint32_t totals[N_SCLR];
};

struct Sclr_Rates_t {
  uint32_t rates[N_SCLR];
};

/// Scalers
struct Sclr_t {
  uint32_t mmtm;
  char fName[4];      ///< bank name
  uint32_t fType;     ///< type of data (see midas.h TID_xxx)
  uint32_t fDataSize;  
  uint32_t scaler[24];
};


//User data masks
//FW VERSION WORD
#define MAJREV_MASK                         0x000000FF  //bits 0 to 7 inclusive
#define MINREV_MASK                         0x00003F00  //bits 8 to 13 inclusive
#define MODTYPE_MASK                        0x03FFC000  //bits 14 to 25 inclusive
#define BOARDID_MASK                        0xFC000000  //bits 26 to 31 inclusive

struct User_Data_t {
  //WORD 1 
  uint8_t fw_majrev;
  uint8_t fw_minrev;
  uint16_t modtype;
  uint8_t boardid;
  //WORD 2
  uint32_t user_extra;
};


//CAEN 2015 
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
  double timestamp;        // timestamp
  uint16_t Ns;               // number of samples in waveform
  uint16_t Ifast;            // short integral
  uint16_t Islow;            // long integral
  uint8_t board;             // board number
  uint8_t channel;           // channel number
  double TOF;                // Time of Flight in ns
  double TOF_Corr;           // Corrected Time of Flight in ns
  double Eslow;              // Calibrated slow integral
  double Efast;              // Calibrated fast integral
  uint16_t ID;               // Extracted ID using DANCE Map
  uint8_t Valid;             // Valid Flag
  uint8_t IsGamma;           // Gamma Flag
  uint8_t IsAlpha;           // Alpha Flag
  uint8_t InvalidReason;    // Reason the entry is invalid
  int pileup_detected;        //Pileup detection from integral ratios

} DEVT_BANK;


//Output of stage 0 bin
typedef struct {
  uint16_t Ifast;            // short integral
  uint16_t Islow;            // long integral
  double timestamp;                // Time-Of-Flight in ns
  uint8_t ID;                // ID from DANCE Map 0 to 161 are dance //241 to 244 are beam monitors //200 is t0
} DEVT_STAGE1;


// DANCE event
typedef struct{
  double En[162];                 //Neutron energy from TOF
  double En_corr[162];            //Neutron energy from TOF corrected for moderator function
  double timestamp[162];     //Timestamp of the crystal
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
  int pileup_detected;
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

// Background Beam Monitor event
typedef struct{
  double En;                 //Neutron energy from TOF
  double En_corr;            //Neutron energy from TOF corrected for moderator function
  double tof;                //Neutron Time-Of-Flight 
  double tof_corr;           //Neutron Time-Of-Flight corrected for moderator function
  uint16_t Ifast;            //Uncalibrated fast integral  
  uint16_t Islow;            //Uncalibrated slow integral 
  uint16_t Valid;            //Valid event flag
} Bkg_Event;


//Input parameters 
typedef struct{
  double Crystal_Blocking_Time;
  double DEvent_Blocking_Time;
  double Coincidence_Window;
  double Energy_Threshold; //MeV
  //Bools
  bool Read_Binary;
  bool Write_Binary;
  bool Read_Simulation;
  bool HAVE_Threshold;
  bool FitTimeDev;
  bool SingleSubrun;
  //Strings
  std::string DataFormat;
  std::string Simulation_File_Name;
  //QGated
  bool QGatedSpectra;
  int NQGates;
  double QGates[20];  //This gives 10 pairs
  //Isomers
  bool IsomerSpectra;
  int NIsomers;
  double IsomerPromptQGates[20];     //10 pairs
  int IsomerPromptMclGates[20];   //10 pairs
  double IsomerPromptTOFGates[20];   //10 pairs
  double IsomerDelayedQGates[20];    //10 pairs
  int IsomerDelayedMclGates[20]; //10 pairs
  double IsomerDelayedTOFGates[20];  //10 pairs
  bool JMOD_Background;
  int RunNumber;
  int SubRunNumber;
  int NumSubRun;
  //Ways to evaluate efficiency from TOF
  bool Evaluate_DeadTime;
  double Artificial_TOF;
  std::string DetectorLoad_FileName;
  std::string DetectorLoad_HistName;
  int Long_Gate;
  bool Use_Firmware_FineTime;  
  int Analysis_Stage;

  //Unpacker variables
  double Buffer_Depth;



} Input_Parameters;


//Analysis parameters
typedef struct{
  double last_timestamp[256];
  double last_Islow[256];
  double last_Eslow[256];
  double last_Efast[256];
  uint16_t last_Alpha[256];
  uint16_t last_Gamma[256];
  uint8_t last_InvalidReason[256];
  double last_valid_timestamp[256];
  double last_valid_Islow[256];
  double last_valid_Eslow[256];

  uint32_t entries_unpacked;            //Entries that have been unpacked
  uint32_t entries_awaiting_timesort;   //Entries in the devt array waiting for timesort
  uint32_t entries_written_to_binary;   //Entries written to binary
  uint32_t entries_processed;           //Entries that have been through eventbuilding (valid or invalid)
  uint32_t entries_invalid;             //Entries invalid before analysis
  uint32_t entries_built;               //Entries built into events
  uint32_t events_built;                //Events built

  uint32_t entries_analyzed;
  uint32_t DANCE_entries_analyzed;
  uint32_t T0_entries_analyzed;
  uint32_t He3_entries_analyzed;
  uint32_t Li6_entries_analyzed;
  uint32_t U235_entries_analyzed;
  uint32_t Bkg_entries_analyzed;
  uint32_t Unknown_entries;

  uint32_t events_analyzed;
  uint32_t DANCE_events_analyzed;
  uint32_t T0_events_analyzed;
  uint32_t He3_events_analyzed;
  uint32_t Li6_events_analyzed;
  uint32_t Bkg_events_analyzed;
  uint32_t U235_events_analyzed;

  double max_buffer_utilization;

  bool first_sort;
  bool event_building_active;  //this says whether or not we are event building yet
  double smallest_timestamp;
  double largest_timestamp;
  double largest_subrun_timestamp;
  double last_subrun_timestamp;
  double wf_integral;

} Analysis_Parameters;
  



#endif
