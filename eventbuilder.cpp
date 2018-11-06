//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  eventbuilder.cpp       *// 
//*  Last Edit: 11/06/18    *//  
//***************************//

//File includes
#include "analyzer.h"
#include "eventbuilder.h"

//#define Eventbuilder_Verbose

int build_events(deque<DEVT_BANK> &datadeque, bool &event_building_active, Input_Parameters input_params, uint32_t &Events_Sent) {
  
#ifdef Eventbuilder_Verbose
  cout<<"About to event build deque size: "<<datadeque.size()<<endl;
#endif
  
  vector<DEVT_BANK> eventvector;   //Vector to store events for analysis
  
  //Eventbuild
  while(true) {
		  
    //check to see if the buffer is longer than the length specificed in global.h
    if((datadeque[datadeque.size()-1].TOF - datadeque[0].TOF) >= (double)1000000000.0*input_params.Buffer_Depth) {
		    
      //we have started to build
      event_building_active=true;
		    
      //clear the event vector
      eventvector.clear();
		    
      //first timestamp in the deque
      double first_entry_time = datadeque[0].TOF;  //start of the event in time
      eventvector.push_back(datadeque[0]); //put the first event in the events vector
      datadeque.pop_front();  //remove the first entry in the deque
      Events_Sent++;
      
      bool event_build =true;  //bool to do eventbuilding     
		    
      while(event_build && datadeque.size()>0) {
	if(datadeque[0].TOF < (first_entry_time + input_params.Coincidence_Window)) {
	  eventvector.push_back(datadeque[0]); //put the first event in the events vector
	  datadeque.pop_front();  //remove the first entry in the deque     
	  Events_Sent++;
	}
	else {
	  event_build = false;
	  break;
	}
      }
			
      if(eventvector.size()>0) {
#ifdef Eventbuilder_Verbose
	cout<<"Processing Event with Size: "<<eventvector.size()<<"  " <<datadeque.size()<<" Entries in the deque"<<endl;
#endif
	//Send it to the analyzer
	Analyze_Data(eventvector, input_params);
      }
      if(datadeque.size()==0) {
	event_build=false;
	break;
      }
      
    }
    else {
#ifdef Eventbuilder_Verbose
      cout<<"event build complete: "<<datadeque.size()<<endl;
#endif
      break;
    }
  }


  return 0;

}
