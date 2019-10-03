//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  sort_functions.cpp     *// 
//*  Last Edit: 11/06/18    *//  
//***************************//

//File includes
#include "sort_functions.h"
#include "global.h"


// To heapify a subtree rooted with node i which is
// an index in arr[]. n is size of heap
void heapify(DEVT_BANK arr[], int n, int i) {
  int largest = i;  // Initialize largest as root
  int l = 2*i + 1;  // left = 2*i + 1
  int r = 2*i + 2;  // right = 2*i + 2
  
  // If left child is larger than root
  if (l < n && arr[l].TOF > arr[largest].TOF)
    largest = l;
  
  // If right child is larger than largest so far
  if (r < n && arr[r].TOF > arr[largest].TOF)
    largest = r;
  
  // If largest is not root
  if (largest != i) {
    std::swap(arr[i], arr[largest]);
    
    // Recursively heapify the affected sub-tree
    heapify(arr, n, largest);
  }
}

// main function to do heap sort
void heapSort(DEVT_BANK arr[], int n) {
  // Build heap (rearrange array)
  for (int i = n / 2 - 1; i >= 0; i--)
    heapify(arr, n, i);
  
  // One by one extract an element from heap
  for (int i=n-1; i>=0; i--) {
    // Move current root to end
    std::swap(arr[0], arr[i]);
    
    // call max heapify on the reduced heap
    heapify(arr, i, 0);
  }
}


int sort_array(DEVT_BANK db_arr[], deque<DEVT_BANK> &datadeque, uint32_t EVTS, Input_Parameters input_params, Analysis_Parameters *analysis_params) {

  
  int EVT_SORT=EVTS;
  int deque_size=datadeque.size();
  
#ifdef EventSort_Verbose
  cout<<endl<<"total events at start: "<<EVTS+datadeque.size()<<" deque size: "<<datadeque.size()<<endl;
#endif
  
  //the first time through we want to sort and push everything from the first "block" onto the buffer
  if(!analysis_params->first_sort) {
    
    //check to see if there are any timestamps that originate before the first one in the buffer.  If so event building is broken and the analysis is wrong
    if(analysis_params->event_building_active) {
      if(analysis_params->smallest_timestamp < datadeque[0].TOF) {
	cout<<RED<<"WARNING THE SMALLEST TIMESTAMP IS LOWER THAN THE SMALLEST ONE IN THE DEQUE!!"<<endl;
	cout<<"smallest: "<<analysis_params->smallest_timestamp<<"  smallest in deque: "<<datadeque[0].TOF<<" largest in deque: "<<datadeque[datadeque.size()-1].TOF<<endl;
	cout<<"Deque Depth: "<<(datadeque[datadeque.size()-1].TOF - datadeque[0].TOF)/(1.0e9)<<" seconds"<<endl;
	cout<<"Make the deque: "<<(datadeque[0].TOF-analysis_params->smallest_timestamp)/(1.0e9)<<" seconds deeper"<<endl;
	cout<<"Exiting"<<RESET<<endl;
	ofstream failfile;
	failfile.open("Failed_Analysis.txt", ios::out | ios::app);
	failfile << "Run: "<<input_params.RunNumber<<" Failed due to insufficient buffer depth...  Add: "<<(datadeque[0].TOF-analysis_params->smallest_timestamp)/(1.0e9)<<" seconds\n";
	failfile.close();
        return -1;
      }   
    }
    
    //location where the smallest timestamp sits in the sorted data
    int first_index=0;
    
    //find where the smallest time stamp sits in the already time sorted data
    for(int k=datadeque.size()-1; k>=0; k--) {
      if(analysis_params->smallest_timestamp >= datadeque[k].TOF) {
	first_index=k-1;
	break;
      }
    }
#ifdef EventSort_Verbose
    cout<<"start at "<<first_index<<" of "<<datadeque.size()<<endl;
#endif
    
    //place everything after that onto the unsorted array
    for(uint k=first_index; k<datadeque.size(); k++) {
#ifdef EventSort_Verbose
      cout<<"EVT_SORT: "<<EVT_SORT<<"  k: "<<k<<endl;
#endif
      db_arr[EVT_SORT]=datadeque[k];
      EVT_SORT++;

#ifdef CheckBufferDepth
      if(EVT_SORT > MaxDEVTArrSize) {
	DANCE_Error("Unpacker","Size of the buffer has exceeded MaxDEVTArrSize.  Change MaxDEVTArrSiz to a higher value");
	return -1;
      }
      else {
	double temp = (1.0*EVTS/(1.0*MaxDEVTArrSize));
	if(temp>analysis_params->max_buffer_utilization) {
	  analysis_params->max_buffer_utilization = temp;
	}
      }
#endif	
    }
    
    //remove the ones put onto the array so we dont double things
    for(int k=0; k<(deque_size-first_index); k++) {
      datadeque.pop_back();
    }
  }
  
#ifdef EventSort_Verbose
  cout<<"Deque Size after reduction: "<<datadeque.size()<<"  About to sort "<<EVT_SORT<<" events"<<endl;
#endif
  
  //the first sort is over
  analysis_params->first_sort=false;
  
  //sort the unsorted data
  heapSort(db_arr, EVT_SORT);
  
  //push the now sorted data onto the sorted buffer
  for(int j=0; j<EVT_SORT; j++) {
    //Update and push it onto the deque if valid
    datadeque.push_back(db_arr[j]);
  }
  
#ifdef CheckTheDeque
  cout<<"Checking deque"<<endl;
  for(int k=0; k<(int)datadeque.size()-1; k++) {
    if(datadeque[k+1].timestamp < datadeque[k].timestamp) {
      cout<<"problem with entry "<<k<<endl;
      return -1;
    }
#ifdef EventSort_Verbose
    cout<<k<<"  "<<datadeque[k].timestamp<<endl;
#endif
  }
#endif

  return 0;

}


