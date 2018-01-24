//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  sort_functions.cpp     *// 
//*  Last Edit: 01/23/18    *//  
//***************************//

//File includes
#include "sort_functions.h"

//C/C++ includes
#include <iostream>

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

/* STAGE 0 SORTING */
