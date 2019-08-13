/******************************************************************************/
/* File:             imputation.cpp                                           */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains general functions for imputation.                                 */
/*                                                                            */
/******************************************************************************/

#include "imputeDepth.h"

int selectPointCoords(int* totalCoords, int nTotalCoords, int d,
                      int** points, int* nPoints, int*** coords, int** nCoords,
                      int** coordsRaw){
  // Sort the total missing coordinates
  int* tmpTotalCoords = new int[nTotalCoords];
  memcpy(tmpTotalCoords, totalCoords, nTotalCoords * sizeof(int));
  quick_sort(tmpTotalCoords, 0, nTotalCoords - 1);
  // Obtain the number of points containing missing coordinates
  *nPoints = 0;
  int lastPointNumber = -1;
  for (int i = 0; i < nTotalCoords; i++){
    int curPointNumber = tmpTotalCoords[i] / d;
    if (curPointNumber > lastPointNumber){
      lastPointNumber = curPointNumber;
      (*nPoints)++;
    }
  }
  // Allocate the structures
  *points = new int[*nPoints]; // indices of points containing missingness
  *coords = new int*[*nPoints]; // double-indexed array of missing coordinates
  *coordsRaw = new int[nTotalCoords]; // ... and its raw memory
  *nCoords = new int[*nPoints]; // number of missing coordinates for each point
  // Fill the structures
  int curIndex = 0; // current local (among missing) index of a point
  // Fill the first point
  (*points)[0] = tmpTotalCoords[0] / d;
  (*coords)[0] = (*coordsRaw);
  (*nCoords)[0] = 0;
  // Fill all the indices
  for (int i = 0; i < nTotalCoords; i++){
    int curPointNumber = tmpTotalCoords[i] / d;
    // Update local index if needed
    if (curPointNumber > (*points)[curIndex]){
      curIndex++;
      (*points)[curIndex] = curPointNumber;
      (*coords)[curIndex] = &(*coordsRaw)[i];
      (*nCoords)[curIndex] = 0;
    }
    (*nCoords)[curIndex]++; // increase number of coordinates for current point
    (*coordsRaw)[i] = tmpTotalCoords[i] % d; // set current coordinate
  }
  // DEBUG: Output the arrays
  // for (int i = 0; i < *nPoints; i++){
  //   cout << (*points)[i] << " (" << (*nCoords)[i] << "): ";
  //     for (int j = 0; j < (*nCoords)[i]; j++){
  //       cout << (*coords)[i][j] << " ";
  //       }
  //   cout << endl;
  // }
  // DEBUG (end)
  // Release memory
  delete[] tmpTotalCoords;
  return 0;
}

void imputeAllZonoid(double** xx, int n, int d,
                     double** mm, int m,
                     int* points, int nPoints, int** coords, int* nCoords,
                     double* imputed){
  double* tmpPoint = new double[d]; // a contemporary point for imputation
  // Loop through points with missing coordinates
  int totImpIndex = 0; // global index in the output array 'imputed'
  for (int i = 0; i < nPoints; i++){
    LinearConditionalZonoid(mm[points[i]], d, xx, n, coords[i], nCoords[i],
                            tmpPoint); // optimize
    // Copy the missing coordinates into the output array 'imputed'
    for (int j = 0; j < nCoords[i]; j++){
      imputed[totImpIndex++] = tmpPoint[coords[i][j]];
    }
  }
  // Release memory
  delete[] tmpPoint;
}
