#include "imputeDepth.h"

double LocalDepth(double* x, double** xx, int n, int d,
                  double locPar, DEPTHNOTION notion){
  // Create doubled sample
  double* doubledSmplRaw = new double[2 * n * d];
  double** doubledSmpl = asMatrix(doubledSmplRaw, 2 * n, d);
  memcpy(doubledSmplRaw, xx[0], n * d * sizeof(double));
  for (int i = 0; i < n; i++){
    for (int j = 0; j < d; j++){
      doubledSmpl[n + i][j] = 2 * x[j] - xx[i][j];
    }
  }
  // Calculate depths
  SortIndex* depths = new SortIndex[n];
  double* tmpPoint = new double[d];
  for (int i = 0; i < n; i++){
    // Prepare the point
    memcpy(tmpPoint, doubledSmpl[i], d * sizeof(double));
    // Calculate its depth
    double tmpDepth = -1;
    switch(notion){
    case HALFSPACE: tmpDepth = HD_Rec(tmpPoint, doubledSmpl, 2 * n, d);
      break;
    case ZONOID: tmpDepth = ZonoidDepth(tmpPoint, doubledSmpl, 2 * n, d);
      break;
    case MAHALANOBIS: MahalanobisDepth(doubledSmpl, &tmpPoint, d, 2 * n, 1, 1,
                                       &tmpDepth);
      break;
    }
    // Assign depth value
    depths[i].index = i;
    //if (i < n){
    depths[i].value = tmpDepth;
    //}else{
    //  depths[i].value = tmpDepth + 1.1;
    //}
  }
  // Sort depths
  quick_sort(depths, 0, n - 1);
  // Create shortened sample
  int nPoints = n * locPar;
  for (int i = 0; i < nPoints; i++){
    memcpy(doubledSmpl[i], xx[depths[i].index], d * sizeof(double));
  }
  // Calculate the depth
  double tmpDepth = -1;
  switch(notion){
  case HALFSPACE: tmpDepth = HD_Rec(x, doubledSmpl, nPoints, d);
    break;
  case ZONOID: tmpDepth = ZonoidDepth(x, doubledSmpl, nPoints, d);
    break;
  case MAHALANOBIS: MahalanobisDepth(doubledSmpl, &x, d, nPoints, 1, 1,
                                     &tmpDepth);
    break;
  }
  // Release memory
  delete[] tmpPoint;
  delete[] depths;
  delete[] doubledSmpl;
  delete[] doubledSmplRaw;
  return tmpDepth;
}
