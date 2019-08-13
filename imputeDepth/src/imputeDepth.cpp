/******************************************************************************/
/* File:             imputeDepth.cpp                                          */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains entry C++-functions for R-package imputeDepth.                    */
/*                                                                            */
/******************************************************************************/

#include "imputeDepth.h"

#define EOF (-1)

#ifdef __cplusplus
extern "C" {
#endif

/* Export functions --------------------------------------------------------- */
void MedianQs(double* x, int* n, double* med){
  med[0] = quick_select(x, n[0]);
}

void MedianHfsp(double* X, int* n, int* d, int* r, double* median){
  double** x = asMatrix(X, *n, *d);
  median_hfsp(x, *n, *d, *r, median);
  delete[] x;
}

void MedianProj(double* X, int* n, int* d, int* r, double* median){
  double** x = asMatrix(X, *n, *d);
  median_proj(x, *n, *d, *r, median);
  delete[] x;
}

void AnSvd(double* X, int* n, int* d, double* s, double* u, double* vt){
  an_svd(X, *n, *d, s, u, vt);
}

void GetW(double* X, int* n, int* d, double* mu, double* w){
  double** xx = asMatrix(X, *n, *d);
  get_w(xx, *n, *d, mu, w);
  delete[] xx;
}

void SvdDepthProj(double* X, int* n, int* d, int* r, int* typeCenter,
                  double* s, double* vt, double* e, double* c){
  double** x = asMatrix(X, *n, *d);
  svd_depth_proj(x, *n, *d, *r, *typeCenter, s, vt, e, c);
  delete[] x;
}

void SvdDepthHfsp(double* X, int* n, int* d, int* r,
                  int* typeDepth, int* typeCenter,
                  double* s, double* vt, double* e, double* c){
  double** x = asMatrix(X, *n, *d);
  svd_depth_hfsp(x, *n, *d, *r, *typeDepth, *typeCenter, s, vt, e, c);
  delete[] x;
}

void DepthsProj(double* X, int* n, int* d, int* r, double* depths){
  double** x = asMatrix(X, *n, *d);
  depths_proj(x, *n, *d, *r, depths);
  delete[] x;
}

void DepthsHfsp(double* X, int* n, int* d, int* r, double* depths){
  double** x = asMatrix(X, *n, *d);
  depthsn_hfsp(x, *n, *d, x, *n, *r, depths);
  delete[] x;
}

void DepthsSptl(double* X, int* n, int* d, int* sd, double* depths){
  double** x = asMatrix(X, *n, *d);
  depthsn_sptl(x, *n, *d, x, *n, *sd != 0, depths);
  delete[] x;
}

void BinSearchlRoutine(double* x, int* n, double* val, int* pos){
  *pos = bin_searchl_routine(x, *n, *val);
}

void BinSearchrRoutine(double* x, int* n, double* val, int* pos){
  *pos = bin_searchr_routine(x, *n, *val);
}

void KthElement(double* x, int* n, int* k, double* kthel){
  double* xCopy = new double[*n];
  memcpy(xCopy, x, *n * sizeof(double));
  *kthel = kth_smallest(xCopy, *n, *k);
  delete[] xCopy;
}

void ImputeOnceZonoid(double* X, int* n, int* d, int* coords, int* nCoords,
                      double* input, double* output){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Optimize
  LinearConditionalZonoid(input, *d, xx, *n, coords, *nCoords, output);
  // Release memory
  delete[] xx;
}

void ImputeOnceHalfspace(double* X, int* n, int* d, int* r,
                         int* coords, int* nCoords,
                         double* input, double* output){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Calculate projections
  double* dirsRaw;
  double* prjsRaw;
  get_prjs(xx, *n, *d, *r, &dirsRaw, &prjsRaw);
  double** dirs = asMatrix(dirsRaw, *r, *d);
  double** prjs = asMatrix(prjsRaw, *r, *n);
  for (int i = 0; i < *r; i++){
    quick_sort(prjs[i], 0, *n - 1);
  }
  // Optimize
  NelderMeadConditionalHfsp(depth_hfsp, input, *d, *r, *n, dirs, prjs,
                            coords, *nCoords, xx, output);
  // Release memory
  delete[] dirs;
  delete[] dirsRaw;
  delete[] prjs;
  delete[] prjsRaw;
  delete[] xx;
}

void ImputeOnceProjection(double* X, int* n, int* d, int* r,
                          int* coords, int* nCoords,
                          double* input, double* output){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Calculate projections
  double* dirsRaw;
  double* medians;
  double* MADs;
  get_MaMADs(xx, *n, *d, *r, &dirsRaw, &medians, &MADs);
  double** dirs = asMatrix(dirsRaw, *r, *d);
  // Optimize
  NelderMeadConditionalProj(depth_proj, input, *d, *r, dirs, medians, MADs,
                            coords, *nCoords, xx, output);
  // Release memory
  delete[] dirs;
  delete[] medians;
  delete[] MADs;
  delete[] xx;
}

void ImputeZonoid(double* X, int* n, int* d,
                  double* M, int* m, int* totalCoords, int* nTotalCoords,
                  double* imputed){
  // Double-index X and M
  double** xx = asMatrix(X, *n, *d);
  double** mm = asMatrix(M, *m, *d);
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Optimize
  imputeAllZonoid(xx, *n, *d, mm, *m, points, nPoints, coords, nCoords,
                  imputed);
  // Release memory
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] mm;
}

void ImputeZonoidIter(double* X, int* n, int* d,
                  int* totalCoords, int* nTotalCoords,
                  double* imputed){
  // Copy and double-index X
  double* tmpX = new double[*n * *d];
  memcpy(tmpX, X, *n * *d * sizeof(double));
  double** xx = asMatrix(tmpX, *n, *d);
  // Working space for imputed coordinates
  double* tmpImputed = new double[*nTotalCoords];
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Iterative imputation loop
  cout << "Zonoid iters: ";
  for (int i = 0; i < 25; i++){
    // Optimize
    imputeAllZonoid(xx, *n, *d, xx, *n, points, nPoints, coords, nCoords,
                    tmpImputed);
    // Calculate distance to the previous values
    double maxDist = 0;
    for (int j = 0; j < *nTotalCoords; j++){
      double curDist = fabs(tmpX[totalCoords[j]] - tmpImputed[j]);
      if (curDist > maxDist){
        maxDist = curDist;
      }
    }
    cout << i + 1 << "(" << maxDist << ") ";
    //Replace
    for (int j = 0; j < *nTotalCoords; j++){
      tmpX[totalCoords[j]] = tmpImputed[j];
    }
    // Break if converged
    if (maxDist < 0.001){
      break;
    }
  }
  cout << "." << endl;
  // Copy to the output
  memcpy(imputed, tmpImputed, *nTotalCoords * sizeof(double));
  // Release memory
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpX;
  delete[] tmpImputed;
}

void ImputeProjectionIter(double* X, int* n, int* d, int* r,
                          int* totalCoords, int* nTotalCoords,
                          double* imputed){
  // Copy and double-index X
  double* tmpX = new double[*n * *d];
  memcpy(tmpX, X, *n * *d * sizeof(double));
  double** xx = asMatrix(tmpX, *n, *d);
  // Working space for imputed coordinates
  double* tmpImputed = new double[*nTotalCoords];
  double* tmpPoint = new double[*d]; // currently imputed point
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  double* depths = new double[nPoints];
  // Generate random directions
  double* dirsRaw;
  double* medians;
  double* MADs;
  get_MaMADs(xx, *n, *d, *r, &dirsRaw, &medians, &MADs);
  double** dirs = asMatrix(dirsRaw, *r, *d);
  // Iterative imputation loop
  cout << "Projection iters: ";
  int numTries = 1;
  double** tmpPoints = new double*[numTries];
  for (int j = 0; j < numTries; j++){
    tmpPoints[j] = new double[*d];
  }
  double tmpBestDepth = -1.;
  int numTmps = 0;
  double prevMaxDist = -1.;
  for (int i = 0; i < *nTotalCoords; i++){
    tmpImputed[i] = tmpX[totalCoords[i]];
  }
  for (int i = 0; i < 30; i++){
    // Calculate medians and MADs
    // if (i < 10 || i % 10 == 0){
    if (i == 5){
    // if (0 == 1){
      delete[] dirsRaw;
      delete[] medians;
      delete[] MADs;
      delete[] dirs;
      get_MaMADs(xx, *n, *d, *r, &dirsRaw, &medians, &MADs);
      dirs = asMatrix(dirsRaw, *r, *d);
    }else{
      get_mads(xx, *n, *d, *r, dirs, medians, MADs);
    }
    for (int j = 0; j < nPoints; j++){
      depths[j] = depth_proj(xx[points[j]], *d, *r, dirs, medians, MADs);
      // cout << depths[j] << " ";
    }
    // cout << endl;
    // Optimize and save
    int curIndex = 0;
    for (int j = 0; j < nPoints; j++){
      tmpBestDepth = -1.;
      numTmps = 0;
      double optDepth = -1.;
      for (int k = 0; k < numTries; k++){
        optDepth = NelderMeadConditionalProj(depth_proj, xx[points[j]], *d, *r,
                                dirs, medians, MADs,
                                coords[j], nCoords[j], xx, tmpPoint);
        // DEBUG: output depths
        // cout << "(" << depths[j] << "->" << optDepth << ") ";
        // DEBUG (end)
        if (optDepth - tmpBestDepth >= -eps){
          if (optDepth - tmpBestDepth > eps){
            tmpBestDepth = optDepth;
            numTmps = 0;
          }
          memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
        }
      }
      // cout << "(" << points[j] << "): " << depths[j] << " -> " << optDepth << "; ";
      for (int k = 0; k < nCoords[j]; k++){
        // if (tmpBestDepth - depths[j] > eps){
          tmpImputed[curIndex] = 0;
          for (int l = 0; l < numTmps; l++){
            tmpImputed[curIndex] += tmpPoints[l][coords[j][k]];
          }
          tmpImputed[curIndex] /= (double)numTmps;
        // }
        curIndex++;
      }
      // if (optDepth > depths[j]){
      //   cout << "1 ";
      // }else{
      //   cout << "0 ";
      // }
    }
    // cout << endl;
    // Calculate distance to the previous values
    double maxDist = -2.;
    for (int j = 0; j < *nTotalCoords; j++){
      double curDist = fabs(tmpX[totalCoords[j]] - tmpImputed[j]);
      if (curDist > maxDist){
        maxDist = curDist;
      }
    }
    cout << i + 1 << "(" << maxDist << ") ";
    //Replace
    for (int j = 0; j < *nTotalCoords; j++){
      tmpX[totalCoords[j]] = tmpImputed[j];
    }
    // Break if converged
    // if (maxDist < 0.001 || fabs(maxDist - prevMaxDist) < 0.000001){
    if (maxDist < 0.001){
      break;
    }
    prevMaxDist = maxDist;
  }
  cout << "." << endl;
  // Copy to the output
  for (int i = 0; i < *nTotalCoords; i++){
    imputed[i] = tmpImputed[i];
  }
  // Release memory
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpX;
  delete[] tmpPoint;
  delete[] tmpImputed;
  delete[] dirsRaw;
  delete[] medians;
  delete[] MADs;
  delete[] dirs;
}

void ImputeHalfspaceIter(double* X, int* n, int* d, int* r,
                         int* totalCoords, int* nTotalCoords,
                         double* imputed){
  // Copy and double-index X
  double* tmpX = new double[*n * *d];
  memcpy(tmpX, X, *n * *d * sizeof(double));
  double** xx = asMatrix(tmpX, *n, *d);
  // Working space for imputed coordinates
  double* tmpImputed = new double[*nTotalCoords];
  double* tmpPoint = new double[*d]; // currently imputed point
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  double* depths = new double[nPoints];
  // Generate random directions
  double* dirsRaw;
  double* prjsRaw;
  get_prjs(xx, *n, *d, *r, &dirsRaw, &prjsRaw);
  double** dirs = asMatrix(dirsRaw, *r, *d);
  double** prjs = asMatrix(prjsRaw, *r, *n);
  // Iterative imputation loop
  cout << "Halfspace iters: ";
  int numTries = *d;
  double** tmpPoints = new double*[numTries];
  for (int j = 0; j < numTries; j++){
    tmpPoints[j] = new double[*d];
  }
  double tmpBestDepth = -1.;
  int numTmps = 0;
  double prevMaxDist = -1.;
  for (int i = 0; i < *nTotalCoords; i++){
    tmpImputed[i] = tmpX[totalCoords[i]];
  }
  for (int i = 0; i < 30; i++){
    // Calculate and sort projections
    //if (i < 10 || i % 10 == 0){
    if (i == 5){
    // if (1 == 0){
      delete[] dirsRaw;
      delete[] prjsRaw;
      delete[] dirs;
      delete[] prjs;
      get_prjs(xx, *n, *d, *r, &dirsRaw, &prjsRaw);
      dirs = asMatrix(dirsRaw, *r, *d);
      prjs = asMatrix(prjsRaw, *r, *n);
    }else{
      upd_prjs(xx, *n, *d, *r, dirs, prjs);
    }
    // DEBUG:
    // for (int i = 0; i < *r; i++){
    //   for (int j = 0; j < *d; j++){
    //     cout << prjs[i][j] << " ";
    //   }
    //   cout << endl;
    //}
    // DEBUG (end)
    for(int j = 0; j < *r; j++){
      quick_sort(prjs[j], 0, *n - 1);
    }
    for (int j = 0; j < nPoints; j++){
      depths[j] = depth_hfsp(xx[points[j]], *d, *r, *n, dirs, prjs);
      // cout << depths[j] << " ";
    }
    // cout << endl;
    // Optimize and save
    int curIndex = 0;
    for (int j = 0; j < nPoints; j++){
      tmpBestDepth = -1.;
      numTmps = 0;
      for (int k = 0; k < numTries; k++){
	      double optDepth = NelderMeadConditionalHfsp(depth_hfsp, xx[points[j]], *d, *r, *n,
	                                dirs, prjs,
	                                coords[j], nCoords[j], &xx[k * (*d + 1)], tmpPoint,
									                2., 4., 0.5, 0.5, 0.00001);
        // DEBUG: output depths
        // cout << "(" << depths[j] << "->" << optDepth << ") ";
        // DEBUG (end)
        if (optDepth - tmpBestDepth > -0.5 / (double)(*n)){
          if (optDepth - tmpBestDepth > 0.5 / (double)(*n)){
            tmpBestDepth = optDepth;
            numTmps = 0;
          }
          memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
        }
      }
      for (int k = 0; k < nCoords[j]; k++){
        if (tmpBestDepth - depths[j] > -0.5 / (double)(*n)){
          tmpImputed[curIndex] = 0;
          for (int l = 0; l < numTmps; l++){
            tmpImputed[curIndex] += tmpPoints[l][coords[j][k]];
          }
          tmpImputed[curIndex] /= (double)numTmps;
        }
        curIndex++;
      }
    }
    // Calculate distance to the previous values
    double maxDist = -2.;
    for (int j = 0; j < *nTotalCoords; j++){
      double curDist = fabs(tmpX[totalCoords[j]] - tmpImputed[j]);
      if (curDist > maxDist){
        maxDist = curDist;
      }
    }
    cout << i + 1 << "(" << maxDist << ") ";
    //Replace
    for (int j = 0; j < *nTotalCoords; j++){
      tmpX[totalCoords[j]] = tmpImputed[j];
    }
    // DEBUG: Output of the imputed values
    // for (int j = 0; j < *nTotalCoords; j++){
    //   cout << tmpX[totalCoords[j]] << " ";
    // }
    // cout << endl;
    // DEBUG (end)
    // Break if converged
    // if (maxDist < 0.001 || fabs(maxDist - prevMaxDist) < 0.000001){
    if (maxDist < 0.001){
      break;
    }
    prevMaxDist = maxDist;
  }
  cout << "." << endl;
  // Copy to the output
  for (int i = 0; i < *nTotalCoords; i++){
    imputed[i] = tmpImputed[i];
  }
  // Release memory
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpX;
  delete[] tmpPoint;
  delete[] tmpImputed;
  delete[] dirsRaw;
  delete[] dirs;
  delete[] prjsRaw;
  delete[] prjs;
  for (int i = 0; i < numTries; i++){
	  delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
}

void ImputeHalfspaceIterEx(double* X, int* n, int* d,
                           int* totalCoords, int* nTotalCoords,
                           double* imputed){
  // Copy and double-index X
  double* tmpX = new double[*n * *d];
  memcpy(tmpX, X, *n * *d * sizeof(double));
  double** xx = asMatrix(tmpX, *n, *d);
  // Working space for imputed coordinates
  double* tmpImputed = new double[*nTotalCoords];
  double* tmpPoint = new double[*d]; // currently imputed point
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  double* depths = new double[nPoints];
  // Iterative imputation loop
  cout << "Halfspace iters: ";
  int numTries = *d;
  double prevMaxDist = -1.;
  for (int i = 0; i < *nTotalCoords; i++){
    tmpImputed[i] = tmpX[totalCoords[i]];
  }
  double** tmpPoints = new double*[numTries];
  for (int j = 0; j < numTries; j++){
    tmpPoints[j] = new double[*d];
  }
  double tmpBestDepth = -1.;
  for (int i = 0; i < 25; i++){
    for (int j = 0; j < nPoints; j++){
      depths[j] = HD_Rec(xx[points[j]], xx, *n, *d);
    }
    // Optimize and save
    int curIndex = 0;
    for (int j = 0; j < nPoints; j++){
      tmpBestDepth = -1.;
      int numTmps = 0;
      for (int k = 0; k < numTries; k++){
        double optDepth = NelderMeadConditionalHfspEx(HD_Rec, xx[points[j]], xx, *n, *d,
                                                      coords[j], nCoords[j], &xx[k * (*d + 1)], tmpPoint,
                                                      2., 4., 0.5, 0.5, 0.00001);
        if (optDepth - tmpBestDepth > -0.5 / (double)(*n)){
          if (optDepth - tmpBestDepth > 0.5 / (double)(*n)){
            tmpBestDepth = optDepth;
            numTmps = 0;
          }
          memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
        }
      }

      // double optDepth = NelderMeadConditionalHfspEx(HD_Rec, xx[points[j]], xx, *n, *d,
      //                                        coords[j], nCoords[j], xx, tmpPoint);
      // cout << "(" << points[j] << "): " << depths[j] << " -> " << optDepth << "; ";
      for (int k = 0; k < nCoords[j]; k++){
        // if (optDepth - depths[j] > 0.5 / (double)(*n)){
        if (tmpBestDepth - depths[j] > -0.5 / (double)(*n)){
          tmpImputed[curIndex] = 0;
          for (int l = 0; l < numTmps; l++){
            tmpImputed[curIndex] += tmpPoints[l][coords[j][k]];
          }
          tmpImputed[curIndex] /= (double)numTmps;
        }
        curIndex++;
      }
        // cout << "1 ";
      // }else{
        // cout << "0 ";
      // }
    }
    // cout << endl;
    // Calculate distance to the previous values
    double maxDist = -2.;
    for (int j = 0; j < *nTotalCoords; j++){
      double curDist = fabs(tmpX[totalCoords[j]] - tmpImputed[j]);
      if (curDist > maxDist){
        maxDist = curDist;
      }
    }
    cout << i + 1 << "(" << maxDist << ") ";
    //Replace
    for (int j = 0; j < *nTotalCoords; j++){
      tmpX[totalCoords[j]] = tmpImputed[j];
    }
    // DEBUG: Output of the imputed values
    // for (int j = 0; j < *nTotalCoords; j++){
    //   cout << tmpX[totalCoords[j]] << " ";
    //}
    // cout << endl;
    // DEBUG (end)
    // Break if converged
    // if (maxDist < 0.001 || fabs(maxDist - prevMaxDist) < 0.000001){
    if (maxDist < 0.001){
      break;
    }
    prevMaxDist = maxDist;
  }
  cout << "." << endl;
  // Copy to the output
  for (int i = 0; i < *nTotalCoords; i++){
    imputed[i] = tmpImputed[i];
  }
  // Release memory
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpX;
  delete[] tmpPoint;
  delete[] tmpImputed;
  delete[] depths;
  for (int i = 0; i < numTries; i++){
    delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
}

void ImputeSpatialIter(double* X, int* n, int* d,
                   int* totalCoords, int* nTotalCoords,
                   double* imputed){
  // Copy and double-index X
  double* tmpX = new double[*n * *d];
  memcpy(tmpX, X, *n * *d * sizeof(double));
  double** xx = asMatrix(tmpX, *n, *d);
  // Working space for imputed coordinates
  double* tmpImputed = new double[*nTotalCoords];
  double* tmpPoint = new double[*d]; // currently imputed point
  double* tmpMu = new double[*d]; // current mean
  double* tmpW = new double[*d * *d]; // current whitening matrix
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  double* depths = new double[nPoints];
  // Iterative imputation loop
  cout << "Spatial iters: ";
  int numTries = *d * 2;
  double prevMaxDist = -1.;
  for (int i = 0; i < *nTotalCoords; i++){
    tmpImputed[i] = tmpX[totalCoords[i]];
  }
  double** tmpPoints = new double*[numTries];
  for (int j = 0; j < numTries; j++){
    tmpPoints[j] = new double[*d];
  }
  double tmpBestDepth = -1.;
  for (int i = 0; i < 25; i++){
    get_w(xx, *n, *d, tmpMu, tmpW);
    for (int j = 0; j < nPoints; j++){
      depths[j] = depth_sptl(xx[points[j]], xx, *n, *d, true, tmpW);
    }
    // Optimize and save
    int curIndex = 0;
    for (int j = 0; j < nPoints; j++){
      tmpBestDepth = -1.;
      int numTmps = 0;
      for (int k = 0; k < numTries; k++){
        double optDepth = NelderMeadConditionalSptl(depth_sptl, xx[points[j]], xx, *n, *d, true, tmpW,
                                                    coords[j], nCoords[j], &xx[k * (*d + 1)], tmpPoint,
                                                    2., 4., 0.5, 0.5, 0.00001);
        if (optDepth - tmpBestDepth >= 0){
          if (optDepth - tmpBestDepth > 0){
            tmpBestDepth = optDepth;
            numTmps = 0;
          }
          memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
        }
      }
      // cout << "(" << points[j] << "): " << depths[j] << " -> " << optDepth << "; ";
      for (int k = 0; k < nCoords[j]; k++){
        // if (optDepth - depths[j] > 0.5 / (double)(*n)){
        // if (tmpBestDepth - depths[j] > 0){
          tmpImputed[curIndex] = 0;
          for (int l = 0; l < numTmps; l++){
            tmpImputed[curIndex] += tmpPoints[l][coords[j][k]];
          }
          tmpImputed[curIndex] /= (double)numTmps;
        // }
        curIndex++;
        //tmpImputed[curIndex++] = tmpPoints[0][coords[j][k]];
      }
      // cout << "1 ";
      // }else{
      // cout << "0 ";
      // }
    }
    // cout << endl;
    // Calculate distance to the previous values
    double maxDist = -2.;
    for (int j = 0; j < *nTotalCoords; j++){
      double curDist = fabs(tmpX[totalCoords[j]] - tmpImputed[j]);
      if (curDist > maxDist){
        maxDist = curDist;
      }
    }
    cout << i + 1 << "(" << maxDist << ") ";
    //Replace
    for (int j = 0; j < *nTotalCoords; j++){
      tmpX[totalCoords[j]] = tmpImputed[j];
    }
    // DEBUG: Output of the imputed values
    // for (int j = 0; j < *nTotalCoords; j++){
    //   cout << tmpX[totalCoords[j]] << " ";
    //}
    // cout << endl;
    // DEBUG (end)
    // Break if converged
    // if (maxDist < 0.001 || fabs(maxDist - prevMaxDist) < 0.000001){
    if (maxDist < 0.001){
      break;
    }
    prevMaxDist = maxDist;
  }
  cout << "." << endl;
  // Copy to the output
  for (int i = 0; i < *nTotalCoords; i++){
    imputed[i] = tmpImputed[i];
  }
  // Release memory
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpX;
  delete[] tmpPoint;
  delete[] tmpMu;
  delete[] tmpW;
  delete[] tmpImputed;
  delete[] depths;
  for (int i = 0; i < numTries; i++){
    delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
}

void ImputeZonoidIter2(double* X, int* n, int* d,
                       int* totalCoords, int* nTotalCoords,
                       double* imputed){
  // Copy and double-index X
  double* tmpX = new double[*n * *d];
  memcpy(tmpX, X, *n * *d * sizeof(double));
  double** xx = asMatrix(tmpX, *n, *d);
  // Working space for imputed coordinates
  double* tmpImputed = new double[*nTotalCoords];
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Iterative imputation loop
  cout << "Zonoid iters: ";
  for (int i = 0; i < 25; i++){
    // Optimize
    LinearConditionalZonoidAll(xx, *n, *d, points, nPoints, coords, nCoords,
                               tmpImputed);
    // Calculate distance to the previous values
    double maxDist = 0;
    for (int j = 0; j < *nTotalCoords; j++){
      double curDist = fabs(tmpX[totalCoords[j]] - tmpImputed[j]);
      if (curDist > maxDist){
        maxDist = curDist;
      }
    }
    cout << i + 1 << "(" << maxDist << ") ";
    //Replace
    for (int j = 0; j < *nTotalCoords; j++){
      tmpX[totalCoords[j]] = tmpImputed[j];
    }
    // Break if converged
    if (maxDist < 0.001){
      break;
    }
  }
  cout << "." << endl;
  // Copy to the output
  memcpy(imputed, tmpImputed, *nTotalCoords * sizeof(double));
  // Release memory
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpX;
  delete[] tmpImputed;
}

void ImputeZonoidOneIter(double* X, int* n, int* d,
                         int* totalCoords, int* nTotalCoords,
                         double* imputed){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  // Form indexed point structures
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Imputation step
  // LinearConditionalZonoidAll(xx, *n, *d, points, nPoints, coords, nCoords,
  //                            imputed);

  imputeAllZonoid(xx, *n, *d, xx, *n, points, nPoints, coords, nCoords,
                  imputed);

  // Release memory
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
}

void ImputeHalfspaceExOneIter(double* X, int* n, int* d,
                              int* totalCoords, int* nTotalCoords,
                              double* imputed){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Working space for imputed coordinates
  double* tmpPoint = new double[*d]; // currently imputed point
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Calculate current depths of the points to be imputed
  double* depths = new double[nPoints];
  for (int i = 0; i < nPoints; i++){
    depths[i] = HD_Rec(xx[points[i]], xx, *n, *d);
  }
  // Storage for multiple optimization runs
  int numTries = *d;
  double** tmpPoints = new double*[numTries];
  for (int i = 0; i < numTries; i++){
    tmpPoints[i] = new double[*d];
  }
  double tmpBestDepth = -1.; // best achieved depth for a point
  // Optimize and save
  int curIndex = 0; // index in the global arrays of points
  for (int i = 0; i < nPoints; i++){ // for each point to be imputed
	  // Run optimization several times
	  tmpBestDepth = -1.;
    int numTmps = 0; // number of points with best depth
    for (int j = 0; j < numTries; j++){
      double optDepth = NelderMeadConditionalHfspEx(HD_Rec, xx[points[i]], xx, *n, *d,
                                                    coords[i], nCoords[i], &xx[j * (*d + 1)], tmpPoint,
                                                    2., 4., 0.5, 0.5, 0.00001);
      // Collect the points with the highest depth
      if (optDepth - tmpBestDepth > -0.5 / (double)(*n)){
        if (optDepth - tmpBestDepth > 0.5 / (double)(*n)){
          tmpBestDepth = optDepth;
          numTmps = 0;
        }
        memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
      }
    }
	  // Imputed coordinates to be imputed if depth of the point is not worse
    for (int j = 0; j < nCoords[i]; j++){
      if (tmpBestDepth - depths[i] > -0.5 / (double)(*n)){
        imputed[curIndex] = 0;
        for (int k = 0; k < numTmps; k++){
          imputed[curIndex] += tmpPoints[k][coords[i][j]];
        }
        imputed[curIndex] /= (double)numTmps;
      }else{
        imputed[curIndex] = xx[points[i]][coords[i][j]];
      }
      curIndex++;
    }
  }
  // Release memory
  for (int i = 0; i < numTries; i++){
	  delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpPoint;
  delete[] depths;
}

void ImputeSpatialOneIter(double* X, int* n, int* d, double* W,
                          int* totalCoords, int* nTotalCoords,
                          double* imputed){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Working space for imputed coordinates
  double* tmpPoint = new double[*d]; // currently imputed point
  double* tmpMu = new double[*d]; // current mean
  double* tmpW = new double[*d * *d]; // current whitening matrix
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Calculate whitening matrix and current depths of the points to be imputed
  if (*W == 0){
    tmpMu = new double[*d];
    tmpW = new double[*d * *d];
    get_w(xx, *n, *d, tmpMu, tmpW);
  }else{
    tmpW = W;
  }
  double* depths = new double[nPoints];
  for (int i = 0; i < nPoints; i++){
    depths[i] = depth_sptl(xx[points[i]], xx, *n, *d, true, tmpW);
  }
  // Storage for multiple optimization runs
  int numTries = *d * 2;
  double** tmpPoints = new double*[numTries];
  for (int i = 0; i < numTries; i++){
    tmpPoints[i] = new double[*d];
  }
  double tmpBestDepth = -1.; // best achieved depth for a point
  // Optimize and save
  int curIndex = 0; // index in the global arrays of points
  for (int i = 0; i < nPoints; i++){ // for each point to be imputed
    // Run optimization several times
    tmpBestDepth = -1.;
    int numTmps = 0; // number of points with best depth
    for (int j = 0; j < numTries; j++){
      double optDepth = NelderMeadConditionalSptl(depth_sptl, xx[points[i]], xx, *n, *d, true, tmpW,
                                                  coords[i], nCoords[i], &xx[j * (*d + 1)], tmpPoint,
                                                  2., 4., 0.5, 0.5, 0.00001);
      // Collect the points with the highest depth
      if (optDepth - tmpBestDepth >= 0){
        if (optDepth - tmpBestDepth > 0){
          tmpBestDepth = optDepth;
          numTmps = 0;
        }
        memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
      }
    }
    // Imputed coordinates to be imputed if depth of the point is not worse
    for (int j = 0; j < nCoords[i]; j++){
      // if (tmpBestDepth - depths[i] > -0.5 / (double)(*n)){
        imputed[curIndex] = 0;
        for (int k = 0; k < numTmps; k++){
          imputed[curIndex] += tmpPoints[k][coords[i][j]];
        }
        imputed[curIndex] /= (double)numTmps;
      // }else{
      //   imputed[curIndex] = xx[points[i]][coords[i][j]];
      // }
      curIndex++;
    }
  }
  // Release memory
  for (int i = 0; i < numTries; i++){
    delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpPoint;
  delete[] depths;
  if (W == 0){
    delete[] tmpMu;
    delete[] tmpW;
  }
}

void ImputeHalfspaceOneIter(double* X, int* n, int* d,
                            int* totalCoords, int* nTotalCoords,
                            int* oldDirs, int* nDirs,
                            double* dirsRaw, double* imputed){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Working space for imputed coordinates
  double* tmpPoint = new double[*d]; // currently imputed point
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Double-index random directions and projections
  double** dirs = asMatrix(dirsRaw, *nDirs, *d);
  double* prjsRaw = new double[*nDirs * *n];
  double** prjs = asMatrix(prjsRaw, *nDirs, *n);
  // Obtain random (directions if needed) and projections
  if (!(*oldDirs)){
    get_prjs2(xx, *n, *d, *nDirs, dirsRaw, prjsRaw);
  }else{
    upd_prjs(xx, *n, *d, *nDirs, dirs, prjs);
  }
  // Sort random projections
  for (int i = 0; i < *nDirs; i++){
    quick_sort(prjs[i], 0, *n - 1);
  }
  // Calculate current depths of the points to be imputed
  double* depths = new double[nPoints];
  for (int i = 0; i < nPoints; i++){
    depths[i] = depth_hfsp(xx[points[i]], *d, *nDirs, *n, dirs, prjs);
  }
  // Storage for multiple optimization runs
  int numTries = *d;
  double** tmpPoints = new double*[numTries];
  for (int i = 0; i < numTries; i++){
    tmpPoints[i] = new double[*d];
  }
  double tmpBestDepth = -1.; // best achieved depth for a point
  // Optimize and save
  int curIndex = 0; // index in the global arrays of points
  for (int i = 0; i < nPoints; i++){ // for each point to be imputed
    // Run optimization several times
    tmpBestDepth = -1.;
    int numTmps = 0; // number of points with best depth
    for (int j = 0; j < numTries; j++){
      double optDepth = NelderMeadConditionalHfsp(depth_hfsp, xx[points[i]], *d, *nDirs, *n,
                                                  dirs, prjs,
                                                  coords[i], nCoords[i], &xx[j * (*d + 1)], tmpPoint,
                                                  2., 4., 0.5, 0.5, 0.00001);
      // Collect the points with the highest depth
      if (optDepth - tmpBestDepth > -0.5 / (double)(*n)){
        if (optDepth - tmpBestDepth > 0.5 / (double)(*n)){
          tmpBestDepth = optDepth;
          numTmps = 0;
        }
        memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
      }
    }
    // Imputed coordinates to be imputed if depth of the point is not worse
    for (int j = 0; j < nCoords[i]; j++){
      if (tmpBestDepth - depths[i] > -0.5 / (double)(*n)){
        imputed[curIndex] = 0;
        for (int k = 0; k < numTmps; k++){
          imputed[curIndex] += tmpPoints[k][coords[i][j]];
        }
        imputed[curIndex] /= (double)numTmps;
      }else{
        imputed[curIndex] = xx[points[i]][coords[i][j]];
      }
      curIndex++;
    }
  }
  // Release memory
  for (int i = 0; i < numTries; i++){
    delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpPoint;
  delete[] depths;
  delete[] dirs;
  delete[] prjs;
  delete[] prjsRaw;
}

void ImputeHalfExtrOneIter(double* X, int* n, int* d,
                           int* totalCoords, int* nTotalCoords,
                           int* oldDirs, int* nDirs, int* exact, int* kPar,
                           double* dirsRaw, double* imputed){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Working space for imputed coordinates
  double* tmpPoint = new double[*d]; // currently imputed point
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Double-index random directions and projections
  double** dirs = asMatrix(dirsRaw, *nDirs, *d);
  double* prjsRaw = new double[*nDirs * *n];
  double** prjs = asMatrix(prjsRaw, *nDirs, *n);
  // Obtain random (directions if needed) and projections
  if (!(*oldDirs)){
    get_prjs2(xx, *n, *d, *nDirs, dirsRaw, prjsRaw);
  }else{
    upd_prjs(xx, *n, *d, *nDirs, dirs, prjs);
  }
  // Calculate absolute values of the sample points
  double* absXvals = new double[*n];
  for (int i = 0; i < *n; i++){
    absXvals[i] = 0;
    for (int j = 0; j < *d; j++){
      absXvals[i] += pow(xx[i][j], 2);
    }
    absXvals[i] = sqrt(absXvals[i]);
  }
  quick_sort(absXvals, 0, *n - 1);
  // Calculate moments
  double m1 = 0;
  double m2 = 0;
  for (int i = 0; i < *kPar; i++){
    m1 += log(absXvals[*n - i - 1]) - log(absXvals[*n - *kPar - 1]);
    m2 += pow(log(absXvals[*n - i - 1]) - log(absXvals[*n - *kPar - 1]), 2);
  }
  m1 /= (double)(*kPar);
  m2 /= (double)(*kPar);
  delete[] absXvals;
  // Sort random projections
  for (int i = 0; i < *nDirs; i++){
    quick_sort(prjs[i], 0, *n - 1);
  }
  // Calculate current depths of the points to be imputed
  double* depths = new double[nPoints];
  for (int i = 0; i < nPoints; i++){
    depths[i] = depth_hfspExtr(xx[points[i]], xx, *d, *nDirs, *n, dirs, prjs,
                               *exact, *kPar, m1, m2);
    //depths[i] = depth_hfsp(xx[points[i]], *d, *nDirs, *n, dirs, prjs);
  }
  // Storage for multiple optimization runs
  int numTries = *d * 2;
  double** tmpPoints = new double*[numTries];
  for (int i = 0; i < numTries; i++){
    tmpPoints[i] = new double[*d];
  }
  double tmpBestDepth = -1.; // best achieved depth for a point
  // Optimize and save
  int curIndex = 0; // index in the global arrays of points
  for (int i = 0; i < nPoints; i++){ // for each point to be imputed
    // Run optimization several times
    tmpBestDepth = -1.;
    int numTmps = 0; // number of points with best depth
    for (int j = 0; j < numTries; j++){
      double optDepth = NelderMeadCondHfspExtr(depth_hfspExtr, xx[points[i]],
                                               xx, *d, *nDirs, *n, dirs, prjs,
                                               *exact, *kPar, m1, m2,
                                               coords[i], nCoords[i],
                                               &xx[j * (*d + 1)], tmpPoint,
                                               2., 4., 0.5, 0.5, 0.00001);
      //double optDepth = NelderMeadConditionalHfsp(depth_hfsp, xx[points[i]],
      //                                         *d, *nDirs, *n,
      //                                         dirs, prjs,
      //                                         coords[i], nCoords[i],
      //                                                           &xx[j * (*d + 1)], tmpPoint,
      //                                                           2., 4., 0.5, 0.5, 0.00001);
      // Collect the points with the highest depth
      if (optDepth - tmpBestDepth > -0.5 / (double)(*n)){
        if (optDepth - tmpBestDepth > 0.5 / (double)(*n)){
          tmpBestDepth = optDepth;
          numTmps = 0;
        }
        memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
      }
    }
    // Imputed coordinates to be imputed if depth of the point is not worse
    for (int j = 0; j < nCoords[i]; j++){
      if (tmpBestDepth - depths[i] > -0.5 / (double)(*n)){
        imputed[curIndex] = 0;
        for (int k = 0; k < numTmps; k++){
          imputed[curIndex] += tmpPoints[k][coords[i][j]];
        }
        imputed[curIndex] /= (double)numTmps;
      }else{
        imputed[curIndex] = xx[points[i]][coords[i][j]];
      }
      curIndex++;
    }
  }
  // Release memory
  for (int i = 0; i < numTries; i++){
    delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpPoint;
  delete[] depths;
  delete[] dirs;
  delete[] prjs;
  delete[] prjsRaw;
}

void ImputeProjectionOneIter(double* X, int* n, int* d,
                             int* totalCoords, int* nTotalCoords,
                             int* oldDirs, int* nDirs,
                             double* dirsRaw, double* imputed){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Working space for imputed coordinates
  double* tmpPoint = new double[*d]; // currently imputed point
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Double-index random directions and projections
  double** dirs = asMatrix(dirsRaw, *nDirs, *d);
  double* medians = new double[*nDirs];
  double* MADs = new double[*nDirs];
  // Obtain random (directions if needed) and projections
  if (!(*oldDirs)){
    get_MaMADs2(xx, *n, *d, *nDirs, dirsRaw, medians, MADs);
  }else{
    get_mads(xx, *n, *d, *nDirs, dirs, medians, MADs);
  }
  // Calculate current depths of the points to be imputed
  double* depths = new double[nPoints];
  for (int i = 0; i < nPoints; i++){
    depths[i] = depth_proj(xx[points[i]], *d, *nDirs, dirs, medians, MADs);
  }
  // Storage for multiple optimization runs
  int numTries = 1;
  double** tmpPoints = new double*[numTries];
  for (int i = 0; i < numTries; i++){
    tmpPoints[i] = new double[*d];
  }
  double tmpBestDepth = -1.; // best achieved depth for a point
  // Optimize and save
  int curIndex = 0; // index in the global arrays of points
  for (int i = 0; i < nPoints; i++){ // for each point to be imputed
    // Run optimization several times
    tmpBestDepth = -1.;
    int numTmps = 0; // number of points with best depth
    for (int j = 0; j < numTries; j++){
      double optDepth = NelderMeadConditionalProj(depth_proj, xx[points[i]], *d, *nDirs,
                                                  dirs, medians, MADs,
                                                  coords[i], nCoords[i], xx, tmpPoint);
      // Collect the points with the highest depth
      if (optDepth - tmpBestDepth >= -eps){
        if (optDepth - tmpBestDepth > eps){
          tmpBestDepth = optDepth;
          numTmps = 0;
        }
        memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
      }
    }
    // Imputed coordinates to be imputed if depth of the point is not worse
    for (int j = 0; j < nCoords[i]; j++){
      // if (tmpBestDepth - depths[i] > -eps){
        imputed[curIndex] = 0;
        for (int k = 0; k < numTmps; k++){
          imputed[curIndex] += tmpPoints[k][coords[i][j]];
        }
        imputed[curIndex] /= (double)numTmps;
      // }else{
      //   imputed[curIndex] = xx[points[i]][coords[i][j]];
      // }
      curIndex++;
    }
  }
  // Release memory
  for (int i = 0; i < numTries; i++){
    delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpPoint;
  delete[] depths;
  delete[] dirs;
  delete[] medians;
  delete[] MADs;
}

void ImputeSimplicialExOneIter(double* X, int* n, int* d,
                               int* totalCoords, int* nTotalCoords,
                               double* imputed){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Working space for imputed coordinates
  double* tmpPoint = new double[*d]; // currently imputed point
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Calculate current depths of the points to be imputed
  double* depths = new double[nPoints];
  for (int i = 0; i < nPoints; i++){
    depths[i] = HD_Rec(xx[points[i]], xx, *n, *d);
  }
  // Storage for multiple optimization runs
  int numTries = *d;
  double** tmpPoints = new double*[numTries];
  for (int i = 0; i < numTries; i++){
    tmpPoints[i] = new double[*d];
  }
  double tmpBestDepth = -1.; // best achieved depth for a point
  // Optimize and save
  int curIndex = 0; // index in the global arrays of points
  for (int i = 0; i < nPoints; i++){ // for each point to be imputed
    // Run optimization several times
    tmpBestDepth = -1.;
    int numTmps = 0; // number of points with best depth
    for (int j = 0; j < numTries; j++){
      double optDepth = NelderMeadConditionalHfspEx(SimplicialDepthEx, xx[points[i]], xx, *n, *d,
                                                    coords[i], nCoords[i], &xx[j * (*d + 1)], tmpPoint,
                                                    2., 4., 0.5, 0.5, 0.00001);
      // Collect the points with the highest depth
      if (optDepth - tmpBestDepth > -0.5 / (double)(*n)){
        if (optDepth - tmpBestDepth > 0.5 / (double)(*n)){
          tmpBestDepth = optDepth;
          numTmps = 0;
        }
        memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
      }
    }
    // Imputed coordinates to be imputed if depth of the point is not worse
    for (int j = 0; j < nCoords[i]; j++){
      if (tmpBestDepth - depths[i] > -0.5 / (double)(*n)){
        imputed[curIndex] = 0;
        for (int k = 0; k < numTmps; k++){
          imputed[curIndex] += tmpPoints[k][coords[i][j]];
        }
        imputed[curIndex] /= (double)numTmps;
      }else{
        imputed[curIndex] = xx[points[i]][coords[i][j]];
      }
      curIndex++;
    }
  }
  // Release memory
  for (int i = 0; i < numTries; i++){
    delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpPoint;
  delete[] depths;
}


void ImputeLocalOneIter(double* X, int* n, int* d, double* locPar, int* notion,
                        int* totalCoords, int* nTotalCoords,
                        double* imputed){
  // Double-index X
  double** xx = asMatrix(X, *n, *d);
  // Working space for imputed coordinates
  double* tmpPoint = new double[*d]; // currently imputed point
  // Double-index missing coordinates
  int* points; // indices of points containing missingness
  int nPoints; // ... and their number
  int** coords; // double-indexed array of missing coordinates
  int* coordsRaw; // ... and its raw memory
  int* nCoords; // number of missing coordinates for each point
  selectPointCoords(totalCoords, *nTotalCoords, *d,
                    &points, &nPoints, &coords, &nCoords, &coordsRaw);
  // Calculate current depths of the points to be imputed
  double* depths = new double[nPoints];
  for (int i = 0; i < nPoints; i++){
    depths[i] = LocalDepth(xx[points[i]], xx, *n, *d, *locPar, (DEPTHNOTION)*notion);
  }
  // Storage for multiple optimization runs
  int numTries = 1;
  double** tmpPoints = new double*[numTries];
  for (int i = 0; i < numTries; i++){
    tmpPoints[i] = new double[*d];
  }
  double tmpBestDepth = -1.; // best achieved depth for a point
  // Optimize and save
  int curIndex = 0; // index in the global arrays of points
  for (int i = 0; i < nPoints; i++){ // for each point to be imputed
    // Run optimization several times
    tmpBestDepth = -1.;
    int numTmps = 0; // number of points with best depth
    for (int j = 0; j < numTries; j++){
      double optDepth = NelderMeadConditionalLocal(LocalDepth, xx[points[i]], xx, *n, *d, *locPar, (DEPTHNOTION)*notion,
                                                   coords[i], nCoords[i], &xx[j * (*d + 1)], tmpPoint,
                                                   2., 4., 0.5, 0.5, 0.00001);
      // Collect the points with the highest depth
      if (optDepth - tmpBestDepth > -0.5 / (double)(*n)){
        if (optDepth - tmpBestDepth > 0.5 / (double)(*n)){
          tmpBestDepth = optDepth;
          numTmps = 0;
        }
        memcpy(tmpPoints[numTmps++], tmpPoint, *d * sizeof(double));
      }
    }
    // Imputed coordinates to be imputed if depth of the point is not worse
    for (int j = 0; j < nCoords[i]; j++){
      if (tmpBestDepth - depths[i] > -0.5 / (double)(*n)){
        imputed[curIndex] = 0;
        for (int k = 0; k < numTmps; k++){
          imputed[curIndex] += tmpPoints[k][coords[i][j]];
        }
        imputed[curIndex] /= (double)numTmps;
      }else{
        imputed[curIndex] = xx[points[i]][coords[i][j]];
      }
      curIndex++;
    }
  }
  // Release memory
  for (int i = 0; i < numTries; i++){
    delete[] tmpPoints[i];
  }
  delete[] tmpPoints;
  delete[] points;
  delete[] coords;
  delete[] coordsRaw;
  delete[] nCoords;
  delete[] xx;
  delete[] tmpPoint;
  delete[] depths;
}

double LocalDepth(double* x, double* xx, int* n, int* d, int* nx,
                  double* locPar, int* notion, double* depths){
  double* XRaw = new double[*n * *d];
  memcpy(XRaw, xx, *n * *d * sizeof(double));
  double** X = asMatrix(XRaw, *n, *d);
  for (int i = 0; i < *nx; i++){
    depths[i] = LocalDepth(x + *d * i, X, *n, *d, *locPar, (DEPTHNOTION)*notion);
  }
  delete[] X;
  delete[] XRaw;
}

/* Export functions (end) --------------------------------------------------- */

#ifdef __cplusplus
}
#endif
