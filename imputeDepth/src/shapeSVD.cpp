/******************************************************************************/
/* File:             shapeSVD.cpp                                             */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains functions for the depth-based singular value decomposition of the */
/* shape matrix.                                                              */
/*                                                                            */
/******************************************************************************/

#include "imputeDepth.h"

void an_svd(double* X, int n, int d, double* s, double* u, double* vt){
  // Initialization
  const char* ju = "S";
  double* xvals = new double[n * d];
  memcpy(xvals, X, n * d * sizeof(double));
  double tmp;
  int lwork = -1;
  int* iwork = new int[8 * n];
  int info;
  // Ask for the size of the working array
  F77_CALL(dgesdd)(ju, &n, &d, xvals, &n, s, u, &n, vt, &d, &tmp, &lwork, iwork, &info);
  if (info != 0){
    error("error code %d from Lapack routine '%s'", info, "dgesdd");
  }
  // Allocate the working array
  lwork = (int)(tmp + 0.5);
  double* work = new double[lwork];
  // Main call to FORTRAN
  F77_CALL(dgesdd)(ju, &n, &d, xvals, &n, s, u, &n, vt, &d, work, &lwork, iwork, &info);
  if (info != 0)
  {
    error("error code %d from Lapack routine '%s'", info, "dgesdd");
  }
  // Release memory
  delete[] xvals;
  delete[] iwork;
  delete[] work;
}

void svd_depth_proj(double** x, int n, int d, int r, int typeCenter,
                    double* s, double* vt, double* scale, double* center){
  // Calculate depths
  double* depths = new double[n];
  depths_proj(x, n, d, r, depths);
  // Find center
  if ((typeCenter > 0) && (typeCenter < n + 1)){
    // Calculate depth-weighted trimmed mean
    // a) Sort points w.r.t. depths
    SortIndex* points = new SortIndex[n];
    for (int i = 0; i < n; i++){
      points[i].value = depths[i];
      points[i].index = i;
    }
    quick_sort(points, 0, n - 1);
    // b) Initialize
    double weightsSum = 0;
    for (int i = 0; i < d; i++){
      center[i] = 0;
    }
    // c) Weight 'centerPortion' most central points by depth
    for (int i = 0; i < typeCenter; i++){
      weightsSum += points[i].value;
      for (int j = 0; j < d; j++){
        center[j] += x[points[i].index][j] * points[i].value;
      }
    }
    for (int i = 0; i < d; i++){
      center[i] /= weightsSum;
    }
    // Release memory
    delete[] points;
  }else{
    // Calculate depth median
    median_proj(x, n, d, r, center);
  }
  // Depth-Tyler X
  double* XTylerTRaw = new double[n * d];
  double maxDepth = 1.;
  for (int i = 0; i < n; i++){ // For each point in X
    // Get its distance from center
    double curLength = 0;
    for (int j = 0; j < d; j++){
      curLength += pow(x[i][j] - center[j], 2);
    }
    curLength = sqrt(curLength);
    // Center, normalize, and scale by depth-based distance
    for (int j = 0; j < d; j++){
      XTylerTRaw[j * n + i] = (x[i][j] - center[j]) /
      curLength * (maxDepth - depths[i]);
    }
  }
  //Perform SVD
  double* tmpU = new double[n * d];
  an_svd(XTylerTRaw, n, d, s, tmpU, vt);
  // Calculate shape components
  // a) Copy and double-index vt
  double* tmpVRaw = new double[d * d];
  double** tmpV = new double*[d];
  for (int i = 0; i < d; i++){
    tmpV[i] = &tmpVRaw[i * d];
    for (int j = 0; j < d; j++){
      tmpV[i][j] = vt[j * d + i];
    }
  }
  // b) Obtain MADs
  double* pca_medians = new double[d];
  double* pca_MADs = new double[d];
  get_mads(x, n, d, d, tmpV, pca_medians, pca_MADs);
  // c) Sort MADs
  SortIndex* mads = new SortIndex[d];
  for (int i = 0; i < d; i++){
    mads[i].value = pca_MADs[i];
    mads[i].index = i;
  }
  quick_sort(mads, 0, d - 1);
  // c) Normalize the last (d - 1) shape components w.r.t. to the first one
  //    and reorder the vt
  scale[0] = mads[0].value;
  s[0] = 1.;
  for (int i = 0; i < d; i++){
    if (i > 0){
      s[i] = mads[i].value / mads[0].value;
    }
    for (int j = 0; j < d; j++){
      vt[j * d + i] = tmpV[mads[i].index][j];
    }
  }
  // Release memory
  delete[] depths;
  delete[] XTylerTRaw;
  delete[] tmpU;
  delete[] tmpV;
  delete[] tmpVRaw;
  delete[] pca_medians;
  delete[] pca_MADs;
  delete[] mads;
}

void svd_depth_hfsp(double** x, int n, int d, int r,
                    int typeDepth, int typeCenter,
                    double* s, double* vt, double* scale, double* center){
  // Calculate depths
  double* depths = new double[n];
  if (typeDepth == 0){
    depthsn_hfsp(x, n, d, x, n, r, depths);
  }else{
    for (int i = 0; i < n; i++){
      depths[i] = HD_Rec(x[i], x, n, d);
    }
  }
  // Find center
  if ((typeCenter > 0) && (typeCenter < n + 1)){
    // Calculate depth-weighted trimmed mean
    // a) Sort points w.r.t. depths
    SortIndex* points = new SortIndex[n];
    for (int i = 0; i < n; i++){
      points[i].value = depths[i];
      points[i].index = i;
    }
    quick_sort(points, 0, n - 1);
    // b) Initialize
    double weightsSum = 0;
    for (int i = 0; i < d; i++){
      center[i] = 0;
    }
    // c) Weight 'centerPortion' most central points by depth
    for (int i = 0; i < typeCenter; i++){
      weightsSum += points[i].value;
      for (int j = 0; j < d; j++){
        center[j] += x[points[i].index][j] * points[i].value;
      }
    }
    for (int i = 0; i < d; i++){
      center[i] /= weightsSum;
    }
    // Release memory
    delete[] points;
  }else{
    // Calculate depth median
    median_hfsp(x, n, d, r, center);
  }
  // Depth-Tyler X
  double* XTylerTRaw = new double[n * d];
  double maxDepth = 0.5;
  for (int i = 0; i < n; i++){ // For each point in X
    // Get its distance from center
    double curLength = 0;
    for (int j = 0; j < d; j++){
      curLength += pow(x[i][j] - center[j], 2);
    }
    curLength = sqrt(curLength);
    // Center, normalize, and scale by depth-based distance
    for (int j = 0; j < d; j++){
      XTylerTRaw[j * n + i] = (x[i][j] - center[j]) /
      curLength * (maxDepth - depths[i]);
    }
  }
  //Perform SVD
  double* tmpU = new double[n * d];
  an_svd(XTylerTRaw, n, d, s, tmpU, vt);
  // Calculate shape components
  // a) Copy and double-index vt
  double* tmpVRaw = new double[d * d];
  double** tmpV = new double*[d];
  for (int i = 0; i < d; i++){
    tmpV[i] = &tmpVRaw[i * d];
    for (int j = 0; j < d; j++){
      tmpV[i][j] = vt[j * d + i];
    }
  }
  // b) Obtain MADs
  double* pca_medians = new double[d];
  double* pca_MADs = new double[d];
  get_mads(x, n, d, d, tmpV, pca_medians, pca_MADs);
  // c) Sort MADs
  SortIndex* mads = new SortIndex[d];
  for (int i = 0; i < d; i++){
    mads[i].value = pca_MADs[i];
    mads[i].index = i;
  }
  quick_sort(mads, 0, d - 1);
  // c) Normalize the last (d - 1) shape components w.r.t. to the first one
  //    and reorder the vt
  scale[0] = mads[0].value;
  s[0] = 1.;
  for (int i = 0; i < d; i++){
    if (i > 0){
      s[i] = mads[i].value / mads[0].value;
    }
    for (int j = 0; j < d; j++){
      vt[j * d + i] = tmpV[mads[i].index][j];
    }
  }
  // Release memory
  delete[] depths;
  delete[] XTylerTRaw;
  delete[] tmpU;
  delete[] tmpV;
  delete[] tmpVRaw;
  delete[] pca_medians;
  delete[] pca_MADs;
  delete[] mads;
}

void get_w(double** xx, int n, int d, double* mu, double* wRaw){
  // Calculate the mean
  for (int i = 0; i < d; i++){
    mu[i] = 0;
    for (int j = 0; j < n; j++){
      mu[i] += xx[j][i];
    }
    mu[i] /= (double)n;
  }
  // Demean and transpose
  double* X = new double[n * d];
  for (int i = 0; i < n; i++){
    for (int j = 0; j < d; j++){
      X[j * n + i] = xx[i][j] - mu[j];
    }
  }
  // Run svd
  double* tmpU = new double[n * d];
  double* s = new double[d];
  double* vt = new double[d * d];
  an_svd(X, n, d, s, tmpU, vt);
  // Gather the "wRaw"
  for (int i = 0; i < d; i++){
    s[i] = sqrt(n - 1) / s[i];
    for (int j = 0; j < d; j++){
      wRaw[i * d + j] = vt[j * d + i] * s[i];
    }
  }
  // Release memory
  delete[] X;
  delete[] tmpU;
  delete[] s;
  delete[] vt;
}

/* Routines ----------------------------------------------------------------- */
void svd_routine(char* jobu, double* x, int* pN, int *pD, double* s, double* u,
                double* vt){
  int n, p, info = 0;
  n = pN[0];
  p = pD[0];

  /* work on a copy of x  */
  double *xvals;
  xvals = new double[n * p];
  memcpy(xvals, x, n * p * sizeof(double));

  int ldu = n;
  int ldvt = p;
  double tmp;
  int* iwork = new int[8 * (n < p ? n : p)];

  /* ask for optimal size of work array */
  const char *ju = "S";
  int lwork = -1;

  F77_CALL(dgesdd)(ju, &n, &p, xvals, &n, s,
           u, &ldu, vt, &ldvt,
           &tmp, &lwork, iwork, &info);
  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dgesdd");
  lwork = (int) (tmp + 0.5);
  double* work = new double[lwork];
  F77_CALL(dgesdd)(ju, &n, &p, xvals, &n, s,
           u, &ldu, vt, &ldvt,
           work, &lwork, iwork, &info);
  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dgesdd");
}
