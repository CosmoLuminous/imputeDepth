/******************************************************************************/
/* File:             depthCalc.cpp                                            */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains functions for depth calculation.                                  */
/*                                                                            */
/******************************************************************************/

#include "imputeDepth.h"

void depths_proj2(double** x, int n, int d, int r, double* depths,
                  double** pDirsRaw, double** pMedians, double** pMADs){
  // Initialize arrays
  *pDirsRaw = new double[r * d];
  *pMedians = new double[r];
  *pMADs = new double[r];
  double* proj = new double[n]; // projection of X onto 'curDir'
  double* projTmp = new double[n]; // working copy of proj. of X onto 'curDir'
  int* rndIndices = new int[d + 1]; // random indices between 0 and (n - 1)
  double* curDir;
  for (int i = 0; i < n; i++){
    depths[i] = -1.; // these are outlyingnesses
  }
  // Calculate outlyingnesses by maximizing over random directions
  for (int i = 0; i < r; i++){
    // 1) Generate a random direction
    // Pick (d + 1) points from X
    for (int j = 0; j < d + 1; j++){
      rndIndices[j] = rand() / (double)RAND_MAX * n;
    }
    // Average from 'rndIndices[1]' to 'rndIndices[d + 1]'
    curDir = &(*pDirsRaw)[i * d];
    for (int j = 0; j < d; j++){curDir[j] = 0;}
    for (int j = 1; j < d + 1; j++){
      for (int k = 0; k < d; k++){
        curDir[k] += x[rndIndices[j]][k];
      }
    }
    // Substract 'rndIndices[0]'
    for (int k = 0; k < d; k++){
      curDir[k] /= (double)d;
      curDir[k] -= x[rndIndices[0]][k];
    }
    // 2) Project X
    for (int j = 0; j < n; j++){
      proj[j] = 0;
      for (int k = 0; k < d; k++){
        proj[j] += x[j][k] * curDir[k];
        projTmp[j] = proj[j]; // a copy for 'quick_select'
      }
    }
    // 3) Calculate median and MAD
    double curMedian = quick_select(projTmp, n);
    (*pMedians)[i] = curMedian;
    for (int j = 0; j < n; j++){
      proj[j] = fabs(proj[j] - curMedian);
      projTmp[j] = proj[j]; // a copy for 'quick_select'
    }
    double curMAD = quick_select(projTmp, n);
    (*pMADs)[i] = curMAD;
    // 4) Minimize outlyingnesses
    for (int j = 0; j < n; j ++){
      depths[j] = fmax(depths[j], proj[j] / curMAD);
    }
  }
  // Return depths
  for (int i = 0; i < n; i++){
    depths[i] = 1./(1. + depths[i]);
  }
  // Release memory
  delete[] proj;
  delete[] projTmp;
  delete[] rndIndices;
}

void depths_proj(double** x, int n, int d, int r, double* depths){
  // Initialize arrays
  double* proj = new double[n]; // projection of X onto 'curDir'
  double* projTmp = new double[n]; // working copy of proj. of X onto 'curDir'
  int* rndIndices = new int[d + 1]; // random indices between 0 and (n - 1)
  double* curDir = new double[d]; // current direction
  for (int i = 0; i < n; i++){
    depths[i] = -1.; // these are outlyingnesses
  }
  // Calculate outlyingnesses by maximizing over random directions
  for (int i = 0; i < r; i++){
    // 1) Generate a random direction
    // Pick (d + 1) points from X
    for (int j = 0; j < d + 1; j++){
      rndIndices[j] = rand() / (double)RAND_MAX * n;
    }
    // Average from 'rndIndices[1]' to 'rndIndices[d + 1]'
    for (int j = 0; j < d; j++){curDir[j] = 0;}
    for (int j = 1; j < d + 1; j++){
      for (int k = 0; k < d; k++){
        curDir[k] += x[rndIndices[j]][k];
      }
    }
    // Substract 'rndIndices[0]'
    for (int k = 0; k < d; k++){
      curDir[k] /= (double)d;
      curDir[k] -= x[rndIndices[0]][k];
    }
    // 2) Project x
    for (int j = 0; j < n; j++){
      proj[j] = 0;
      for (int k = 0; k < d; k++){
        proj[j] += x[j][k] * curDir[k];
        projTmp[j] = proj[j]; // a copy for 'quick_select'
      }
    }
    // 3) Calculate median and MAD
    double curMedian = quick_select(projTmp, n);
    for (int j = 0; j < n; j++){
      projTmp[j] = fabs(proj[j] - curMedian);
    }
    double curMAD = quick_select(projTmp, n);
    // 4) Minimize outlyingnesses
    for (int j = 0; j < n; j ++){
      depths[j] = fmax(depths[j], proj[j] / curMAD);
    }
  }
  // Return depths
  for (int i = 0; i < n; i++){
    depths[i] = 1./(1. + depths[i]);
  }
  // Release memory
  delete[] proj;
  delete[] projTmp;
  delete[] rndIndices;
  delete[] curDir;
}

void depthsn_hfsp(double** x, int n, int d, double** data, int m,
                  int r, double* depths){
  // Calculate projections
  double* dirsRaw;
  double* prjsRaw;
  get_prjs(data, m, d, r, &dirsRaw, &prjsRaw);
  double** dirs = asMatrix(dirsRaw, r, d);
  double** prjs = asMatrix(prjsRaw, r, m);
  for (int i = 0; i < r; i++){
    quick_sort(prjs[i], 0, m - 1);
  }
  // Loop through points
  for (int i = 0; i < n; i++){
    depths[i] = depth_hfsp(x[i], d, r, m, dirs, prjs);
  }
  // Release memory
  delete[] dirs;
  delete[] dirsRaw;
  delete[] prjs;
  delete[] prjsRaw;
}

void depthsn_sptl(double** x, int n, int d, double** data, int m,
                  bool sd, double* depths){
  double* lambda = new double[d * d];
  // Calculate the whitening matrix if necessary
  if (sd){
    double* mu = new double[d];
    get_w(data, m, d, mu, lambda);
    delete[] mu;
  }
  // Loop through points
  for (int i = 0; i < n; i++){
    depths[i] = depth_sptl(x[i], data, m, d, sd, lambda);
  }
  // Release memory
  delete[] lambda;
}

// 'prjs' should be sorted
double depth_hfsp(double* x, int d, int r, int n, double** dirs, double** prjs){
  // Initialization
  int depth = n + 1;
  double curPrj = 0;
  // Loop through directions
  for (int i = 0; i < r; i++){
    // Project point
    curPrj = 0;
    for (int j = 0; j < d; j++){
      curPrj += x[j] * dirs[i][j];
    }
    // Get univariate depth from left
    int depthFromLeft = bin_searchl_routine(prjs[i], n, curPrj);
    if (depthFromLeft < depth){
      depth = depthFromLeft;
    }
    // Get univarite depth from right
    int depthFromRight = bin_searchr_routine(prjs[i], n, curPrj);
    if (depthFromRight < depth){
      depth = depthFromRight;
    }
  }
  // Return depth
  return depth / (double)n;
}

// If the depth is smaller than 'kPar/n' at least in one direction
// then anyway it is at most 'kPar/n', and thus one applies the extreme
// extention
double depth_hfspExtr(double* x, double** xx, int d, int r, int n,
                      double** dirs, double** prjs,
                      bool exact, int kPar, double m1, double m2){
  // Initialization
  int depth = n + 1;
  double pDepth = 1.1;
  double curPrj = 0;
  // Flag indicating whether extreme depth shold be computed
  bool isRandomTukey = true;
  // First calculate exact depth if demanded
  if (exact){
    isRandomTukey = false;
    int intDepth = HD_Rec(x, xx, n, d) * n + 1 / (double)n / 2;
    if (intDepth >= kPar){ // if deep enough
      return intDepth / (double)n;
    }else{
      pDepth = intDepth / (double)n;
    }
  }
  // Loop through directions
  for (int i = 0; i < r; i++){
    // Project point
    curPrj = 0;
    for (int j = 0; j < d; j++){
      curPrj += x[j] * dirs[i][j];
    }
    // The border of the "good" Tukey region
    double b = prjs[i][n - kPar - 1];
    if (curPrj > b){ // if extreme
      // Compute extreme values depth
      double gammaMinus = 1 - 0.5 / (1 - pow(m1, 2) / m2);
      double gamma = m1 + gammaMinus;
      double a = b * m1 * (1 - gammaMinus);
      double pr = fmax(0, 1 + gamma * (curPrj - b) / a);
      pr = pow(pr, -1/gamma) * kPar / n;
      if (pr < pDepth){
        pDepth = pr;
        isRandomTukey = false;
      }
    }
    if (isRandomTukey){
      // Get univariate depth from left
      int depthFromLeft = bin_searchl_routine(prjs[i], n, curPrj);
      if (depthFromLeft < depth){
        depth = depthFromLeft;
      }
      // Get univarite depth from right
      int depthFromRight = bin_searchr_routine(prjs[i], n, curPrj);
      if (depthFromRight < depth){
        depth = depthFromRight;
      }
    }
  }
  // Return depth
  if (isRandomTukey){
    return depth / (double)n;
  }else{
    return pDepth;
  }
}

double depth_proj(double* x, int d, int r, double** dirs, double* medians,
                  double* MADs){
  double outlyingness = -1.;
  // Maximize through directions
  for (int i = 0; i < r; i++){
    // Obtain current projection
    double curProj = 0;
    for (int j = 0; j < d; j++){
      curProj += x[j] * dirs[i][j];
    }
    double curOutlyingness = fabs(curProj - medians[i]) / MADs[i];
    // Update outlyingness
    if (curOutlyingness > outlyingness){
      outlyingness = curOutlyingness;
    }
  }
  return 1. / (1. + outlyingness);
}

double depth_sptl(double* x, double** xx, int n, int d, bool sd, double* w){
  // Get whitening matrix if needed
  double* lambda = 0;
  if (sd){
    if(w){
      lambda = w;
    }else{
      double* mu = new double[d];
      lambda = new double[d * d];
      get_w(xx, n, d, mu, lambda);
      delete[] mu;
    }
  }
  // Calculate the depth
  double* sum = new double[d];
  for (int i = 0; i < d; i++){
    sum[i] = 0;
  }
  double* tmpDir = new double[d];
  double* tmpWork = new double[d];
  // Go through all points
  for (int i = 0; i < n; i++){
    // Calculate the difference
    for (int j = 0; j < d; j++){
      tmpWork[j] = x[j] - xx[i][j];
    }
    // Whiten if necessary
    for (int j = 0; j < d; j++){
      if (sd){
        tmpDir[j] = 0;
        for (int k = 0; k < d; k++){
          tmpDir[j] += lambda[j * d + k] * tmpWork[k];
        }
      }else{
        tmpDir[j] = tmpWork[j];
      }
    }
    // Normalize and add
    double norm = 0;
    for (int j = 0; j < d; j++){
      norm += pow(tmpDir[j], 2);
    }
    norm = sqrt(norm);
    for (int j = 0; j < d; j++){
      if (norm > eps){
        sum[j] += tmpDir[j] / norm;
      }
    }
  }
  // Get norm and depth
  double norm = 0;
  for (int i = 0; i < d; i++){
    norm += pow(sum[i] / (double)n, 2);
  }
  norm = sqrt(norm);
  // Release memory
  if (sd && !w){
    delete[] lambda;
  }
  delete[] sum;
  delete[] tmpDir;
  delete[] tmpWork;
  // Return depth
  return (double)1 - norm;
}

void get_mads(double** x, int n, int d, int r, double** dirs, double* medians,
              double* MADs){
  // Initialization
  double* proj = new double[n];
  double* tmpProj = new double[n];
  // For each direction:
  for (int i = 0; i < r; i++){
    // Project
    for (int j = 0; j < n; j++){
      proj[j] = 0;
      for (int k = 0; k < d; k++){
        proj[j] += dirs[i][k] * x[j][k];
        tmpProj[j] = proj[j];
      }
    }
    // Obtain the median
    double curMedian = quick_select(tmpProj, n);
    medians[i] = curMedian;
    // Calculate the MAD
    for (int j = 0; j < n; j++){
      proj[j] = fabs(proj[j] - curMedian);
    }
    MADs[i] = quick_select(proj, n);
  }
  // Release memory
  delete[] proj;
  delete[] tmpProj;
}

void get_prjs(double** x, int n, int d, int r,
              double** pDirsRaw, double** pPrjsRaw){
  // Initialize arrays
  *pDirsRaw = new double[r * d]; // random directions
  *pPrjsRaw = new double[r * n]; // corresponding random projections
  int* rndIndices = new int[d + 1]; // random indices between 0 and (n - 1)
  double* curDir;
  double* curPrj;
  // Calculate projections on random directions
  for (int i = 0; i < r; i++){
    // 1) Generate a random direction
    curDir = &(*pDirsRaw)[i * d];
    // Pick (d + 1) points from x
    for (int j = 0; j < d + 1; j++){
      rndIndices[j] = rand() / (double)RAND_MAX * n;
    }
    // Average from 'rndIndices[1]' to 'rndIndices[d + 1]'
    for (int j = 0; j < d; j++){curDir[j] = 0;}
    for (int j = 1; j < d + 1; j++){
      for (int k = 0; k < d; k++){
        curDir[k] += x[rndIndices[j]][k];
      }
    }
    // Substract 'rndIndices[0]'
    for (int k = 0; k < d; k++){
      curDir[k] /= (double)d;
      curDir[k] -= x[rndIndices[0]][k];
    }
    // 2) Project x
    curPrj = &(*pPrjsRaw)[i * n];
    for (int j = 0; j < n; j++){
      curPrj[j] = 0;
      for (int k = 0; k < d; k++){
        curPrj[j] += x[j][k] * curDir[k];
      }
    }
  }
  // Release memory
  delete[] rndIndices;
}

void get_prjs2(double** x, int n, int d, int r,
               double* dirsRaw, double* prjsRaw){
  // Initialize arrays
  int* rndIndices = new int[d + 1]; // random indices between 0 and (n - 1)
  double* curDir;
  double* curPrj;
  // Calculate projections on random directions
  for (int i = 0; i < r; i++){
    // 1) Generate a random direction
    curDir = &dirsRaw[i * d];
    // Pick (d + 1) points from x
    for (int j = 0; j < d + 1; j++){
      rndIndices[j] = rand() / (double)RAND_MAX * n;
    }
    // Average from 'rndIndices[1]' to 'rndIndices[d + 1]'
    for (int j = 0; j < d; j++){curDir[j] = 0;}
    for (int j = 1; j < d + 1; j++){
      for (int k = 0; k < d; k++){
        curDir[k] += x[rndIndices[j]][k];
      }
    }
    // Substract 'rndIndices[0]'
    for (int k = 0; k < d; k++){
      curDir[k] /= (double)d;
      curDir[k] -= x[rndIndices[0]][k];
    }
    // 2) Project x
    curPrj = &prjsRaw[i * n];
    for (int j = 0; j < n; j++){
      curPrj[j] = 0;
      for (int k = 0; k < d; k++){
        curPrj[j] += x[j][k] * curDir[k];
      }
    }
  }
  // Release memory
  delete[] rndIndices;
}

void get_prjs_norm(double** x, int n, int d, int r,
                   double** pDirsRaw, double** pPrjsRaw){
  // Initialize arrays
  *pDirsRaw = new double[r * d]; // random directions
  *pPrjsRaw = new double[r * n]; // corresponding random projections
  int* rndIndices = new int[d + 1]; // random indices between 0 and (n - 1)
  double* curDir;
  double* curPrj;
  // Calculate projections on random directions
  for (int i = 0; i < r; i++){
    // 1) Generate a random direction
    curDir = &(*pDirsRaw)[i * d];
    bool found = false;
    while (!found){
      // Pick (d + 1) points from x
      for (int j = 0; j < d + 1; j++){
        rndIndices[j] = rand() / (double)RAND_MAX * n;
      }
      // Average from 'rndIndices[1]' to 'rndIndices[d + 1]'
      for (int j = 0; j < d; j++){curDir[j] = 0;}
      for (int j = 1; j < d + 1; j++){
        for (int k = 0; k < d; k++){
          curDir[k] += x[rndIndices[j]][k];
        }
      }
      // Substract 'rndIndices[0]'
      for (int j = 0; j < d; j++){
        curDir[j] /= (double)d;
        curDir[j] -= x[rndIndices[0]][j];
      }
      // Check the direction
      for (int j = 0; j < d; j++){
        if (fabs(curDir[j]) > eps){
          found = true;
          break;
        }
      }
    }
    // 2) Normalize the random direction
    double tmpNorm = 0;
    for (int j = 0; j < d; j++){
      tmpNorm += pow(curDir[j], 2);
    }
    tmpNorm = sqrt(tmpNorm);
    for (int j = 0; j < d; j++){
      curDir[j] /= tmpNorm;
    }
    // 3) Project x
    curPrj = &(*pPrjsRaw)[i * n];
    for (int j = 0; j < n; j++){
      curPrj[j] = 0;
      for (int k = 0; k < d; k++){
        curPrj[j] += x[j][k] * curDir[k];
      }
    }
  }
  // Release memory
  delete[] rndIndices;
}

void upd_prjs(double** x, int n, int d, int r,
              double** dirs, double** prjs){
  // Calculate projections on given random directions
  for (int i = 0; i < r; i++){
    // Project x
    for (int j = 0; j < n; j++){
      prjs[i][j] = 0;
      for (int k = 0; k < d; k++){
        prjs[i][j] += dirs[i][k] * x[j][k];
      }
    }
  }
}

void get_MaMADs(double** x, int n, int d, int r,
                double** pDirsRaw, double** pMedians, double** pMADs){
  // Initialize arrays
  *pDirsRaw = new double[r * d]; // random directions
  *pMedians = new double[r]; // medians
  *pMADs = new double[r]; // MADs
  double* proj = new double[n]; // projection of X onto 'curDir'
  double* projTmp = new double[n]; // working copy of proj. of X onto 'curDir'
  int* rndIndices = new int[d + 1]; // random indices between 0 and (n - 1)
  double* curDir;
  double* curPrj;
  // Calculate projections on random directions
  for (int i = 0; i < r; i++){
    // 1) Generate a random direction
    curDir = &(*pDirsRaw)[i * d];
    bool found = false;
    while (!found){
      // Pick (d + 1) points from x
      for (int j = 0; j < d + 1; j++){
        rndIndices[j] = rand() / (double)RAND_MAX * n;
      }
      // Average from 'rndIndices[1]' to 'rndIndices[d + 1]'
      for (int j = 0; j < d; j++){curDir[j] = 0;}
      for (int j = 1; j < d + 1; j++){
        for (int k = 0; k < d; k++){
          curDir[k] += x[rndIndices[j]][k];
        }
      }
      // Substract 'rndIndices[0]'
      for (int j = 0; j < d; j++){
        curDir[j] /= (double)d;
        curDir[j] -= x[rndIndices[0]][j];
        // Save direction
        (*pDirsRaw)[i * d + j] = curDir[j];
      }
      // Check the direction
      for (int j = 0; j < d; j++){
        if (fabs(curDir[j]) > eps){
          found = true;
          break;
        }
      }
    }
    // 2) Project x
    for (int j = 0; j < n; j++){
      proj[j] = 0;
      for (int k = 0; k < d; k++){
        proj[j] += x[j][k] * curDir[k];
        projTmp[j] = proj[j]; // a copy for 'quick_select'
      }
    }
    // 3) Calculate median and MAD
    (*pMedians)[i] = quick_select(projTmp, n);
    for (int j = 0; j < n; j++){
      projTmp[j] = fabs(proj[j] - (*pMedians)[i]);
    }
    (*pMADs)[i] = quick_select(projTmp, n);
  }
  // Release memory
  delete[] proj;
  delete[] projTmp;
  delete[] rndIndices;
}

void get_MaMADs2(double** x, int n, int d, int r,
                 double* dirsRaw, double* medians, double* MADs){
  // Initialize arrays
  double* proj = new double[n]; // projection of X onto 'curDir'
  double* projTmp = new double[n]; // working copy of proj. of X onto 'curDir'
  int* rndIndices = new int[d + 1]; // random indices between 0 and (n - 1)
  double* curDir;
  double* curPrj;
  // Calculate projections on random directions
  for (int i = 0; i < r; i++){
    // 1) Generate a random direction
    curDir = &dirsRaw[i * d];
    bool found = false;
    while (!found){
      // Pick (d + 1) points from x
      for (int j = 0; j < d + 1; j++){
        rndIndices[j] = rand() / (double)RAND_MAX * n;
      }
      // Average from 'rndIndices[1]' to 'rndIndices[d + 1]'
      for (int j = 0; j < d; j++){curDir[j] = 0;}
      for (int j = 1; j < d + 1; j++){
        for (int k = 0; k < d; k++){
          curDir[k] += x[rndIndices[j]][k];
        }
      }
      // Substract 'rndIndices[0]'
      for (int j = 0; j < d; j++){
        curDir[j] /= (double)d;
        curDir[j] -= x[rndIndices[0]][j];
        // Save direction
        dirsRaw[i * d + j] = curDir[j];
      }
      // Check the direction
      for (int j = 0; j < d; j++){
        if (fabs(curDir[j]) > eps){
          found = true;
          break;
        }
      }
    }
    // 2) Project x
    for (int j = 0; j < n; j++){
      proj[j] = 0;
      for (int k = 0; k < d; k++){
        proj[j] += x[j][k] * curDir[k];
        projTmp[j] = proj[j]; // a copy for 'quick_select'
      }
    }
    // 3) Calculate median and MAD
    medians[i] = quick_select(projTmp, n);
    for (int j = 0; j < n; j++){
      projTmp[j] = fabs(proj[j] - medians[i]);
    }
    MADs[i] = quick_select(projTmp, n);
  }
  // Release memory
  delete[] proj;
  delete[] projTmp;
  delete[] rndIndices;
}

void median_hfsp(double** x, int n, int d, int r, double* median){
  // Calculate projections
  double* dirsRaw;
  double* prjsRaw;
  get_prjs(x, n, d, r, &dirsRaw, &prjsRaw);
  double** dirs = asMatrix(dirsRaw, r, d);
  double** prjs = asMatrix(prjsRaw, r, n);
  for (int i = 0; i < r; i++){
    quick_sort(prjs[i], 0, n - 1);
  }
  // Optimize
  double* x_o = new double[d];
  int* coords = new int[d];
  for (int i = 0; i < d; i++){coords[i] = i;}
  NelderMeadConditionalHfsp(depth_hfsp, x_o, d, r, n, dirs, prjs,
                            coords, d, x, median);
  // Release memory
  delete[] dirs;
  delete[] dirsRaw;
  delete[] prjs;
  delete[] prjsRaw;
  delete[] x_o;
  delete[] coords;
}

void median_proj(double** x, int n, int d, int r, double* median){
  // Calculate projections
  double* dirsRaw;
  double* medians;
  double* MADs;
  get_MaMADs(x, n, d, r, &dirsRaw, &medians, &MADs);
  double** dirs = asMatrix(dirsRaw, r, d);
  // Optimize
  double* x_o = new double[d];
  int* coords = new int[d];
  for (int i = 0; i < d; i++){coords[i] = i;}
  NelderMeadConditionalProj(depth_proj, x_o, d, r, dirs, medians, MADs,
                            coords, d, x, median);
  // Release memory
  delete[] dirs;
  delete[] medians;
  delete[] MADs;
  delete[] x_o;
  delete[] coords;
}

// Taken from the R-package ddalpha
// by Pokotylo, Mozharovskyi, Dyckerhoff, Nagy
double SimplicialDepthEx(double* x, TDMatrix X, int n, int d){

  double* b = new double[d + 1]; b[d] = 1;
  double* z = new double[d + 1];
  int* counters = new int[d + 1];
  TDMatrix A = newM(d + 1, d + 1);
  unsigned long long div0 = choose(n, d + 1);

//  for (int obs = 0; obs < nx; obs++){
    unsigned long long theCounter = 0;
    unsigned long long numSimplicesChecked = 0;

    for (int i = 0; i < d; i++){ counters[i] = i; }counters[d] = d - 1;
    while (counters[0] != n - (d + 1)){
      int i = d;
      while (i > 0 && counters[i] == n - (d + 1) + i){ i--; }
      counters[i]++; int j = i + 1;
      while (j < d + 1){ counters[j] = counters[j - 1] + 1; j++; }

      for (int j = 0; j < d; j++){
        for (int k = 0; k < d + 1; k++){
          A[j][k] = X[counters[k]][j];
        }
      }
      for (int k = 0; k < d + 1; k++){
        A[d][k] = 1;
      }
      memcpy(b, x, d*sizeof(double)); b[d] = 1;
      if (solveUnique(A, b, z, d + 1)){
        bool isInside = true;
        for (int j = 0; j < d + 1; j++){
          if (z[j] < 0){ isInside = false; break; }
        }
        if (isInside){ theCounter++; }
      }
      (numSimplicesChecked) ++;
    }
    bool sc = numSimplicesChecked == div0;
    delete[] b;
    delete[] z;
    delete[] counters;
    deleteM(A);
    double depth = (double)theCounter / div0;
//    depths[obs] = depth;
    return depth;
//  }

//  delete[] b;
//  delete[] z;
//  delete[] counters;
//  deleteM(A);
}
