/******************************************************************************/
/* File:             condOpt.cpp                                              */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains functions for conditional depth maximization.                     */
/*                                                                            */
/******************************************************************************/

#include "imputeDepth.h"

double NelderMeadConditionalSptl(
    double (*fun)(double* x, double** xx, int n, int d, bool sd, double* w),
    double* x, double** xx, int n, int d, bool sd, double* w,
    int* coords, int nCoords, double** start, double* output,
    double alpha, double gamma, double rho, double sigma, double eps){
  // Allocate memory
  Feval* fevals = new Feval[nCoords + 1]; // function evaluations
  for (int i = 0; i < nCoords + 1; i++){
    fevals[i].arg = new double[nCoords];
    for (int j = 0; j < nCoords; j++){
      fevals[i].arg[j] = 0;
    }
  }
  double* x_o = new double[nCoords]; // centroid
  double* x_r = new double[nCoords]; double f_r; // reflected point
  double* x_e = new double[nCoords]; double f_e; // expanded point
  double* x_c = new double[nCoords]; double f_c; // contracted point
  double* x_full = new double[d]; // complete point to evaluate
  memcpy(x_full, x, d * sizeof(double)); // the condition (template) point
  // Generate initial simplex
  if (start == 0){
    for (int i = 0; i < nCoords; i++){
      fevals[i].arg[i] = 1;
    }
  }else{
    for (int i = 0; i < nCoords + 1; i++){
      for (int j = 0; j < nCoords; j++){
        fevals[i].arg[j] = start[i][coords[j]];
      }
    }
  }
  // Perform initial evaluation
  for (int i = 0; i < nCoords; i++){
    // Complete the point
    for (int j = 0; j < nCoords; j++){
      x_full[coords[j]] = fevals[i].arg[j];
    }
    // Evaluate it
    fevals[i].val = -(*fun)(x_full, xx, n, d, sd, w);
  }
  // Main optimization cycle
  int iter = 0;
  bool converged = false;
  while(!converged){
    iter++;
    // cout << "Iteration " << iter << ": ";
    // Evaluate the new "worst" point
    for (int i = 0; i < nCoords; i++){
      x_full[coords[i]] = fevals[nCoords].arg[i];
    }
    fevals[nCoords].val = -(*fun)(x_full, xx, n, d, sd, w);
    // Order the values
    quick_sort(fevals, 0, nCoords);
    // DEBUG -------------------------------------------------------------------
    // Print current simplex
    //cout << "Current simplex consists of vertices:" << endl;
    //for (int i = 0; i < d + 1; i++){
    //  cout << fevals[i].val << ":\t";
    //  for (int j = 0; j < d; j++){
    //    cout << fevals[i].arg[j] << " ";
    //  }
    //  cout << endl;
    //}
    // DEBUG (end) -------------------------------------------------------------
    // Check convergence criterium
    double maxDist = 0;
    for (int i = 0; i < nCoords; i++){
      double tmpVal = fabs(fevals[0].arg[i] - fevals[nCoords].arg[i]);
      if (tmpVal > maxDist){
        maxDist = tmpVal;
      }
    }
    if (maxDist < eps){
      converged = true;
      break;
    }
    // Calculate the centroid of nCoords best points
    for (int i = 0; i < nCoords; i++){
      x_o[i] = 0;
      for (int j = 0; j < nCoords; j++){
        x_o[i] += fevals[j].arg[i];
      }
      x_o[i] /= (double)nCoords;
    }
    // Calculate and evaluate reflected point
    for (int i = 0; i < nCoords; i++){
      x_r[i] = x_o[i] + alpha * (x_o[i] - fevals[nCoords].arg[i]);
      x_full[coords[i]] = x_r[i];
    }
    f_r = -(*fun)(x_full, xx, n, d, sd, w);
    // Choose what to do
    if ((fevals[0].val <= f_r) && (f_r < fevals[nCoords - 1].val)){
      // Reflection
      // cout << "Reflection" << endl;
      for (int i = 0; i < nCoords; i++){
        fevals[nCoords].arg[i] = x_r[i];
      }
    }else{
      if (f_r < fevals[0].val){
        // Calculate and evaluate expanded point
        for (int i = 0; i < nCoords; i++){
          x_e[i] = x_r[i] + gamma * (x_r[i] - x_o[i]);
          x_full[coords[i]] = x_e[i];
        }
        f_e = -(*fun)(x_full, xx, n, d, sd, w);
        if (f_e < f_r){
          // Expansion
          // cout << "Expansion" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_e[i];
          }
        }else{
          // Still (just) reflection
          // cout << "Reflection" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_r[i];
          }
        }
      }else{
        // Calculate and evaluate contracted point
        for (int i = 0; i < nCoords; i++){
          x_c[i] = x_o[i] + rho * (fevals[nCoords].arg[i] - x_o[i]);
          x_full[coords[i]] = x_c[i];
        }
        f_c = -(*fun)(x_full, xx, n, d, sd, w);
        if (f_c < fevals[nCoords].val){
          // Contraction
          // cout << "Contraction" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_c[i];
          }
        }else{
          // Reduction
          // cout << "Reduction" << endl;
          for (int i = 1; i < nCoords + 1; i++){
            for (int j = 0; j < nCoords; j++){
              fevals[i].arg[j] = fevals[0].arg[j] +
                sigma * (fevals[i].arg[j] - fevals[0].arg[j]);
            }
          }
          // Partial evaluation
          for (int i = 1; i < nCoords; i++){
            // Complete the point
            for (int j = 0; j < nCoords; j++){
              x_full[coords[j]] = fevals[i].arg[j];
            }
            // Evaluate it
            fevals[i].val = -(*fun)(x_full, xx, n, d, sd, w);
          }
        }
      }
    }
  }
  // Extract results
  if (converged){
    memcpy(output, x, d * sizeof(double));
    for (int i = 0; i < nCoords; i++){
      output[coords[i]] = fevals[0].arg[i];
    }
  }
  double res = -fevals[0].val;
  // Release memory
  for (int i = 0; i < nCoords + 1; i++){
    delete[] fevals[i].arg;
  }
  delete[] fevals;
  delete[] x_o;
  delete[] x_r;
  delete[] x_e;
  delete[] x_c;
  delete[] x_full;
  return res;
}

double NelderMeadConditionalHfspEx(
    double (*fun)(double* z, double** xx, int n, int d),
    double* x, double** xx, int n, int d,
    int* coords, int nCoords, double** start, double* output,
    double alpha, double gamma, double rho, double sigma, double eps){
  // Allocate memory
  Feval* fevals = new Feval[nCoords + 1]; // function evaluations
  for (int i = 0; i < nCoords + 1; i++){
    fevals[i].arg = new double[nCoords];
    for (int j = 0; j < nCoords; j++){
      fevals[i].arg[j] = 0;
    }
  }
  double* x_o = new double[nCoords]; // centroid
  double* x_r = new double[nCoords]; double f_r; // reflected point
  double* x_e = new double[nCoords]; double f_e; // expanded point
  double* x_c = new double[nCoords]; double f_c; // contracted point
  double* x_full = new double[d]; // complete point to evaluate
  memcpy(x_full, x, d * sizeof(double)); // the condition (template) point
  // Generate initial simplex
  if (start == 0){
    for (int i = 0; i < nCoords; i++){
      fevals[i].arg[i] = 1;
    }
  }else{
    for (int i = 0; i < nCoords + 1; i++){
      for (int j = 0; j < nCoords; j++){
        fevals[i].arg[j] = start[i][coords[j]];
      }
    }
  }
  // Perform initial evaluation
  for (int i = 0; i < nCoords; i++){
    // Complete the point
    for (int j = 0; j < nCoords; j++){
      x_full[coords[j]] = fevals[i].arg[j];
    }
    // Evaluate it
    fevals[i].val = -(*fun)(x_full, xx, n, d);
  }
  // Main optimization cycle
  int iter = 0;
  bool converged = false;
  while(!converged){
    iter++;
    // cout << "Iteration " << iter << ": ";
    // Evaluate the new "worst" point
    for (int i = 0; i < nCoords; i++){
      x_full[coords[i]] = fevals[nCoords].arg[i];
    }
    fevals[nCoords].val = -(*fun)(x_full, xx, n, d);
    // Order the values
    quick_sort(fevals, 0, nCoords);
    // DEBUG -------------------------------------------------------------------
    // Print current simplex
    //cout << "Current simplex consists of vertices:" << endl;
    //for (int i = 0; i < d + 1; i++){
    //  cout << fevals[i].val << ":\t";
    //  for (int j = 0; j < d; j++){
    //    cout << fevals[i].arg[j] << " ";
    //  }
    //  cout << endl;
    //}
    // DEBUG (end) -------------------------------------------------------------
    // Check convergence criterium
    double maxDist = 0;
    for (int i = 0; i < nCoords; i++){
      double tmpVal = fabs(fevals[0].arg[i] - fevals[nCoords].arg[i]);
      if (tmpVal > maxDist){
        maxDist = tmpVal;
      }
    }
    if (maxDist < eps){
      converged = true;
      break;
    }
    // Calculate the centroid of nCoords best points
    for (int i = 0; i < nCoords; i++){
      x_o[i] = 0;
      for (int j = 0; j < nCoords; j++){
        x_o[i] += fevals[j].arg[i];
      }
      x_o[i] /= (double)nCoords;
    }
    // Calculate and evaluate reflected point
    for (int i = 0; i < nCoords; i++){
      x_r[i] = x_o[i] + alpha * (x_o[i] - fevals[nCoords].arg[i]);
      x_full[coords[i]] = x_r[i];
    }
    f_r = -(*fun)(x_full, xx, n, d);
    // Choose what to do
    if ((fevals[0].val <= f_r) && (f_r < fevals[nCoords - 1].val)){
      // Reflection
      // cout << "Reflection" << endl;
      for (int i = 0; i < nCoords; i++){
        fevals[nCoords].arg[i] = x_r[i];
      }
    }else{
      if (f_r < fevals[0].val){
        // Calculate and evaluate expanded point
        for (int i = 0; i < nCoords; i++){
          x_e[i] = x_r[i] + gamma * (x_r[i] - x_o[i]);
          x_full[coords[i]] = x_e[i];
        }
        f_e = -(*fun)(x_full, xx, n, d);
        if (f_e < f_r){
          // Expansion
          // cout << "Expansion" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_e[i];
          }
        }else{
          // Still (just) reflection
          // cout << "Reflection" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_r[i];
          }
        }
      }else{
        // Calculate and evaluate contracted point
        for (int i = 0; i < nCoords; i++){
          x_c[i] = x_o[i] + rho * (fevals[nCoords].arg[i] - x_o[i]);
          x_full[coords[i]] = x_c[i];
        }
        f_c = -(*fun)(x_full, xx, n, d);
        if (f_c < fevals[nCoords].val){
          // Contraction
          // cout << "Contraction" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_c[i];
          }
        }else{
          // Reduction
          // cout << "Reduction" << endl;
          for (int i = 1; i < nCoords + 1; i++){
            for (int j = 0; j < nCoords; j++){
              fevals[i].arg[j] = fevals[0].arg[j] +
                sigma * (fevals[i].arg[j] - fevals[0].arg[j]);
            }
          }
          // Partial evaluation
          for (int i = 1; i < nCoords; i++){
            // Complete the point
            for (int j = 0; j < nCoords; j++){
              x_full[coords[j]] = fevals[i].arg[j];
            }
            // Evaluate it
            fevals[i].val = -(*fun)(x_full, xx, n, d);
          }
        }
      }
    }
  }
  // Extract results
  if (converged){
    memcpy(output, x, d * sizeof(double));
    for (int i = 0; i < nCoords; i++){
      output[coords[i]] = fevals[0].arg[i];
    }
  }
  double res = -fevals[0].val;
  // Release memory
  for (int i = 0; i < nCoords + 1; i++){
    delete[] fevals[i].arg;
  }
  delete[] fevals;
  delete[] x_o;
  delete[] x_r;
  delete[] x_e;
  delete[] x_c;
  delete[] x_full;
  return res;
}

double NelderMeadConditionalHfsp(
  double (*fun)(double* x, int d, int r, int n, double** dirs, double** prjs),
  double* x, int d, int r, int n, double** dirs, double** prjs,
  int* coords, int nCoords, double** start, double* output,
  double alpha, double gamma, double rho, double sigma, double eps){
  // Allocate memory
  Feval* fevals = new Feval[nCoords + 1]; // function evaluations
  for (int i = 0; i < nCoords + 1; i++){
    fevals[i].arg = new double[nCoords];
    for (int j = 0; j < nCoords; j++){
      fevals[i].arg[j] = 0;
    }
  }
  double* x_o = new double[nCoords]; // centroid
  double* x_r = new double[nCoords]; double f_r; // reflected point
  double* x_e = new double[nCoords]; double f_e; // expanded point
  double* x_c = new double[nCoords]; double f_c; // contracted point
  double* x_full = new double[d]; // complete point to evaluate
  memcpy(x_full, x, d * sizeof(double)); // the condition (template) point
  // Generate initial simplex
  if (start == 0){
    for (int i = 0; i < nCoords; i++){
      fevals[i].arg[i] = 1;
    }
  }else{
    for (int i = 0; i < nCoords + 1; i++){
      for (int j = 0; j < nCoords; j++){
        fevals[i].arg[j] = start[i][coords[j]];
      }
    }
  }
  // Perform initial evaluation
  for (int i = 0; i < nCoords; i++){
    // Complete the point
    for (int j = 0; j < nCoords; j++){
      x_full[coords[j]] = fevals[i].arg[j];
    }
    // Evaluate it
    fevals[i].val = -(*fun)(x_full, d, r, n, dirs, prjs);
  }
  // Main optimization cycle
  int iter = 0;
  bool converged = false;
  while(!converged){
    iter++;
    // cout << "Iteration " << iter << ": ";
    // Evaluate the new "worst" point
    for (int i = 0; i < nCoords; i++){
      x_full[coords[i]] = fevals[nCoords].arg[i];
    }
    fevals[nCoords].val = -(*fun)(x_full, d, r, n, dirs, prjs);
    // Order the values
    quick_sort(fevals, 0, nCoords);
    // DEBUG -------------------------------------------------------------------
    // Print current simplex
    // cout << "Current simplex consists of vertices:" << endl;
    // for (int i = 0; i < nCoords + 1; i++){
    //   cout << fevals[i].val << ":\t";
    //   for (int j = 0; j < nCoords; j++){
    //     cout << fevals[i].arg[j] << " ";
    //   }
    //   cout << endl;
    // }
    // DEBUG (end) -------------------------------------------------------------
    // Check convergence criterium
    double maxDist = 0;
    for (int i = 0; i < nCoords; i++){
      double tmpVal = fabs(fevals[0].arg[i] - fevals[nCoords].arg[i]);
      if (tmpVal > maxDist){
        maxDist = tmpVal;
      }
    }
    if (maxDist < eps){
      converged = true;
      break;
    }
    // Calculate the centroid of nCoords best points
    for (int i = 0; i < nCoords; i++){
      x_o[i] = 0;
      for (int j = 0; j < nCoords; j++){
        x_o[i] += fevals[j].arg[i];
      }
      x_o[i] /= (double)nCoords;
    }
    // Calculate and evaluate reflected point
    for (int i = 0; i < nCoords; i++){
      x_r[i] = x_o[i] + alpha * (x_o[i] - fevals[nCoords].arg[i]);
      x_full[coords[i]] = x_r[i];
    }
    f_r = -(*fun)(x_full, d, r, n, dirs, prjs);
    // Choose what to do
    if ((fevals[0].val <= f_r) && (f_r < fevals[nCoords - 1].val)){
      // Reflection
      // cout << "Reflection" << endl;
      for (int i = 0; i < nCoords; i++){
        fevals[nCoords].arg[i] = x_r[i];
      }
    }else{
      if (f_r < fevals[0].val){
        // Calculate and evaluate expanded point
        for (int i = 0; i < nCoords; i++){
          x_e[i] = x_r[i] + gamma * (x_r[i] - x_o[i]);
          x_full[coords[i]] = x_e[i];
        }
        f_e = -(*fun)(x_full, d, r, n, dirs, prjs);
        if (f_e < f_r){
          // Expansion
          // cout << "Expansion" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_e[i];
          }
        }else{
          // Still (just) reflection
          // cout << "Reflection" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_r[i];
          }
        }
      }else{
        // Calculate and evaluate contracted point
        for (int i = 0; i < nCoords; i++){
          x_c[i] = x_o[i] + rho * (fevals[nCoords].arg[i] - x_o[i]);
          x_full[coords[i]] = x_c[i];
        }
        f_c = -(*fun)(x_full, d, r, n, dirs, prjs);
        if (f_c < fevals[nCoords].val){
          // Contraction
          // cout << "Contraction" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_c[i];
          }
        }else{
          // Reduction
          // cout << "Reduction" << endl;
          for (int i = 1; i < nCoords + 1; i++){
            for (int j = 0; j < nCoords; j++){
              fevals[i].arg[j] = fevals[0].arg[j] +
                sigma * (fevals[i].arg[j] - fevals[0].arg[j]);
            }
          }
          // Partial evaluation
          for (int i = 1; i < nCoords; i++){
            // Complete the point
            for (int j = 0; j < nCoords; j++){
              x_full[coords[j]] = fevals[i].arg[j];
            }
            // Evaluate it
            fevals[i].val = -(*fun)(x_full, d, r, n, dirs, prjs);
          }
        }
      }
    }
  }
  // Extract results
  if (converged){
    memcpy(output, x, d * sizeof(double));
    for (int i = 0; i < nCoords; i++){
      output[coords[i]] = fevals[0].arg[i];
    }
  }
  double res = -fevals[0].val;
  // Release memory
  for (int i = 0; i < nCoords + 1; i++){
    delete[] fevals[i].arg;
  }
  delete[] fevals;
  delete[] x_o;
  delete[] x_r;
  delete[] x_e;
  delete[] x_c;
  delete[] x_full;
  return res;
}

double NelderMeadCondHfspExtr(
    double (*fun)(double* x, double** xx, int d, int r, int n,
            double** dirs, double** prjs,
            bool exact, int kPar, double m1, double m2),
    double* x, double** xx, int d, int r, int n, double** dirs, double** prjs,
    bool exact, int kPar, double m1, double m2,
    int* coords, int nCoords, double** start, double* output,
    double alpha, double gamma, double rho, double sigma, double eps){
  // Allocate memory
  Feval* fevals = new Feval[nCoords + 1]; // function evaluations
  for (int i = 0; i < nCoords + 1; i++){
    fevals[i].arg = new double[nCoords];
    for (int j = 0; j < nCoords; j++){
      fevals[i].arg[j] = 0;
    }
  }
  double* x_o = new double[nCoords]; // centroid
  double* x_r = new double[nCoords]; double f_r; // reflected point
  double* x_e = new double[nCoords]; double f_e; // expanded point
  double* x_c = new double[nCoords]; double f_c; // contracted point
  double* x_full = new double[d]; // complete point to evaluate
  memcpy(x_full, x, d * sizeof(double)); // the condition (template) point
  // Generate initial simplex
  if (start == 0){
    for (int i = 0; i < nCoords; i++){
      fevals[i].arg[i] = 1;
    }
  }else{
    for (int i = 0; i < nCoords + 1; i++){
      for (int j = 0; j < nCoords; j++){
        fevals[i].arg[j] = start[i][coords[j]];
      }
    }
  }
  // Perform initial evaluation
  for (int i = 0; i < nCoords; i++){
    // Complete the point
    for (int j = 0; j < nCoords; j++){
      x_full[coords[j]] = fevals[i].arg[j];
    }
    // Evaluate it
    fevals[i].val = -(*fun)(x_full, xx, d, r, n, dirs, prjs,
                      exact, kPar, m1, m2);
  }
  // Main optimization cycle
  int iter = 0;
  bool converged = false;
  while(!converged){
    iter++;
    // cout << "Iteration " << iter << ": ";
    // Evaluate the new "worst" point
    for (int i = 0; i < nCoords; i++){
      x_full[coords[i]] = fevals[nCoords].arg[i];
    }
    fevals[nCoords].val = -(*fun)(x_full, xx, d, r, n, dirs, prjs,
                            exact, kPar, m1, m2);
    // Order the values
    quick_sort(fevals, 0, nCoords);
    // DEBUG -------------------------------------------------------------------
    // Print current simplex
    // cout << "Current simplex consists of vertices:" << endl;
    // for (int i = 0; i < nCoords + 1; i++){
    //   cout << fevals[i].val << ":\t";
    //   for (int j = 0; j < nCoords; j++){
    //     cout << fevals[i].arg[j] << " ";
    //   }
    //   cout << endl;
    // }
    // DEBUG (end) -------------------------------------------------------------
    // Check convergence criterium
    double maxDist = 0;
    for (int i = 0; i < nCoords; i++){
      double tmpVal = fabs(fevals[0].arg[i] - fevals[nCoords].arg[i]);
      if (tmpVal > maxDist){
        maxDist = tmpVal;
      }
    }
    if (maxDist < eps){
      converged = true;
      break;
    }
    // Calculate the centroid of nCoords best points
    for (int i = 0; i < nCoords; i++){
      x_o[i] = 0;
      for (int j = 0; j < nCoords; j++){
        x_o[i] += fevals[j].arg[i];
      }
      x_o[i] /= (double)nCoords;
    }
    // Calculate and evaluate reflected point
    for (int i = 0; i < nCoords; i++){
      x_r[i] = x_o[i] + alpha * (x_o[i] - fevals[nCoords].arg[i]);
      x_full[coords[i]] = x_r[i];
    }
    f_r = -(*fun)(x_full, xx, d, r, n, dirs, prjs,
            exact, kPar, m1, m2);
    // Choose what to do
    if ((fevals[0].val <= f_r) && (f_r < fevals[nCoords - 1].val)){
      // Reflection
      // cout << "Reflection" << endl;
      for (int i = 0; i < nCoords; i++){
        fevals[nCoords].arg[i] = x_r[i];
      }
    }else{
      if (f_r < fevals[0].val){
        // Calculate and evaluate expanded point
        for (int i = 0; i < nCoords; i++){
          x_e[i] = x_r[i] + gamma * (x_r[i] - x_o[i]);
          x_full[coords[i]] = x_e[i];
        }
        f_e = -(*fun)(x_full, xx, d, r, n, dirs, prjs,
                exact, kPar, m1, m2);
        if (f_e < f_r){
          // Expansion
          // cout << "Expansion" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_e[i];
          }
        }else{
          // Still (just) reflection
          // cout << "Reflection" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_r[i];
          }
        }
      }else{
        // Calculate and evaluate contracted point
        for (int i = 0; i < nCoords; i++){
          x_c[i] = x_o[i] + rho * (fevals[nCoords].arg[i] - x_o[i]);
          x_full[coords[i]] = x_c[i];
        }
        f_c = -(*fun)(x_full, xx, d, r, n, dirs, prjs,
                exact, kPar, m1, m2);
        if (f_c < fevals[nCoords].val){
          // Contraction
          // cout << "Contraction" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_c[i];
          }
        }else{
          // Reduction
          // cout << "Reduction" << endl;
          for (int i = 1; i < nCoords + 1; i++){
            for (int j = 0; j < nCoords; j++){
              fevals[i].arg[j] = fevals[0].arg[j] +
                sigma * (fevals[i].arg[j] - fevals[0].arg[j]);
            }
          }
          // Partial evaluation
          for (int i = 1; i < nCoords; i++){
            // Complete the point
            for (int j = 0; j < nCoords; j++){
              x_full[coords[j]] = fevals[i].arg[j];
            }
            // Evaluate it
            fevals[i].val = -(*fun)(x_full, xx, d, r, n, dirs, prjs,
                              exact, kPar, m1, m2);
          }
        }
      }
    }
  }
  // Extract results
  if (converged){
    memcpy(output, x, d * sizeof(double));
    for (int i = 0; i < nCoords; i++){
      output[coords[i]] = fevals[0].arg[i];
    }
  }
  double res = -fevals[0].val;
  // Release memory
  for (int i = 0; i < nCoords + 1; i++){
    delete[] fevals[i].arg;
  }
  delete[] fevals;
  delete[] x_o;
  delete[] x_r;
  delete[] x_e;
  delete[] x_c;
  delete[] x_full;
  return res;
}

double NelderMeadConditionalProj(
    double (*fun)(double* x, int d, int r, double** dirs, double* medians,
            double* MADs),
    double* x, int d, int r, double** dirs, double* medians, double* MADs,
    int* coords, int nCoords, double** start, double* output,
    double alpha, double gamma, double rho, double sigma, double eps){
  // Allocate memory
  Feval* fevals = new Feval[nCoords + 1]; // function evaluations
  for (int i = 0; i < nCoords + 1; i++){
    fevals[i].arg = new double[nCoords];
    for (int j = 0; j < nCoords; j++){
      fevals[i].arg[j] = 0;
    }
  }
  double* x_o = new double[nCoords]; // centroid
  double* x_r = new double[nCoords]; double f_r; // reflected point
  double* x_e = new double[nCoords]; double f_e; // expanded point
  double* x_c = new double[nCoords]; double f_c; // contracted point
  double* x_full = new double[d]; // complete point to evaluate
  memcpy(x_full, x, d * sizeof(double)); // the condition (template) point
  // Generate initial simplex
  if (start == 0){
    for (int i = 0; i < nCoords; i++){
      fevals[i].arg[i] = 1;
    }
  }else{
    for (int i = 0; i < nCoords + 1; i++){
      for (int j = 0; j < nCoords; j++){
        fevals[i].arg[j] = start[i][coords[j]];
      }
    }
  }
  // Perform initial evaluation
  for (int i = 0; i < nCoords; i++){
    // Complete the point
    for (int j = 0; j < nCoords; j++){
      x_full[coords[j]] = fevals[i].arg[j];
    }
    // Evaluate it
    fevals[i].val = -(*fun)(x_full, d, r, dirs, medians, MADs);
  }
  // Main optimization cycle
  int iter = 0;
  bool converged = false;
  while(!converged && iter < 1000){
    iter++;
    // cout << "Iteration " << iter << ": ";
    // Evaluate the new "worst" point
    for (int i = 0; i < nCoords; i++){
      x_full[coords[i]] = fevals[nCoords].arg[i];
    }
    fevals[nCoords].val = -(*fun)(x_full, d, r, dirs, medians, MADs);
    // for (int i = 0; i < nCoords + 1; i++){
    //   cout << fevals[i].val << " ";
    // }
    // cout << endl;
    // Order the values
    quick_sort(fevals, 0, nCoords);
    // Check convergence criterium
    double maxDist = 0;
    for (int i = 0; i < nCoords; i++){
      double tmpVal = fabs(fevals[0].arg[i] - fevals[nCoords].arg[i]);
      if (tmpVal > maxDist){
        maxDist = tmpVal;
      }
    }
    if (maxDist < eps){
      converged = true;
      break;
    }
    // Calculate the centroid
    for (int i = 0; i < nCoords; i++){
      x_o[i] = 0;
      for (int j = 0; j < nCoords; j++){
        x_o[i] += fevals[j].arg[i];
      }
      x_o[i] /= (double)nCoords;
    }
    // Calculate and evaluate reflected point
    for (int i = 0; i < nCoords; i++){
      x_r[i] = x_o[i] + alpha * (x_o[i] - fevals[nCoords].arg[i]);
      x_full[coords[i]] = x_r[i];
    }
    f_r = -(*fun)(x_full, d, r, dirs, medians, MADs);
    // Choose what to do
    if ((fevals[0].val <= f_r) && (f_r < fevals[nCoords - 1].val)){
      // Reflection
      // cout << "Reflection" << endl;
      for (int i = 0; i < nCoords; i++){
        fevals[nCoords].arg[i] = x_r[i];
      }
    }else{
      if (f_r < fevals[0].val){
        // Calculate and evaluate expanded point
        for (int i = 0; i < nCoords; i++){
          x_e[i] = x_r[i] + gamma * (x_r[i] - x_o[i]);
          x_full[coords[i]] = x_e[i];
        }
        f_e = -(*fun)(x_full, d, r, dirs, medians, MADs);
        if (f_e < f_r){
          // Expansion
          // cout << "Expansion" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_e[i];
          }
        }else{
          // Still (just) reflection
          // cout << "Reflection" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_r[i];
          }
        }
      }else{
        // Calculate and evaluate contracted point
        for (int i = 0; i < nCoords; i++){
          x_c[i] = x_o[i] + rho * (fevals[nCoords].arg[i] - x_o[i]);
          x_full[coords[i]] = x_c[i];
        }
        f_c = -(*fun)(x_full, d, r, dirs, medians, MADs);
        if (f_c < fevals[nCoords].val){
          // Contraction
          // cout << "Contraction" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_c[i];
          }
        }else{
          // Reduction
          // cout << "Reduction" << endl;
          for (int i = 1; i < nCoords + 1; i++){
            for (int j = 0; j < nCoords; j++){
              fevals[i].arg[j] = fevals[0].arg[j] +
                sigma * (fevals[i].arg[j] - fevals[0].arg[j]);
            }
          }
          // Partial evaluation
          for (int i = 1; i < nCoords; i++){
            // Complete the point
            for (int j = 0; j < nCoords; j++){
              x_full[coords[j]] = fevals[i].arg[j];
            }
            // Evaluate it
            fevals[i].val = -(*fun)(x_full, d, r, dirs, medians, MADs);
          }
        }
      }
    }
  }
  // Extract results
  if (converged){
    memcpy(output, x, d * sizeof(double));
    for (int i = 0; i < nCoords; i++){
      output[coords[i]] = fevals[0].arg[i];
    }
  }
  double res = -fevals[0].val;
  // Release memory
  for (int i = 0; i < nCoords + 1; i++){
    delete[] fevals[i].arg;
  }
  delete[] fevals;
  delete[] x_o;
  delete[] x_r;
  delete[] x_e;
  delete[] x_c;
  delete[] x_full;
  return res;
}

int LinearConditionalZonoid(double* x, int d, double** data, int n,
                 int* coords, int nCoords, double* output){
  // If nothing to do
  if (nCoords == 0){
    memcpy(output, x, d * sizeof(double));
    return 0;
  }
  // If the trivial case, without constraints, return the average
  if (nCoords == d){
    for (int i = 0; i < d; i++){
      output[i] = 0;
      for (int j = 0; j < n; j++){
        output[i] += data[j][i];
      }
      output[i] /= (double)n;
    }
    return 0;
  }
  // The expected case (0 < nCoords < d): run linea programming
  // Fix coordinates
  int nFixedCoords = d - nCoords;
  int* fixedCoords = new int[nFixedCoords];
  int indexFixedCoords = 0;
  int indexCoords = 0;
  for (int i = 0; i < d; i++){
    if (indexCoords < nCoords && coords[indexCoords] == i){
      indexCoords++;
    }else{
      fixedCoords[indexFixedCoords++] = i;
    }
  }
  // Create LP-structures
  glp_prob *lp;
  lp = glp_create_prob();
  glp_term_out(GLP_OFF);
  // Add rows
  int nrow = nFixedCoords + 1 + n;
  glp_add_rows(lp, nrow);
  for (int i = 0; i < nFixedCoords; i++){
    glp_set_row_bnds(lp, i + 1, GLP_FX, x[fixedCoords[i]], x[fixedCoords[i]]);
  }
  glp_set_row_bnds(lp, nFixedCoords + 1, GLP_FX, 1., 1.);
  for (int i = 0; i < n; i++){
    glp_set_row_bnds(lp, (nFixedCoords + 1) + i + 1, GLP_LO, 0., 0.);
  }
  // Add columns
  int ncol = n + 1;
  glp_add_cols(lp, ncol);
  glp_set_col_bnds(lp, 1, GLP_LO, 0., 0.);
  glp_set_obj_coef(lp, 1, 1.);
  for (int i = 1; i < n + 1; i++){
    glp_set_col_bnds(lp, i + 1, GLP_LO, 0., 0.);
    glp_set_obj_coef(lp, i + 1, 0.);
  }
  // Set constraint matrix
  int ne = nrow * ncol;
  int *rowIndices = new int[ne];
  int *colIndices = new int[ne];
  double *values = new double[ne];
  // Make all entries to zeros first and fill in indices
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      rowIndices[i * ncol + j] = i + 1;
      colIndices[i * ncol + j] = j + 1;
      values[i * ncol + j] = 0.0;
    }
  }
  // Now fill the nonzero entries
  // The data constraint(s)
  for (int i = 0; i < n; i++){
    for (int j = 0; j < nFixedCoords; j++){
      values[j * (n + 1) + (i + 1)] = data[i][fixedCoords[j]];
    }
  }
  // The multiplier constraint
  for (int i = 1; i < n + 1; i++){
    values[nFixedCoords * (n + 1) + i] = 1.;
  }
  // The minimization constraints
  for (int i = 0; i < n; i++){
    values[(nFixedCoords + 1) * (n + 1) + i * (n + 1)] = 1.;
    values[(nFixedCoords + 1) * (n + 1) + i * (n + 1) + (i + 1)] = -1.;
  }
  // Set the direction of optimization
  glp_set_obj_dir(lp, GLP_MIN);
  // DEBUG: output the constraint matrix
  // for (int i = 0; i < nFixedCoords + 1 + n; i++){
  //   for (int j = 0; j < n + 1; j++){
  //     cout << values[i * (n + 1) + j] << " ";
  //   }
  //   cout << endl;
  // }
  // DEBUG (end)
  // Load constraint matrix
  glp_load_matrix(lp, ne, &rowIndices[-1], &colIndices[-1], &values[-1]);
  // Execute linear solver
  glp_simplex(lp, NULL);
  // Collect the point
  if (glp_get_status(lp) == GLP_OPT){
    memcpy(output, x, d * sizeof(double));
    for (int i = 0; i < nCoords; i++){
      output[coords[i]] = 0;
      for (int j = 0; j < n; j++){
        output[coords[i]] += glp_get_col_prim(lp, j + 2) * data[j][coords[i]];
      }
    }
  }
  // Release memory
  delete[] fixedCoords;
  delete[] rowIndices;
  delete[] colIndices;
  delete[] values;
  glp_delete_prob(lp);
  return 0;
}

int LinearConditionalZonoidAll(double** xx, int n, int d,
                               int* points, int nPoints,
                               int** coords, int* nCoords,
                               double* outputs){
  // Create a list of coordinate structures
  Coords* fullCoords = new Coords[nPoints];
  int curTotalIndex = 0; // index in the total "output"
  for (int i = 0; i < nPoints; i++){
    // Fill existing fields (coordinates to impute)
    fullCoords[i].pointIndex = points[i];
    fullCoords[i].totalIndex = curTotalIndex;
    fullCoords[i].d = d;
    fullCoords[i].nCoords = nCoords[i];
    fullCoords[i].coords = new int[nCoords[i]];
    memcpy(fullCoords[i].coords, coords[i], nCoords[i] * sizeof(int));
    fullCoords[i].fixedCoords = new int[d - nCoords[i]];
    // Fill nonexisting fields (present coordinates) by exclusion
    int indexFixedCoords = 0;
    int indexCoords = 0;
    for (int j = 0; j < d; j++){
      if (indexCoords < nCoords[i] && coords[i][indexCoords] == j){
        indexCoords++;
      }else{
        fullCoords[i].fixedCoords[indexFixedCoords++] = j;
      }
    }
    curTotalIndex += nCoords[i];
  }
  // Sort coordinate structures
  quick_sort(fullCoords, 0, nPoints - 1);
  // DEBUG: The sorted coordinates
  // for (int i = 0; i < nPoints; i++){
  //   for (int j = 0; j < fullCoords[i].d - fullCoords[i].nCoords; j++){
  //     cout << fullCoords[i].fixedCoords[j] << " ";
  //   }
  //   cout << endl;
  // }
  // Loop through points
  int lpNFixed = 0; // number of constraints of the existing LP-problem
  // The instance of the LP-problem
  glp_prob *lp;
  int *rowIndices;
  int *colIndices;
  double *values;
  int ne;
  for (int i = 0; i < nPoints; i++){
    // Check the current needed number of constraints
    int nFixedCoords = fullCoords[i].d - fullCoords[i].nCoords;
    if (nFixedCoords == d){
      continue;
    }
    if (nFixedCoords == 0){
      // If the trivial case, without constraints, return the average
      for (int j = 0; j < d; j++){
        outputs[fullCoords[i].totalIndex + j] = 0;
        for (int k = 0; k < n; k++){
          outputs[fullCoords[i].totalIndex + j] += xx[k][j];
        }
        outputs[fullCoords[i].totalIndex + j] /= (double)n;
      }
      continue;
    }
    // The expected case (0 < nCoords[i] < d): run linear programming
    if (nFixedCoords > lpNFixed){ // if we need more constraints than we have
      // Recreate the problem from scratch
      if (lpNFixed > 0){
        delete[] rowIndices;
        delete[] colIndices;
        delete[] values;
        glp_delete_prob(lp);
      }
      // Create LP-structure
      lp = glp_create_prob();
      glp_term_out(GLP_OFF);
      // Add rows
      int nrow = nFixedCoords + 1 + n;
      glp_add_rows(lp, nrow);
      // Right-side constraint(s) of the current point
      for (int j = 0; j < nFixedCoords; j++){
        glp_set_row_bnds(lp, j + 1, GLP_FX,
                         xx[fullCoords[i].pointIndex]
                           [fullCoords[i].fixedCoords[j]],
                         xx[fullCoords[i].pointIndex]
                           [fullCoords[i].fixedCoords[j]]);
      }
      // Right-side weighted-mean constraint
      glp_set_row_bnds(lp, nFixedCoords + 1, GLP_FX, 1., 1.);
      // Right-side minimization constraint
      for (int j = 0; j < n; j++){
        glp_set_row_bnds(lp, (nFixedCoords + 1) + j + 1, GLP_LO, 0., 0.);
      }
      // Add columns
      int ncol = n + 1;
      glp_add_cols(lp, ncol);
      glp_set_col_bnds(lp, 1, GLP_LO, 0., 0.);
      glp_set_obj_coef(lp, 1, 1.);
      for (int j = 1; j < n + 1; j++){
        glp_set_col_bnds(lp, j + 1, GLP_LO, 0., 0.);
        glp_set_obj_coef(lp, j + 1, 0.);
      }
      // Set constraint matrix
      ne = nrow * ncol;
      rowIndices = new int[ne];
      colIndices = new int[ne];
      values = new double[ne];
      // Make all entries to zeros first and fill in indices
      for (int j = 0; j < nrow; j++){
        for (int k = 0; k < ncol; k++){
          rowIndices[j * ncol + k] = j + 1;
          colIndices[j * ncol + k] = k + 1;
          values[j * ncol + k] = 0.0;
        }
      }
      // Now fill the nonzero entries
      // The data constraint(s)
      for (int j = 0; j < n; j++){
        for (int k = 0; k < nFixedCoords; k++){
          values[k * (n + 1) + (j + 1)] = xx[j][fullCoords[i].fixedCoords[k]];
        }
      }
      // The weighted-mean multiplier constraint
      for (int j = 1; j < n + 1; j++){
        values[nFixedCoords * (n + 1) + j] = 1.;
      }
      // The minimization constraints
      for (int j = 0; j < n; j++){
        values[(nFixedCoords + 1) * (n + 1) + j * (n + 1)] = 1.;
        values[(nFixedCoords + 1) * (n + 1) + j * (n + 1) + (j + 1)] = -1.;
      }
      // Set the direction of optimization
      glp_set_obj_dir(lp, GLP_MIN);
      // Update the number of constraints in the existing problem
      lpNFixed = nFixedCoords;
      // Load constraint matrix
      glp_load_matrix(lp, ne, &rowIndices[-1], &colIndices[-1], &values[-1]);
    }else{ // ... we need exactly as many constraints as we have
      if (Compare(fullCoords[i - 1], fullCoords[i])){ // if constraints differ
        // Change the values on both left and right side
        for (int j = 0; j < n; j++){
          for (int k = 0; k < nFixedCoords; k++){
            values[k * (n + 1) + (j + 1)] = xx[j][fullCoords[i].fixedCoords[k]];
          }
        }
        // Right-side constraint(s) of the current point
        for (int j = 0; j < nFixedCoords; j++){
          glp_set_row_bnds(lp, j + 1, GLP_FX,
                           xx[fullCoords[i].pointIndex]
                             [fullCoords[i].fixedCoords[j]],
                           xx[fullCoords[i].pointIndex]
                             [fullCoords[i].fixedCoords[j]]);
        }
        // Load constraint matrix
        glp_load_matrix(lp, ne, &rowIndices[-1], &colIndices[-1], &values[-1]);
      }else{ // ... exactly the same set of constraints
        // Change the right side only
        for (int j = 0; j < nFixedCoords; j++){
          glp_set_row_bnds(lp, j + 1, GLP_FX,
                           xx[fullCoords[i].pointIndex]
                             [fullCoords[i].fixedCoords[j]],
                           xx[fullCoords[i].pointIndex]
                             [fullCoords[i].fixedCoords[j]]);
        }
      }
    }
    // Execute linear solver
    glp_simplex(lp, NULL);
    // Collect the point coordinates
    if (glp_get_status(lp) == GLP_OPT){
      for (int j = 0; j < fullCoords[i].nCoords; j++){
        outputs[fullCoords[i].totalIndex + j] = 0;
        for (int k = 0; k < n; k++){
          outputs[fullCoords[i].totalIndex + j] += glp_get_col_prim(lp, k + 2) *
            xx[k][fullCoords[i].coords[j]];
        }
      }
    }
  }
  // Release memory
  for (int i = 0; i < nPoints; i++){
    delete[] fullCoords[i].coords;
    delete[] fullCoords[i].fixedCoords;
  }
  delete[] fullCoords;
  delete[] rowIndices;
  delete[] colIndices;
  delete[] values;
  glp_delete_prob(lp);
}

double NelderMeadConditionalLocal(
    double (*fun)(double* z, double** xx, int n, int d,
            double locPar, DEPTHNOTION notion),
    double* x, double** xx, int n, int d, double locPar, DEPTHNOTION notion,
    int* coords, int nCoords, double** start, double* output,
    double alpha, double gamma, double rho, double sigma, double eps){
  // Allocate memory
  Feval* fevals = new Feval[nCoords + 1]; // function evaluations
  for (int i = 0; i < nCoords + 1; i++){
    fevals[i].arg = new double[nCoords];
    for (int j = 0; j < nCoords; j++){
      fevals[i].arg[j] = 0;
    }
  }
  double* x_o = new double[nCoords]; // centroid
  double* x_r = new double[nCoords]; double f_r; // reflected point
  double* x_e = new double[nCoords]; double f_e; // expanded point
  double* x_c = new double[nCoords]; double f_c; // contracted point
  double* x_full = new double[d]; // complete point to evaluate
  memcpy(x_full, x, d * sizeof(double)); // the condition (template) point
  // Generate initial simplex
  if (start == 0){
    for (int i = 0; i < nCoords; i++){
      fevals[i].arg[i] = 1;
    }
  }else{
    for (int i = 0; i < nCoords + 1; i++){
      for (int j = 0; j < nCoords; j++){
        fevals[i].arg[j] = start[i][coords[j]];
      }
    }
  }
  // Perform initial evaluation
  for (int i = 0; i < nCoords; i++){
    // Complete the point
    for (int j = 0; j < nCoords; j++){
      x_full[coords[j]] = fevals[i].arg[j];
    }
    // Evaluate it
    fevals[i].val = -(*fun)(x_full, xx, n, d, locPar, notion);
  }
  // Main optimization cycle
  int iter = 0;
  bool converged = false;
  while(!converged){
    iter++;
    // cout << "Iteration " << iter << ": ";
    // Evaluate the new "worst" point
    for (int i = 0; i < nCoords; i++){
      x_full[coords[i]] = fevals[nCoords].arg[i];
    }
    fevals[nCoords].val = -(*fun)(x_full, xx, n, d, locPar, notion);
    // Order the values
    quick_sort(fevals, 0, nCoords);
    // DEBUG -------------------------------------------------------------------
    // Print current simplex
    //cout << "Current simplex consists of vertices:" << endl;
    //for (int i = 0; i < d + 1; i++){
    //  cout << fevals[i].val << ":\t";
    //  for (int j = 0; j < d; j++){
    //    cout << fevals[i].arg[j] << " ";
    //  }
    //  cout << endl;
    //}
    // DEBUG (end) -------------------------------------------------------------
    // Check convergence criterium
    double maxDist = 0;
    for (int i = 0; i < nCoords; i++){
      double tmpVal = fabs(fevals[0].arg[i] - fevals[nCoords].arg[i]);
      if (tmpVal > maxDist){
        maxDist = tmpVal;
      }
    }
    if (maxDist < eps){
      converged = true;
      break;
    }
    // Calculate the centroid of nCoords best points
    for (int i = 0; i < nCoords; i++){
      x_o[i] = 0;
      for (int j = 0; j < nCoords; j++){
        x_o[i] += fevals[j].arg[i];
      }
      x_o[i] /= (double)nCoords;
    }
    // Calculate and evaluate reflected point
    for (int i = 0; i < nCoords; i++){
      x_r[i] = x_o[i] + alpha * (x_o[i] - fevals[nCoords].arg[i]);
      x_full[coords[i]] = x_r[i];
    }
    f_r = -(*fun)(x_full, xx, n, d, locPar, notion);
    // Choose what to do
    if ((fevals[0].val <= f_r) && (f_r < fevals[nCoords - 1].val)){
      // Reflection
      // cout << "Reflection" << endl;
      for (int i = 0; i < nCoords; i++){
        fevals[nCoords].arg[i] = x_r[i];
      }
    }else{
      if (f_r < fevals[0].val){
        // Calculate and evaluate expanded point
        for (int i = 0; i < nCoords; i++){
          x_e[i] = x_r[i] + gamma * (x_r[i] - x_o[i]);
          x_full[coords[i]] = x_e[i];
        }
        f_e = -(*fun)(x_full, xx, n, d, locPar, notion);
        if (f_e < f_r){
          // Expansion
          // cout << "Expansion" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_e[i];
          }
        }else{
          // Still (just) reflection
          // cout << "Reflection" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_r[i];
          }
        }
      }else{
        // Calculate and evaluate contracted point
        for (int i = 0; i < nCoords; i++){
          x_c[i] = x_o[i] + rho * (fevals[nCoords].arg[i] - x_o[i]);
          x_full[coords[i]] = x_c[i];
        }
        f_c = -(*fun)(x_full, xx, n, d, locPar, notion);
        if (f_c < fevals[nCoords].val){
          // Contraction
          // cout << "Contraction" << endl;
          for (int i = 0; i < nCoords; i++){
            fevals[nCoords].arg[i] = x_c[i];
          }
        }else{
          // Reduction
          // cout << "Reduction" << endl;
          for (int i = 1; i < nCoords + 1; i++){
            for (int j = 0; j < nCoords; j++){
              fevals[i].arg[j] = fevals[0].arg[j] +
                sigma * (fevals[i].arg[j] - fevals[0].arg[j]);
            }
          }
          // Partial evaluation
          for (int i = 1; i < nCoords; i++){
            // Complete the point
            for (int j = 0; j < nCoords; j++){
              x_full[coords[j]] = fevals[i].arg[j];
            }
            // Evaluate it
            fevals[i].val = -(*fun)(x_full, xx, n, d, locPar, notion);
          }
        }
      }
    }
  }
  // Extract results
  if (converged){
    memcpy(output, x, d * sizeof(double));
    for (int i = 0; i < nCoords; i++){
      output[coords[i]] = fevals[0].arg[i];
    }
  }
  double res = -fevals[0].val;
  // Release memory
  for (int i = 0; i < nCoords + 1; i++){
    delete[] fevals[i].arg;
  }
  delete[] fevals;
  delete[] x_o;
  delete[] x_r;
  delete[] x_e;
  delete[] x_c;
  delete[] x_full;
  return res;
}
