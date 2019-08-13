/******************************************************************************/
/* File:             depthCalc.h                                              */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains declarations of functions for depth calculation.                  */
/*                                                                            */
/******************************************************************************/

#include<vector>

void depths_proj2(double** x, int n, int d, int r, double* depths,
                 double* dirsRaw, double* medians, double* MADs);
void depths_proj(double** x, int n, int d, int r, double* depths);
void depthsn_hfsp(double** x, int n, int d, double** data, int m, int r,
                  double* depths);
void depthsn_sptl(double** x, int n, int d, double** data, int m, bool sd,
                  double* depths);
double depth_hfsp(double* x, int d, int r, int n, double** dirs, double** prjs);
double depth_proj(double* x, int d, int r, double** dirs, double* medians,
                  double* MADs);
double depth_sptl(double* x, double** xx, int n, int d, bool sd, double* w);
double depth_hfsp(double* x, int d, int r, int n, double** dirs, double** prjs);
double depth_hfspExtr(double* x, double** xx, int d, int r, int n,
                      double** dirs, double** prjs,
                      bool exact, int kPar, double m1, double m2);
void get_prjs(double** x, int n, int d, int r,
              double** pDirsRaw, double** pPrjsRaw);
void get_prjs_norm(double** x, int n, int d, int r,
                   double** pDirsRaw, double** pPrjsRaw);
void get_prjs2(double** x, int n, int d, int r,
               double* dirsRaw, double* prjsRaw);
void upd_prjs(double** x, int n, int d, int r,
              double** dirs, double** prjs);
void get_MaMADs(double** x, int n, int d, int r,
                double** pDirsRaw, double** pMedians, double** pMADs);
void get_MaMADs2(double** x, int n, int d, int r,
                 double* dirsRaw, double* medians, double* MADs);
void get_mads(double** x, int n, int d, int r, double** dirs, double* medians,
              double* MADs);
void median_hfsp(double** x, int n, int d, int r, double* median);
void median_proj(double** x, int n, int d, int r, double* median);
double SimplicialDepthEx(double* x, TDMatrix X, int n, int d);
