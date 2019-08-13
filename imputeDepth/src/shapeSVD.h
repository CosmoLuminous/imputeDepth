/******************************************************************************/
/* File:             shapeSVD.h                                               */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains declarations of functions singular value decomposition of the     */
/* shape matrix.                                                              */
/*                                                                            */
/******************************************************************************/

void an_svd(double* X, int n, int d, double* s, double* u, double* vt);
void svd_depth_proj(double** x, int n, int d, int r, int typeCenter,
                    double* s, double* vt, double* scale, double* center);
void svd_depth_hfsp(double** x, int n, int d, int r,
                    int typeDepth, int typeCenter,
                    double* s, double* vt, double* scale, double* center);
void get_w(double** xx, int n, int d, double* mu, double* wRaw);
