/******************************************************************************/
/* File:             imputation.h                                             */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains declarations of general functions for imputation.                 */
/*                                                                            */
/******************************************************************************/

int selectPointCoords(int* totalCoords, int nTotalCoords, int d,
                      int** points, int* nPoints, int*** coords, int** nCoords,
                      int** coordsRaw);
void imputeAllZonoid(double** xx, int n, int d,
                     double** mm, int m,
                     int* points, int nPoints, int** coords, int* nCoords,
                     double* imputed);
