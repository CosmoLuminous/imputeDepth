/******************************************************************************/
/* File:             condOpt.h                                                */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains declarations of functions for conditional depth maximization.     */
/*                                                                            */
/******************************************************************************/

double NelderMeadConditionalSptl(
    double (*fun)(double* x, double** xx, int n, int d, bool sd, double* w),
    double* x, double** xx, int n, int d, bool sd, double* w,
    int* coords, int nCoords, double** start, double* output,
    double alpha = 1., double gamma = 2., double rho = 0.5, double sigma = 0.5,
    double eps = 0.00001);
double NelderMeadConditionalHfsp(
    double (*fun)(double* x, int d, int r, int n, double** dirs, double** prjs),
    double* x, int d, int r, int n, double** dirs, double** prjs,
    int* coords, int nCoords, double** start, double* output,
    double alpha = 1., double gamma = 2., double rho = 0.5, double sigma = 0.5,
    double eps = 0.00001);
double NelderMeadCondHfspExtr(
    double (*fun)(double* x, double** xx, int d, int r, int n,
            double** dirs, double** prjs,
            bool exact, int kPar, double m1, double m2),
            double* x, double** xx, int d, int r, int n,
            double** dirs, double** prjs,
            bool exact, int kPar, double m1, double m2,
            int* coords, int nCoords, double** start, double* output,
            double alpha = 1., double gamma = 2., double rho = 0.5,
            double sigma = 0.5, double eps = 0.00001);
double NelderMeadConditionalHfspEx(
    double (*fun)(double* z, double** xx, int n, int d),
    double* x, double** xx, int n, int d,
    int* coords, int nCoords, double** start, double* output,
    double alpha = 1., double gamma = 2., double rho = 0.5, double sigma = 0.5,
    double eps = 0.00001);
double NelderMeadConditionalProj(
    double (*fun)(double* x, int d, int r, double** dirs, double* medians,
            double* MADs),
    double* x, int d, int r, double** dirs, double* medians, double* MADs,
    int* coords, int nCoords, double** start, double* output,
    double alpha = 1., double gamma = 2., double rho = 0.5, double sigma = 0.5,
    double eps = 0.00001);
int LinearConditionalZonoid(double* x, int d, double** data, int n,
                            int* coords, int nCoords, double* output);
int LinearConditionalZonoidAll(double** data, int n, int d,
                               int* points, int nPoints,
                               int** coords, int* nCoords,
                               double* outputs);
double NelderMeadHfspEx(
    double (*fun)(double* data, int n, int d,
            double* z, int* coords, int nCoords),
            double* data, int n, int d, int* coords, int nCoords,
            double* start, double* output,
            double alpha = 1., double gamma = 2., double rho = 0.5,
            double sigma = 0.5, double eps = 0.001);
double NelderMeadConditionalLocal(
    double (*fun)(double* z, double** xx, int n, int d,
            double locPar, DEPTHNOTION notion),
            double* x, double** xx, int n, int d, double locPar, DEPTHNOTION notion,
            int* coords, int nCoords, double** start, double* output,
            double alpha = 1., double gamma = 2., double rho = 0.5,
            double sigma = 0.5, double eps = 0.001);
