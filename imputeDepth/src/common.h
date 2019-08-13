/******************************************************************************/
/* File:             common.h                                                 */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains declarations of common structures and search functions for        */
/* R-package imputeDepth.                                                     */
/*                                                                            */
/******************************************************************************/

enum DEPTHNOTION {
  NONE = 0,
  HALFSPACE = 1,
  ZONOID = 2,
  MAHALANOBIS = 3
};

typedef double** TDMatrix;
typedef vector<double> TPoint;
typedef vector<vector<double> > TMatrix;
typedef vector<vector<int> > TIntMatrix;
typedef vector<int> TVariables;

TDMatrix asMatrix(double* arr, int n, int d);

struct SortRec {
  double v;
  TPoint* p;
  SortRec(double v = 0, TPoint* p = NULL) {
    this->v = v;
    this->p = p;
  }
};

struct SortIndex{
  double value;
  int index;
};

struct Feval{
  double* arg;
  double val;
};

struct Coords{
  int pointIndex;
  int totalIndex;
  int d;
  int nCoords;
  int* coords;
  int* fixedCoords;
};

template<typename T>
int Compare(T &p1, T &p2) {
  return p1 < p2;
}

template<typename T>
void Swap(T *p1, T *p2) {
  T pTmp = *p1;
  *p1 = *p2;
  *p2 = pTmp;
};

int Compare(SortIndex &si1, SortIndex &si2);

void Swap(SortIndex *si1, SortIndex *si2);

int Compare(Feval &fe1, Feval &fe2);

void Swap(Feval *fe1, Feval *fe2);

int Compare(Coords &co1, Coords &co2);

void Swap(Coords *co1, Coords *co2);

#define elem_type double

const double eps = 10e-8; // precision constant

elem_type quick_select(elem_type arr[], int n);
elem_type kth_smallest(elem_type a[], int n, int k);

/* -------------------------------------------------------------------------- */
/* quickSort from http://www.proggen.org/doku.php?id=algo:quicksort           */
/* (modified, templated)                                                      */
/* -------------------------------------------------------------------------- */
template<typename T>
void quick_sort(T *values, int left, int right, int(*cmp)(T& x, T& y),
                void(*swap)(T* x, T* y)){
  int i = left, j = right; // Z?hlindizes
  // Pivot-Element (Array-Mitte) bestimmen
  T pivot = values[(left + right) >> 1];
  // Solange Paare suchen und tauschen, bis sich die Z?hlindizes "getroffen"
  // haben
  do{
    // Paar finden, dass getauscht werden muss
    while (cmp(values[i], pivot)){++i;}
    while (cmp(pivot, values[j])){--j;}
    // Wenn sich die Z?hlindizes noch nicht "getroffen" haben, dann
    // tauschen und weiterz?hlen. Sollten die Z?hlindizes gleich sein, dann
    // z?hle einfach weiter und unterbrich die Schleife
    if (i < j){swap(&values[i], &values[j]);++i;--j;
    }else{if (i == j){++i;--j;break;}}
  }while (i <= j);
  // Wenn die Teillisten mehr als ein Element enthalten, dann wende quickSort
  // auf sie an
  if (left < j){quick_sort(values, left, j, cmp, swap);}
  if (i < right){quick_sort(values, i, right, cmp, swap);}
}
template<typename T>
void quick_sort(T *values, int left, int right){
  quick_sort(values, left, right, Compare, Swap);
}

int bin_searchl_routine(double* x, int n, double val);
int bin_searchr_routine(double* x, int n, double val);

double* means(TDMatrix X, int n, int d);
TDMatrix cov(TDMatrix X, int n, int d);
double** newM(int n, int d);
void deleteM(TDMatrix X);
bool solveUnique(TDMatrix A, double* b, double* x, int d);
unsigned long long choose(unsigned long long n, unsigned long long k);
