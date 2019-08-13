/******************************************************************************/
/* File:             common.cpp                                               */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* Contains common structures and search functions for R-package imputeDepth. */
/*                                                                            */
/******************************************************************************/

#include "imputeDepth.h"

TDMatrix asMatrix(double* arr, int n, int d){
  TDMatrix mat = new double*[n];
  for (int i = 0; i < n; i++)
    mat[i] = arr + i*d;
  return mat;
}

int Compare(SortIndex &si1, SortIndex &si2){
  return si2.value < si1.value;
}

void Swap(SortIndex *si1, SortIndex *si2){
  SortIndex siTmp = *si1;
  *si1 = *si2;
  *si2 = siTmp;
}

int Compare(Feval &fe1, Feval &fe2){
  return fe1.val < fe2.val;
}

void Swap(Feval *fe1, Feval *fe2){
  Feval feTmp = *fe1;
  *fe1 = *fe2;
  *fe2 = feTmp;
}

int Compare(Coords &co1, Coords &co2){
  if (co1.d == co2.d){
    if (co1.d - co1.nCoords == co2.d - co2.nCoords){
      for (int i = 0; i < co1.d - co1.nCoords - 1; i++){
        if (co1.fixedCoords[i] != co2.fixedCoords[i]){
          return co1.fixedCoords[i] < co2.fixedCoords[i];
        }
      }
      return co1.fixedCoords[co1.d - co1.nCoords - 1] <
        co2.fixedCoords[co2.d - co2.nCoords - 1];
    }else{
      return co1.d - co1.nCoords < co2.d - co2.nCoords;
    }
  }else{
    return co1.d < co2.d;
  }
}

void Swap(Coords *co1, Coords *co2){
  Coords coTmp = *co1;
  *co1 = *co2;
  *co2 = coTmp;
}

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 *
 *  This code has been downloaded from http://ndevilla.free.fr/median/median/
 */
#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }
elem_type quick_select(elem_type arr[], int n) {
  int low, high ;
  int median;
  int middle, ll, hh;

  low = 0 ; high = n-1 ; median = (low + high) / 2;
  for (;;) {
    if (high <= low) /* One element only */
      return arr[median] ;

    if (high == low + 1) {  /* Two elements only */
      if (arr[low] > arr[high])
      ELEM_SWAP(arr[low], arr[high]) ;
      return arr[median] ;
    }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do ll++; while (arr[low] > arr[ll]) ;
      do hh--; while (arr[hh]  > arr[low]) ;

      if (hh < ll)
        break;

      ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
      low = ll;
    if (hh >= median)
      high = hh - 1;
  }
}

/*---------------------------------------------------------------------------
Function :   kth_smallest()
In       :   array of elements, # of elements in the array, rank k
Out      :   one element
Job      :   find the kth smallest element in the array
Notice   :   use the median() macro defined below to get the median.

Reference:

Author: Wirth, Niklaus
Title: Algorithms + data structures = programs
Publisher: Englewood Cliffs: Prentice-Hall, 1976
Physical description: 366 p.
Series: Prentice-Hall Series in Automatic Computation

This code has been downloaded from http://ndevilla.free.fr/median/median/,
and more precisely: http://ndevilla.free.fr/median/median/src/wirth.c
---------------------------------------------------------------------------*/
elem_type kth_smallest(elem_type a[], int n, int k)
{
  int i,j,l,m ;
  elem_type x ;

  l=0 ; m=n-1 ;
  while (l<m) {
    x=a[k] ;
    i=l ;
    j=m ;
    do {
      while (a[i]<x) i++ ;
      while (x<a[j]) j-- ;
      if (i<=j) {
        ELEM_SWAP(a[i],a[j]) ;
        i++ ; j-- ;
      }
    } while (i<=j) ;
    if (j<k) l=i ;
    if (k<i) m=j ;
  }
  return a[k] ;
}

#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

#undef ELEM_SWAP

int bin_searchl_routine(double* x, int n, double val){
  int min = 0;
  int max = n + 1;
  int mid = (min + max) / 2;
  while((x[mid - 1] > val) || (val >= x[mid])){
    mid = (min + max) / 2;
    if (mid == 0){
      return 0;
    }
    if (mid == n){
      if ((x[0] == val) && (x[mid - 1] == val)){
        // return -1; // coincide with all the points
        return n;
      }else{
        return n;
      }
    }
    if (x[mid - 1] > val){
      max = mid;
    }else{
      min = mid;
    }
  }
  return mid;
}

int bin_searchr_routine(double* x, int n, double val){
  int min = 0;
  int max = n + 1;
  int mid = (min + max) / 2;
  while((x[mid - 1] >= val) || (val > x[mid])){
    mid = (min + max) / 2;
    if (mid == 0){
      if ((x[0] == val) && (x[mid - 1] == val)){
        // return -1; // coincide with all the points
        return n;
      }else{
        return n;
      }
    }
    if (mid == n){
      return 0;
    }
    if (x[mid - 1] >= val){
      max = mid;
    }else{
      min = mid;
    }
  }
  return n - mid;
}

double* means(TDMatrix X, int n, int d) {
  double* ms = new double[d];
  for (unsigned i = 0; i < d; i++) {
    ms[i] = 0.0;
    for (unsigned j = 0; j < n; j++)
      ms[i] += X[j][i];
    ms[i] /= n;
  }
  return ms;
}

double** newM(int n, int d){
  double* a = new double[n*d];
  return asMatrix(a, n, d);
}

TDMatrix cov(TDMatrix X, int n, int d) {
  double* means = new double[d];
  double* dev = new double[d];
  // zeroing TDMatrix
  TDMatrix covX = newM(d, d);
  for (unsigned k = 0; k < d; k++)
    for (unsigned j = 0; j < d; j++)
      covX[k][j] = 0;
  // means
  for (unsigned i = 0; i < d; i++) {
    means[i] = 0.0;
    for (unsigned j = 0; j < n; j++)
      means[i] += X[j][i];
    means[i] /= n;
  }
  for (unsigned i = 0; i < n; i++) {
    // deviations
    for (unsigned k = 0; k < d; k++) {
      dev[k] = X[i][k] - means[k];
    }
    // add to cov
    for (unsigned k = 0; k < d; k++) {
      for (unsigned j = 0; j < d; j++) {
        covX[k][j] += dev[k] * dev[j];
      }
    }
  }
  //scale
  for (unsigned i = 0; i < d; i++) {
    for (unsigned j = 0; j < d; j++) {
      covX[i][j] /= n - 1;
    }
  }
  delete[] means;
  delete[] dev;
  return covX;
}

void deleteM(TDMatrix X){
  delete[] X[0];
  delete[] X;
}

unsigned long long choose(unsigned long long n, unsigned long long k){
  unsigned long long r = n--; unsigned long long d = 2;
  while (d <= k){ r *= n--; r /= d++; }
  return r;
}

/* -------------------------------------------------------------------------- */
/* By Rainer Dyckerhoff, modified by Pavlo Mozharovskyi                       */
/* Solves a uniquely solvable system of linear equations                      */
/* -------------------------------------------------------------------------- */
bool solveUnique(TDMatrix A, double* b, double* x, int d){
  int imax, jmax;
  int* colp = new int[d];
  double amax;
  for (int k = 0; k < d - 1; k++) {
    imax = k;
    amax = abs(A[k][k]);
    colp[k] = k;
    // Spaltenmaximum finden
    for (int i = k + 1; i < d; i++) {
      if (abs(A[i][k]) > amax) {
        amax = abs(A[i][k]);
        imax = i;
      }
    }
    // Spaltenmaximum gleich null => complete pivoting
    if (amax < eps) {
      for (int j = k + 1; j < d; j++) {
        for (int i = k; i < d; i++) {
          if (abs(A[i][j]) > amax) {
            amax = abs(A[i][j]);
            imax = i;
            jmax = j;
          }
        }
      }
      if (amax < eps) {
        delete[] colp;
        return false;
      }
      // Spaltentausch
      for (int i = 0; i < d; i++) {
        double tmp = A[i][k];
        A[i][k] = A[i][jmax];
        A[i][jmax] = tmp;
      }
      colp[k] = jmax;
    }
    // Zeilentausch
    if (imax != k) {
      for (int j = k; j < d; j++) {
        double tmp = A[k][j];
        A[k][j] = A[imax][j];
        A[imax][j] = tmp;
      }
      double tmp = b[k];
      b[k] = b[imax];
      b[imax] = tmp;
    }
    // Elimination
    for (int i = k + 1; i < d; i++) {
      double factor = A[i][k] / A[k][k];
      for (int j = k + 1; j < d; j++){
        A[i][j] -= factor * A[k][j];
      }
      b[i] -= factor * b[k];
    }
  }
  // R?cksubstituition
  colp[d - 1] = d - 1;
  for (int k = d - 1; k >= 0; k--) {
    x[k] = b[k] / A[k][k];
    for (int i = k - 1; i >= 0; i--) b[i] -= x[k] * A[i][k];
  }
  // Spaltenvertauschungen r?ckg?ngig machen
  for (int k = d - 1; k >= 0; k--) {
    if (colp[k] != k) {
      double temp = x[k];
      x[k] = x[colp[k]];
      x[colp[k]] = temp;
    }
  }
  delete[] colp;
  return true;
}
