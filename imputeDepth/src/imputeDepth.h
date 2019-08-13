/******************************************************************************/
/* File:             imputeDepth.h                                            */
/* Created by:       Pavlo Mozharovskyi                                       */
/* Last revised:     23.12.2016                                               */
/*                                                                            */
/* The general header for R-package imputeDepth.                              */
/*                                                                            */
/******************************************************************************/

#pragma once

#include <math.h>
#include <cstdlib>
#include <string.h>
#include <vector>
#include <algorithm>
#include <R_ext/Lapack.h>

using namespace std;

#ifndef _MSC_VER
#include <Rcpp.h>
using namespace Rcpp;
#endif

#include "common.h"
#include "depthCalc.h"
#include "condOpt.h"
#include "shapeSVD.h"
#include "HD.h"
#include "ZonoidDepth.h"
#include "LocalDepth.h"
#include "Mahalanobis.h"
#include "imputation.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/random.hpp>
#include <boost/numeric/ublas/matrix.hpp>
extern "C" {
  #include "glpk.h"
}

#include <iostream>
