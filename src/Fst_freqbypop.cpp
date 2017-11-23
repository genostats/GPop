// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"

using namespace Rcpp;
using namespace RcppParallel;

uint8_t NN[256] = {
4, 3, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};

uint8_t NNN[256] = {
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
2, 3, 2, 2, 3, 4, 3, 3, 2, 3, 2, 2, 2, 3, 2, 2, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0};


List countbypop(matrix4 & A, IntegerVector pop, int beg, int end) {

  int p=end-beg+1;
  int r=max(pop);
  NumericMatrix NA(p,r);
  NumericMatrix N0(p,r);
  NumericMatrix N1(p,r);
  NumericMatrix N2(p,r);
  
  for(int i = 1; i <= r; i++) {
    // création vecteur de masques pop
    uint8_t * pop_ = new uint8_t[A.true_ncol];
    std::fill(pop_, pop_+A.true_ncol, 0);  
    int nbo = 0;
    for(size_t j = 0; j < A.true_ncol; j++) {
      if(pop(4*j)!=i) { pop_[j] |= 3; nbo++; }
      if(4*j+1 < A.ncol && pop(4*j+1)!=i) { pop_[j] |= 12;  nbo++; }
      if(4*j+2 < A.ncol && pop(4*j+2)!=i) { pop_[j] |= 48;  nbo++; }
      if(4*j+3 < A.ncol && pop(4*j+3)!=i) { pop_[j] |= 192; nbo++; }
    }

	for(int k = beg; k <= end; k++) {
      for(size_t j = 0; j < A.true_ncol; j++) {
        uint8_t d = A.data[k-beg][j];
        d |= pop_[j];
        N0(k-beg,i-1) += NN[d];
        NA(k-beg,i-1) += NN[255-d];
        N1(k-beg,i-1) += NNN[d];
        N2(k-beg,i-1) += NNN[255-d];
      }
	  // le masque met les AUTRES à NA et ne touche pas à la bordure
      NA(k-beg,i-1) -= (4*A.true_ncol - A.ncol) + nbo;
    }
  }

  List L;

  L["N0"] = N0;
  L["N1"] = N1;
  L["N2"] = N2;
  L["NA"] = NA;

  return L;
}

//[[Rcpp::export]]
List countbypop(XPtr<matrix4> p_A,  IntegerVector pop, int beg, int end) {
  return countbypop(*p_A, pop, beg, end);
}

RcppExport SEXP gg_Fst_countbypop(SEXP p_ASEXP, SEXP popSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pop(popSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    __result = Rcpp::wrap(countbypop(p_A, pop, beg, end));
    return __result;
END_RCPP
}


