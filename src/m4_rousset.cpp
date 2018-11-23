// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include "matrix4.h"
#include "m4_kinship_type.h"

using namespace Rcpp;
using namespace RcppParallel;


// ************* [on ne symmétrise pas] ***********

struct paraRouss : public Worker {
  // input and others
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol;
  const size_t sizeK;
  double het_het;
  int share0_het_het;
  int share1_het_het;
  int share2_het_het;

  // output
  Ktype * K;
  Ktype * N;
  Ktype * Share;
  Ktype * SShare;
  Ktype * SSShare;
  
  // constructeurs
  paraRouss(uint8_t ** data, const size_t ncol, const size_t true_ncol, double het_het, int share0_het_het, int share1_het_het, int share2_het_het) : data(data),
           ncol(ncol), true_ncol(true_ncol), het_het(het_het),
           share0_het_het(share0_het_het), share1_het_het(share1_het_het), share2_het_het(share2_het_het),
           sizeK((4*true_ncol)*(4*true_ncol+1)/2) { 
          K = new Ktype[sizeK];  // K is padded to a multiple of 4...
          std::fill(K, K+sizeK, 0);
          N = new Ktype[sizeK];  
          std::fill(N, N+sizeK, 0);
          Share = new Ktype[sizeK];  
          std::fill(Share, Share+sizeK, 0);
          SShare = new Ktype[sizeK];  
          std::fill(SShare, SShare+sizeK, 0);
          SSShare = new Ktype[sizeK];  
          std::fill(SSShare, SSShare+sizeK, 0);
        }
  paraRouss(paraRouss & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), het_het(Q.het_het),
                                    share0_het_het(Q.share0_het_het), share1_het_het(Q.share1_het_het), share2_het_het(Q.share2_het_het),
                                    sizeK(Q.sizeK) {
          K = new Ktype[sizeK];  
          std::fill(K, K+sizeK, 0);
          N = new Ktype[sizeK];  
          std::fill(N, N+sizeK, 0); 
          Share = new Ktype[sizeK];  
          std::fill(Share, Share+sizeK, 0);
          SShare = new Ktype[sizeK];  
          std::fill(SShare, SShare+sizeK, 0);
          SSShare = new Ktype[sizeK];  
          std::fill(SSShare, SSShare+sizeK, 0);
        }

  // destructeur
  ~paraRouss() { 
          delete [] K; 
          delete [] N; 
          delete [] Share; 
          delete [] SShare; 
          delete [] SSShare; 
  }

  // worker !
  void operator()(size_t beg, size_t end) {
    Ktype H[16];
    Ktype HH[16];
    Ktype S[16];
    Ktype SS[16];
    Ktype SSS[16];
    H[0] = 1;        HH[0] = 1;    S[0] = 0;        SS[0] = 0;    SSS[0] = 1;
    H[1] = 0.5;      HH[1] = 1;    S[1] = 0;        SS[1] = 1;    SSS[1] = 0;
    H[2] = 0;        HH[2] = 1;    S[2] = 1;        SS[2] = 0;    SSS[2] = 0;
    H[3] = 0;        HH[3] = 0;    S[3] = 0;        SS[3] = 0;    SSS[3] = 0;
    H[4] = 0.5;      HH[4] = 1;    S[4] = 0;        SS[4] = 1;    SSS[4] = 0;
    H[5] = het_het;  HH[5] = 1;    S[5] = share0_het_het;        SS[5] = share1_het_het;    SSS[5] = share2_het_het;
    H[6] = 0.5;      HH[6] = 1;    S[6] = 0;        SS[6] = 1;    SSS[6] = 0;
    H[7] = 0;        HH[7] = 0;    S[7] = 0;        SS[7] = 0;    SSS[7] = 0;
    H[8] = 0;        HH[8] = 1;    S[8] = 1;        SS[8] = 0;    SSS[8] = 0;
    H[9] = 0.5;      HH[9] = 1;    S[9] = 0;        SS[9] = 1;    SSS[9] = 0;
    H[10] = 1;       HH[10] = 1;   S[10] = 0;       SS[10] = 0;   SSS[10] = 1;
    H[11] = 0;       HH[11] = 0;   S[11] = 0;       SS[11] = 0;   SSS[11] = 0;
    H[12] = 0;       HH[12] = 0;   S[12] = 0;       SS[12] = 0;   SSS[12] = 0;
    H[13] = 0;       HH[13] = 0;   S[13] = 0;       SS[13] = 0;   SSS[13] = 0;
    H[14] = 0;       HH[14] = 0;   S[14] = 0;       SS[14] = 0;   SSS[14] = 0;
    H[15] = 0;       HH[15] = 0;   S[15] = 0;       SS[15] = 0;   SSS[15] = 0;

    for(size_t i = beg; i < end; i++) {
 
      uint8_t * dd = data[i];
      
      size_t k = 0;
      for(size_t j1 = 0; j1 < true_ncol; j1++) {
        uint8_t x1 = dd[j1];
        for(unsigned int ss1 = 0; (ss1 < 4); ss1++) {
          for(size_t j2 = 0; j2 < j1; j2++) {
            uint8_t x2 = dd[j2];
	    N[k] += HH[ (x1&3)*4 + (x2&3) ];
        K[k] += H[ (x1&3)*4 + (x2&3) ];
        Share[k] += S[ (x1&3)*4 + (x2&3) ];
        SShare[k] += SS[ (x1&3)*4 + (x2&3) ];
            SSShare[k++] += SSS[ (x1&3)*4 + (x2&3) ];
	    
	    N[k] += HH[ (x1&3)*4 + ((x2>>2)&3) ];
        K[k] += H[ (x1&3)*4 + ((x2>>2)&3) ];
        Share[k] += S[ (x1&3)*4 + ((x2>>2)&3) ];
        SShare[k] += SS[ (x1&3)*4 + ((x2>>2)&3) ];
            SSShare[k++] += SSS[ (x1&3)*4 + ((x2>>2)&3) ];
	    
        N[k] += HH[ (x1&3)*4 + ((x2>>4)&3) ];
        K[k] += H[ (x1&3)*4 + ((x2>>4)&3) ];
        Share[k] += S[ (x1&3)*4 + ((x2>>4)&3) ];
        SShare[k] += SS[ (x1&3)*4 + ((x2>>4)&3) ];
            SSShare[k++] += SSS[ (x1&3)*4 + ((x2>>4)&3) ];
	    
        N[k] += HH[ (x1&3)*4 + ((x2>>6)&3) ];
        K[k] += H[ (x1&3)*4 + ((x2>>6)&3) ];
        Share[k] += S[ (x1&3)*4 + ((x2>>6)&3) ];
        SShare[k] += SS[ (x1&3)*4 + ((x2>>6)&3) ];
            SSShare[k++] += SSS[ (x1&3)*4 + ((x2>>6)&3) ];
          }
          size_t j2 = j1;
          uint8_t x2 = dd[j2];
          for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
	    /* if (ss2 == ss1) {
	      if (x1&3==3) {k++;} else {
                N[k] += 1;
                K[k] += 1;
                SSShare[k++] += 1;
        }
            } else {*/
              N[k] += HH[ (x1&3)*4 + (x2&3) ];
	          K[k] += H[ (x1&3)*4 + (x2&3) ];
            Share[k] += S[ (x1&3)*4 + (x2&3) ];
            SShare[k] += SS[ (x1&3)*4 + (x2&3) ];
            SSShare[k++] += SSS[ (x1&3)*4 + (x2&3) ];
	    //}
            x2 >>= 2;
          } 
          x1 >>= 2; 
        }
      } 
    }
  }

  // recoller
  void join(const paraRouss & Q) {
    std::transform(K, K + sizeK, Q.K, K, std::plus<Ktype>());
    // autrement dit : K += Q.K;
    std::transform(N, N + sizeK, Q.N, N, std::plus<Ktype>());
    // autrement dit : N += Q.N;
    std::transform(Share, Share + sizeK, Q.Share, Share, std::plus<Ktype>());
    std::transform(SShare, SShare + sizeK, Q.SShare, SShare, std::plus<Ktype>());
    std::transform(SSShare, SSShare + sizeK, Q.SSShare, SSShare, std::plus<Ktype>());
  }


};


List Rousset(XPtr<matrix4> p_A, LogicalVector snps, double het_het, int share0_het_het, int share1_het_het, int share2_het_het, int chunk) {
  int nb_snps = sum(snps);

  if(snps.length() != p_A->nrow) 
    stop("Dimensions mismatch");

  uint8_t ** data = new uint8_t * [nb_snps];
  size_t k = 0;
  for(size_t i = 0; i < p_A->nrow; i++) {
    if(snps[i]) data[k++] = p_A->data[i];
  }

  paraRouss X(data, p_A->ncol, p_A->true_ncol, het_het, share0_het_het, share1_het_het, share2_het_het);
  parallelReduce(0, nb_snps, X, chunk);

  delete [] data;

  NumericMatrix Y(p_A->ncol,p_A->ncol);
  NumericMatrix S0(p_A->ncol,p_A->ncol);
  NumericMatrix S1(p_A->ncol,p_A->ncol);
  NumericMatrix S2(p_A->ncol,p_A->ncol);
  NumericMatrix NN(p_A->ncol,p_A->ncol);
  double Q = 0;
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y(j,i) = (double) X.K[k]/X.N[k];
      Y(i,j) = (double) X.K[k]/X.N[k];
      NN(j,i) = (double) X.N[k];
      NN(i,j) = (double) X.N[k];
      S0(j,i) = (double) X.Share[k]/X.N[k];
      S0(i,j) = (double) X.Share[k]/X.N[k];
      S1(j,i) = (double) X.SShare[k]/X.N[k];
      S1(i,j) = (double) X.SShare[k]/X.N[k];
      S2(j,i) = (double) X.SSShare[k]/X.N[k];
      S2(i,j) = (double) X.SSShare[k]/X.N[k];
      if (i!=j & X.N[k]!=0) Q += (double) X.K[k]/X.N[k];
      k++;
    }
  }
  
  Q /= p_A->ncol*(p_A->ncol-1)/2;
  NumericMatrix Z;
  Z = (Y-Q)/(1-Q);

  List R;
  R["Rousset"] = Z;
  R["F"] = Y;
  R["Qm"] = Q;
  R["Nsnp"] = NN;
  R["Share0"] = S0;
  R["Share1"] = S1;
  R["Share2"] = S2;
  return R;
}


RcppExport SEXP gg_Rousset(SEXP p_ASEXP, SEXP snpsSEXP, SEXP het_hetSEXP, SEXP share0_het_hetSEXP, SEXP share1_het_hetSEXP, SEXP share2_het_hetSEXP, SEXP chunkSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
	Rcpp::traits::input_parameter< LogicalVector >::type snps(snpsSEXP);
        Rcpp::traits::input_parameter< const double >::type het_het(het_hetSEXP);
        Rcpp::traits::input_parameter< int >::type chunk(chunkSEXP );
        Rcpp::traits::input_parameter< int >::type share0_het_het(share0_het_hetSEXP );
        Rcpp::traits::input_parameter< int >::type share1_het_het(share1_het_hetSEXP );
        Rcpp::traits::input_parameter< int >::type share2_het_het(share2_het_hetSEXP );
        List __result = Rousset(p_A, snps, het_het, share0_het_het, share1_het_het, share2_het_het, chunk);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

