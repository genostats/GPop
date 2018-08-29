#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gaston/matrix4.h"

using namespace Rcpp;

/*** 
  Dans le format bed : 00 et 11 pour les deux homoz 
                       10 pour hétéroz
                       01 pour manquant

  Donc il faut une conversion 0 -> 0 ou 00 -> 00
                              1 -> 3    01 -> 11
                              2 -> 1    10 -> 01
                              3 -> 2    11 -> 10 

  on peut le faire simultanément sur les 8 génotypes contenus
  dans un uint8_t
***/

uint8_t bedc[256] = {
  0,   3,   1,   2,  12,  15,  13,  14,   4,   7,   5,   6,   8,  11,   9,  10, 
 48,  51,  49,  50,  60,  63,  61,  62,  52,  55,  53,  54,  56,  59,  57,  58, 
 16,  19,  17,  18,  28,  31,  29,  30,  20,  23,  21,  22,  24,  27,  25,  26, 
 32,  35,  33,  34,  44,  47,  45,  46,  36,  39,  37,  38,  40,  43,  41,  42, 
192, 195, 193, 194, 204, 207, 205, 206, 196, 199, 197, 198, 200, 203, 201, 202, 
240, 243, 241, 242, 252, 255, 253, 254, 244, 247, 245, 246, 248, 251, 249, 250, 
208, 211, 209, 210, 220, 223, 221, 222, 212, 215, 213, 214, 216, 219, 217, 218, 
224, 227, 225, 226, 236, 239, 237, 238, 228, 231, 229, 230, 232, 235, 233, 234, 
 64,  67,  65,  66,  76,  79,  77,  78,  68,  71,  69,  70,  72,  75,  73,  74,  
112, 115, 113, 114, 124, 127, 125, 126, 116, 119, 117, 118, 120, 123, 121, 122, 
 80,  83,  81,  82,  92,  95,  93,  94,  84,  87,  85,  86,  88,  91,  89,  90,  
 96,  99,  97,  98, 108, 111, 109, 110, 100, 103, 101, 102, 104, 107, 105, 106, 
128, 131, 129, 130, 140, 143, 141, 142, 132, 135, 133, 134, 136, 139, 137, 138, 
176, 179, 177, 178, 188, 191, 189, 190, 180, 183, 181, 182, 184, 187, 185, 186, 
144, 147, 145, 146, 156, 159, 157, 158, 148, 151, 149, 150, 152, 155, 153, 154, 
160, 163, 161, 162, 172, 175, 173, 174, 164, 167, 165, 166, 168, 171, 169, 170};


// bed magic numbers : 108 27 1
// [[Rcpp::export]]
XPtr<matrix4> read_bed_file_part(CharacterVector filename, int n_ind, int beg, int end) {
  int n_snp = end - beg + 1;
  std::ifstream file(filename[0], std::ifstream::binary);
  if(!file.is_open()) {
    Rf_error("Cannot open file");
  }
  uint8_t m1, m2, m3;
  file.read(reinterpret_cast<char *>(&m1), 1);
  file.read(reinterpret_cast<char *>(&m2), 1);
  file.read(reinterpret_cast<char *>(&m3), 1);
  if(m1 != 108 || m2 != 27) {
    Rf_error("Not a bed file");
  }
  if(m3 != 1) {
    Rf_error("Not a bed file in SNP major mode");
  }

  XPtr<matrix4> p_A(new matrix4(n_snp, n_ind));
  uint8_t b;
  uint8_t bordermask;
  switch(4*p_A->true_ncol - n_ind) {
    case 0:
      bordermask = 0;
      break;
    case 1:
      bordermask = 192;
      break;
    case 2: 
      bordermask = 240;
      break;
    case 3:
      bordermask = 252;
      break;
    default:
      Rf_error("Some shit hit the fan very hard");
  }
  
  // fast forward to beg...
  file.seekg((beg-1)*p_A->true_ncol, file.cur);
  for(int i = 0; i < n_snp; i++) {
    for(int j = 0; j < p_A->true_ncol; j++) {
      file.read(reinterpret_cast<char *> (&b),1);
      p_A->data[i][j] = bedc[ (int) b ];
    }
    // être sûr d'être bordé de 11 -> NA
    p_A->data[i][p_A->true_ncol - 1] |= bordermask;
  }
  
  file.close();
  return p_A;
}

RcppExport SEXP read_bed_file_part(SEXP filenameSEXP, SEXP n_indSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type n_ind(n_indSEXP);
    Rcpp::traits::input_parameter< int >::type beg(begSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(read_bed_file_part(filename, n_ind, beg, end));
    return rcpp_result_gen;
END_RCPP
}


