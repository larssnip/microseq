#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


// [[Rcpp::export]]
SEXP extractSeq(SEXP Gseq, SEXP Left, SEXP Right, SEXP Strand) {
  std::string gseq                = as<std::string>(Gseq);
  std::vector< int > left         = as<std::vector< int > >(Left);
  std::vector< int > right        = as<std::vector< int > >(Right);
  std::vector< int > strand       = as<std::vector< int > >(Strand);
  
  std::vector< std::string > sequences(left.size());
  
  for(unsigned i=0; i<sequences.size() ; i++){
    if(left[i] < right[i]){
      sequences[i] = gseq.substr(left[i]-1, right[i]-left[i]+1);
    } else {
      sequences[i] = gseq.substr(left[i]-1, gseq.length());
      sequences[i].append(gseq.substr(0, right[i]));
    }
  }
  return wrap(sequences);
}
