#include <Rcpp.h>
using namespace Rcpp;

int switchNUC(int s){
  s = toupper(s);
  switch(s){
  case 'A' :
    return('T');
  case 'T' :
    return('A');
  case 'U' :
    return('S');
  case 'C' :
    return('G');
  case 'G' :
    return('C');
    
  case 'R' :
    return('Y');
  case 'Y' :
    return('R');
  case 'K' :
    return('M');
  case 'M' :
    return('K');
    
  case 'D' :
    return('H');
  case 'H' :
    return('D');
  case 'B' :
    return('V');
  case 'V' :
    return('B');
  default : 
    return(s);
  }
}

// [[Rcpp::export]]
CharacterVector revComp(CharacterVector Seq, bool rev) {
  std::vector<std::string> seq = as<std::vector<std::string> >(Seq);
  for(unsigned int i=0; i<seq.size(); ++i){
    std::transform(seq[i].begin(), seq[i].end(),seq[i].begin(), switchNUC);
    if(rev){
      std::reverse(seq[i].begin(), seq[i].end());
    }
  }
  return(wrap(seq));
}

