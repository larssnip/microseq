#include <Rcpp.h>
using namespace Rcpp;

// CharacterVector revComp1(CharacterVector Seq) {
//   std::vector<std::string> seq = as<std::vector<std::string> >(Seq);
//   unsigned int N = seq.size();
//   std::vector<std::string> aa = std::vector<std::string>(N);
//   std::string codon;
//   char c1, c2, c3;
//   
//   static const std::string arr[MAX_KEYS]  = { "ABC", "DEF", "GHI", "JKL" , ...};
//   std::set<std::string> myset(&arr[0], &arr[MAX_KEYS]);
//   
//   for(unsigned int i=0; i<N; ++i){
//     unsigned int n = seq[i].length()/3;
//     aa[i].reserve(n);
//     for(unsigned int j=0; j<n; ++j){
//       c1 = (char)seq[i].substr(3*j);
//       c2 = (char)seq[i].substr(3*j+1);
//       c3 = (char)seq[i].substr(3*j+2);
//       
//       codon = seq[i].substr(3*j, 3*j+2);
//       switch(codon.substr(0,0)  ){
//       case 'A' :
//         switch(codon.substr(1,1)){
//         case 'A' :
//           switch(codon.substr(2,2)){
//           case 'A' :
//             return('K');
//           default :
//             return('-');
//           }
//         default :
//           return('-');
//         }
//       default :
//         return('-');
//       }
//     }
//   }
//   
//   return x * 2;
// }


// [[Rcpp::export]]
CharacterVector transl(CharacterVector Seq) {
  std::vector<std::string> seq = as<std::vector<std::string> >(Seq);
  unsigned int N = seq.size();
  std::vector<std::string> aa = std::vector<std::string>(N);
  std::string codon;
  std::string c1, c2, c3;
  for(unsigned int i=0; i<N; ++i){
    unsigned int n = seq[i].length()/3;
    aa[i].reserve(n);
    for(unsigned int j=0; j<n; ++j){
      c1 = seq[i].substr(3*j, 1);
      c2 = seq[i].substr(3*j+1, 1);
      c3 = seq[i].substr(3*j+2, 1);
      if(c1 == "A"){
        if(c2 == "A"){
          if(c3 == "T" || c3 == "C" || c3 == "Y"){
            aa[i].replace(j,1,"N");
          } else {
            if(c3 == "A" || c3 == "G" || c3 == "R"){
              aa[i].replace(j,1,"K");
            } else {
              aa[i].replace(j,1,"-");
            }
          }
        } else {
          if(c2 == "C"){
            aa[i].replace(j,1,"T");
          } else {
            if(c2 == "G"){
              if(c3 == "T" || c3 == "C" || c3 == "Y"){
                aa[i].replace(j,1,"S");
              } else {
                if(c3 == "A" || c3 == "G" || c3 == "R"){
                  aa[i].replace(j,1,"R");
                } else {
                  aa[i].replace(j,1,"-");
                }
              }
            } else { // c2 == T
              if(c3 == "T" || c3 == "C" || c3 == "A" || c3 == "Y" || c3 == "M" || c3 == "W" || c3 == "H"){
                aa[i].replace(j,1,"I");
              } else {
                if(c3 == "G"){
                  aa[i].replace(j,1,"M");
                } else {
                  aa[i].replace(j,1,"-");
                }
              }
            }
          }
        }
      } else {
        if(c1 == "C"){
          if(c2 == "A"){
            if(c3 == "T" || c3 == "C" || c3 == "Y"){
              aa[i].replace(j,1,"H");
            } else {
              if(c3 == "A" || c3 == "G" || c3 == "R"){
                aa[i].replace(j,1,"Q");
              } else {
                aa[i].replace(j,1,"-");
              }
            }
          } else {
            if(c2 == "C"){
              aa[i].replace(j,1,"P");
            } else {
              if(c2 == "G"){
                aa[i].replace(j,1,"R");
              } else { // T
                aa[i].replace(j,1,"L");
              }
            }
          }
        } else {
          if(c1 == "G"){
            if(c2 == "A"){
              if(c3 == "T" || c3 == "C" || c3 == "Y"){
                aa[i].replace(j,1,"D");
              } else {
                if(c3 == "A" || c3 == "G" || c3 == "R"){
                  aa[i].replace(j,1,"E");
                } else {
                  aa[i].replace(j,1,"-");
                }
              }
            } else {
              if(c2 == "C"){
                aa[i].replace(j,1,"A");
              } else {
                if(c2 == "G"){
                  aa[i].replace(j,1,"G");
                } else { // T
                  aa[i].replace(j,1,"V");
                }
              }
            }
          } else { // c1 == T
            if(c2 == "A"){
              if(c3 == "T" || c3 == "C" || c3 == "Y"){
                aa[i].replace(j,1,"Y");
              } else {
                if(c3 == "A" || c3 == "G" || c3 == "R"){
                  aa[i].replace(j,1,"*");
                } else {
                  aa[i].replace(j,1,"-");
                }
              }
            } else {
              if(c2 == "C"){
                aa[i].replace(j,1,"S");
              } else {
                if(c2 == "G"){
                  if(c3 == "T" || c3 == "C" || c3 == "Y"){
                    aa[i].replace(j,1,"C");
                  } else {
                    if(c3 == "A"){
                      aa[i].replace(j,1,"*");
                    } else {
                      if(c3 == "G"){
                        aa[i].replace(j,1,"W");
                      } else {
                        aa[i].replace(j,1,"-");
                      }
                    }
                  }
                } else { // c2 == T
                  if(c3 == "T" || c3 == "C" || c3 == "Y"){
                    aa[i].replace(j,1,"F");
                  } else {
                    if(c3 == "A" || c3 == "G" || c3 == "R"){
                      aa[i].replace(j,1,"L");
                    } else {
                      aa[i].replace(j,1,"-");
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(wrap(aa));
}
