#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
CharacterVector transl(CharacterVector Seq, int trans_tab) {
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
      if(c1=="G" && c2=="C"){
        aa[i].replace(j,1,"A");
      } else if(c1=="C" && c2=="G"){
        aa[i].replace(j,1,"R");
      } else if((c1=="A" || c1=="M") && (c2=="G") && (c3=="A" || c3=="G" || c3=="R")){
        aa[i].replace(j,1,"R");
      } else if(c1=="A" && c2=="A" && (c3=="T" || c3=="C" || c3=="Y")){
        aa[i].replace(j,1,"N");
      } else if(c1=="G" && c2=="A" && (c3=="T" || c3=="C" || c3=="Y")){
        aa[i].replace(j,1,"D");
      } else if(c1=="T" && c2=="G" && (c3=="T" || c3=="C" || c3=="Y")){
        aa[i].replace(j,1,"C");
      } else if(c1=="C" && c2=="A" && (c3=="A" || c3=="G" || c3=="R")){
        aa[i].replace(j,1,"Q");
      } else if(c1=="G" && c2=="A" && (c3=="A" || c3=="G" || c3=="R")){
        aa[i].replace(j,1,"E");
      } else if(c1=="G" && c2=="G"){
        aa[i].replace(j,1,"G");
      } else if(c1=="C" && c2=="A" && (c3=="T" || c3=="C" || c3=="Y")){
        aa[i].replace(j,1,"H");
      } else if(c1=="A" && c2=="T" && (c3=="T" || c3=="C" || c3=="A" || c3=="H")){
        aa[i].replace(j,1,"I");
      } else if(c1=="A" && c2=="T" && c3=="G"){
        aa[i].replace(j,1,"M");
      } else if(c1=="C" && c2=="T"){
        aa[i].replace(j,1,"L");
      } else if((c1=="T" || c1=="Y") && c2=="T" && (c3=="A" || c3=="G" || c3=="R")){
        aa[i].replace(j,1,"L");
      } else if(c1=="A" && c2=="A" && (c3=="A" || c3=="G" || c3=="R")){
        aa[i].replace(j,1,"K");
      } else if(c1=="T" && c2=="T" && (c3=="T" || c3=="C" || c3=="Y")){
        aa[i].replace(j,1,"F");
      } else if(c1=="C" && c2=="C"){
        aa[i].replace(j,1,"P");
      } else if(c1=="T" && c2=="C"){
        aa[i].replace(j,1,"S");
      } else if(c1=="A" && c2=="G" && (c3=="T" || c3=="C" || c3=="Y")){
        aa[i].replace(j,1,"S");
      } else if(c1=="A" && c2=="C"){
        aa[i].replace(j,1,"T");
      } else if(c1=="T" && c2=="G" && c3=="G"){
        aa[i].replace(j,1,"W");
      } else if(c1=="T" && c2=="A" && (c3=="T" || c3=="C" || c3=="Y")){
        aa[i].replace(j,1,"Y");
      } else if(c1=="G" && c2=="T"){
        aa[i].replace(j,1,"V");
      } else if(c1=="T" && c2=="A" && (c3=="A" || c3=="G" || c3=="R")){
        aa[i].replace(j,1,"*");
      } else if(c1=="T" && (c2=="A" || c2=="G" || c2=="R") && c3=="A"){
        if(trans_tab!=11 && c1=="T" && c2=="G" && c3=="A"){
          aa[i].replace(j,1,"W");
        } else {
          aa[i].replace(j,1,"*");
        }
      } else {
        aa[i].replace(j,1,"X");
      }
    }
  }
  return(wrap(aa));
}

// [[Rcpp::export]]
CharacterVector translCodon(CharacterVector Seq) {
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
      if(c1=="A" && c2=="A"){
        if(c3=="A"){
          aa[i].replace(j,1,"A");
        } else if(c3=="C"){
          aa[i].replace(j,1,"B");
        } else if(c3=="G"){
          aa[i].replace(j,1,"C");
        } else if(c3=="T"){
          aa[i].replace(j,1,"D");
        }
      } else if(c1=="A" && c2=="C"){
        if(c3=="A"){
          aa[i].replace(j,1,"E");
        } else if(c3=="C"){
          aa[i].replace(j,1,"F");
        } else if(c3=="G"){
          aa[i].replace(j,1,"G");
        } else if(c3=="T"){
          aa[i].replace(j,1,"H");
        }
      } else if(c1=="A" && c2=="G"){
        if(c3=="A"){
          aa[i].replace(j,1,"I");
        } else if(c3=="C"){
          aa[i].replace(j,1,"J");
        } else if(c3=="G"){
          aa[i].replace(j,1,"K");
        } else if(c3=="T"){
          aa[i].replace(j,1,"L");
        }
      } else if(c1=="A" && c2=="T"){
        if(c3=="A"){
          aa[i].replace(j,1,"M");
        } else if(c3=="C"){
          aa[i].replace(j,1,"N");
        } else if(c3=="G"){
          aa[i].replace(j,1,"O");
        } else if(c3=="T"){
          aa[i].replace(j,1,"P");
        }
      } else if(c1=="C" && c2=="A"){
        if(c3=="A"){
          aa[i].replace(j,1,"Q");
        } else if(c3=="C"){
          aa[i].replace(j,1,"R");
        } else if(c3=="G"){
          aa[i].replace(j,1,"S");
        } else if(c3=="T"){
          aa[i].replace(j,1,"T");
        }
      } else if(c1=="C" && c2=="C"){
        if(c3=="A"){
          aa[i].replace(j,1,"U");
        } else if(c3=="C"){
          aa[i].replace(j,1,"V");
        } else if(c3=="G"){
          aa[i].replace(j,1,"W");
        } else if(c3=="T"){
          aa[i].replace(j,1,"X");
        }
      } else if(c1=="C" && c2=="G"){
        if(c3=="A"){
          aa[i].replace(j,1,"Y");
        } else if(c3=="C"){
          aa[i].replace(j,1,"Z");
        } else if(c3=="G"){
          aa[i].replace(j,1,"a");
        } else if(c3=="T"){
          aa[i].replace(j,1,"b");
        }
      } else if(c1=="C" && c2=="T"){
        if(c3=="A"){
          aa[i].replace(j,1,"c");
        } else if(c3=="C"){
          aa[i].replace(j,1,"d");
        } else if(c3=="G"){
          aa[i].replace(j,1,"e");
        } else if(c3=="T"){
          aa[i].replace(j,1,"f");
        }
      } else if(c1=="G" && c2=="A"){
        if(c3=="A"){
          aa[i].replace(j,1,"g");
        } else if(c3=="C"){
          aa[i].replace(j,1,"h");
        } else if(c3=="G"){
          aa[i].replace(j,1,"i");
        } else if(c3=="T"){
          aa[i].replace(j,1,"j");
        }
      } else if(c1=="G" && c2=="C"){
        if(c3=="A"){
          aa[i].replace(j,1,"k");
        } else if(c3=="C"){
          aa[i].replace(j,1,"l");
        } else if(c3=="G"){
          aa[i].replace(j,1,"m");
        } else if(c3=="T"){
          aa[i].replace(j,1,"n");
        }
      } else if(c1=="G" && c2=="G"){
        if(c3=="A"){
          aa[i].replace(j,1,"o");
        } else if(c3=="C"){
          aa[i].replace(j,1,"p");
        } else if(c3=="G"){
          aa[i].replace(j,1,"q");
        } else if(c3=="T"){
          aa[i].replace(j,1,"r");
        }
      } else if(c1=="G" && c2=="T"){
        if(c3=="A"){
          aa[i].replace(j,1,"s");
        } else if(c3=="C"){
          aa[i].replace(j,1,"t");
        } else if(c3=="G"){
          aa[i].replace(j,1,"u");
        } else if(c3=="T"){
          aa[i].replace(j,1,"v");
        }
      } else if(c1=="T" && c2=="A"){
        if(c3=="A"){
          aa[i].replace(j,1,"w");
        } else if(c3=="C"){
          aa[i].replace(j,1,"x");
        } else if(c3=="G"){
          aa[i].replace(j,1,"y");
        } else if(c3=="T"){
          aa[i].replace(j,1,"z");
        }
      } else if(c1=="T" && c2=="C"){
        if(c3=="A"){
          aa[i].replace(j,1,"0");
        } else if(c3=="C"){
          aa[i].replace(j,1,"1");
        } else if(c3=="G"){
          aa[i].replace(j,1,"2");
        } else if(c3=="T"){
          aa[i].replace(j,1,"3");
        }
      } else if(c1=="T" && c2=="G"){
        if(c3=="A"){
          aa[i].replace(j,1,"4");
        } else if(c3=="C"){
          aa[i].replace(j,1,"5");
        } else if(c3=="G"){
          aa[i].replace(j,1,"6");
        } else if(c3=="T"){
          aa[i].replace(j,1,"7");
        }
      } else if(c1=="T" && c2=="T"){
        if(c3=="A"){
          aa[i].replace(j,1,"8");
        } else if(c3=="C"){
          aa[i].replace(j,1,"9");
        } else if(c3=="G"){
          aa[i].replace(j,1,"-");
        } else if(c3=="T"){
          aa[i].replace(j,1,"+");
        }
      }
    }
  }
  return(wrap(aa));
}
