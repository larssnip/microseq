#include <fstream>
#include <sstream>
#include <string>
#include <Rcpp.h>
using namespace Rcpp;

std::istream& safeGetline(std::istream& is, std::string& t)
{
  t.clear();
  
  // Source: https://stackoverflow.com
  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.
  
  std::istream::sentry se(is, true);
  std::streambuf* sb = is.rdbuf();
  
  for(;;) {
    int c = sb->sbumpc();
    switch (c) {
    case '\n':
      return is;
    case '\r':
      if(sb->sgetc() == '\n')
        sb->sbumpc();
      return is;
    case EOF:
      // Also handle the case when the last line has no line ending
      if(t.empty())
        is.setstate(std::ios::eofbit);
      return is;
    default:
      t += (char)c;
    }
  }
}

// [[Rcpp::export]]
Rcpp::List read_fastq(std::string path) { // for standard fastq, not Sanger
  std::vector<std::string> header, sequence, quality;
  std::ifstream t(path.c_str());
  std::stringstream ss;
  ss << t.rdbuf();
  std::string to;

  while(safeGetline(ss,to) && !ss.eof()){
    header.push_back(to);
    safeGetline(ss,to);
    sequence.push_back(to);
    safeGetline(ss,to); // +
    safeGetline(ss,to); // Quality
    quality.push_back(to);
  }
  // return Rcpp::DataFrame::create( Rcpp::Named("Header")= header, Rcpp::Named("Sequence") = sequence);
  return Rcpp::List::create(Rcpp::Named("Header") = header, Rcpp::Named("Sequence") = sequence, Rcpp::Named("Quality") = quality);
}

// [[Rcpp::export]]
Rcpp::List read_fastq_Sanger(std::string path) { // for fastq with Sanger compatibility
  std::vector<std::string> header, sequence, quality;
  std::ifstream t(path.c_str());
  std::stringstream ss;
  ss << t.rdbuf();
  std::string to;
  unsigned int n_lin = 1, n_seq = 0;
  
  // First round
  safeGetline(ss,to);
  header.push_back(to);
  safeGetline(ss,to);
  sequence.push_back(to);
  while(safeGetline(ss,to)){
    if(to == "+" || (to.substr(0,1) == "+" && to.substr(1, std::string::npos) == header[n_seq].substr(1,std::string::npos))){
      break;
    } else { // Count number of lines used for current sequence
      sequence[n_seq].append(to);
      n_lin++;
    }
  }

  while(safeGetline(ss,to)){
    // Read number of lines used for current sequence (quality)
    quality.push_back(to);
    for(unsigned int i = 1; i<n_lin; ++i){
      safeGetline(ss,to); 
      quality[n_seq].append(to);
    }
    n_lin = 1;

    // Store header and sequence
    safeGetline(ss,to);
    if(to == "" || ss.eof()){
      break;
    } else {
      n_seq++;
    }
    header.push_back(to);
    safeGetline(ss,to);
    sequence.push_back(to);
    while(safeGetline(ss,to)){
      if(to == "+" || (to.substr(0,1) == "+" && to.substr(1, std::string::npos) == header[n_seq].substr(1,std::string::npos))){
        break;
      } else { // Count number of lines used for current sequence
        sequence[n_seq].append(to);
        n_lin++;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("Header") = header, Rcpp::Named("Sequence") = sequence, Rcpp::Named("Quality") = quality);
}

// [[Rcpp::export]]
bool write_fastq(std::vector<std::string> header, std::vector<std::string> sequence, std::vector<std::string> quality, std::string path) { // for standard FASTA format
  std::ofstream t(path.c_str());
  for(unsigned int j=0; j<header.size();++j){
    t << header[j] << std::endl;
    t << sequence[j] << std::endl;
    t << '+' << std::endl;
    t << quality[j] << std::endl;
  }
  t.close();
  return(true);
}


// [[Rcpp::export]]
Rcpp::List read_fasta(std::string path) { // for standard FASTA format
  std::vector<std::string> header, sequence;
  std::ifstream t(path.c_str());
  std::stringstream ss;
  ss << t.rdbuf();
  std::string to;
  unsigned int n_seq = 0;
  
  safeGetline(ss,to);
  header.push_back(to.substr(1));
  safeGetline(ss,to);
  sequence.push_back(to);
  while(safeGetline(ss,to)){
    if(to.substr(0,1) == ">"){
      ++n_seq;
      header.push_back(to.substr(1));
      safeGetline(ss,to);
      sequence.push_back(to);
    } else {
      sequence[n_seq].append(to);
    }
  }
  return Rcpp::List::create(Rcpp::Named("Header") = header, Rcpp::Named("Sequence") = sequence);
}

// [[Rcpp::export]]
bool write_fasta(std::vector<std::string> header, std::vector<std::string> sequence, std::string path, int width) { // for standard FASTA format
  std::ofstream t(path.c_str());
  for(unsigned int j=0; j<header.size();++j){
    t << ">" << header[j] << std::endl;
    for (unsigned i = 0; i < sequence[j].length(); i += width) {
      t << sequence[j].substr(i, width) << std::endl;
    }
  }
  t.close();
  return(true);
}
