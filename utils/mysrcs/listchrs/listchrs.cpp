#include "../myinc/macro.h"
#include <locale>
#include <sstream>

int main(int argc, char *argv[])
{ 
  
  if (argc == 1) {
    cout << " List all contig/scaffold/chr and their sizes " << endl;
    fprintf(stderr, "\nUsage: %s <reads.fq/fa>  \n", argv[0]);   
    return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file 1 !! \n");
    return 1;
  }

  int isfq=fasttype(argv[1]);
  int err=0;
  int saveinfo=1;
  int readseq=0;
  if(!isfq){
    err=readfasta(argv[1],saveinfo);
  }else{
    err=readfastq(argv[1],saveinfo);
  }
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }

  cout << endl << " Please pick a Chromosome from: " << endl;
  for(int i=0; i<rname.size(); i++){
    if(i < 30){
      std::stringstream ss;
      ss.imbue(std::locale(""));
      ss << "  " << rname[i] << " lenght=" <<  std::fixed << rlen[i] << " bp";
      cout << ss.str() << endl;
    }else if(i==31){
      cout << endl << " Too many Contig/Chr, only printing 30." <<endl;
    }
  }
  cout << endl;
  
  return 0;
}

