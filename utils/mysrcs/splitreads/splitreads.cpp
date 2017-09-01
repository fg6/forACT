
#include "../myinc/macro.h"

static  string myname;

int main(int argc, char *argv[])
{ 

  if (argc < 3) {
   fprintf(stderr, "Shred each read/contig in read/contig of length <length> \n Usage: %s <reads.fasta> <length> \n", argv[0]);
   
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file  !! \n");
    return 1;
  }
  int len=0;
  len=atoi(argv[2]);  
 
  int write=1; // fasta file
  string seqfile = argv[1];

  // output file name
  size_t pos = 0;
  string token;
  string delimiter = "/";
  myname = seqfile;
  while ((pos = myname.find(delimiter)) != std::string::npos) {
    token = myname.substr(0, pos);
    myname.erase(0, pos + delimiter.length());
  }
  myname = "shred" + to_string(len) + "_" + myname; 
  if(write)myfile.open(myname.c_str());
 
  //read&write
  int isfq=fasttype(argv[1]);
  int err=0;
  int saveinfo=1;
  int readseq=1;
  if(!isfq){
    err=readfasta(argv[1],saveinfo,"",readseq);
  }else{
    err=readfastq(argv[1],saveinfo,"",readseq);
  }
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }

  char out[5]={">"};

  for (unsigned i=0; i < rseq.size(); i++) { // loop over scaffold/contigs 
    string s = rseq[i];
    string name = rname[i];
  
    for (unsigned i = 0; i < s.length(); i += len) {
      string thisname = name + "_" + to_string(i+1); //+ " split_in_" + to_string(len);
      myfile <<  out[0] << thisname <<  endl;
      myfile << s.substr(i, len) << endl;
    }
  }

  if(write) myfile.close();
  
  return 0;
}

