
#include "../myinc/macro.h"

static vector<int> chrs;

int readals(char* file);

int main(int argc, char *argv[])
{ 
  int select=0;
  
  if (argc < 2) {
    cout << " Write each contig/scaffold in a single fasta/fastq file (same as original if not specified)" << endl;
    cout << " If there's an alignment file in input, only contig/scaffolds in the alignment file are written out " << endl;
    fprintf(stderr, "\n Usage: %s <reads.fq/fa>  [al.txt] [fasta/fastq] \n", argv[0]);
    return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file 1 !! \n");
    return 1;
  }
  gzclose(argv[1]);
  if(argc > 2){ 
    select=1;
    if((fp = gzopen(argv[2],"r")) == NULL){ 
      printf("ERROR main:: missing input file 2 !! \n");
      return 1;
    }
    gzclose(argv[2]); 
  }

  string otype="same";
  if(argc > 3) otype=argv[3];     // write fasta from fastq

  // fasta file 
  string seqfile = argv[1];
  size_t ns=0;
  string myname=seqfile;
  ns=myname.find_last_of("/");
  if(ns!=std::string::npos) { 
    myname=myname.substr(ns+1,myname.size());
  }
  

  int isfq=fasttype(argv[1]);
  if(otype=="fastq" && isfq) otype="same";

  int err=0;
  int saveinfo=1;
  int readseq=1;
  if(!isfq){
    err=readfasta(argv[1],saveinfo,"",readseq); //"same" for writing fasta
  }else{
    err=readfastq(argv[1],saveinfo,otype,readseq); //"same" for fastq, "fasta" for fasta
  }
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }

  for(int i=0; i<rname.size(); i++){
    string name=rname[i];
    seqmap[name] = i;
  }


  // al file
  if(select){
    chrs.resize(rname.size(),0);
    readals(argv[2]);
  }else{
    for(int i=0; i<rname.size(); i++){
      chrs.push_back(1);
    }
  }

  char out[5]={">"};
  for(int i=0; i<rname.size(); i++){
    string seq = rseq[i];
    string name = rname[i];
     
    string thisname= name + "_"+ myname;
    
    if(chrs[i]){
      myfile.open(thisname.c_str());
      myfile <<  out[0] << name << endl << seq <<endl;
      myfile.close();
    }
  }
  rname.clear();
  rlen.clear();
  rseq.clear();


  return 0;
}


int readals(char* file){
  std::ifstream infile(file);
  string line;

  int ii=0;
  while(getline(infile,line)){
        std::stringstream ss(line);
        string ctg, chr;
        int ctgi, ctgf, chri,chrf;
        vector<string> more(7);
        
        ss >> ctg >> chr >> more[0] >> more[1]  >> more[2]  >> more[3] >> ctgi >> ctgf 
	   >> chri >> chrf >> more[4] >> more[5] >> more[6];
	
	if(seqmap.count(ctg)){
	  chrs[seqmap[ctg]]=1;
	  if(pri)cout << ctg << " " << seqmap[ctg]<< endl;
	}else
	  cout << " ERROR!! " << ctg << " is not in the fasta!!" << endl;        
	ii++;
  }
}
