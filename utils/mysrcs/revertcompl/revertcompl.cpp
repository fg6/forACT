#include "../myinc/macro.h"

vector<int> alforward;
vector<int> alall;
vector<int> forward;

int readals(char* file);

int main(int argc, char *argv[])
{ 
  
  if (argc == 1) {
   cout << " Revert complement contigs/scaffolds which have opposite direction to reference according to alignment"<<endl;
   fprintf(stderr, "Usage: %s <reads.fq/fa>  <alignment> \n", argv[0]);
   
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file 1 !! \n");
    return 1;
  }
  gzclose(argv[1]);
  if((fp = gzopen(argv[2],"r")) == NULL){ 
  printf("ERROR main:: missing input file 2 !! \n");
       return 1;
  }
   gzclose(argv[2]);


  // fasta file - contigs_sel - minl
  string seqfile = argv[1];
  string alfile = argv[2];

  // output file name
  string myname=myrename(seqfile,"forw");
    
  //draft assembly
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


  for(int i=0; i<rname.size(); i++){
    string name=rname[i];
    seqmap[name] = i;
  }
  

  readals(argv[2]);
  char out[5]={">"};
  myfile.open(myname.c_str());

  int forw=0;
  int comp=0;
  int notfound=0;
  for(int i=0; i<rname.size(); i++){
    string seq = rseq[i];
    string name = rname[i];
    string forseq;

    if(0)cout << name << " " << alall[i] << " " << alforward[i] << endl;

    if(alall[i]>0){
      if( alall[i]-alforward[i] >  alall[i]*0.5){
	if(0)cout << "  complement! " << endl;
	forseq=comple(seq);
	comp++;
      }else{
	if(0)cout << "  forward! " << endl;
	forseq=seq;
	forw++;
      }

    }else{
      forseq=seq;
      notfound++;
    }

    myfile <<  out[0] << name << endl << forseq <<endl;
    forseq.clear(); 
  }
  myfile.close();
  rname.clear();
  rlen.clear();
  rseq.clear();

  cout << " I found " << forw << " forward scaffolds and " << comp << " complement scaffolds " 
       << notfound << " scaffolds were not found in the alignments "
       << endl; 

  return 0;
}



int readals(char* file){
  std::ifstream infile(file);
  string line;

  alforward.resize(rname.size(),0);
  alall.resize(rname.size(),0);
  int ii=-1;
  while(getline(infile,line)){
        std::stringstream ss(line);
	ii++;
	
        string ctg, chr;
	int ctgi, ctgf, chri,chrf;
	vector<string> more(7);

	
        ss >> ctg >> chr >> more[0] >> more[1]  >> more[2]  >> more[3] >> ctgi >> ctgf >> chri >> chrf >> more[4] >> more[5] >> more[6];
	
	string newctg=ctg;
	string sy ("_");   
	std::size_t  found = newctg.find_last_of(sy);
		

	if (found!=std::string::npos){
	  string nn=newctg.erase(found);
	}else{
	  cout << " warning! initial position not found! " << ctg << " line " << ii << endl;
	  newctg.clear();   
	}


	int thisctg=seqmap[newctg];
	if (more[6]=="*")continue;
	//cout << more[6] << " " << newctg << " " << thisctg << endl;

	if(more[6] =="F") alforward[thisctg]++;
	alall[thisctg]++;

  }
}
