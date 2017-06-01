
#include "/nfs/users/nfs_f/fg6/ana/cpp/myinclude/readfaq.h"
#include "/nfs/users/nfs_f/fg6/ana/cpp/myinclude/macro.h"

//static int pri=0;

#include <map>
template<typename Map>
void print_map(Map& m)
{
   std::cout << '{';
   for(auto& p: m)
        std::cout << p.first << ':' << p.second << ' ';
   std::cout << "}\n";
}
static  std::map<string, int> seqmap;
static vector<int> chrs;

int readals(char* file);

int main(int argc, char *argv[])
{ 
  int select=0;
  
  if (argc == 1) {
   fprintf(stderr, "Usage: %s <reads.fq/fa>  \n", argv[0]);
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


  // fasta file 
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
 
  fp = gzopen(argv[1],"r");
  readseqs(1);
  for(int i=0; i<rname.size(); i++){
    string name=rname[i];
    seqmap[name] = i;
  }
  gzclose(argv[1]);

  if(select){
    chrs.resize(rname.size(),0);
    readals(argv[2]);
  }else{
    for(int i=0; i<rname.size(); i++){
      chrs.push_back(1);
    }
  }

  char out[5]={">"};
  string thisname= "selctg_"+ myname;
  myfile.open(thisname.c_str());

  for(int i=0; i<rname.size(); i++){
    string seq = rseq[i];
    string name = rname[i];
     
    
    if(chrs[i]){
      myfile <<  out[0] << name << endl << seq <<endl;
    }
  }
  myfile.close();
      
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
        string ctg;
        
        ss >> ctg; 
	
	if(seqmap.count(ctg)){
	  chrs[seqmap[ctg]]=1;
	  if(pri)cout << ctg << " " << seqmap[ctg]<< endl;
	  }else
	  cout << " ERROR!! " << ctg << " is not in the fasta!!" << endl;

        
        ii++;
  }
}
