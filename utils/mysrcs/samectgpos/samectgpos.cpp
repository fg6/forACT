
#include "../myinc/macro.h"

static  vector<int> seqpos;
static  vector<string> seqname;

static  std::ofstream myals;
static  string myname;
static  string myname1;


int readals(char* file);

int main(int argc, char *argv[])
{ 
  
  if (argc == 1) {
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
  size_t pos = 0;
  string token;
  string delimiter = "/";

  myname = alfile;
  while ((pos = myname.find(delimiter)) != std::string::npos) {
    token = myname.substr(0, pos);
    myname.erase(0, pos + delimiter.length());
  }
  myname= "ctgpos_"+ myname;
   
  //draft assembly
  int isfq=fasttype(argv[1]);
  int saveinfo=1;
  fp = gzopen(argv[1],"r");
  if(!isfq){
    readfasta(argv[1],saveinfo);
  }
  for(int i=0; i<rname.size(); i++){
    string name=rname[i];
    string names=rname[i];
    seqmap[name] = i;

    string sy ("_");   
    std::size_t  found = name.find_last_of(sy);
    if (found!=std::string::npos){
      string nn=name.erase(0,found+1);
      seqpos.push_back(to_int(nn.c_str())-1);
      seqname.push_back(names.erase(found));
    }else{
      cout << " warning! initial position not found! " << rname[i] << endl;
      name.clear();   
    }
  }
  rname.clear();
  rlen.clear();
  gzclose(argv[1]);
  
  myals.open(myname.c_str()); // written in readals
  readals(argv[2]);
  myals.close();
   
  return 0;
}


int readals(char* file){
  std::ifstream infile(file);
  string line;


  string prevctg="";
  string prevchr="";

  while(getline(infile,line)){
        std::stringstream ss(line);
        string ctg, chr;
	int ctgi, ctgf, chri,chrf;
	vector<string> more(7);
        ss >> ctg >> chr >> more[0] >> more[1]  >> more[2]  >> more[3] >> ctgi >> ctgf >> chri >> chrf >> more[4] >> more[5]>> more[6];

	int thisctg=seqmap[ctg];

	int fctgi=ctgi+seqpos[thisctg];
	int fctgf=ctgf+seqpos[thisctg];
	
	//cout <<  seqname[thisctg] << " " <<  ctgi << "\t" <<  ctgf << " " 
	//   <<  fctgi << "\t" <<  fctgf << endl;

	myals <<  seqname[thisctg] << "\t" <<  chr << "\t" <<  more[0] << "\t" <<  more[1]  <<  "\t" <<  more[2]  << "\t" <<  more[3] << "\t" 
	   <<  fctgi << "\t" <<  fctgf << "\t" <<  chri << "\t" <<  chrf << "\t" <<  more[4] << "\t" <<  more[5] << "\t" << more[6] << endl;
  }
}
