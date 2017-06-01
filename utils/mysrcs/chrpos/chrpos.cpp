#include "/nfs/users/nfs_f/fg6/ana/cpp/myinclude/readnwritefaq.h"
#include "/nfs/users/nfs_f/fg6/ana/cpp/myinclude/macro.h"

static  vector<int> seqpos;
static  vector<int> refpos;

#include <map>
template<typename Map>
void print_map(Map& m)
{
   std::cout << '{';
   for(auto& p: m)
        std::cout << p.first << ':' << p.second << ' ';
   std::cout << "}\n";
}
 
static  std::map<string, int> refmap;
static  std::map<string, int> seqmap;
static  std::ofstream myals;


int readals(char* file);

int main(int argc, char *argv[])
{ 
  
  if (argc == 1) {
    cout << " Change ref and draft -position of alignments according to position of contigs and chrs in fasta files " << endl;
    cout << " (to be compatible with ACT positioning scheme: " << endl 
	 << "     not restart at each chr/contig, but cumulative according to order of contigs/chrs in fasta files)" << endl;
   fprintf(stderr, "Usage: %s <reads.fq/fa>  <reads.fq/fa> <alignment> \n", argv[0]);
   
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
  if((fp = gzopen(argv[3],"r")) == NULL){ 
    printf("ERROR main:: missing input file 3 !! \n");
    return 1;
  }
  gzclose(argv[3]);

  // fasta file - contigs_sel - minl
  string reffile = argv[1];
  string seqfile = argv[2];
  string alfile = argv[3];


  string myname=myrename(alfile,"foract");
  myals.open(myname.c_str()); // written in readals

  
  // find chr/contigs lengths and creating maps
  
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

  refpos.resize(rlen.size());

  refpos[0]=0;
  if(pri)cout << "0 " << refpos[0] << " size is " << rlen[0] << endl;
  for(int i=0; i<rname.size(); i++){
    string name=rname[i];
    refmap[name] = i;   // map name to chr order
    
    if(i>0){
      refpos[i]=std::accumulate(rlen.begin(), rlen.begin()+i, 0); 
      if(pri)cout << i << " " << refpos[i] << " size is " << rlen[i] << endl;
    }
  }
  if(pri)print_map(refmap);
  rname.clear();
  rlen.clear();
  
  isfq=fasttype(argv[2]);
  err=0;
  saveinfo=1;
  readseq=0;
  if(!isfq){
    err=readfasta(argv[2],saveinfo);
  }else{
    err=readfastq(argv[2],saveinfo);
  }
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }

  seqpos.resize(rlen.size());
  seqpos[0]=0;
  if(pri)cout << "0 " << seqpos[0] << " size is " << rlen[0] << endl;
  for(int i=0; i<rname.size(); i++){
    string name=rname[i];
    seqmap[name] = i;
    if(i>0){
      seqpos[i]=std::accumulate(rlen.begin(), rlen.begin()+i, 0); 
      if(pri)cout << i << " " << seqpos[i] << " size is " << rlen[i] << endl;
    }
  }
  if(pri)print_map(seqmap);
  rname.clear();
  rlen.clear();

  readals(argv[3]);

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
	vector<string> more(6);
        ss >> ctg >> chr >> more[0] >> more[1]  >> more[2]  >> more[3] >> ctgi >> ctgf >> chri >> chrf >> more[4] >> more[5];

	int fctgi=ctgi+seqpos[seqmap[ctg]];
	int fctgf=ctgf+seqpos[seqmap[ctg]];
	int ri=refmap[chr];

	int fchri=chri+refpos[refmap[chr]];
	int fchrf=chrf+refpos[refmap[chr]];

	//cout <<  ctg << "\t" <<  chr << "\t" <<  more[0] << "\t" <<  more[1]  <<  "\t" <<  more[2]  << "\t" <<  more[3] << "\t" 
	//   <<  fctgi << "\t" <<  fctgf << "\t" <<  fchri << "\t" <<  fchrf << "\t" <<  more[4] << "\t" <<  more[5] << endl;
	myals <<  ctg << "\t" <<  chr << "\t" <<  more[0] << "\t" <<  more[1]  <<  "\t" <<  more[2]  << "\t" <<  more[3] << "\t" 
	     <<  fctgi << "\t" <<  fctgf << "\t" <<  fchri << "\t" <<  fchrf << "\t" <<  more[4] << "\t" <<  more[5] << endl;
  }
}
