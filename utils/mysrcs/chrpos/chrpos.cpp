#include "../myinc/macro.h"

static  vector<long int> seqpos;
static  vector<long int> refpos;

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
  gzclose(fp);
  if((fp = gzopen(argv[2],"r")) == NULL){ 
    printf("ERROR main:: missing input file 2 !! \n");
    return 1;
  }
  gzclose(fp);
  if((fp = gzopen(argv[3],"r")) == NULL){ 
    printf("ERROR main:: missing input file 3 !! \n");
    return 1;
  }
  gzclose(fp);

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
  int w10=0;
  for(long int i=0; i<rname.size(); i++){
    string name=rname[i];
    refmap[name] = i;   // map name to chr order
    
    if(i>0){
      //refpos[i]=std::accumulate(rlen.begin(), rlen.begin()+i, 0); 
      // accumulate does not handle long int??
      if(0)cout << endl;
      for (long int kk=0; kk<i; kk++){
	refpos[i]+=rlen[kk];
	if(0)cout << i << " " << kk << " " << refpos[i] << " " << rlen[kk] << endl;
      }
      
      if(((refpos[i]<0 && w10<10) || (w10<10)) && 0)
	cout << "    final " << i << " " << refpos[i] << " size is " << rlen[i] << endl;
      w10++;
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
  for(long int i=0; i<rname.size(); i++){
    string name=rname[i];
    seqmap[name] = i;
    if(i>0){
      
      //seqpos[i]=std::accumulate(rlen.begin(), rlen.begin()+i, 0); 
      // accumulate does not handle long int??
      for (long int kk=0; kk<i; kk++){
	seqpos[i]+=rlen[kk];
      }

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
  int pi=0;
  while(getline(infile,line)){
        std::stringstream ss(line);
        string ctg, chr;
	long int ctgi, ctgf, chri,chrf;
	vector<string> more(6);
        ss >> ctg >> chr >> more[0] >> more[1]  >> more[2]  >> more[3] >> ctgi >> ctgf >> chri >> chrf >> more[4] >> more[5];

	long int fctgi=ctgi+seqpos[seqmap[ctg]];
	long int fctgf=ctgf+seqpos[seqmap[ctg]];
	long int ri=refmap[chr];

	long int fchri=chri+refpos[refmap[chr]];
	long int fchrf=chrf+refpos[refmap[chr]];

	if(fctgi < 0 ){ //&& chr == "chr21"){

	  if(pi<0){
	    cout << chri << " " << refpos[refmap[chr]] << " " 
		 << refmap[chr] << " " << fchri
		 << endl;
	    if(0)cout <<  ctg << "\t" <<  chr << "\t" <<  more[0] << "\t" <<  more[1]  <<  "\t" <<  more[2]  
	       << "\t" <<  more[3] << "\t"  <<  fctgi << "\t" <<  fctgf << "\t" <<  fchri << "\t"
	       <<  fchrf << "\t" <<  more[4] << "\t" <<  more[5] << endl;
	  }
	  pi++;
	}

	if(fctgi < 0 || fctgf < 0 || fchri<0 || fchrf<0)
	  cout << " Error, possible overflow?? " <<  fctgi << "\t" <<  fctgf << "\t" <<  fchri << "\t" <<  fchrf << endl;

	myals <<  ctg << "\t" <<  chr << "\t" <<  more[0] << "\t" <<  more[1]  <<  "\t" <<  more[2]  << "\t" <<  more[3] << "\t" 
	      <<  fctgi << "\t" <<  fctgf << "\t" <<  fchri << "\t" <<  fchrf << "\t" <<  more[4] << "\t" <<  more[5] << endl;
  }
}
