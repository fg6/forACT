#include "../myinc/macro.h"  //misf.h"
#include <ctime>
#include <locale> // for imbue

bool mysort (int i,int j) { return (i>j); }

struct ALIGNMENTS {
  std::string refname,draftname,strand;
  long int refstart,refstop,draftstart,draftstop;
  long int reflength,draftlength;
  int numbases, score;
  float identity;

  bool operator() (ALIGNMENTS i, ALIGNMENTS j) { 
    return( i.numbases > j.numbases ); 
  }
} alignments;
static  std::vector<ALIGNMENTS> myaligns;

int Read(FILE *namef);
int CheckIfAllChrs(void);
int CheckChr(std::string chrname, int nchr);
int CheckCtg(std::string ctgname);
int readals(char* file);

static char nameout[100];
static  FILE *outfile;
static  int nseq=0;
static int printsome=10;
static int onecontig=0;
static int misassembled=0;
static std::vector<std::string> misass;
static std::vector<std::string> chrs;
static std::vector<std::string> ctgs;
static std::vector<std::string> tocheck;
static std::vector<std::string> singlectg90;
static std::vector<std::string> singlectg80;
static std::vector<std::string> singlectg70;
static std::ofstream ofile;
static std::ofstream o2file;

// minimum alignment to declare single contig per chromosome
static int covtosingle=90; // percentage of chromosome covered 

// minimum alignment to multimple chromosome to declare possible misassembly
static int mismin=3; //at least 2% of the contig is aligned to a chromosome other than the main one

// maximum distance on reference between 2 pieces of a contig to declare countiguity 
static int maxdist=2000; // max number of bases between the end position on ref of one piece of 
                         // contig and the start position of next piece
// how much of each chr is covered
static std::vector<float> chrcov;
static std::vector<std::vector<float>> chrcovperc;
static std::vector<float> chrcovn;
static std::vector<std::string> chrcovname;
static  std::vector<std::string> allchrs;


// new
static std::map<string, std::pair<int,int> > mchrmapped;    // one el per chr: 0 or 1 if it is not (0) or is (1) mapped READY
static std::map<string, std::pair<int,int> > mctgmapped;    // one el per ctg: 0 or 1 if it is not (0) or is (1) mapped READY

static std::map<string, vector<std::tuple<string,int,int,float>> > mctgs; // one ele per ctg: a vector of chrs mapped to the ctg READY
static std::map<string, vector<std::tuple<string,int,int,float>>> mchrs;  // one ele per chr: a vector of ctgs mapped to the ctg


static std::vector<std::pair<string,int> > fastachrs;  //vector of chr names, length READY
static std::vector<std::pair<string,int> > fastactgs;  // vector of ctg names, length  READY
static  vector<string> unmappedchrs;  //vector of unmapped chrs READY
static std::ofstream misfile;
static int nmisctg=0;
static int nmis=0;
static long int  refsize=0;
static  vector<string> misassembledctgs;  //vector of unmapped chrs READY


struct nALIGNMENTS {
  std::string chr,ctg,strand;
  long int chri,chrf,ctgf,ctgi;
  int albases;
  float id;

  bool operator() (nALIGNMENTS i, nALIGNMENTS j) { 
    return( i.albases > j.albases ); 
  }
} nalignments;

static int CheckOneCtg(std::vector<nALIGNMENTS> ctgals);




int main(int argc, char **argv)
{
  
  if(argc < 2){
    printf("Usage: %s  <ref> <draft> <alignment_file>\n", argv[1]);
    exit(1);
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

  if(0)cout << "\nParameters: " 
      << "\n Minimum alignment to declare single contig per chromosome: " << covtosingle
      << "% of chromosome\n Minimum alignment to multimple chromosome to declare possible misassembly: " << mismin
      << "% of contig \n Maximum distance on reference between 2 pieces of a contig to declare countiguity: " 
      << maxdist << " bases " 
      << std::endl;
   

  //read reference chrs
  int err=0; int saveinfo=1; 
  err=readfasta(argv[1],saveinfo);
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }
  
  for(int c=0; c<rname.size(); c++){
    mchrmapped[rname[c]] = std::make_pair(rlen[c],0);
    fastachrs.push_back(std::make_pair(rname[c],rlen[c]));
    refsize+=rlen[c];
  }
  rname.clear();
  rlen.clear();
  //read draft ctgs
  err=0;  
  err=readfasta(argv[2],saveinfo);
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }
  for(int c=0; c<rname.size(); c++){
    mctgmapped[rname[c]]=std::make_pair(rlen[c],0);
    fastactgs.push_back(std::make_pair(rname[c],rlen[c]));
  }
  rname.clear();
  rlen.clear();


  //cout << " Found " << fastachrs.size() << " chromosomes and " << fastactgs.size() << " contigs " << endl;

  myaligns.reserve(10000);
  readals(argv[3]);

  cout << " Possibly " << nmisctg
       << " contig/s is/are misassembled " 
       << " with a total of possibly " << nmis << " break points "
       << endl;

  /* for(int c=0; c<misassembledctgs.size();c++){
    string thisctg=misassembledctgs[c];
    cout <<  " Contig " << thisctg << " is aligned to:" << endl;
    for ( const auto &p :  mctgs[misassembledctgs[c]] ){
      cout  << "    chromosome " << std::get<0>(p) 
	    << " :  Total length aligned=" << std::get<1>(p) 
	    << "  (of which " << std::get<2>(p) << " unmapped)"
	    << std::fixed << std::setprecision(2)
	    << " covers about " <<  std::get<1>(p)*100./std::get<0>(mchrmapped[ std::get<0>(p)]) << "% of the chr" 
	    << std::fixed << std::setprecision(1)
	    << " with an average identity " <<  std::get<3>(p) << "%"
	    << endl;
    }
    }*/

  return 0;
}


int readals(char* file){
  std::ifstream infile(file);
  string line;


  string prevctg="";
  string prevchr="";
  int nn=-1;
  int newctg=0;
  int nmappedctgs=1; // not counting nn=0 so adding it here
  long int refcov=0;

  std::vector<nALIGNMENTS> ctgals;

  myfile.open("ctg_report.txt");
  misfile.open("misassembly_report.txt");
  while(getline(infile,line)){
    std::stringstream ss(line);
    nn++; // count alignments read
    string ctg, chr;
    int ctgi, ctgf, chri,chrf,albases, kk;
    float id;
    
    ss >> ctg >> chr >> id >> albases  >> kk  >> kk >> ctgi >> ctgf >> chri >> chrf >> kk >> kk;
    if(0) cout <<  ctg << "\t" <<  chr << "\t" << id << "\t" <<  albases  <<  "\t" 
	       <<  ctgi << "\t" <<  ctgf << "\t" <<  chri << "\t" <<  chrf << endl;
    

    refcov+=std::abs(chrf-chri);


    //this chr is mapped, add in mchrmapped
    if(!mchrmapped.count(chr))cout << " Error chr missing in map! " << chr << endl;
    std::get<1>(mchrmapped[chr])=1;
    
    //this ctg is mapped, add in mctgmapped  // can reduce time significantly if only checked once
    if(!mctgmapped.count(ctg))cout << " Error ctg missing in map! " << ctg << endl;
    std::get<1>(mctgmapped[ctg])=1;

    nALIGNMENTS thisal;     
    thisal.chr = chr;
    thisal.ctg = ctg;
    thisal.chri =  chri;
    thisal.chrf = chrf;
    thisal.ctgi = ctgi;
    thisal.ctgf =  ctgf; 
    if(ctgf-ctgi < 0)
      thisal.strand = "C" ;
    else
      thisal.strand = "F" ;
    thisal.albases = albases;
    thisal.id = id;
     
    // check if new ctg:
    if(nn>0 && prevctg != ctg) newctg=1;
    else newctg=0;
    
    
    if(newctg || infile.eof()){
      nmappedctgs++;
      CheckOneCtg(ctgals);

      // reset all vectors
      ctgals.clear(); ctgals.resize(0);

      // fill first al
      ctgals.push_back(thisal);
    }else{
      // keep filling vectors
      ctgals.push_back(thisal);
    }
   
    thisal = nALIGNMENTS(); // re-initialize thisal
    prevctg=ctg;  

  }// end loop 
  //atch last contig
  nmappedctgs++;
  CheckOneCtg(ctgals);
 
  myfile.close();   
  misfile.close();   

  int nmappedchrs=0;
  for ( const auto &p : mchrmapped ){ 
    if(std::get<1>(p.second)==1) nmappedchrs++;
    else unmappedchrs.push_back(p.first);
  }

  cout << " Found " << nn << " alignments, " 
       << nmappedctgs << " contigs, and " << nmappedchrs << " chromosomes " 
       << endl;
  cout << std::fixed << std::setprecision(1) 
       << " Reference coverage: " <<  refcov*100./refsize << "%"
       << " (" << refcov << " bases)"
       << endl;

  //for(int cc=0; cc<unmappedchrs.size(); cc++)
  //cout << unmappedchrs[cc]<<endl;

  return 0;
}

// ******************************************** //
int CheckOneCtg(std::vector<nALIGNMENTS> ctgals)
// ******************************************** //
{
  int sum=0;
  string thisctg=ctgals[0].ctg;
  std::map<string, std::tuple<int,int,float,int> > mappedtochr; 
  vector<std::tuple<string,int,int,float> > thisctglinks;


  // add per each covered chr: min_pos, max_pos, are there holes in coverage? 
  //                           min_ctg_pos, max_ctg_pos, are there unmapped region of ctg?

  std::for_each(ctgals.begin(),ctgals.end(), [&] (nALIGNMENTS const& a) {
      int naligned=std::abs(a.ctgf-a.ctgi);
      int numapped=naligned-a.albases;

      sum+=std::abs(a.ctgf-a.ctgi); // total number of alignments
          
      if(!mappedtochr.count(a.chr)){
	int ii=1;
	mappedtochr[a.chr]=std::make_tuple(naligned,numapped,a.id,1);
      }else{
	int sumtoals=std::get<0>(mappedtochr[a.chr])+naligned;
	int sumtoun=std::get<1>(mappedtochr[a.chr])+numapped;
	int sumid=std::get<2>(mappedtochr[a.chr])+a.id;
	int nid=std::get<3>(mappedtochr[a.chr])+1;
	mappedtochr[a.chr]=std::make_tuple(sumtoals,sumtoun,sumid,nid);
      }
    });

  int tsum=0;
  int usum=0;
  for ( const auto &p : mappedtochr ){ 
    tsum+=std::get<0>(p.second);  //length aligned
    usum+=std::get<1>(p.second);  // unmapped
    string chrname=p.first;
    int g0=std::get<0>(p.second);//length aligned
    int g1=std::get<1>(p.second); // unmapped

    float g2=std::get<2>(p.second);
    int g3=std::get<3>(p.second);
    float avgid=g2/g3;
    thisctglinks.push_back(std::make_tuple(chrname,g0,g1,avgid));
  }
  if(!mctgs.count(thisctg)){
    mctgs[thisctg]=thisctglinks;
  }
  
  /*float cov=(fabs(sum)*100.0)/aligns[0].draftlength;
  if(cov>mismin){
    cout << "  About " << std::fixed << std::setprecision(1)  
	 << cov << "% of "  << ctgname
	 << " aligns to " << thisname << std::endl;
    sigals++;
  }*/
  cout.imbue(std::locale(""));


  if(thisctglinks.size()==1){
    string thischr=std::get<0>(thisctglinks[0]);

    myfile << " "  << std::fixed << std::setprecision(1)  << thisctg << " :  Significant alignment only to chromosome "
	 << thischr
	 << " :  Total length aligned=" << std::get<1>(thisctglinks[0]) 
	 << "  (of which " << std::get<2>(thisctglinks[0]) << " unmapped)"
	 << std::fixed << std::setprecision(2)
	 << " covers about " <<  std::get<1>(thisctglinks[0])*100./std::get<0>(mchrmapped[thischr])<< "% of the chr"  //<< " " << std::get<0>(mchrmapped[thischr])
	 << std::fixed << std::setprecision(1)
	 << " with an average identity " <<  std::get<3>(thisctglinks[0])<< "%"
	 << endl;
  }else{
    myfile << " " << thisctg 
	 << " :  Multiple chromosome alignments: " << endl;
    misfile << " " << thisctg 
	 << " :  Multiple chromosome alignments: " << endl;
    nmisctg++;
    misassembledctgs.push_back(thisctg);

    int nbreaks=-1;
    for ( const auto &p : thisctglinks ){
      nbreaks++;
      myfile  << std::fixed << std::setprecision(1)  
	      << "                  chromosome " << std::get<0>(p) 
	      << " :  Total length aligned=" << std::get<1>(p) 
	      << "  (of which " << std::get<2>(p) << " unmapped)"
	      << std::fixed << std::setprecision(2)
	      << " covers about " <<  std::get<1>(p)*100./std::get<0>(mchrmapped[ std::get<0>(p)]) << "% of the chr"//<< " " << std::get<0>(mchrmapped[ std::get<0>(p)]) 
	      << std::fixed << std::setprecision(1)
	      << " with an average identity " <<  std::get<3>(p) << "%"
	      << endl;
      
      misfile  << std::fixed << std::setprecision(1)  
	       << "                  chromosome " << std::get<0>(p) 
	       << " :  Total length aligned=" << std::get<1>(p) 
	       << "  (of which " << std::get<2>(p) << " unmapped)"
	       << std::fixed << std::setprecision(2)
	       << " covers about " <<  std::get<1>(p)*100./std::get<0>(mchrmapped[ std::get<0>(p)]) << "% of the chr" 
	       << std::fixed << std::setprecision(1)
	       << " with an average identity " <<  std::get<3>(p) << "%"
	       << endl;

   }
    nmis+=nbreaks;
  }

 return(0);
}




int restofmain(){
  CheckIfAllChrs();
 
 if(1 == 0){
   pri=0;
   for(int i=0; i<allchrs.size();i++){  // using allchr so always have same order as in the ref
     if(std::find(chrs.begin(), chrs.end(), allchrs[i]) != chrs.end()) {
       CheckChr(allchrs[i],i);
     }
   }

   /*  out << "Chr cov ";
   for(int i=0;i<chrcov.size();i++){
     out<< std::setprecision(2) << " (" << 
       chrcov[i] << "," << chrcovn[i] << ")" << " " ;
     out1 << std::setprecision(2) << chrcov[i] <<  " " << chrcovn[i]; 
     std::sort(chrcovperc[i].begin(), chrcovperc[i].end(),mysort);
     for(int j=0;j<chrcovn[i];j++)out1 << " " << chrcovperc[i][j];
     out1 << std::endl;
   }
   out << std::endl;


   std::cout << out.str();
   out.clear();
   out.str("");
   o2file << out1.str();
   out1.clear();
   out1.str("");

   out << "Chromosomes coverage: " << std::endl;
   for(int i=0;i<chrcov.size();i++)
     out<< std::setprecision(2) << " (" << chrcovname[i] << "," << 
       chrcov[i] << "," << chrcovn[i] << ")" << " " ;
   out << std::endl << std::endl;


   ofile << out.str();
   out.clear();
   out.str("");

   out << "The following " << onecontig  << " chromosomes are in 1 contig (at least " 
       << covtosingle << "% of chromosome covered by a single contig):" 
       << std::endl;
   for(int i=0;i<singlectg90.size();i++)
     out<< singlectg90[i] << " ";
   out << std::endl;

   out << "The following chromosomes are in 1 contig (at least "
       << "80% of chromosome covered by a single contig):"
       << std::endl;
   for(int i=0;i<singlectg80.size();i++)
     out<< singlectg80[i] << " ";
   out << std::endl;

   out << "The following chromosomes are in 1 contig (at least "
       << "70% of chromosome covered by a single contig):"
       << std::endl;
   for(int i=0;i<singlectg70.size();i++)
     out<< singlectg70[i] << " ";
   out << std::endl;


   out << "\nThere are possibly " << tocheck.size() 
       << " misassembled contig ...checking... " << std::endl;
   
   ofile << out.str();
   out.clear();
   out.str("");


   tocheck.erase( Unique( tocheck.begin(), tocheck.end() ), tocheck.end() );

   for(int i=0; i<tocheck.size();i++){
     CheckCtg(tocheck[i]);
   }

   if(misassembled){
     out << "\n... Possibly " << misassembled 
	 << " contigs are misassembled (more than "
	 << mismin << "% of them align to more than 1 chromosome): ";
     out1 << "Mis " << misassembled;
     for(int m=0; m<misass.size();m++){ out << misass[m] ; if(m<misass.size()-1) out << ", ";}
   }else{
     out << "\n... If something is misassembled is very small (<" << mismin << "% of contig) " << std::endl;
     out1 << "Mis " << misassembled;
   }
   std::cout << out1.str() << std::endl;
   ofile << out.str() << std::endl << std::endl;

 }else{
   std::cout<< "\nERROR:: Alignment file not properly read!! \n" << std::endl;
   exit(1);
 }

 ofile.close();
   */

  

 return(0);
  
 } /* end of the main */
}





// ****************************** //
int CheckCtg(std::string ctgname)
// ****************************** //
{
  std::stringstream out;
  std::stringstream out1;

  std::vector<ALIGNMENTS> aligns;
  std::copy_if(myaligns.begin(), myaligns.end(), std::back_inserter(aligns),
	       [ctgname](const ALIGNMENTS& a) {return a.draftname==ctgname;});

  std::vector<std::string> tctgs;
  std::transform(aligns.begin(), aligns.end(), std::back_inserter(tctgs),
		 [](ALIGNMENTS const& x) { return x.refname;});
  tctgs.erase( Unique( tctgs.begin(), tctgs.end() ), tctgs.end() );
  
  if(no)out << "\n" << ctgname 
	    << " length " << aligns[0].draftlength
	    << " is aligned to " << tctgs.size() << " chromosomes " 
	    << std::endl;
 
  if(tctgs.size()==1)return(0);

  
  int sigals=0;
  for(int i=0; i<tctgs.size();i++){

    std::string thisname=tctgs[i];

    std::vector<ALIGNMENTS> taligns;
    std::copy_if(aligns.begin(), aligns.end(), std::back_inserter(taligns),
                 [thisname](const ALIGNMENTS& a) {return a.refname==thisname;});
    
    int sum=0;
    std::for_each(taligns.begin(),taligns.end(), [&] (ALIGNMENTS const& a) {
	sum+=(a.draftstop-a.draftstart);
      });
    
    float cov=(fabs(sum)*100.0)/aligns[0].draftlength;
    if(cov>mismin){
      out << "  About " << std::fixed << std::setprecision(1)  
	  << cov << "% of "  << ctgname
	  << " aligns to " << thisname << std::endl;
      sigals++;
    }
  }
  if(sigals==1){
    out1 << " " << ctgname << " Significant alignment only to one chromosome, the rest of alignments <" 
	 << mismin << "% of contig " 
	 << std::endl;
  }else{
   out1 << ctgname 
	 << std::endl;
    misassembled++;
    misass.push_back(ctgname);
  }


  if(pri)std::cout << out1.str() << out.str();
  ofile << out1.str() << out.str();
  
  out.clear();
  out.str("");
  out1.clear();
  out1.str("");

 return(0);
}


// ****************************** //
int CheckIfAllChrs(void)
// ****************************** //
{
  FILE *namef;
  std::stringstream out;
  std::stringstream out1;
  char line[2000]={0};
  namef = fopen("chrlist.txt","r");

  std::vector<std::string> tallchrs;
  while (fgets(line,2000,namef))   
    {
      std::istringstream split(line);
      for(std::string each; getline(split, each, '\n'); 
	  tallchrs.push_back(each));
    }
  
  allchrs=tallchrs;
  if(chrs.size() != tallchrs.size()){
    out << " !! There is/are " << tallchrs.size()-chrs.size() 
	<< " Chrs not aligned at all:  " ;
   
    for(int i; i<chrs.size();i++)
      tallchrs.erase(std::remove(tallchrs.begin(), tallchrs.end(), chrs[i]), tallchrs.end());
    for(int i; i<tallchrs.size();i++)
      out << tallchrs[i] << " " ;
    out << std::endl<< std::endl;

  }else{
    out << " All Chrs aligned " << std::endl<< std::endl;
  }
  out1 << "Chrs aligned " << chrs.size()<< " Not aligned " 
       << allchrs.size()-chrs.size() << std::endl;

  std::cout << out1.str();
  ofile << out.str();
  out.clear();
  out.str("");
  return(0);
}
// ****************************** //
int CheckChr(std::string chrname, int nchr)
// ****************************** //
{	       
  std::vector<ALIGNMENTS> aligns;
  std::copy_if(myaligns.begin(), myaligns.end(), std::back_inserter(aligns),
	       [chrname](const ALIGNMENTS& a) {return a.refname==chrname;});

  std::sort(aligns.begin(), aligns.end(), alignments);

  std::vector<std::string> tctgs;
  std::transform(aligns.begin(), aligns.end(), std::back_inserter(tctgs),
		 [](ALIGNMENTS const& x) { return x.draftname;});
  tctgs.erase( Unique( tctgs.begin(), tctgs.end() ), tctgs.end() );

  float thischrcov=0;   
  int thischrcovn=0;  
  std::vector<float> thisperc;
  for(int j=0; j<tctgs.size();j++){  // loop over contigs aligned to this chromosome

    std::stringstream out;
    if(j==0) out << "\n************************************************\n"
		 << "Chromosome " << chrname 
		 << " (length " << aligns[0].reflength
		 << "): there are " << tctgs.size() 
		 << " contigs aligned to it. "
		 <<  std::endl <<  std::endl;
    
    std::string thisname=tctgs[j];
    
    std::vector<ALIGNMENTS> taligns;
    std::copy_if(aligns.begin(), aligns.end(), std::back_inserter(taligns),
		 [thisname](const ALIGNMENTS& a) {return a.draftname==thisname;});
 
    

    int first=std::min(taligns[0].refstart,taligns[0].refstop);
    int last=std::max(taligns[taligns.size()-1].refstart,taligns[taligns.size()-1].refstop);   
    int forw=0, rev=0, totbases=0, cont=0;
    float id=0.;
    int prevstop;
    
 
    for(int i=0; i<taligns.size();i++){      
      if(fabs(taligns[i].refstart-prevstop)<maxdist)cont++;
      totbases+= fabs(taligns[i].draftstop-taligns[i].draftstart);
      id+= taligns[i].identity;
      if(taligns[i].strand=="F") forw++; else rev++; 
      prevstop=taligns[i].refstop;
    }

    id /= taligns.size();
 
    
    float cov=(fabs(totbases)*100.0)/aligns[0].reflength;
    int F=-1;
    if(forw==taligns.size()) F=1;
    else if(rev==taligns.size()) F=0;
    

   
    if(cov >= 2){
      thisperc.push_back(cov);
      thischrcov+=cov;
      thischrcovn++;   

      if(cov >= covtosingle){
	onecontig++;
	out << " " << thisname 
	    << "  Single-contig chromosome! "  << std::fixed << std::setprecision(1)
	    << cov << "% coverage!, " <<  totbases << " bases aligned, identity=" 
	    << id << "%"
	    << std::endl; 
	singlectg90.push_back(chrname);
      } else{ 
	out  << " " << thisname  << ", " << std::fixed << std::setprecision(1)
	     << cov << "% coverage "  <<  totbases << " bases aligned, identity=" 
	    << id << "%" << std::endl;      
      }
     
    if(cov >= 80){
        out << " " << thisname
            << "  Single-contig chromosome! "  << std::fixed << std::setprecision(1)
            << "80% coverage!, " <<  totbases << " bases aligned, identity="
            << id << "%"
            << std::endl;
        singlectg80.push_back(chrname);
      } 
    if(cov >= 70){
        out << " " << thisname
            << "  Single-contig chromosome! "  << std::fixed << std::setprecision(1)
            << "70% coverage!, " <<  totbases << " bases aligned, identity="
            << id << "%"
            << std::endl;
        singlectg70.push_back(chrname);
      }

      out << "  " << taligns.size() << " pieces and " << totbases
	  << " bases aligned to this chromosome " 
	  << "\n  with " << cont << " contigous pieces, " ;
	  
      
      if(forw==taligns.size())  out << " in Forward direction " << std::endl;
      else if(rev==taligns.size())  out << " in Reverse direction " << std::endl;
      else  out << " with " << forw << " pieces in Forward direction and " 
		<< rev << " in Reverse direction " << std::endl;
      
      if(taligns[0].draftlength > taligns[0].reflength
	 || (totbases*100.0)/taligns[0].draftlength < covtosingle){ // possible misassembly
	tocheck.push_back(thisname);
      }
      
      if(pri)std::cout << out.str() ;
      ofile <<  out.str() ;
      out.clear();
      out.str("");
   }else{
      out  << " " << thisname << " aligns to about "  << std::fixed << std::setprecision(1) 
	     <<  cov << "% of the chromosome, " <<  totbases << " bases aligned " << std::endl;
      
      if(pri)std::cout << out.str();
      ofile <<  out.str();
      out.clear();
      out.str("");
    }
    
  }// end loop on contigs
  
  
  chrcovperc.push_back(thisperc);


  chrcov.push_back(thischrcov);
  chrcovn.push_back(thischrcovn);	
  chrcovname.push_back(chrname);	
  std::stringstream out;
  out  << std::endl << " This chromosome covered for " << std::setprecision(2) << thischrcov << "% by "
       << thischrcovn << " contigs" << std::endl << std::endl;
  if(pri)std::cout << out.str();
  ofile <<  out.str();
  out.clear();
  out.str("");


  return(0);
}

// ******************* //
int Read(FILE *namef)
// ******************* //
{
  nseq=-1;
  char line[2000]={0};
  std::stringstream out;

  
  while (fgets(line,2000,namef))   
    {
      nseq++;
      if(nseq<printsome)pri=1;
      else pri=0;
      
      std::vector<std::string> words;
      std::istringstream split(line);
      for(std::string each; getline(split, each, ' '); 
	  words.push_back(each));
      
      
      ALIGNMENTS thisal;
      
      thisal.score= to_int(words[0]);
      thisal.refname = words[1];
      thisal.draftname = words[2];
      thisal.refstart =  to_int(words[3]);
      thisal.refstop = to_int(words[4]);
      thisal.draftstart =  to_int(words[5]);
      thisal.draftstop =  to_int(words[6]);
      thisal.strand = words[7];
      thisal.numbases =  to_int(words[8]);
      thisal.identity = to_int(words[9]);
      thisal.reflength   =  to_int(words[10]);
      thisal.draftlength =  to_int(words[11]);
      
      myaligns.push_back(thisal);
  
      if(no) std::cout << nseq << " " << words[1] 
		       << " " << words[3] << std:: endl;
     
      thisal = ALIGNMENTS(); // re-initialize thisal
    }
  fclose(namef);

  

  std::transform(myaligns.begin(), myaligns.end(), std::back_inserter(chrs),
               [](ALIGNMENTS const& x) { return x.refname; });
  std::transform(myaligns.begin(), myaligns.end(), std::back_inserter(ctgs),
               [](ALIGNMENTS const& x) { return x.draftname; });
  
  chrs.erase( Unique( chrs.begin(), chrs.end() ), chrs.end() );
  ctgs.erase( Unique( ctgs.begin(), ctgs.end() ), ctgs.end() );
  
  out << "\nTotal alignments: " <<  nseq << " shreded pieces " << std:: endl;
  out << "Reference has " << chrs.size() << " chromosomes " << std:: endl;
  out << "Query has " << ctgs.size() << " contigs\n" << std:: endl;
  

  //std::cout << out.str();
  ofile <<  out.str();
  out.clear();
  out.str("");

  return(0);
}

