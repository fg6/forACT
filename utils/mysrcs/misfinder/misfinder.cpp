#include "../myinc/macro.h"  //misf.h"
#include <ctime>
#include <locale> // for imbue



static std::map<string, std::pair<long int,int> > mchrmapped;    // one el per chr: 0 or 1 if it is not (0) or is (1) mapped READY
static std::map<string, std::pair<long int,int> > mctgmapped;    // one el per ctg: 0 or 1 if it is not (0) or is (1) mapped READY
//
static std::map<string, vector<std::tuple<string,long int,long int,float>> > mctgs; // one ele per ctg: a vector of chrs mapped to the ctg READY
static std::map<string, vector<std::tuple<string,long int,long int,float, long int, long int>>> mchrs;  // one ele per chr: a vector of ctgs mapped to the ctg
//
static std::vector<std::pair<string,int> > fastachrs;  //vector of chr names, length READY
static std::vector<std::pair<string,int> > fastactgs;  // vector of ctg names, length  READY
static  vector<string> unmappedchrs;  //vector of unmapped chrs READY
static int nmisctg=0;
static int nmis=0;
static long int  refsize=0;
static  vector<string> misassembledctgs;  //vector of unmapped chrs READY   Filled NOT used anywhere?
static std::map<string, long int > chrcovmap;    // one el per chr: 0 or 1 if it is not (0) or is (1) mapped READY

struct nALIGNMENTS {
  std::string chr,ctg,strand;
  long int chri,chrf,ctgf,ctgi;
  int albases;
  float id;

  bool operator() (nALIGNMENTS i, nALIGNMENTS j) { 
    return( i.albases > j.albases ); 
  }
} nalignments;

static int readals(char* file);
static int CheckOneCtg(std::vector<nALIGNMENTS> ctgals);
static int CheckChrs(void);

static std::ofstream misfile;
static std::ofstream chrfile;
static std::ofstream chrdetails;


// **** not in use **** //

//static int onecontig=0;
//static std::vector<std::string> chrs;
//static std::vector<std::string> ctgs;
//static std::vector<std::string> tocheck;
//static std::vector<std::string> singlectg90;
//static std::vector<std::string> singlectg80;
//static std::vector<std::string> singlectg70;
//static std::ofstream ofile;
// minimum alignment to declare single contig per chromosome
//static int covtosingle=90; // percentage of chromosome covered 
// minimum alignment to multimple chromosome to declare possible misassembly
//static int mismin=3; //at least 2% of the contig is aligned to a chromosome other than the main one
// maximum distance on reference between 2 pieces of a contig to declare countiguity 
//static int maxdist=2000; // max number of bases between the end position on ref of one piece of 
                         // contig and the start position of next piece
// how much of each chr is covered
//static std::vector<float> chrcov;
//static std::vector<std::vector<float>> chrcovperc;
//static std::vector<float> chrcovn;
//static std::vector<std::string> chrcovname;
//static  std::vector<std::string> allchrs;

//int CheckIfAllChrs(void);
//int CheckChr(std::string chrname, int nchr);
//int CheckCtg(std::string ctgname);


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

  //  ****************************** //
  //  **** READ REFERENCE CHRS ***** //
  //  ****************************** //
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
  //  ****************************** //
  //  ****** READ DRAFT CTGS ******* //
  //  ****************************** //
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

  
  //  ****************************** //
  //  ****** READ ALIGNMENTS ******* //
  //  ****************************** //
  readals(argv[3]);


  cout << " Possibly " << nmisctg
       << " contig/s is/are misassembled " 
       << " with a total of possibly " << nmis << " break points "
       << endl;

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
    
    if(!chrcovmap.count(chr)) chrcovmap[chr]=std::abs(chrf-chri);
    else chrcovmap[chr]+=std::abs(chrf-chri);

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

  CheckChrs();


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

 
  


  return 0;
}

// ******************************************** //
int CheckOneCtg(std::vector<nALIGNMENTS> ctgals)
// ******************************************** //
{
  int sum=0;
  string thisctg=ctgals[0].ctg;
  std::map<string, std::tuple<long int,long int,float,int, long int, long int> > mappedtochr; 
  vector<std::tuple<string,long int,long int,float> > thisctglinks;


  // add per each covered chr: min_pos, max_pos, are there holes in coverage? 
  //                           min_ctg_pos, max_ctg_pos, are there unmapped region of ctg?

  std::for_each(ctgals.begin(),ctgals.end(), [&] (nALIGNMENTS const& a) {
      long int naligned=std::abs(a.ctgf-a.ctgi);
      long int numapped=naligned-a.albases;

      sum+=std::abs(a.ctgf-a.ctgi); // total number of alignments
  
      if(!mappedtochr.count(a.chr)){
	int ii=1;
	long int minpos=std::min(a.chrf,a.chri);
	long int maxpos=std::max(a.chrf,a.chri);
	mappedtochr[a.chr]=std::make_tuple(naligned,numapped,a.id,1, minpos, maxpos);
      }else{
	long int sumtoals=std::get<0>(mappedtochr[a.chr])+naligned;
	long int sumtoun=std::get<1>(mappedtochr[a.chr])+numapped;
	int sumid=std::get<2>(mappedtochr[a.chr])+a.id;
	int nid=std::get<3>(mappedtochr[a.chr])+1;
	long int minpos=std::min(std::get<4>(mappedtochr[a.chr]),std::min(a.chrf,a.chri));
	long int maxpos=std::max(std::get<5>(mappedtochr[a.chr]),std::max(a.chrf,a.chri));
	
	mappedtochr[a.chr]=std::make_tuple(sumtoals,sumtoun,sumid,nid, minpos, maxpos);
      }
    });

  int tsum=0;
  int usum=0;
  for ( const auto &p : mappedtochr ){ 
    tsum+=std::get<0>(p.second);  //length aligned
    usum+=std::get<1>(p.second);  // unmapped
    string chrname=p.first;
    long int g0=std::get<0>(p.second);//length aligned
    long int g1=std::get<1>(p.second); // unmapped

    float g2=std::get<2>(p.second);
    long int g3=std::get<3>(p.second);
    float avgid=g2/g3;
    
    long int minpos=std::get<4>(p.second);
    long int maxpos=std::get<5>(p.second);
    thisctglinks.push_back(std::make_tuple(chrname,g0,g1,avgid));
    
    mchrs[chrname].push_back(std::make_tuple(thisctg,g0,g1,avgid,minpos,maxpos));
    
  }
  if(!mctgs.count(thisctg)){
    mctgs[thisctg]=thisctglinks;
  }

  cout.imbue(std::locale(""));


  if(thisctglinks.size()==1){
    string thischr=std::get<0>(thisctglinks[0]);

    myfile << " "  << std::fixed << std::setprecision(1)  << thisctg << " :  Map significantly only to chromosome "
	   << thischr
	   << ":  contig length= " << std::get<0>(mctgmapped[thisctg])
	   << ", Total length aligned=" << std::get<1>(thisctglinks[0]) 
	   << "  (of which " << std::get<2>(thisctglinks[0]) << " unmapped)"
	   << std::fixed << std::setprecision(2)
	   << " covers ~" <<  std::get<1>(thisctglinks[0])*100./std::get<0>(mchrmapped[thischr])<< "% of the chr"  
	   << " from ~" << std::get<4>(mappedtochr[thischr]) << " to ~" << std::get<5>(mappedtochr[thischr])
	   << std::fixed << std::setprecision(1)
	   << " with an average identity " <<  std::get<3>(thisctglinks[0])<< "%"
	   << endl;
  }else{
    myfile << " " << thisctg 
	 << " :  Multiple chromosome mapping: " << endl;
    misfile << " " << thisctg 
	 << " :  Multiple chromosome mapping: ";
    nmisctg++;
    misassembledctgs.push_back(thisctg);

    //vector<string> links;
    for ( const auto &p : thisctglinks ){
      misfile <<std::get<0>(p) << ", ";
	//links.push_back(std::get<0>(p));
    }
    misfile <<endl;

    int nbreaks=-1;
    for ( const auto &p : thisctglinks ){
      string thischr=std::get<0>(p);
      nbreaks++;
      myfile  << std::fixed << std::setprecision(1)  
	      << "                  chromosome " << std::get<0>(p) 
	      << ":  contig length= " << std::get<0>(mctgmapped[thisctg])
	      << ",  Total length aligned=" << std::get<1>(p) 
	      << "  (of which " << std::get<2>(p) << " unmapped)"
	      << std::fixed << std::setprecision(2)
	      << " covers ~" <<  std::get<1>(p)*100./std::get<0>(mchrmapped[ std::get<0>(p)]) << "% of the chr"//<< " " << std::get<0>(mchrmapped[ std::get<0>(p)]) 
	      << " from ~" << std::get<4>(mappedtochr[thischr]) << " to ~" << std::get<5>(mappedtochr[thischr])
	      << std::fixed << std::setprecision(1)
	      << " with an average identity " <<  std::get<3>(p) << "%"
	      << endl;
      
      misfile  << std::fixed << std::setprecision(1)  
	       << "                  chromosome " << std::get<0>(p) 
	       << ":  contig length= " << std::get<0>(mctgmapped[thisctg])
	       << ",  Total length aligned=" << std::get<1>(p) 
	       << "  (of which " << std::get<2>(p) << " unmapped)"
	       << std::fixed << std::setprecision(2)
	       << " covers ~" <<  std::get<1>(p)*100./std::get<0>(mchrmapped[ std::get<0>(p)]) << "% of the chr" 
	       << " from ~" << std::get<4>(mappedtochr[thischr]) << " to ~" << std::get<5>(mappedtochr[thischr])
	       << std::fixed << std::setprecision(1)
	       << " with an average identity " <<  std::get<3>(p) << "%"
	       << endl;

   }
    nmis+=nbreaks;
  }


 return(0);
}




// ****************************** //
int CheckChrs(void)
// ****************************** //
{

  chrfile.open("chromosomes_report.txt");
  chrdetails.open("chromosomes_details.txt");

  for (  const auto &p : fastachrs ){
    string thischr=std::get<0>(p);
    
    float thiscov=0;
    float structcov=0;
    int ctglinked=0;
    double maxctg=0.;
    long int minpos=std::get<0>(mchrmapped[thischr]); //chr length
    long int maxpos=0;
    chrdetails<< std::fixed << std::setprecision(1)
	      << thischr << " num_of_ctg:" <<  mchrs[thischr].size()<< endl;
    for (int c=0; c< mchrs[thischr].size(); c++){
      // struc coverage for details need to account for fas away mapped areas, see seabass_22 for instance, fAnaTes1_78
      chrdetails << "\t\t" << std::fixed << std::setprecision(2)
		 << std::get<0>(mchrs[thischr][c])   // contig 
		 << "\tmin " << std::get<4>(mchrs[thischr][c]) << "\tmax "<< std::get<5>(mchrs[thischr][c])
		 << "\tbig-structure-coverage:" 
		 << std::abs(std::get<4>(mchrs[thischr][c])- std::get<5>(mchrs[thischr][c]))*100./std::get<0>(mchrmapped[thischr])<< "%"
		 << "\tactual-coverage:" << std::get<1>(mchrs[thischr][c])*100./std::get<0>(mchrmapped[thischr]) << "%"   // contig length mapped here
		 << endl;
      thiscov+=std::get<1>(mchrs[thischr][c])*100.;
      structcov+=std::abs(std::get<4>(mchrs[thischr][c])- std::get<5>(mchrs[thischr][c]))*100./std::get<0>(mchrmapped[thischr]);
      ctglinked++;
      maxctg=std::max(maxctg, std::get<1>(mchrs[thischr][c])*100./std::get<0>(mchrmapped[thischr]));

      minpos=std::min(minpos,std::min(std::get<4>(mchrs[thischr][c]),std::get<5>(mchrs[thischr][c])));
      maxpos=std::max(maxpos,std::max(std::get<4>(mchrs[thischr][c]),std::get<5>(mchrs[thischr][c])));
      
    }
    chrdetails<<endl;
    thiscov/=std::get<0>(mchrmapped[thischr]);
    
    chrfile << std::fixed << std::setprecision(1)
	 << thischr << "\tmapped contigs: " << ctglinked
	 << "\tlongest " << maxctg  << "%" 
	    << "\tbig-structure-coverage:" 
	    <<  std::abs(minpos-maxpos)*100./std::get<0>(mchrmapped[thischr])
	    << "\tactual-coverage ~" << thiscov << "%" << endl;
  }
  chrfile.close();   
  chrdetails.close();   

  return(0);
}

// ****************************** //
//int CheckIfAllChrs(void)
// ****************************** //
/*{
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
*/





