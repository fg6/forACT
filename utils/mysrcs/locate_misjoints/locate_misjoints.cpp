#include "../myinc/macro.h"  //misf.h"
#include <ctime>
#include <locale> // for imbue
#include <typeinfo>
static std::map<string, std::pair<long int,int> > mchrmapped;    // one el per chr: 0 or 1 if it is not (0) or is (1) mapped READY
static std::map<string, std::pair<long int,int> > mctgmapped;    // one el per ctg: 0 or 1 if it is not (0) or is (1) mapped READY
//
static std::map<string, vector<std::tuple<string,long int,long int,float>> > mctgs; // one ele per ctg: a vector of chrs mapped to the ctg READY
static std::map<string, vector<std::tuple<string,long int,long int,float, long int, long int>>> mchrs;  // one ele per chr: a vector of ctgs mapped to the ctg
//
static std::vector<std::pair<string,int> > fastachrs;  //vector of chr names, length READY
static std::vector<std::tuple<string,int,string> > fastactgs;  // vector of ctg names, length  READY
static std::map<string, string> ctg_seqs;    // one el per ctg: 0 or 1 if it is not (0) or is (1) mapped READY

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

 
  //bool operator() (nALIGNMENT& i, nALIGNMENTS& j) { 
  // return(std::min(i.ctgi,i.ctgf) < std::min(j.ctgi,j.ctgf)); 
  //}

} nalignments;

static int readals(char* file);
static int CheckOneCtg(std::vector<nALIGNMENTS> ctgals);
static int CheckChrs(void);
static int LocateMisJ(std::vector<nALIGNMENTS> ctgals);
static int AnalyseCtg(std::vector<nALIGNMENTS> ctgals);

static std::ofstream misfile;
static std::ofstream chrfile;
static std::ofstream chrdetails;


// possible misjoint selection:
static  float min_len_perc = 0.25; // perc min_lenght to consider as possible misjoint wrt major al (reduce noise and small repeats)
static  long int min_len = 500000; // min_lenght to consider as possible misjoint
static  long int min_len_for_max=100000; // min_lenght for a major al


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
  int err=0; int saveinfo=1; int readseq=1;
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
  err=readfasta(argv[2],saveinfo,"",readseq);
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }
  for(int c=0; c<rname.size(); c++){
    mctgmapped[rname[c]]=std::make_pair(rlen[c],0);
    fastactgs.push_back(std::make_tuple(rname[c],rlen[c],rseq[c]));
    ctg_seqs[rname[c]]=rseq[c];
  }
  rname.clear();
  rlen.clear();

  
  //  ****************************** //
  //  ****** READ ALIGNMENTS ******* //
  //  ****************************** //
  readals(argv[3]);


  if(0)cout << " Possibly " << nmisctg
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
  int nmappedctgs=0; // ?? not counting nn=0 so adding it here
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
    

    // how much of the ref and this chr is covered //
    refcov+=std::abs(chrf-chri);
    if(!chrcovmap.count(chr)) chrcovmap[chr]=std::abs(chrf-chri);
    else chrcovmap[chr]+=std::abs(chrf-chri);

    //this chr is mapped, add in mchrmapped
    if(!mchrmapped.count(chr))cout << " Error chr missing in map! " << chr << endl;
    std::get<1>(mchrmapped[chr])=1;
    
    //this ctg is mapped, add in mctgmapped  // can reduce time significantly if only checked once
    if(!mctgmapped.count(ctg))cout << " Error ctg missing in map! " << ctg << endl;
    std::get<1>(mctgmapped[ctg])=1;


    // save alignment in struct
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
      // save info from previous contig, before moving on
      nmappedctgs++;
      //CheckOneCtg(ctgals);
      AnalyseCtg(ctgals); 

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

  //catch last contig
  nmappedctgs++;

  //CheckOneCtg(ctgals); 
  AnalyseCtg(ctgals); 
  myfile.close();   
  misfile.close();   

  /* pause for now:
     CheckChrs();


  int nmappedchrs=0;
  for ( const auto &p : mchrmapped ){ 
    if(std::get<1>(p.second)==1) nmappedchrs++;
    else unmappedchrs.push_back(p.first);
    }*/
  int nmappedchrs=0;

  cout << "\n Found " << nn+1 << " alignments, "   // nn is from 0
       << nmappedctgs << " contigs, and " 
       << nmappedchrs << " chromosomes " 
       << endl;


  if(0)cout << std::fixed << std::setprecision(1) 
       << " Reference coverage: " <<  refcov*100./refsize << "%"
       << " (" << refcov << " bases)"  
       << endl;

  return 0;
}

// ******************************************** //
int AnalyseCtg(std::vector<nALIGNMENTS> ctgals)
// ******************************************** //
{
  //
  // in ctgals struct:   chr, ctg, chri, chrf, ctgi, ctgf, strand, albases, id
  //

  int sum=0;
  string thisctg=ctgals[0].ctg;
  std::map<string, long int> len_mappedtochr; 
  std::map<string, long int> thisctglinks;
  vector<string> linkedchr;


  if(1) cout << "\n Analising contig " << thisctg << endl;

  std::for_each(ctgals.begin(),ctgals.end(), [&] (nALIGNMENTS const& a) {
      long int naligned=std::abs(a.ctgf-a.ctgi);
      if(!len_mappedtochr.count(a.chr)){
	len_mappedtochr[a.chr] = naligned;
      }else{
	len_mappedtochr[a.chr] += naligned;
      }
    });
 
  // find major alignment
  std::pair<string, int> major_al = *std::max_element(len_mappedtochr.begin(), len_mappedtochr.end(),
				[](const std::pair<string, long int>& p1, const std::pair<string, long int>& p2) {
				  return p1.second < p2.second; });
  
  long int max_al = major_al.second; 
  string major_chr= major_al.first ; 
  

  // ************* Select out ************* //
  // not considering too small alignments for locating misjoints
  if(max_al < min_len_for_max)  return 0;
  
  
  int major_als=0;
  for ( const auto &p : len_mappedtochr ){ 
    long int al_len=p.second; //length aligned
    
    if (al_len*1./max_al >= min_len_perc || al_len > min_len){
      major_als++;
      string chrname=p.first;
      if(0)cout << "   Major alignment " << chrname 
	   << " al lenght: "<<  al_len ;
      if(chrname == major_chr ) if(0)cout << " Major chr " << endl;
      else if(0)cout  << endl;
      thisctglinks[chrname]=al_len;
      linkedchr.push_back(chrname);
    }
  }

  if(0)cout << "\n   This ctg aligns mostly to " <<  major_als << " chr(s)" << endl;

  // ************* Select out ************* //
  // only a major mapping
  if(major_als == 1 )  return 0;


  // Remove from ctgals all alignments to chrs others than majors
  ctgals.erase(std::remove_if(ctgals.begin(), ctgals.end(), [&](nALIGNMENTS const& x){
	return !(std::find(linkedchr.begin(), linkedchr.end(), x.chr) != linkedchr.end()); }), ctgals.end());

  // all remaining alignments:
  int ff=0;
  std::for_each(ctgals.begin(),ctgals.end(), [&] (nALIGNMENTS const& a) {
      ff++;
    });


  LocateMisJ(ctgals);
 
  return(0);
}

// ******************************************** //
int LocateMisJ(std::vector<nALIGNMENTS> ctgals)
// ******************************************** //
{
  // sort alignment by position in contig
  std::stable_sort(ctgals.begin(), ctgals.end(), [](const nALIGNMENTS &a, const nALIGNMENTS &b){
      return ( std::min(a.ctgi,a.ctgf) < std::min(b.ctgi,b.ctgf) );
    });

  string ctg=ctgals[0].ctg;
  vector<long int> checksort;
  vector<string> chrs;
  vector<string> blocks;
  string previous_chr=ctgals[0].chr;

  long int previous_pos;
  int started=0;
  long int last_pos;

 
  std::stringstream ss;
  ss.imbue(std::locale(""));
  ss << "   Possible misjoint locations:" << endl;
  
  // only works if sorted by ctg position:
  std::for_each(ctgals.begin(),ctgals.end(), [&] (nALIGNMENTS const& a) {
      checksort.push_back(std::min(a.ctgi,a.ctgf));
      blocks.push_back(a.chr);  
      
      long int this_pos=std::min(a.ctgi,a.ctgf);

      if(a.chr != previous_chr){
	string thisseq = ctg_seqs[a.ctg].substr(previous_pos,this_pos-previous_pos);
	std::size_t start_misjoint = thisseq.find_last_of("NNN") + previous_pos + 1; 
	//ss << endl << previous_pos << " " << this_pos << " string len: " << thisseq.size() << endl;

	std::size_t previous_misjoint_end;   
	if(started){
	  previous_misjoint_end = thisseq.find_first_of("NNN") + previous_pos + 1; 
	  ss << " and ends at " << previous_misjoint_end 
	     << endl;
	}

	ss << "       break between positions " 
	   << previous_pos << " and " << this_pos
	   << ":   Possibly block starts at " << start_misjoint;

	started=this_pos;
      }    


      last_pos=this_pos;
    
      previous_chr=a.chr;
      previous_pos=std::max(a.ctgi,a.ctgf);

    });
  if(0)cout << "   is ctgals sorted?" << std::is_sorted(checksort.begin(),checksort.end()) << endl;


  // unique works only if elements are sorted (same element close to each other,
  // as ctgals is sorted by ctg position, unique gives the blocks !
  blocks.erase( unique(blocks.begin(), blocks.end() ), blocks.end() );
  // to get only unique chrs, need to sort first:
  chrs=blocks;
  sort(chrs.begin(),chrs.end());
  chrs.erase( unique( chrs.begin(), chrs.end() ), chrs.end() );

  cout << "   This sequence aligns mostly to " << chrs.size() << " chromosomes "
       << " and it is divided in about " << blocks.size() 
       << " blocks, ordered as follows: " << endl << "      ";
  for(string n : blocks) { std::cout << " one_block_mapped_to " << n << ", "; }

  cout << endl <<  ss.str() << endl;




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
	 << " :  Multiple chromosome mapping: " << endl;
    nmisctg++;
    misassembledctgs.push_back(thisctg);

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
		 << "\tactual-coverage:" << std::get<1>(mchrs[thischr][c])*100./std::get<0>(mchrmapped[thischr]) << "%" // contig length mapped here
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
