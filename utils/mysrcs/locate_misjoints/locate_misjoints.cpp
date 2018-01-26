#include "../myinc/macro.h"  //misf.h"
#include <ctime>
#include <locale> // for imbue
#include <typeinfo>


static std::map<string, long int> ref_seqs;    // one el per ctg: 0 or 1 if it is not (0) or is (1) mapped READY
static std::map<string, int> chrs_found;    // one el per ctg: 0 or 1 if it is not (0) or is (1) mapped READY
static std::map<string, string> ctg_seqs;    // one el per ctg: 0 or 1 if it is not (0) or is (1) mapped READY
static std::ofstream major_file;
static std::ofstream bedfile;

static vector<string> chrs;

struct nALIGNMENTS {
  std::string chr,ctg,strand;
  long int chri,chrf,ctgf,ctgi;
  int albases;
  float id;
} nalignments;

static int ReadAls(char* file);
static int LocateMisJ(std::vector<nALIGNMENTS> ctgals, std::pair<string, int> major_al);
static int AnalyseCtg(std::vector<nALIGNMENTS> ctgals);
static int checknewblock(long int pchri,long int nchri, long int pchrf, long int nchrf);
static int nbreak=0;
static int verbose=1;

// possible misjoint selection:
static  float min_len_perc = 100.; // perc min_lenght to consider as possible misjoint wrt major al (reduce noise and small repeats)
static  long int min_len = 100000; // min_lenght to consider as possible misjoint
static  long int min_len_for_max=100000; // min_lenght for a major al
static  long int newblock=5000000; // only for extreme cases (also accounts for the fact that als not sorted by reference)
//static  long int newblock=30000000; // only for extreme cases (also accounts for the fact that als not sorted by reference)
                                      // with this option plot is fine (one line per block) but number of breakings overestimated, as some blocks should merge bc ref-pos not sorted
                                       // take number of breaks from forACT
//static  long int newblock=50000000; 

int main(int argc, char **argv)
{
  
  if(argc < 3){
    printf("Usage: %s  <draft> <alignment_file>  <chr.list> <newblock> <min_len_max> <min_len> \n", argv[1]);
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
  if (argc >= 5)  newblock = atol(argv[4]);
  if (argc >= 6)  min_len_for_max = atol(argv[5]);
  if (argc >= 7)  min_len = atol(argv[6]);


  if(1)cout << endl<< " Selected filters: " << endl
       << "  Min length for longest alignment to a single chr: " << min_len_for_max << " bp" <<endl
       << "  Min length for shorter alignments to supplementary chrs: " 
       << min_len << " bp" << endl;

  

  int err=0; int saveinfo=1; int readseq=1;
 
  //  ****************************** //
  //  ****** READ REF CTGS ******* //
  //  ****************************** //
  //readseqlist(argv[3]);
  err=0;  
  err=readfasta(argv[3],saveinfo);
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }
  for(int c=0; c<rname.size(); c++){
    ref_seqs[rname[c]]=rlen[c];
    chrs_found[rname[c]]=0;
    //chrs.push_back(rname[c]);
  }
  chrs=rname;
  rname.clear();
  rlen.clear();
 
  //  ****************************** //
  //  ****** READ DRAFT CTGS ******* //
  //  ****************************** //
  err=0;  
  err=readfasta(argv[1],saveinfo,"",readseq);
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }
  for(int c=0; c<rname.size(); c++){
    ctg_seqs[rname[c]]=rseq[c];
  }
  rname.clear();
  rseq.clear();
  rlen.clear();

  
  //  ****************************** //
  //  ****** READ ALIGNMENTS ******* //
  //  ****************************** //
  
  myfile.open("alfile.txt");
  bedfile.open("bedfile.txt");
  major_file.open("majorfiles.txt");
  ReadAls(argv[2]);
  major_file.close();
  myfile.close();
  bedfile.close();


  for(int c=0; c<chrs.size(); c++){
    if(!chrs_found[chrs[c]] )
      cout << "Warning: missed chr "
	   << chrs[c] << endl;
  }

  std::stringstream tt; 
  if(0){
    if(nbreak)
       tt << endl << "\n **************  REPORT SUMMARY  *****************" << endl 
	 << "  ****** " << nbreak << " possible breaking points found ******" 
	 << endl <<   " *************************************************\n  "<< endl;
    else
      tt << endl << " No breaking points found " << endl;
    cout << tt.str() <<  endl ;
  }


  return 0;
}


int ReadAls(char* file){
  std::ifstream infile(file);
  string line;


  string prevctg="";
  string prevchr="";
  int nn=-1;
  int newctg=0;
  int nmappedctgs=0; // ?? not counting nn=0 so adding it here
  long int refcov=0;
  std::vector<nALIGNMENTS> ctgals;

  while(getline(infile,line)){
    std::stringstream ss(line);
    nn++; // count alignments read
    string ctg, chr;
    int ctgi, ctgf, chri,chrf,albases, kk;
    float id;
    
    ss >> ctg >> chr >> id >> albases  >> kk  >> kk  >> ctgi >> ctgf >> chri >> chrf >> kk >> kk;

    // how much of the ref and this chr is covered //
    refcov+=std::abs(chrf-chri);
 
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
     
    //check if new ctg:
    if(nn>0 && prevctg != ctg) newctg=1;
    else newctg=0;
    
    
    if(newctg || infile.eof()){
      // save info from previous contig, before moving on
      nmappedctgs++;
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


  AnalyseCtg(ctgals); 

  
  return 0;
}

// ******************************************** //
int AnalyseCtg(std::vector<nALIGNMENTS> ctgals)
// ******************************************** //
{
  string thisctg=ctgals[0].ctg;
  std::map<string, long int> len_mappedtochr; 
  std::map<string, long int> thisctglinks;
  vector<string> linkedchr;
  
  
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

  //if(thisctg.find("ynteny") != std::string::npos)  // only synteny groups, not unplaced // not valid in every case, do it later
  //cout << "writing majorfile " << endl;
  major_file<< thisctg << " " <<  ctg_seqs[thisctg].size() << " " << major_chr 
	    << " " << ref_seqs[major_chr] <<  endl;
  chrs_found[major_chr]=1; // to check that all chrs are mapped


 
  // ************* Select out ************* //
  // not considering too small alignments for locating misjoints

  if(max_al < min_len_for_max)  return 0;
  
  int major_als=0;
  for ( const auto &p : len_mappedtochr ){ 
    long int al_len=p.second; //length aligned
    
 
    if(al_len > min_len){
      major_als++;
      string chrname=p.first;
      
      thisctglinks[chrname]=al_len;
      linkedchr.push_back(chrname);
    }else{
      if(0) 
	cout << thisctg << " maj: " << major_chr << " " << p.first 
	     << " too small " << al_len << " min is " << min_len << endl;
    }
  }


  // ************* Select out ************* //
  // only a major mapping: no inter-chr misjoint
  if(major_als == 1 )  {
    return 0;
  }

 
  // Remove from ctgals all alignments to chrs others than majors
  ctgals.erase(std::remove_if(ctgals.begin(), ctgals.end(), [&](nALIGNMENTS const& x){
	return !(std::find(linkedchr.begin(), linkedchr.end(), x.chr) != linkedchr.end()); }), ctgals.end());

  LocateMisJ(ctgals, major_al);
  //LocateMisJ_plot(ctgals, major_al);

  return(0);
}

// ******************************************** //
int LocateMisJ(std::vector<nALIGNMENTS> ctgals, std::pair<string, int> major_al)
// ******************************************** //
{
  // sort alignment by position in contig
  std::stable_sort(ctgals.begin(), ctgals.end(), [](const nALIGNMENTS &a, const nALIGNMENTS &b){
      return ( std::min(a.ctgi,a.ctgf) < std::min(b.ctgi,b.ctgf) );
    });

  string ctg=ctgals[0].ctg;
  vector<string> chrs;
  vector<string> als_per_block;
  vector<long int> bases_per_block;
  vector<std::pair<long int, long int> > block_length;
  vector<string> blocks;
  vector<std::pair<long int, long int> > block_bp; //block breaking points
  vector<string> strand_per_block;

  long int max_al = major_al.second; 
  string major_chr= major_al.first ; 

    
  long int previous_pos;
  long int last_pos;

  std::stringstream tt;
  tt.imbue(std::locale(""));
  std::stringstream ss;
  tt << endl << " Synteny Group:  " << ctg << " major chr " <<  major_chr << endl;
  
  // Start set to 0
  string previous_chr=ctgals[0].chr;
  long int previous_chri=ctgals[0].chri;
  long int previous_chrf=ctgals[0].chrf;
  int ntt=1; // counter for summing aligned bases to each separate block
  bases_per_block.push_back(std::abs(ctgals[0].ctgf-ctgals[0].ctgi));// start aligned bases for first block
  long int min_pos=std::min(ctgals[0].ctgf,ctgals[0].ctgi);
  long int max_pos=std::max(ctgals[0].ctgf,ctgals[0].ctgi);
  block_length.push_back(std::make_pair(min_pos,max_pos));
  string thisseq = ctg_seqs[ctgals[0].ctg].substr(0,min_pos);
  block_bp.push_back(std::make_pair(thisseq.find_last_of("NNN") + 1, 0));
  long int this_pos=0;

  long int ref_min=std::min(ctgals[0].chrf,ctgals[0].chri);
  long int ref_max=std::max(ctgals[0].chrf,ctgals[0].chri);
  vector<std::pair<long int, long int> > block_length_inref;
  block_length_inref.push_back(std::make_pair(ref_min,ref_max));

  vector<string> chrblocks;
  chrblocks.push_back(ctgals[0].chr);


  float meanid=ctgals[0].id;
  int nals_block=1;
  vector<float> percid;

  int atleast_onelarge=0;
  int ccount=0;
  int plus_strand=0;
  int nstrand=0;
  std::for_each(ctgals.begin(),ctgals.end(), [&] (nALIGNMENTS const& a) {   // only works if sorted by ctg position:
 
      long int naligned=std::abs(a.ctgf-a.ctgi);
      this_pos=std::min(a.ctgi,a.ctgf);
      max_pos=std::max(a.ctgi,a.ctgf);
      
      ref_min=std::min(a.chrf,a.chri);
      ref_max=std::max(a.chrf,a.chri);
      als_per_block.push_back(a.chr);  
    
      
      // New Block
      if(a.chr != previous_chr
	 || ((a.chr == previous_chr) && a.chr != major_chr
	     && (checknewblock(previous_chri,a.chri,previous_chrf,a.chrf)))){  
	
	if(bases_per_block[ntt-1] > min_len 
	   && previous_chr != major_chr)atleast_onelarge++;

	// Previous 
	if ( plus_strand*1./nstrand > 0.8 ) { strand_per_block.push_back("+"); }
	else if ( plus_strand*1./nstrand < 0.2 ) { strand_per_block.push_back("-"); }
	else {
	  strand_per_block.push_back(".");
	  //cout << plus_strand << " " << nstrand << " " << plus_strand*1./nstrand << endl;
	}
       

	//if(a.chr != previous_chr) nbreak+=1;

	bases_per_block.push_back(naligned); // start aligned bases for new block
	block_length.push_back(std::make_pair(this_pos,max_pos));  
	block_length_inref.push_back(std::make_pair(ref_min,ref_max));
	chrblocks.push_back(a.chr);
	percid.push_back(meanid/nals_block);

	// get relevant part of sequence
	thisseq = ctg_seqs[a.ctg].substr(previous_pos,this_pos-previous_pos);
       	std::get<1>(block_bp[ntt-1])=thisseq.find_first_of("NNN") + previous_pos + 1;   // end of previous block at NNN  (-1, bc ntt starts from 1)
	ntt++;// new block
	
	// new block start
	block_bp.push_back(std::make_pair(thisseq.find_last_of("NNN") + previous_pos + 1, 0));
 
	// First al in block
	if ( a.strand == "F")
	  plus_strand=1;
	nstrand=1;

      }else{  // Continuing previous Block
	meanid+=a.id;
	nals_block++;

	bases_per_block[ntt-1]+=naligned; // sum aligned bases to present block
	block_length[ntt-1]=std::make_pair(std::min(this_pos,std::get<0>(block_length[ntt-1])),std::max(max_pos,std::get<1>(block_length[ntt-1])));  
	block_length_inref[ntt-1]=std::make_pair(std::min(ref_min,std::get<0>(block_length_inref[ntt-1])),std::max(ref_max,std::get<1>(block_length_inref[ntt-1])));  

	if ( a.strand == "F")
	  plus_strand+=1;
	nstrand+=1;

      }
        
      last_pos=this_pos;
      previous_chr=a.chr;
      previous_chri=a.chri;
      previous_chrf=a.chrf;
      previous_pos=std::max(a.ctgi,a.ctgf);
 
      ccount++;
    });
  // last end: 
  thisseq =ctg_seqs[ctgals[ctgals.size()-1].ctg].substr(previous_pos,this_pos-previous_pos);
  std::get<1>(block_bp[ntt-1])=thisseq.find_first_of("NNN") + previous_pos + 1; // end of previous block at NNN  (-1, bc ntt starts from 1)
  percid.push_back(meanid/nals_block);
  if ( plus_strand*1./nstrand > 0.8 ) strand_per_block.push_back("+");
  else if ( plus_strand*1./nstrand < 0.2 ) strand_per_block.push_back("-");
  else strand_per_block.push_back("NA");
  
  // previous block
  if(bases_per_block[ntt-1] > min_len 
     && previous_chr != major_chr)atleast_onelarge++;
 
  // do not write if there is not at least a large block to a not major chr
  if( ! atleast_onelarge )
    return 0;


  // unique works only if elements are sorted (same element close to each other,
  // as ctgals is sorted by ctg position, unique gives the blocks !
  //blocks=als_per_block;
  //blocks.erase( unique(blocks.begin(), blocks.end() ), blocks.end() );
  
  blocks=chrblocks;
  // to get only unique chrs, need to sort first:
  chrs=blocks;
  sort(chrs.begin(),chrs.end());
  chrs.erase( unique( chrs.begin(), chrs.end() ), chrs.end() );

  vector<int> block_len;
  for(string b : blocks) {
    block_len.push_back(std::count(als_per_block.begin(), als_per_block.end(), b));
  } 
  
  if(major_chr != "chrY") tt << "   This sequence aligns mostly to " << chrs.size() << " chromosomes "
       << " and it is divided in about " << blocks.size() 
       << " blocks, ordered as follows: " << endl;

  int ii=0;

  // only blocks > min_len:
  for(string n : blocks) {   
    if(bases_per_block[ii]>min_len){

      
      if(ii>0 && n != "chrY" && major_chr != "chrY") tt << endl;

      if(n != "chrY" && major_chr != "chrY")
	tt << "    Block " << ii+1 << " maps to " << n 
	   << ":  " << bases_per_block[ii] << " bp mapped, first pb: "
	   << std::get<0>(block_length[ii]) << ", last bp: " 
	   << std::get<1>(block_length[ii])  << "  covers " 
	   << std::get<1>(block_length[ii])-std::get<0>(block_length[ii]) 
	   << " bp, sub-scaff : " << std::get<0>(block_bp[ii]) << " - " << std::get<1>(block_bp[ii])
	   << " in ref: " << std::get<0>(block_length_inref[ii])  << " " <<  std::get<1>(block_length_inref[ii])
	   << " percid " << percid[ii];
      
      if(n != major_chr && n != "chrY" && major_chr != "chrY"){
	nbreak+=2;
	ss << major_chr << " " <<  ctg  //<< "synteny-group_" << major_chr.substr(major_chr.find_last_of("_")+1)
	   << " " <<  std::get<0>(block_length[ii])  << " " << std::get<1>(block_length[ii]) << " "
	   << n << " " <<  std::get<0>(block_length_inref[ii])  << " " <<  std::get<1>(block_length_inref[ii])
	   << endl;
	
	//bedfile << ctg << " " << std::get<0>(block_bp[ii]) << " " << std::get<1>(block_bp[ii])
	bedfile << ctg << " " << std::get<0>(block_length[ii]) << " " << std::get<1>(block_length[ii])
		<< " maps_to_" << n << " 150 " << strand_per_block[ii] <<  endl;

      }
    }
    ii++;
  }
  if(verbose)cout << tt.str() <<  endl ;

  myfile << ss.str();
  return 0;
}

// ******************************************** //
int LocateMisJ_plot(std::vector<nALIGNMENTS> ctgals, std::pair<string, int> major_al)
// ******************************************** //
{
  // sort alignment by position in contig
  std::stable_sort(ctgals.begin(), ctgals.end(), [](const nALIGNMENTS &a, const nALIGNMENTS &b){
      return ( std::min(a.ctgi,a.ctgf) < std::min(b.ctgi,b.ctgf) );
    });

  string ctg=ctgals[0].ctg;
  vector<string> chrs;
  vector<string> als_per_block;
  vector<long int> bases_per_block;
  vector<std::pair<long int, long int> > block_length;
  vector<string> blocks;
  vector<std::pair<long int, long int> > block_bp; //block breaking points


  long int max_al = major_al.second; 
  string major_chr= major_al.first ; 

    
  long int previous_pos;
  long int last_pos;

  std::stringstream tt;
  tt.imbue(std::locale(""));
  std::stringstream ss;

  tt << endl << " Synteny Group:  " << ctg << " major chr " <<  major_chr << endl;
  
  if(0) cout << " Synteny Group:  " << ctg << endl;
  


  // Start set to 0
  string previous_chr=ctgals[0].chr;
  long int previous_chri=ctgals[0].chri;
  long int previous_chrf=ctgals[0].chrf;
  int ntt=1; // counter for summing aligned bases to each separate block
  bases_per_block.push_back(std::abs(ctgals[0].ctgf-ctgals[0].ctgi));// start aligned bases for first block
  long int min_pos=std::min(ctgals[0].ctgf,ctgals[0].ctgi);
  long int max_pos=std::max(ctgals[0].ctgf,ctgals[0].ctgi);
  block_length.push_back(std::make_pair(min_pos,max_pos));
  string thisseq = ctg_seqs[ctgals[0].ctg].substr(0,min_pos);
  block_bp.push_back(std::make_pair(thisseq.find_last_of("NNN") + 1, 0));
  long int this_pos=0;

  long int ref_min=std::min(ctgals[0].chrf,ctgals[0].chri);
  long int ref_max=std::max(ctgals[0].chrf,ctgals[0].chri);
  vector<std::pair<long int, long int> > block_length_inref;
  block_length_inref.push_back(std::make_pair(ref_min,ref_max));

  vector<string> chrblocks;
  chrblocks.push_back(ctgals[0].chr);

  int ccount=0;
  std::for_each(ctgals.begin(),ctgals.end(), [&] (nALIGNMENTS const& a) {   // only works if sorted by ctg position:
      long int naligned=std::abs(a.ctgf-a.ctgi);
      this_pos=std::min(a.ctgi,a.ctgf);
      max_pos=std::max(a.ctgi,a.ctgf);
      
      ref_min=std::min(a.chrf,a.chri);
      ref_max=std::max(a.chrf,a.chri);
      als_per_block.push_back(a.chr);  
    
      if(0)
	cout << " Al:" 
	     << a.chr << " prev: " << previous_chr << " "
	     << a.ctgi << " "<< a.ctgf << " pchri " << previous_chri
	     << " " << a.chri<< " pchrf "<< previous_chrf<< " "<< a.chrf << " " 
	     << checknewblock(previous_chri,a.chri,previous_chrf,a.chrf)
	     << endl;      


      if(a.chr != previous_chr || ((a.chr == previous_chr) && (checknewblock(previous_chri,a.chri,previous_chrf,a.chrf)))){
	
	if(a.chr != previous_chr) nbreak+=1;

	bases_per_block.push_back(naligned); // start aligned bases for new block
	block_length.push_back(std::make_pair(this_pos,max_pos));  
	block_length_inref.push_back(std::make_pair(ref_min,ref_max));
	chrblocks.push_back(a.chr);

	if(0)
	  cout << "     New block " 
	       << a.chr << " prev: " << previous_chr << " "
	       << a.ctgi << " "<< a.ctgf << " pchri " << previous_chri
	       << " " << a.chri<< " pchrf "<< previous_chrf<< " "<< a.chrf << " " 
	       << checknewblock(previous_chri,a.chri,previous_chrf,a.chrf)
	       << endl;
	

	// get relevant part of sequence
	thisseq = ctg_seqs[a.ctg].substr(previous_pos,this_pos-previous_pos);
       	std::get<1>(block_bp[ntt-1])=thisseq.find_first_of("NNN") + previous_pos + 1;   // end of previous block at NNN  (-1, bc ntt starts from 1)
	ntt++;// new block
	
	// new block start
	block_bp.push_back(std::make_pair(thisseq.find_last_of("NNN") + previous_pos + 1, 0));

      }else{ 
	bases_per_block[ntt-1]+=naligned; // sum aligned bases to present block
	block_length[ntt-1]=std::make_pair(std::min(this_pos,std::get<0>(block_length[ntt-1])),std::max(max_pos,std::get<1>(block_length[ntt-1])));  
	block_length_inref[ntt-1]=std::make_pair(std::min(ref_min,std::get<0>(block_length_inref[ntt-1])),std::max(ref_max,std::get<1>(block_length_inref[ntt-1])));  

      }
        
      last_pos=this_pos;
      previous_chr=a.chr;
      previous_chri=a.chri;
      previous_chrf=a.chrf;
      previous_pos=std::max(a.ctgi,a.ctgf);
 
      ccount++;
    });
  // last end: 
  thisseq =ctg_seqs[ctgals[ctgals.size()-1].ctg].substr(previous_pos,this_pos-previous_pos);
  std::get<1>(block_bp[ntt-1])=thisseq.find_first_of("NNN") + previous_pos + 1; // end of previous block at NNN  (-1, bc ntt starts from 1)

  // unique works only if elements are sorted (same element close to each other,
  // as ctgals is sorted by ctg position, unique gives the blocks !
  //blocks=als_per_block;
  //blocks.erase( unique(blocks.begin(), blocks.end() ), blocks.end() );
  
  blocks=chrblocks;
  // to get only unique chrs, need to sort first:
  chrs=blocks;
  sort(chrs.begin(),chrs.end());
  chrs.erase( unique( chrs.begin(), chrs.end() ), chrs.end() );

  vector<int> block_len;
  for(string b : blocks) {
    block_len.push_back(std::count(als_per_block.begin(), als_per_block.end(), b));
  } 



  if(major_chr != "chrY") tt << "   This sequence aligns mostly to " << chrs.size() << " chromosomes "
       << " and it is divided in about " << blocks.size() 
       << " blocks, ordered as follows: " << endl;

  int ii=0;
  for(string n : blocks) { 

    if(ii>0 && n != "chrY" && major_chr != "chrY") tt << endl;
    if(n != "chrY" && major_chr != "chrY")
      tt << "\tBlock " << ii+1 << " maps to " << n 
	 << ":\t" << bases_per_block[ii] << " bp mapped,\tfirst pb: "
	 << std::get<0>(block_length[ii]) << ", last bp: " 
	 << std::get<1>(block_length[ii])  << " \tcovers " 
	 << std::get<1>(block_length[ii])-std::get<0>(block_length[ii]) 
	 << " bp\t\tsub-scaffold edges : " << std::get<0>(block_bp[ii]) << " - " << std::get<1>(block_bp[ii])
	 << " in ref: " << std::get<0>(block_length_inref[ii])  << " " <<  std::get<1>(block_length_inref[ii]);

    if(n != major_chr && 
       n != "chrY" && major_chr != "chrY"){
      ss << ctg   //<< "synteny-group_" << major_chr.substr(major_chr.find_last_of("_")+1)
	 << " " <<  std::get<0>(block_length[ii])  << " " << std::get<1>(block_length[ii]) << " "
	 << n << " " <<  std::get<0>(block_length_inref[ii])  << " " <<  std::get<1>(block_length_inref[ii])
	 << endl;

    }
    ii++;
  }
  
  if(verbose)cout << tt.str() <<  endl ;
  myfile << ss.str();
  return 0;
}
int checknewblock(long int pchri,long int nchri, long int pchrf, long int nchrf)
{
  int nblock=1;
  
  if( abs(pchri-nchri) < newblock ||
      abs(pchrf-nchrf) < newblock ||
      abs(pchri-nchrf) < newblock ||
      abs(pchrf-nchri) < newblock ) 
    nblock=0; // same block

  return nblock;
}
