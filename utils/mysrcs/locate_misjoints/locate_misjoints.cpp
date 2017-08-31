#include "../myinc/macro.h"  //misf.h"
#include <ctime>
#include <locale> // for imbue
#include <typeinfo>


static std::map<string, string> ctg_seqs;    // one el per ctg: 0 or 1 if it is not (0) or is (1) mapped READY
 

struct nALIGNMENTS {
  std::string chr,ctg,strand;
  long int chri,chrf,ctgf,ctgi;
  int albases;
  float id;
} nalignments;

static int ReadAls(char* file);
static int LocateMisJ(std::vector<nALIGNMENTS> ctgals);
static int AnalyseCtg(std::vector<nALIGNMENTS> ctgals);
static int nbreak=0;


// possible misjoint selection:
static  float min_len_perc = 0.25; // perc min_lenght to consider as possible misjoint wrt major al (reduce noise and small repeats)
static  long int min_len = 500000; // min_lenght to consider as possible misjoint
static  long int min_len_for_max=100000; // min_lenght for a major al


int main(int argc, char **argv)
{
  
  if(argc < 2){
    printf("Usage: %s  <draft> <alignment_file>  <perc_min_len> <min_len> <min_len_max> \n", argv[1]);
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
  if (argc >= 4)  min_len_for_max = atof(argv[3]);
  if (argc >= 5)  min_len_perc = atof(argv[4]);
  if (argc >= 6)  min_len = atof(argv[5]);

  cout << endl<< " Selected filters: " << endl
       << "  Min length for longest alignment to a single chr: " << min_len_for_max << " bp" <<endl
       << "  Min length for shorter alignments to supplementary chrs: " 
       << min_len_perc*100 << "% of longest or " << min_len << " bp" << endl;


  int err=0; int saveinfo=1; int readseq=1;
 
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
  rlen.clear();

  
  //  ****************************** //
  //  ****** READ ALIGNMENTS ******* //
  //  ****************************** //
  myfile.open("misjoints_details.txt");
  ReadAls(argv[2]);
  myfile.close();   

  
  if(nbreak)
    cout << endl << "\n **************  REPORT SUMMARY  *****************" << endl 
      //<< endl 
	 << "  ****** " << nbreak << " possible breaking points found ******" 
      // << endl
	 << endl <<   " *************************************************\n  "<< endl;
  else
    cout << endl << " No breaking points found " << endl;

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
    
    ss >> ctg >> chr >> id >> albases  >> kk  >> kk >> ctgi >> ctgf >> chri >> chrf >> kk >> kk;
    if(0) cout <<  ctg << "\t" <<  chr << "\t" << id << "\t" <<  albases  <<  "\t" 
	       <<  ctgi << "\t" <<  ctgf << "\t" <<  chri << "\t" <<  chrf << endl;
    

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
     
    // check if new ctg:
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
  

  // ************* Select out ************* //
  // not considering too small alignments for locating misjoints
  if(max_al < min_len_for_max)  return 0;
  
  
  int major_als=0;
  for ( const auto &p : len_mappedtochr ){ 
    long int al_len=p.second; //length aligned
    
    if(0)cout << al_len << " " << max_al << " " << al_len*1./max_al << " " << min_len_perc <<endl;
    if (al_len*1./max_al >= min_len_perc || al_len > min_len){
      major_als++;
      string chrname=p.first;
      if(chrname == major_chr ) if(0)cout << " Major chr " << endl;
      else if(0)cout  << endl;
      thisctglinks[chrname]=al_len;
      linkedchr.push_back(chrname);
    }
  }

  // ************* Select out ************* //
  // only a major mapping: no inter-chr misjoint
  if(major_als == 1 )  return 0;


  // Remove from ctgals all alignments to chrs others than majors
  ctgals.erase(std::remove_if(ctgals.begin(), ctgals.end(), [&](nALIGNMENTS const& x){
	return !(std::find(linkedchr.begin(), linkedchr.end(), x.chr) != linkedchr.end()); }), ctgals.end());


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
  vector<string> chrs;
  vector<string> als_per_block;
  vector<long int> bases_per_block;
  vector<std::pair<long int, long int> > block_length;
  vector<string> blocks;
  vector<std::pair<long int, long int> > block_bp; //block breaking points

     
  long int previous_pos;
  long int last_pos;

  std::stringstream tt;
  tt.imbue(std::locale(""));

  tt << endl << " Synteny Group:  " << ctg << endl;
  
  // Start set to 0
  string previous_chr=ctgals[0].chr;
  int ntt=1; // counter for summing aligned bases to each separate block
  bases_per_block.push_back(std::abs(ctgals[0].ctgf-ctgals[0].ctgi));// start aligned bases for first block
  long int min_pos=std::min(ctgals[0].ctgf,ctgals[0].ctgi);
  long int max_pos=std::max(ctgals[0].ctgf,ctgals[0].ctgi);
  block_length.push_back(std::make_pair(min_pos,max_pos));
  string thisseq = ctg_seqs[ctgals[0].ctg].substr(0,min_pos);
  block_bp.push_back(std::make_pair(thisseq.find_last_of("NNN") + 1, 0));
  long int this_pos=0;

 
  std::for_each(ctgals.begin(),ctgals.end(), [&] (nALIGNMENTS const& a) {   // only works if sorted by ctg position:
      long int naligned=std::abs(a.ctgf-a.ctgi);
      this_pos=std::min(a.ctgi,a.ctgf);
      max_pos=std::max(a.ctgi,a.ctgf);
      
      als_per_block.push_back(a.chr);  
    
      
      if(a.chr != previous_chr){
	bases_per_block.push_back(naligned); // start aligned bases for new block
	block_length.push_back(std::make_pair(this_pos,max_pos));  
		
	// get relevant part of sequence
	thisseq = ctg_seqs[a.ctg].substr(previous_pos,this_pos-previous_pos);
       	std::get<1>(block_bp[ntt-1])=thisseq.find_first_of("NNN") + previous_pos + 1;   // end of previous block at NNN  (-1, bc ntt starts from 1)
	ntt++;// new block
	
	// new block start
	block_bp.push_back(std::make_pair(thisseq.find_last_of("NNN") + previous_pos + 1, 0));


	if(0 && ctg == "fAnaTes1_74"){
	  long int t1=thisseq.find_first_of("NNN");
	  long int t2= thisseq.find_last_of("NNN");
	  cout << t1 << " " << t2 << endl;
	  cout << thisseq.size() << endl;
	  cout << ctg_seqs[a.ctg].substr(thisseq.find_first_of("NNN")-5, 10) << endl;
	}	

      }else{ 
	bases_per_block[ntt-1]+=naligned; // sum aligned bases to present block
	block_length[ntt-1]=std::make_pair(std::min(this_pos,std::get<0>(block_length[ntt-1])),std::max(max_pos,std::get<1>(block_length[ntt-1])));  
	if(0)cout << a.chr << " "<< this_pos << " " << std::get<0>(block_length[ntt-1]) << " " << std::min(this_pos,std::get<0>(block_length[ntt-1])) << endl;
      }
        
      last_pos=this_pos;
      previous_chr=a.chr;
      previous_pos=std::max(a.ctgi,a.ctgf);
 
    });
  // last end: 
  thisseq =ctg_seqs[ctgals[ctgals.size()-1].ctg].substr(previous_pos,this_pos-previous_pos);
  std::get<1>(block_bp[ntt-1])=thisseq.find_first_of("NNN") + previous_pos + 1; // end of previous block at NNN  (-1, bc ntt starts from 1)

  // unique works only if elements are sorted (same element close to each other,
  // as ctgals is sorted by ctg position, unique gives the blocks !
  blocks=als_per_block;
  blocks.erase( unique(blocks.begin(), blocks.end() ), blocks.end() );
  // to get only unique chrs, need to sort first:
  chrs=blocks;
  sort(chrs.begin(),chrs.end());
  chrs.erase( unique( chrs.begin(), chrs.end() ), chrs.end() );

  vector<int> block_len;
  for(string b : blocks) {
    block_len.push_back(std::count(als_per_block.begin(), als_per_block.end(), b));
  } 

  nbreak+=blocks.size()-1;

  tt << "   This sequence aligns mostly to " << chrs.size() << " chromosomes "
       << " and it is divided in about " << blocks.size() 
       << " blocks, ordered as follows: " << endl;


  int ii=0;
  for(string n : blocks) { 
    if(ii>0) tt << endl;
    tt << "\tBlock " << ii+1 << " maps to " << n 
       << ":\t" << bases_per_block[ii] << " bp mapped,\tfirst pb: "
       << std::get<0>(block_length[ii]) << ", last bp: " 
       << std::get<1>(block_length[ii])  << " \tcovers " 
       << std::get<1>(block_length[ii])-std::get<0>(block_length[ii]) 
       << " bp\t\tsub-scaffold edges : " << std::get<0>(block_bp[ii]) << " - " << std::get<1>(block_bp[ii]);
    ii++;
  }

  cout << tt.str() <<  endl ;
  myfile << tt.str() << endl;


  return 0;
}
