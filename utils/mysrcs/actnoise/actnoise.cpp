#include "../myinc/macro.h"

//const int maxals=10000;
static float noise=10; //0.2 for 20% of major
static int minlength=5000; //30K tried for devil // some noise has many als, but all concentrated in a small ref region  
                   // not used anymore
static int minnoise=3000; // no used
//static int maxnoise=10000;  
static int newblock=100000; 
static  float minid=80;  // min id 


static int longestchr=0;
static int yes=1;
static int nchr=0;
static int nctg=0;
static string myname1;

static vector<string> majchr(nctg,""); // major chr al per each ctg
static vector<int> goodctg(nctg,0); // major chr per each ctg
static vector<int> longestmin(nctg,0);
static vector< vector<int> > gblocks(nctg, vector<int>(1,-1));   
static vector< vector<int> > gindex(nctg, vector<int>(1,-1));   

static  vector<int> okblock;
static  vector<int> minblock; 
static  vector<int> maxblock; 
static  vector<int> alblock;
static  vector<int> majorc;    
static string aligner = "smalt";
static  std::map<string, int> ctgsizes;

struct MYALS
{
  vector<string> ctg, chr;
  vector<int> ctgi, ctgf, chri,chrf,len;
  vector<string> more1,more2,more3,more4,more5;
  vector<int> block;
}als;
static vector<MYALS>  myals;
static MYALS empty(string ctg);
static MYALS fillall(string ctg);
static MYALS fillals(vector<string> str, vector<long int> in, MYALS thisals);
struct MYCTGS
{
  string name;
  int position;
  int good;
  string major;
  int chrpos;

  bool operator() (MYCTGS i,  MYCTGS j) { 
    return( (i.chrpos < j.chrpos) || 
            (i.chrpos == j.chrpos && i.position < j.position )  
            ); 
  }
}ctgs;
static vector<MYCTGS>  mycontigs;


int readals(char* file);
int orderals();
int sortnwritectgs();
int checknewblock(int pchri, int nchri,int pctgi, int nctgi);


int main(int argc, char *argv[])
{ 
  int thispri=0;
  pri=0;

  if(noise==0) newblock=longestchr;

  if (argc < 4) {
   fprintf(stderr, "Remove short blocks of isolated alignemnts (noise)\nUsage: %s <reference.fasta>  <draft.fasta>   <alignment>  <noise_level> <minid> <aligner> \n", argv[1]); 
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


  // read reference and draft assembly
  string reffile = argv[1];
  string seqfile = argv[2];
  string alfile = argv[3];
  if (argc >= 5)  noise= to_int(argv[4])*1./10;
  minlength=10000; //maxnoise; // not used anymore, now using ctg length  !!! still using???? !!!
  int tempid=0;
  if (argc >= 6)  tempid= to_int(argv[5]);
  minid=tempid*1.;
  
  if(argc>6) aligner = argv[6];

  string newname="nonoise"+ to_string(noise)+"_minid"+to_string(minid)+"_";
  string myname=myrename(seqfile,newname);
  myname1=myrename(alfile,newname);

  //reference: find longest chr
  int err=0;
  int saveinfo=1;
  int readseq=0;
  err=readfasta(argv[1],saveinfo);
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }
  longestchr= *std::max_element(rlen.begin(),rlen.end());
  nchr=rname.size();
  for(int i=0; i<rname.size(); i++){
    string name=rname[i];
    refmap[name] = i;
  }
  refchrs=rname;
  rname.clear();
  rlen.clear();

  // draft
  err=0;
  saveinfo=1;
  readseq=1;
  err=readfasta(argv[2],saveinfo,"",readseq);
  if(err){
    cout << " Sorry, something went wrong..." << endl;
    return 1;
  }
  for(int i=0; i<rname.size(); i++){
    string name=rname[i];
    seqmap[name] = i;
    ctgsizes[name] = rlen[i];
  }
  seqctgs=rname;
  nctg=rname.size();
  rname.clear();
  rlen.clear();
  if(no)cout << " draft done " << endl;

  majchr.resize(nctg,"");
  goodctg.resize(nctg,0);
  longestmin.resize(nctg,longestchr);
  gblocks.resize(nctg, vector<int>(1,-1));
  gindex.resize(nctg, vector<int>(1,-1));


  // Read alignments
  readals(argv[3]);
  if(no)cout << " als read " << endl;

  //Cut noise and align:
  orderals();
  if(no)cout << " als filtered " << endl;

  // sort ctg positions and write als
  myfile.open(myname1.c_str());
  sortnwritectgs();
  myfile.close();
  if(no)cout << " contigs sorted, als written " << endl;
  
  
  // write new draft fasta with aligned and ordered ctgs 
  //get seqmap: order of ctgs in seq-fasta
  //write new seq-fasta
  //cout << " now writing fasta " << myname << " " << seqctgs.size() << endl;
  myfile.open(myname.c_str());
  char fa[5]={">"};
  for(int c=0; c<seqctgs.size(); c++){
    string newo=newmap[c];
    int newi=seqmap[newo];
    // only write ctgs that are in the alignment: (take out if filtered by noise)
    if(!newo.empty())myfile <<  fa[0] << rname[newi] << endl << rseq[newi]<<endl;
    if(!newo.empty())if(pri)cout <<  fa[0] << rname[newi] << endl << rseq[newi]<<endl;
  }
  rseq.clear();
  rname.clear();
  rlen.clear();

  myfile.close();

  return 0;
}



int readals(char* file){
  std::ifstream infile(file);
  string line;
  string prevctg="";
  int alnum=0;
  MYALS thisals;
  string ctg, chr;

  // create vector of struct with all the ctgs > new function?
  for (int ii=0; ii<seqctgs.size(); ii++){
    MYALS emptystr=empty(seqctgs[ii]);
    myals.push_back(emptystr);
  }

  vector<vector<int> > totals(nctg,vector<int>(nchr,0));
  int ccount=0;
  int tt1=0;
  int tt2=0;
  while(getline(infile,line)){  // loop on alignments
    std::stringstream ss(line);
    int ctgi, ctgf, chri,chrf;
    vector<string> str(7);
    vector<long int> pos(6);
    
    int mapscore;
    string mycigar;

    ss >> str[0] >> str[1]  >> str[2]  >> str[3]   >> str[4]  >> str[5]  
       >> pos[0] >> pos[1] >> pos[2] >> pos[3] >> str[6] >> pos[4] >> mapscore; 
		    
		    
    ctg=str[0];
    chr=str[1];
    
    if(ctg=="fAnaTes1_67" || prevctg=="fAnaTes1_67") pri=0;
    else pri=0;

    if(0) cout <<  str[0] << " " <<  str[1]  << " " <<  str[2]  << " " <<  str[3]   << " " <<  str[4]  << " " <<  str[5]  
	       << " " <<  pos[0] << " " <<  pos[1] << " " <<  pos[2] << " " <<  pos[3] << " " <<  str[6] << " " <<  pos[4] << " " << mapscore << endl;
     if(0) cout << " " << mycigar << endl;
    if(pri) tt1++;

    // Remove alignment with ID < min-ID
    float id;  
    id=to_float(str[2])*1./100;  // id perc from smalt
    
    
    if(id<minid)
      continue;
    
    if(aligner == "smalt") {
      //cout << mapscore << endl;
      if(mapscore<10)
	continue;
    }
    /*else{
      if(mapscore<20)   
	continue;
	}*/
    if(0) cout <<  "   after cuts " << str[0] << " " <<  str[1]  << " " <<  str[2]  << " " <<  str[3]   << " " <<  str[4]  << " " <<  str[5]  
		 << " " <<  pos[0] << " " <<  pos[1] << " " <<  pos[2] << " " <<  pos[3] << " " <<  str[6] << " " <<  pos[4] << endl;
    if(pri) tt2++;

    if(str[0] == prevctg || alnum==0){ // First alignment or same ctg
      thisals=fillals(str,pos,thisals);
      totals[seqmap[ctg]][refmap[chr]]+=(abs(pos[1]-pos[0])); //keep track of how much aligned bases per each chr
    if(pri)ccount++;

    }else{ // New ctg
      
      //Saving info on prev ctg:
      myals[seqmap[prevctg]]=thisals;
      
      if(pri && prevctg!="fAnaTes1_66" && prevctg!="fAnaTes1_10" )
	cout << "  inside loop " << prevctg << " als: " 
	     << myals[seqmap[prevctg]].ctg.size() << " " << ccount<< " " 
	     << "  new ctg is " << ctg << endl;

      // Find to which chr there are more bases aligned ? should this be done after cutting noise??
      vector<int>::iterator max=max_element(totals[seqmap[prevctg]].begin(),totals[seqmap[prevctg]].end());
      int imax= *max_element(totals[seqmap[prevctg]].begin(),totals[seqmap[prevctg]].end());
      int ii = std::distance( totals[seqmap[prevctg]].begin(), std::find(totals[seqmap[prevctg]].begin(), totals[seqmap[prevctg]].end(), imax ) );
      majchr[seqmap[prevctg]]=refchrs[ii];
      if(0) cout << prevctg <<  " Major al against " << majchr[seqmap[prevctg]] << " numb als " << alnum << endl ;

      thisals=MYALS();
      if(myals[seqmap[ctg]].chr[0]!=""){  // if this ctg has already some alignments in myals, add to them (accounts for als file not ordered by ctg name)
	thisals=myals[seqmap[ctg]];
      }
      thisals=fillals(str,pos,thisals);
      totals[seqmap[ctg]][refmap[chr]]+=(abs(pos[1]-pos[0])); //???
    }
    
    prevctg=str[0];
    alnum++;
  } // end reading al file loop

  // catch last alignment:
  myals[seqmap[prevctg]]=thisals;

  // Find to which chr there are more bases aligned ? should this be done after cutting noise??
  vector<int>::iterator max=max_element(totals[seqmap[prevctg]].begin(),totals[seqmap[prevctg]].end());
  int imax= *max_element(totals[seqmap[prevctg]].begin(),totals[seqmap[prevctg]].end());
  int ii = std::distance( totals[seqmap[prevctg]].begin(), std::find(totals[seqmap[prevctg]].begin(), totals[seqmap[prevctg]].end(), imax ) );
  majchr[seqmap[prevctg]]=refchrs[ii];
  if(0)cout << " " << prevctg <<  " Major al against " <<majchr[seqmap[prevctg]]  << " numb als " << alnum << endl << endl;
  


  if(pri)cout << " final numbers: " << tt1 << " " << tt2 << endl;
  if(pri)cout << prevctg << " als: " << myals[seqmap["fAnaTes1_67"]].ctg.size() << " " << ccount<< endl;
  

  pri=0;
  return 0;
}


int orderals()
{

  MYALS thisals;
  MYCTGS thisctg;
  string ctg, chr;
  int thispri=0;
  
  int thisctglenght=0;
  int thisnoise=0;

  // loop over ctgs:
  for (int ictg=0; ictg<myals.size(); ictg++){ // loop over ctgs
     
     if(myals[ictg].chr[0]=="") continue;
     ctg=myals[ictg].ctg[0];
     thisctglenght=ctgsizes[ctg];
     int percnoise=(int)(noise*thisctglenght/100);  // ex. noise=20
     thisnoise=percnoise;   //std::min(percnoise, maxnoise);

     if(ctg == "fAnaTes1_16" || ctg == "fAnaTes1_15") 
       thispri=0;
     else 
       thispri=0;
     if(thispri)cout << "  " << ctg << " " << thisctglenght << " " << noise << " " << thisnoise << endl;
         

     int previ=-1;
     int prevctgi=-1;
     int prevf=-1;
     int block=0;

     vector<int> tempmin;
     vector<int> tempmax;
     int alperblock=0;
     string prevchr="";

     okblock.resize(0);
     //minblockperchr.resize(0,vector<int>(0));
     minblock.resize(0);
     maxblock.resize(0);
     alblock.resize(0);
     majorc.resize(0);
     myals[ictg].block.resize(myals[ictg].ctg.size(),-1);

     // sort aligns for this ctg according to chr and position in chr
     vector<std::tuple<string,int,int> > myindex;
     for(int in=0; in<myals[ictg].ctg.size(); in++)
       myindex.push_back (make_tuple(myals[ictg].chr[in],myals[ictg].chri[in],in)); 
     sort(myindex.begin(), myindex.end(),
	  [](std::tuple<string,int,int>& p1, std::tuple<string,int,int>& p2) {
	    return std::tie(std::get<0>(p1),std::get<1>(p1))  < std::tie(std::get<0>(p2),std::get<1>(p2));
	  });    
     

     if(thispri) cout << endl << " Contig " << ctg << " als: " << myals[ictg].ctg.size() << endl;
     
     int countgb=0;  // count good blocks for this ctg
     int countmchr=0;  // count good block in the major chr for this ctg
     gindex[ictg].resize(myals[ictg].ctg.size());
     int maj=0;
     
     for (int j=0; j<myals[ictg].ctg.size(); j++) // loop over als for this ctg
       {
	 int i=std::get<2>(myindex[j]);
	

	 /*if(ctg == "fAnaTes1_67" && block < 80 && myals[ictg].chr[i]=="seabass_18") 
	   thispri=1;
	 else 
	 thispri=0;*/

	 // these to use no sorted chr list:
	 //i=j;
	 //previn=i-1;

	 gindex[ictg][j]=i;
	 
	 int thisi=myals[ictg].chri[i];
	 int thisf=myals[ictg].chrf[i];      
	 int thisctgi=myals[ictg].ctgi[i];
	 int thisctgf=myals[ictg].ctgf[i];      
	 string thischr=myals[ictg].chr[i];

	 if(0) cout << " is same block? " <<  i << " "<< previ<< " "
			  <<thisi<< " "<< prevctgi<< " "<<thisctgi 
			  << " calc: " <<  abs(previ-thisi) << " "<< abs(prevctgi-thisctgi) 
			  << " "<< newblock << " "
			  << checknewblock(previ,thisi, prevctgi,thisctgi) << endl;

	 if( (previ==-1 && prevf==-1) ||  // first al 
	     (thischr == prevchr 
	      && checknewblock(previ,thisi, prevctgi,thisctgi))){   // define a block from chrs positions
	       
	   // same block
	   alperblock+=abs(myals[ictg].ctgi[i]-myals[ictg].ctgf[i]); 
	   tempmin.push_back(std::min(thisi,thisf));
	   tempmax.push_back(std::max(thisi,thisf));
	   myals[ictg].block[i]=block; // belongs to block #block
	   
	   if(previ==-1 && prevf==-1){ // check if major only at the beginning
	       if(myals[ictg].chr[i] == majchr[ictg]){ //major al
		 maj=1;
		 //if(thispri) cout << " major chr " << endl;
	       }else{
		 maj=0;
		 //if(thispri) cout << " NOT major chr " << endl;
	       }
	     }

	 }else{ 
	   // new block
	   
	   // check previous block
	   int good=0; 
	   
	  
	   if(maj){
	     alblock.push_back(alperblock);
	     minblock.push_back(*min_element(tempmin.begin(),tempmin.end()));
	     maxblock.push_back(*max_element(tempmax.begin(),tempmax.end()));
	   }else{
	     alblock.push_back(0);
	     minblock.push_back(longestchr);
	     maxblock.push_back(0);
	   }

	   int blocklength=(*max_element(tempmax.begin(),tempmax.end()))-(*min_element(tempmin.begin(),tempmin.end()));  // in reference
	   
	   if(alperblock > thisnoise && blocklength >= minlength) // thisnoise) //minlength)  // good block   
	     good++;

	   countgb+=good;
	   goodctg[ictg]+=good;
	   countmchr+=maj;
	   okblock.push_back(good);
	   majorc.push_back(maj);



	   if(thispri){
	     if(good)cout << "  Block  " << block  << " is good: ";
	     else cout << "  Block  " << block  << " is bad: ";
	     if(maj==1) cout << " this is Major chr,   als=";
	     else  cout << " this is NOT Major chr,   als=";
	     cout << alperblock << " for chr " << prevchr <<" " << okblock[block] 
		  << " min is " << (*min_element(tempmin.begin(),tempmin.end()))
		  << " max is " << (*max_element(tempmax.begin(),tempmax.end()))
		  << " lenght is " << blocklength
	       //<<  endl;
	       //cout  
		  << "      alblock size " <<  alblock.size() << endl;
	   }
	   
	   //reset 
	   tempmin.resize(0);
	   tempmax.resize(0);
	   maj=0;
	   
	   // Start new block
	   alperblock=abs(myals[ictg].ctgi[i]-myals[ictg].ctgf[i]); 
	   tempmin.push_back(std::min(thisi,thisf));
	   tempmax.push_back(std::max(thisi,thisf));
	   block++;
	   myals[ictg].block[i]=block; // belongs to block #block
	   if(myals[ictg].chr[i] == majchr[ictg]){ //major al
	     maj=1;
	     //if(thispri) cout << " major chr " << endl;
	   }else{
	     maj=0;
	     // if(thispri) cout << " NOT major chr " << endl;
	   }


	 } // end new block
	 
	 if(no) cout << " Al " << i << " " << ctg << " " << thischr << " " << thisi << " " 
			  << thisf<< " min is " << std::min(thisi,thisf)
			  << " " << thisctgi << " " << thisctgf
			  << endl;
	 
	 
	 prevctgi=thisctgi;
	 previ=thisi;
	 prevf=thisf;
	 prevchr=thischr;
	 
       }// end loop over als for this ctg
	 

     // save last block:

     int blocklength=(*max_element(tempmax.begin(),tempmax.end()))-(*min_element(tempmin.begin(),tempmin.end()));  // in reference
    if(maj){
       alblock.push_back(alperblock);
       minblock.push_back(*min_element(tempmin.begin(),tempmin.end()));
       maxblock.push_back(*max_element(tempmax.begin(),tempmax.end()));
       //minblockperchr[refmap[prevchr]].push_back(*min_element(tempmin.begin(),tempmin.end()));
     }else{
       alblock.push_back(0);
       minblock.push_back(longestchr);
       maxblock.push_back(0);
       //minblockperchr[refmap[prevchr]].push_back(longestchr);
     }
     int good=0;
     if(alperblock > thisnoise && blocklength >= minlength) // thisnoise) //minlength)  // good block   
       good++;
     countgb+=good;
     goodctg[ictg]+=good;
     countmchr+=maj;
     okblock.push_back(good);
     majorc.push_back(maj);

    	 
     if(thispri){
       if(good)cout << "  Block  " << block  << " is good: als= ";
       else cout << "  Block  " << block  << " is bad : als=";
       cout << alperblock << " for chr " << prevchr <<" " << okblock[block] 
	    << " min is " << (*min_element(tempmin.begin(),tempmin.end()))
	    << " max is " << (*max_element(tempmax.begin(),tempmax.end()))
	    << " lenght is " << blocklength
	    <<  endl;
       cout << "   alblock size " <<  alblock.size() << " tempmin size " << tempmin.size()<< endl;
     }
    
     if(thispri)cout << endl << " This ctg had " << countgb << " good blocks and " 
	       << countmchr << " good block in the major chr " << endl;
      
     gblocks[ictg].resize(okblock.size());
     gblocks[ictg]=okblock;
            

     if(goodctg[ictg]<1) continue;  // .. fill vector to sort only if there is at least a single good block

     // Choose major al for this ctg in major chr and min position:
     std::vector<int>::iterator lb=max_element(alblock.begin(),alblock.end());
     int pos=std::distance(std::begin(alblock), lb);
     longestmin[ictg]=minblock[pos];

     if(thispri)cout << " Contig " << ctg << " " << ictg << " major al is against " <<  majchr[seqmap[ctg]]
		 << " longest block is " << *lb 
		 << " global min is " << longestmin[ictg] 
		 << endl;
     
     thisctg.name=ctg;
     thisctg.good=goodctg[ictg];
     thisctg.position=longestmin[ictg];
     thisctg.major=majchr[ictg];
     thisctg.chrpos=refmap[majchr[ictg]];
     mycontigs.push_back(thisctg);
     
   } // end loop over ctgs

   return 0;
}

int sortnwritectgs()
{  
  MYALS thisals;
  int thispri=0;
  
  // Sort and write
  std::sort(mycontigs.begin(), mycontigs.end(), ctgs);  // sort by major chr and min position 
  for(int i=0; i<mycontigs.size(); i++){ // new ctg position for draft assembly
    string name=mycontigs[i].name;
    newmap[i] = name;
    // contig order:
    if(no) cout << i << " " << name << endl;
  }

  for (int ictg=0; ictg<myals.size(); ictg++) {   // loop over ctgs
    if(myals[ictg].chr[0]=="") continue;   // if no alignment 
    thisals= myals[ictg];
    vector<int> gblock=gblocks[ictg];
   
    if(thispri)cout << "Sorting and writing:  ctg " 
		    << myals[ictg].ctg[0] << endl;
    vector<int> myindex=gindex[ictg];


    int nn=0;
    int prevblock=-1;
    for(int j=0; j<thisals.ctg.size();j++) {// loop over al per ctg
      nn++;
      int t=(myindex[j]);
      
      int block=thisals.block[t];
      if(block < 0) cout << " Error, block =-1 " << endl;

      if(thisals.chr[t]=="seabass_0")
	thispri=0;
      else 
	thispri=0;
      
   
      if(thispri && prevblock != block){
	cout << " " << j << " " << t << " " << (myindex[j-1])  << " block " << block << " good " << gblock[block] << endl;
	cout << " ctg " << ictg << " " << thisals.chr[t] << " al " << t << " block " << block
	     << " good? " << gblock[block] << " " <<   thisals.chri[t] << " " <<  thisals.chrf[t] << endl;
      }
      if(gblock[block]==0) continue; 
	   
      myfile <<  thisals.ctg[t] << "\t" <<  thisals.chr[t] << "\t" <<  thisals.more1[t] << "\t" <<  thisals.more2[t]  <<  "\t" 
	     <<  thisals.more3[t]  << "\t" <<  thisals.more4[t] << "\t"
	     <<  thisals.ctgi[t] << "\t" <<  thisals.ctgf[t] << "\t" <<  thisals.chri[t] << "\t" <<  thisals.chrf[t]
	     << "\t" <<  thisals.more5[t] << "\t" <<  thisals.len[t] << endl;
      if(thispri)cout <<  "output " << thisals.ctg[t] << "\t" <<  thisals.chr[t] << "\t" <<  thisals.more1[t] << "\t" <<  thisals.more2[t]  <<  "\t" 
		      <<  thisals.more3[t]  << "\t" <<  thisals.more4[t] << "\t"
		      <<  thisals.ctgi[t] << "\t" <<  thisals.ctgf[t] << "\t" <<  thisals.chri[t] << "\t" <<  thisals.chrf[t]
		      << "\t" <<  thisals.more5[t] << "\t" <<  thisals.len[t] << endl;
    
      prevblock=block;
    }
       
  }
 
  return 0;
}


 int checknewblock(int pchri, int nchri,int pctgi, int nctgi)
{
  int nblock=0;
  
  if( abs(pchri-nchri) < newblock &&
      abs(pctgi-nctgi) < newblock ) nblock=1; // same block

 //(thischr == prevchr && (abs(thisi-previ)<newblock && abs(thisi-prevf)<newblock &&   // or same chr + closer than newblock
	     //  abs(thisf-previ)<newblock && abs(thisf-prevf)<newblock))){

  return nblock;
}

MYALS fillals(vector<string> str, vector<long int> in, MYALS thisals){
  thisals.ctg.push_back(str[0]);
  thisals.chr.push_back(str[1]);
  thisals.more1.push_back(str[2]);
  thisals.more2.push_back(str[3]);
  thisals.more3.push_back(str[4]);
  thisals.more4.push_back(str[5]);
  thisals.more5.push_back(str[6]); 
  thisals.ctgi.push_back(in[0]);
  thisals.ctgf.push_back(in[1]);
  thisals.chri.push_back(in[2]);
  thisals.chrf.push_back(in[3]);
  thisals.len.push_back(in[4]);
  return(thisals);
}

MYALS empty(string ctg){
  MYALS thisals;
  thisals.ctg.push_back(ctg);
  thisals.chr.push_back("");
  thisals.more1.push_back("");
  thisals.more2.push_back("");
  thisals.more3.push_back("");
  thisals.more4.push_back("");
  thisals.more5.push_back(""); 
  thisals.ctgi.push_back(-1);
  thisals.ctgf.push_back(-1);
  thisals.chri.push_back(-1);
  thisals.chrf.push_back(-1);
  thisals.len.push_back(-1);
  thisals.block.push_back(-1);

  return(thisals);
}
