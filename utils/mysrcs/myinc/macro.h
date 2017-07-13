#include <vector>
#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>  //setprecision
#include <algorithm>    // sort, reverse
#include <iomanip>  //setprecision
#include <gzstream.h>
#include <tuple> // C++11, for std::tie
#include <numeric>  // for accumulate
#include <iterator>     // std::iterator, std::input_iterator_tag

using std::cout;
using std::endl;
using std::vector;
using std::string;
static int no=0;
static int pri=0;

static gzFile fp;
static  vector<int> rlen;
static  vector<string> rseq;
static  vector<string> rqual;
static  vector<string> rname;
static  vector<string> rcomment;
static std::ofstream myfile;

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
static  std::map<int, string> newmap;
static vector<string> refchrs;
static vector<string> seqctgs;



// ---------------------------------------- //
string comple(string seq)
// ---------------------------------------- //
{
  // return string complement

  string myalph = "ACGTN";
  string low_myalph = "acgtn";
  string cseq = "";

  if(0)cout << " string size: " << seq.size() << endl;
  reverse(seq.begin(),seq.end());

  for (int i=0; i<seq.size(); i++){
    //if(i>270)cout << i << " " << seq[i] << " " << seq.size() << " " << cseq.size()<< endl;
   if(seq[i] == myalph[0] )
      cseq.append("T");
    else if (seq[i] == myalph[1] )
      cseq.append("G");
    else if (seq[i] == myalph[2] )
      cseq.append("C");
    else if (seq[i] == myalph[3] )
      cseq.append("A");
    else if (seq[i] == myalph[4] )
      cseq.append("N");
    else if(seq[i] == low_myalph[0] )
      cseq.append("t");
    else if (seq[i] == low_myalph[1] )
      cseq.append("g");
    else if (seq[i] == low_myalph[2] )
      cseq.append("c");
    else if (seq[i] == low_myalph[3] )
      cseq.append("a");
    else if (seq[i] == low_myalph[4] )
      cseq.append("n");
    
   /*    if(seq[i] == myalph[0] )
      cseq.insert(i,"T");
    else if (seq[i] == myalph[1] )
      cseq.insert(i,"G");
    else if (seq[i] == myalph[2] )
      cseq.insert(i,"C");
    else if (seq[i] == myalph[3] )
      cseq.insert(i,"A");
    else if (seq[i] == myalph[4] )
    cseq.insert(i,"N");*/
  }
  if(0)cout << " comple: " << cseq.size()  << endl;
  return cseq;
}


std::pair<string, string>  getname(string str){
  size_t ns=0;
  size_t nt=0;
  ns=str.find(" ");
  nt=str.find("\t");

  string name="";
  string comment="";

  if(ns!=std::string::npos) { 
    name=str.substr(1,ns-1);
    comment=str.substr(ns+1,str.size());
  }else if(nt!=std::string::npos) {
    name=str.substr(1,nt-1);
    comment=str.substr(nt+1,str.size());
  }else{
    name=str.erase(0,1);
  }

  return std::make_pair(name,comment);
}

string myrename(string myname, string addthis="", string exte=""){

  size_t pos = 0;
  string token;
  string delimiter = "/";

  while ((pos = myname.find(delimiter)) != std::string::npos) {
    token = myname.substr(0, pos);
    myname.erase(0, pos + delimiter.length());
  }

  pos=0;
  myname= addthis + myname;

  if(exte.size()){  //change file extension
    token=myname;
    pos   = token.find_last_of(".");
    if (pos != std::string::npos) {     
        token.erase(pos,myname.size());
   }
   myname=token+exte;
  }

  return(myname);
}


// ---------------------------------------- //
int fasttype(char* file)
// ---------------------------------------- //
{ 
  char fq[5]={"@"};
  char fa[5]={">"};
  string ttname;
  int isfq=0;
  igzstream infile(file);

  getline(infile,ttname);
  string ftype=ttname.substr(0,1);
  if(ftype==fa) isfq=0;
  else isfq=1;

  return(isfq);
}

// ---------------------------------------- //
int readfastq(char* file, int saveinfo=0, string otype="", int readseq=0, int minlen=0, string selctg="")
// ---------------------------------------- // //otype="same" for writing fastq, fasta for fasta
{ 
  igzstream infile(file);
  char fq[5]={"@"};
  char fa[5]={">"};
  char plus[5]={"+"};
  int nseq=0;

 
  rlen.reserve(100000);
  rname.reserve(100000);
  if(readseq){
    rseq.reserve(100000);
    rqual.reserve(100000);
  }


  string read;
  string lname;
  string lcomment="";   
  string lseq="";
  int seqlen=0;
  int quallen=0;
  string lqual;
  int seqlines=0;
  int ctgfound=0;


  int stop=1;
  while(stop){
    getline(infile,read);
    
    if(read.substr(0,1)==fq){  // name
      nseq++;

      if(nseq>1){ // previous
	//cout << selctg << " " << lname << endl;

	if(seqlen>=minlen){
	  if(selctg.size()==0 || selctg==lname){
	    ctgfound++;

	    if(saveinfo){
	      rname.push_back(lname);
	      if(lcomment.size())rcomment.push_back(lcomment);
	      rlen.push_back(seqlen);
	      if(readseq)rseq.push_back(lseq);
	      if(readseq)rqual.push_back(lqual);
	    }
	
	    if(otype.size()){
	      if(otype!="same"){
		myfile << fa << lname ;
		if(lcomment.size()) myfile << " " << lcomment <<endl;
		else myfile << endl;
		myfile << lseq << endl;
	      }else if(lqual.size()){
		myfile << fq << lname ;
		if(lcomment.size()) myfile << " " << lcomment <<endl;
		else myfile << endl;
		myfile << lseq << endl;
		myfile << "+" << endl << lqual << endl;
	      }else{
		cout << " Error! trying to write a fastq but quality string is empty! " << endl;
		return 1;
	      }
	    }
	  
	    if(quallen != seqlen)
	      cout << " ERROR! seq length different from quality lenght!! " << endl;
	  }
	}
      }
      
      
      std::tie(lname,lcomment) = getname(read); 

      lseq="";
      lqual="";

      seqlen=0;
      seqlines=0;
      quallen=0;
    }else if(read.substr(0,1)==plus){ // + and qual

      for(int ll=0; ll<seqlines; ll++){
	getline(infile,read);
	lqual.append(read);
	quallen+=read.size();
      }
    }else{ // sequence 
      lseq.append(read);
      seqlines++;
      seqlen+=read.size();
    }
 
    // EOF
    if(infile.eof()){ // previous
 
      if(seqlen>=minlen){
	if(selctg.size()==0 || selctg==lname){
	  ctgfound++;

	  if(saveinfo){
	    rname.push_back(lname);
	    if(lcomment.size())rcomment.push_back(lcomment);
	    rlen.push_back(seqlen);
	    if(readseq)rseq.push_back(lseq);
	    if(readseq)rqual.push_back(lqual);
	  }
	  
	  if(otype.size()){
	    if(otype!="same"){
	      myfile << fa << lname ;
	      if(lcomment.size()) myfile << " " << lcomment <<endl;
	      else myfile << endl;
	      myfile << lseq << endl;
	    }else if(lqual.size()){
	      myfile << fq << lname ;
	      if(lcomment.size()) myfile << " " << lcomment <<endl;
	      else myfile << endl;
	      myfile << lseq << endl;
	      myfile << "+" << endl << lqual << endl;
	    }else{
	      cout << " Error! trying to write a fastq but quality string is empty! " << endl;
	      return 1;
	    }
	  }
	  
	  if(quallen != seqlen)
	    cout << " ERROR! seq length different from quality lenght!! " << endl;
	}
      }
      stop=0;
    }


  }//read loop
 

  if(!ctgfound) {
    cout << "Could not find contig/scaffold with requested lenght or name!" << endl;
    return 1;
  }

  return 0;
}


// ---------------------------------------- //
int readfasta(char* file, int saveinfo=0,  string otype="", int readseq=0, int minlen=0, string selctg="")
// ---------------------------------------- //otype="same" for writing 
{ 
  igzstream infile(file);
  char fa[5]={">"};
  int nseq=0;
  int thispri=0;

  rlen.reserve(100000);
  rname.reserve(100000);
  if(readseq)
    rseq.reserve(100000);

  if(thispri)cout << file << " " << saveinfo << " " << otype 
		  << " "<< readseq << " " << minlen << endl;


  string read;
  string lname;
  string lcomment="";   
  string lseq="";
  int seqlen=0;
  int ctgfound=0;

  int stop=1;
  while(stop){
    getline(infile,read);
    
    if(read.substr(0,1)==fa){  // name
      nseq++;
      if(thispri)cout << nseq << endl;

      if(nseq>1){ // previous
	if(seqlen>=minlen){
	  if(selctg.size()==0 || selctg==lname){

	    ctgfound++;
	    if(thispri)cout << " valid contig found " 
			    << nseq << " " << seqlen << endl;

	    if(saveinfo){
	      rname.push_back(lname);
	      rlen.push_back(seqlen);
	      if(readseq)rseq.push_back(lseq);
	      if(lcomment.size())rcomment.push_back(lcomment);
	    }
		  
	    if(otype.size()){
	      if(thispri)cout << " writing in " << myfile << endl; 

	      myfile << fa << lname ;
	      if(lcomment.size()) myfile << " " << lcomment <<endl;
	      else myfile << endl;
	      myfile << lseq << endl;
	    }
	  }
	}
      }
      

      std::tie(lname,lcomment) = getname(read); 


      lseq="";
      seqlen=0;

    }else{ // sequence 
      lseq.append(read);
      seqlen+=read.size();
    }
 
    // EOF
    if(infile.eof()){ // previous
 
      if(seqlen>=minlen){
	if(selctg.size()==0 || selctg==lname){
	  ctgfound++;

	  if(saveinfo){
	    rname.push_back(lname);
	    rlen.push_back(seqlen);
	    if(readseq)rseq.push_back(lseq);
	    if(lcomment.size())rcomment.push_back(lcomment);
	  }      
	  
	  
	  if(otype.size()){
	    myfile << fa << lname ;
	    if(lcomment.size()) myfile << " " << lcomment <<endl;
	    else myfile << endl;
	    myfile << lseq << endl;
	  }
	}
      }
      stop=0;
    }
  }//read loop
 
  if(!ctgfound) {
    cout << "Could not find contig/scaffold with requested lenght or name!" << endl;
    return 1;
  }

  return 0;
}

// numeric to string
template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}
// string of numbers to int 
template <class T1>
inline int to_int (const T1& t)
{
  int ii;
  std::stringstream ss(t);
  ss >> ii;
  return ii;
}
// string of numbers to int 
template <class T2>
inline int to_float (const T2& t)
{
  float ii;
  std::stringstream ss(t);
  ss >> ii;
  float rr=ii*100.;
  return rr;
}




template<class Iterator>
Iterator Unique(Iterator first, Iterator last) 
{ // eliminate repetitive elements from sorted vector
  while (first != last)
    {
      Iterator next(first);
      last = std::remove(++next, last, *first);
      first = next;
    }
  return last;
}

