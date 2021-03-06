#include "../myinc/macro.h"

static int printn=1000; // split in 1K contigs per file
static string seqfile ;
static string aligner="smalt";

int split_by_ctgs(char* file);
int split_by_size(char* file);

int main(int argc, char *argv[])
{ 
  int select=0;
  
  if (argc == 1) {
   fprintf(stderr, "Usage: %s <reads.fq/fa>  <aligner> <splitting_size>\n", argv[0]);
   return 1;
  }	
  if((fp = gzopen(argv[1],"r")) == NULL){ 
    printf("ERROR main:: missing input file 1 !! \n");
    return 1;
  }
  gzclose(fp);
 

  // fasta/q file 
  seqfile = argv[1];
  if(argc>2) 
    aligner=to_string(argv[2]);
  if(argc>3) 
    printn=to_int(argv[3]);
 
  if(0) cout << argc << " " << seqfile << " " <<  aligner << " " << printn << endl;

  if(0)cout << aligner << endl;
  //if(aligner == "smalt"){
  // cout << "ok" << endl;
  //}
  
  if(aligner == "smalt")
    split_by_ctgs(argv[1]);
  else
    split_by_size(argv[1]);

 


  return 0;
}


// ---------------------------------------- //
int split_by_size(char* file)
// ---------------------------------------- //
{ 

  if(0)cout << " splitting by size " << endl;

  igzstream infile(file);
  char fq[5]={"@"};
  char fa[5]={">"};
  char plus[5]={"+"};
  int nseq=0;

 
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
  int ll=0;
  int line=0;
  int lprint=0;

  while(stop){
    
    string myname="";
    string aname="";

    getline(infile,read);
   
    if(read.substr(0,1)==fa){  // name
      nseq++;
      
      if(nseq>1){ // previous
	ctgfound++;
	//cout << " read new contig " << read << " " << lprint << endl;

	
	if(lprint>=printn){ //next file
	  myfile.close();
	  ll++;
	  aname="split"+to_string(ll)+"_";
	  myname=myrename(seqfile,aname);
	  myfile.open(myname);  
	  //cout << "   new file opened " << myname << endl;
	  lprint=0;
	}
	/*else if(lprint==0){ // first ctg
	  aname="split"+to_string(ll)+"_";
	  myname=myrename(seqfile,aname);
	  myfile.open(myname);  
	  cout << "   first file opened " << myname << " " << lname<< endl;
	  }*/

	lprint+=seqlen;	
	//cout << " writing new contig in file " << myname << " " << lname << " " << lprint << endl;
	myfile << lname << endl << lseq << endl;
	
      }
      if(nseq==1){ // first ctg
	aname="split"+to_string(ll)+"_";
	myname=myrename(seqfile,aname);
	myfile.open(myname);  
	//cout << "   first file opened " << myname << " " << lname<< endl;
      }

      lname=read;

      lseq="";
      lqual="";

      seqlen=0;
      seqlines=0;
      quallen=0;
    }else{ // sequence 
      lseq.append(read);
      seqlines++;
      seqlen+=read.size();
    }
 
    // EOF
    if(infile.eof()){ // previous
      ctgfound++;

      if(lprint>=printn){ //next file
	myfile.close();
	ll++;
	aname="split"+to_string(ll)+"_";
	myname=myrename(seqfile,aname);
	myfile.open(myname);  
	//cout << "   new file opened " << myname << endl;
      }
      myfile << lname << endl << lseq << endl;
      myfile.close();

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
int split_by_ctgs(char* file)
// ---------------------------------------- //
{ 

  if(0)cout << "split in printn" << endl;
  igzstream infile(file);
  char fq[5]={"@"};
  char fa[5]={">"};
  char plus[5]={"+"};
  int nseq=0;

 
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
  int ll=0;
  int line=0;
  int lprint=0;
  while(stop){
    
    string myname="";
    string aname="";

    getline(infile,read);
   
    if(read.substr(0,1)==fa){  // name
      nseq++;

      if(nseq>1){ // previous
	ctgfound++;
	
	if(lprint==printn){
	  myfile.close();
	  ll++;
	  aname="split"+to_string(ll)+"_";
	  myname=myrename(seqfile,aname);
	  myfile.open(myname);  
	  lprint=0;
	}else if(lprint==0){
	  aname="split"+to_string(ll)+"_";
	  myname=myrename(seqfile,aname);
	  myfile.open(myname);  
	}
	lprint++;
	
	myfile << fa << lname << endl << lseq << endl;
	
      }
      

      size_t ns=0;
      size_t nt=0;
      ns=read.find(" ");
      nt=read.find("\t");

      if(ns!=std::string::npos) { 
	lname=read.substr(1,ns-1);
      }else if(nt!=std::string::npos) {
	lname=read.substr(1,nt-1);
      }else{
	lname=read.erase(0,1);
      }	

      
      int f1=lname.find("/1");
      if (f1 != string::npos){
	lname.erase(f1,2);
      }
      int f2=lname.find("/2");
      if (f2 != string::npos){
	lname.erase(f2,2);
      }
      

      lseq="";
      lqual="";

      seqlen=0;
      seqlines=0;
      quallen=0;
    }else{ // sequence 
      lseq.append(read);
      seqlines++;
      seqlen+=read.size();
    }
 
    // EOF
    if(infile.eof()){ // previous
      ctgfound++;
      myfile << fa << lname << endl << lseq << endl;

      stop=0;
    }
  }//read loop
  if(!ctgfound) {
    cout << "Could not find contig/scaffold with requested lenght or name!" << endl;
    return 1;
  }
  return 0;
}


