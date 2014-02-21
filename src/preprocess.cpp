#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <string>
#include <vector>
#include <unistd.h>
using namespace std;

/**** user headers ****/
#include "functions.h"
#include "readref.h"
#include "read_CIGAR.h"
//#include "matchreads.h"
#include "color.h"

/**** samtools headers ****/
#include <bam.h>
#include <sam.h>

vector<char> op(100);    // CIGAR operator
vector<int> opnum(100);  // CIGAR operator bases
vector<int> oppos(100);  // positions of operators in REF
vector<int> opposseq(100);  // positions of operators in SEQ
size_t anchor,iclip;


int bamread1(samfile_t *fp_in, bam_iter_t& iter, bam1_t *b, string &bamRegion)
{
  return bamRegion.empty() ? 
    samread(fp_in, b) : 
    bam_iter_read(fp_in->x.bam, iter, b);
}

bool is_keeping_read(string& SEQ, string& CIGAR, int POS1, string& FASTA, int minS, int minq, int minQ)
{
  
  read_POS_CIGAR(POS1,CIGAR,
		 anchor, iclip, 
		 op, opnum,
		 opposseq, oppos );
  
  if ( op[0]=='*' ) return false;
  if ( iclip==string::npos ) return false;
  if ( opnum[iclip] < minS ) return false;
  
  calibrate_cigar(SEQ, POS1, FASTA,
		  anchor, iclip, 
		  op, opnum,
		  opposseq, oppos);
  
  if ( iclip==string::npos ) return false;
  if ( opnum[iclip] < minS ) return false;
  if ( opnum[iclip]*2 > (int)SEQ.size() ) return false;
  
  int Sotherend=0;
  if ( iclip==0 && op.back()=='S' ) Sotherend=opnum.back();
  if ( iclip>0  && op[0]=='S' ) Sotherend=opnum[0];
  if ( Sotherend > SEQ.size()*0.08 ) return false;
  
  //  if ( minQ>0 ) {
  //string CLIPPEDQUAL=QUAL.substr(opposseq[iclip],opnum[iclip]);
  //int bQ=0, bQCount=0;
  //for(i=0;i<CLIPPEDQUAL.size();++i) {
  //  bQ=CLIPPEDQUAL[i]-33;
  //  if ( bQ < minQ ) ++bQCount;
  //}
  //if ( bQCount > CLIPPEDQUAL.size()/3 ) return false;
  //}
  
  string CLIPPEDSEQ=SEQ.substr(opposseq[iclip],opnum[iclip]);
  string REFCLIP=FASTA.substr(oppos[iclip]-1,opnum[iclip]);
  int ndiff=0,Ncount=0;
  for(size_t i=0;i<REFCLIP.size();++i) {
    if ( CLIPPEDSEQ[i]=='N' ) Ncount++;
    if ( REFCLIP[i]==CLIPPEDSEQ[i] ) continue;
    ++ndiff;
  }
  if ( Ncount > minS/2 ) return false;
  if ( ndiff <= 2 ) return false;    
  if ( ndiff <= int(0.25*opnum[iclip]) ) return false;
  
  return true;
}

struct sclip {
  string TYPE;
  string RNAME;
  string SEQM;
  string SEQS;
  string CIGAR;
  int POS;
  int sbeg;
  int send;
  sclip():POS(0),
	  sbeg(0),
	  send(0){};
};


void compact_reads(vector<sclip>& clips)
{
  if ( clips.size()<=1 ) return;
  vector<sclip> kept(0);
  vector<bool> iskept(clips.size(), false);
  
  size_t i,j,k;
  //  for(i=0;i<clips.size();++i) 
  //    cerr << clips[i].sbeg << "\t" << clips[i].send << "\n";

  for(i=0;i<clips.size();++i) {
    //if ( iskept[i] ) continue;
    for(j=i+1;j<clips.size();++j) {
      if ( iskept[j] ) continue;
      int com_beg=max(clips[i].sbeg,clips[j].sbeg);
      int com_end=min(clips[i].send,clips[j].send);
      //  cerr << clips[i].sbeg << "\t" << clips[i].send << "\t"
      //   << clips[j].sbeg << "\t" << clips[j].send << "\t"
      //   << com_beg << "\t" << com_end << endl;
      string com1=clips[i].SEQS.substr(com_beg-clips[i].sbeg, com_end-com_beg+1);
      string com2=clips[j].SEQS.substr(com_beg-clips[j].sbeg, com_end-com_beg+1);
      int ndiff=0;
      for(k=0;k<com1.size();++k) {
	if ( com1[k]!=com2[k] ) ndiff++;
	if ( ndiff>3 ) break;
      }
      if ( ndiff<=3 ) {
	iskept[i]=true;
	iskept[j]=true;
      }
    }
  }
  
  
  for(i=0;i<iskept.size();++i) 
    if ( iskept[i] ) kept.push_back(clips[i]);
  
  clips=kept;

  return;
}


void reduce_reads(string outFile, string TYPE)
{
  
  string output=outFile.substr(0,outFile.size()-2);
  ifstream FIN(outFile.c_str());
  if (!FIN) { cerr << "Can't open " << outFile << endl; exit(0);}
  
  ofstream FOUT(output.c_str());
  if (!FOUT ) { cerr << "Can't write to " << output << endl; exit(0);}
  
  sclip ACP;
  vector<sclip> clips(0);
  
  while ( !FIN.eof() ) {
    string tmps,tmps1;
    getline(FIN,tmps);
    if ( tmps.length()<10 ) continue;
    istringstream iss(tmps);
    
    ACP.TYPE=TYPE;
    iss >> ACP.RNAME 
	>> ACP.POS 
	>> ACP.sbeg 
	>> ACP.send 
	>> ACP.CIGAR;
    if ( TYPE=="MS" ) iss >> ACP.SEQM >> ACP.SEQS ;
    else iss >> ACP.SEQS >> ACP.SEQM ;
    
    if (clips.size()==0) {
      clips.push_back(ACP);
      continue;
    }    
    if ( ACP.RNAME!=clips[0].RNAME  || ACP.sbeg>clips[0].send ) {
      //cerr << clips.size() << endl;
      compact_reads(clips);
      //cerr << "--------------" << endl;
      for(size_t i=0;i<clips.size();++i) {
	FOUT << clips[i].RNAME << "\t" << clips[i].POS << "\t" << clips[i].sbeg << "\t" 
	     << clips[i].send << "\t" <<  clips[i].CIGAR << "\t" ;
	if ( TYPE=="MS" ) FOUT << clips[i].SEQM << "\t" << clips[i].SEQS << "\n";
	else FOUT << clips[i].SEQS << "\t" << clips[i].SEQM << "\n";
      }
      clips.clear();
    }
    clips.push_back(ACP);
    
  }
  FIN.close();
  compact_reads(clips);
  for(size_t i=0;i<clips.size();++i) {
    FOUT << clips[i].RNAME << "\t" << clips[i].POS << "\t" << clips[i].sbeg << "\t" 
	 << clips[i].send << "\t" <<  clips[i].CIGAR << "\t" ;
    if ( TYPE=="MS" ) FOUT << clips[i].SEQM << "\t" << clips[i].SEQS << "\n";
    else FOUT << clips[i].SEQS << "\t" << clips[i].SEQM << "\n";
  }
  clips.clear();
  
}


void preprocess_bam_file(string bamFile, string bamRegion, int minS, int minq, int minQ, string fastaFile)
{
  //  reduce_reads("kirby211442.ms.txt.s", "MS");
  //  reduce_reads("kirby211442.sm.txt.s", "SM");
  //  exit(0);
  
  samfile_t *fp_in = NULL;
  bam1_t *b=NULL;   
  b = bam_init1();
  bam_index_t *bamidx=NULL;
  bam_iter_t iter=0;
  int ref=0, beg=0, end=0x7fffffff;
  
  vector<string> targets(0);
  if ( is_binary(bamFile) ) {
    fp_in = samopen(bamFile.c_str(), "rb", 0);
    bamidx = bam_index_load(bamFile.c_str()); // load BAM index
    for(int i=0;i<fp_in->header->n_targets;++i) 
      targets.push_back( string(fp_in->header->target_name[i]) );
  }
  else return;
  
  if ( bamRegion!="" ) {
    targets.clear();
    targets.push_back(bamRegion);
  }
  
  for(size_t i=0;i<targets.size();++i) cerr << targets[i] << endl;
  //exit(0);
  
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  string host=hostname;
  size_t i=host.find(".");
  if ( i!=string::npos ) host=host.substr(0,i);
  
  string PREFIX=host+to_string(getpid());
  
  string QNAME,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  int FLAG;
  int POS1,mapq=0;
  string FASTA;
  string fastaname="";

  vector<char> op(100);    // CIGAR operator
  vector<int> opnum(100);  // CIGAR operator bases
  vector<int> oppos(100);  // positions of operators in REF
  vector<int> opposseq(100);  // positions of operators in SEQ
  size_t anchor,iclip;
  string  read = "";
  
  for(size_t iRegion=0;iRegion<targets.size();iRegion++) {
    bamRegion=targets[iRegion];
    ref=0;
    beg=0; 
    end=0x7fffffff;
    bam_parse_region(fp_in->header, bamRegion.c_str(), &ref, &beg, &end); 
    iter = bam_iter_query(bamidx, ref, beg, end);
    
    string fileMS=PREFIX+".ms."+bamRegion;
    string fileSM=PREFIX+".sm."+bamRegion;
    ofstream FOUT1(fileMS.c_str());
    ofstream FOUT2(fileSM.c_str());
    
    size_t nRead=0, nSelected=0;
    while( bamread1(fp_in, iter, b, bamRegion) > 0 ) {
      ++nRead;
      if ( nRead%1000000==0 ) cerr << nRead << " reads processed" << endl;
      if ( (int)b->core.qual < minq ) continue;
      
      char *s = bam_format1(fp_in->header, b);
      istringstream iss( s );
      iss >> QNAME
	  >> FLAG
	  >> RNAME 
	  >> POS 
	//>> MAPQ 
	  >> mapq 
	  >> CIGAR 
	  >> MRNM 
	  >> MPOS 
	  >> ISIZE 
	  >> SEQ 
	  >> QUAL 
	  >> OPT;
      free(s);
      
      /* load reference SEQ if necessary */
      if ( RNAME != fastaname ) {
	string RNAMEtmp=RNAME;
	if ( ci_find(RNAME,"chr") != string::npos ) RNAMEtmp=RNAME.substr(3); 
	read_fasta(fastaFile,RNAMEtmp,FASTA);
	fastaname=RNAME;
	cerr << fastaFile << "\t" << fastaname << "\t" << FASTA.size() << endl;
      }
      if ( FASTA.size()<1 ) continue;
      
      if ( QNAME[0]=='@' || QNAME[0]=='#' ) continue;
      if ( CIGAR=="*" ) continue;
      if ( RNAME=="*" ) continue;
      if ( mapq < minq ) continue; 
      
      POS1=atoi(POS.c_str());
      read_POS_CIGAR(POS1,CIGAR,
		     anchor, iclip, 
		     op, opnum,
		     opposseq, oppos );
      
      if ( op[0]=='*' ) continue;
      if ( iclip==string::npos ) continue;
      if ( opnum[iclip] < minS ) continue;
      
      calibrate_cigar(SEQ, POS1, FASTA,
		      anchor, iclip, 
		      op, opnum,
		      opposseq, oppos);
      
      if ( iclip==string::npos ) continue;
      if ( opnum[iclip] < minS ) continue;
      if ( opnum[iclip]*2 > (int)SEQ.size() ) continue;
      
      int Sotherend=0;
      if ( iclip==0 && op.back()=='S' ) Sotherend=opnum.back();
      if ( iclip>0  && op[0]=='S' ) Sotherend=opnum[0];
      if ( Sotherend > SEQ.size()*0.08 ) continue;
      
      // minumum base quality
      if ( minQ>0 ) {
	string CLIPPEDQUAL=QUAL.substr(opposseq[iclip],opnum[iclip]);
	int bQ=0, bQCount=0;
	for(i=0;i<CLIPPEDQUAL.size();++i) {
	  bQ=CLIPPEDQUAL[i]-33;
	  if ( bQ < minQ ) ++bQCount;
	}
	if ( bQCount > (int)CLIPPEDQUAL.size()/3 ) continue;
      }
      
      // check mismatched bases
      string CLIPPEDSEQ=SEQ.substr(opposseq[iclip],opnum[iclip]);
      string REFCLIP=FASTA.substr(oppos[iclip]-1,opnum[iclip]);
      int ndiff=0,Ncount=0;
      for(i=0;i<REFCLIP.size();++i) {
	if ( CLIPPEDSEQ[i]=='N' ) Ncount++;
	if ( REFCLIP[i]==CLIPPEDSEQ[i] ) continue;
	++ndiff;
      }
      if ( Ncount > minS/2 ) continue;
      if ( ndiff <= 2 ) continue;    
      if ( ndiff <= int(0.25*opnum[iclip]) ) continue;
      
      // check mismatched bases for M part
      ndiff=0;
      for(i=0;i<op.size();++i) {
	if ( op[i]!='M' ) continue;
	string SEQM=SEQ.substr(opposseq[i],opnum[i]);
	string REFM=FASTA.substr(oppos[i]-1,opnum[i]);
	for(int j=0;j<opnum[i];++j) if ( SEQM[j]!=REFM[j] ) ++ndiff;
      }
      if ( ndiff>SEQ.size()*0.08 ) continue;
      
      if ( POS1!=oppos[anchor] ) cerr << POS1 << "\t" << oppos[anchor] << endl;
      
      ++nSelected;
      string TYPE = iclip>anchor ? "MS":"SM" ;
      string CIGAR1="";
      for(i=0;i<op.size();++i) CIGAR1+=to_string(opnum[i])+to_string(op[i]);
      if ( TYPE=="MS" ) {
	FOUT1 << RNAME << "\t" 
	  //<< POS1 << "\t" 
	      << oppos[anchor] << "\t"
	      << oppos[iclip] << "\t"
	      << oppos[iclip]+opnum[iclip]-1 << "\t"
	      << CIGAR1 << "\t" 
	      << SEQ.substr(0,opposseq[iclip]) << "\t"
	      << SEQ.substr(opposseq[iclip]) << endl;
      }
      if ( TYPE=="SM" ) {
	FOUT2 << RNAME << "\t" 
	  //<< POS1 << "\t" 
	      << oppos[anchor] << "\t"
	      << oppos[iclip] << "\t"
	      << oppos[iclip]+opnum[iclip]-1 << "\t"
	      << CIGAR1 << "\t" 
	      << SEQ.substr(0,opnum[0]) << "\t"
	      << SEQ.substr(opposseq[1]) << endl;
      }
    }
    FOUT1.close();
    FOUT2.close();
    
    string outFile=fileMS+".s";
    string sortcommand="sort -k4n "+fileMS+" -o "+outFile;
    system( sortcommand.c_str() );
    remove( fileMS.c_str() );
    outFile=fileSM+".s";
    sortcommand="sort -k4n "+fileSM+" -o " + outFile;
    system( sortcommand.c_str() );
    remove( fileSM.c_str() );
    
    outFile=fileMS+".s";
    string TYPE="MS";
    reduce_reads(outFile,TYPE);
    
    outFile=fileSM+".s";
    TYPE="SM";
    reduce_reads(outFile,TYPE);
    
  }
  samclose(fp_in);
  bam_destroy1(b);
  bam_iter_destroy(iter);
  bam_index_destroy(bamidx);
  
  return;
}

