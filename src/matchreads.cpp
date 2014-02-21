#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <complex>
#include <map>
#include <algorithm>
#include <string>
#include <vector>
#include <unistd.h>
#include "color.h"
using namespace std;

/**** user samtools headers ****/
#include "samfunctions.h"

/**** user headers ****/
#include "functions.h"
#include "readref.h"
#include "read_CIGAR.h"
#include "preprocess.h"
#include "filters.h"
#include "matchreads.h"
vector<string> BREAKSEEK_st::TARGET;
map<string, int> BREAKSEEK_st::RNAMEINT;
int BREAKSEEK_st::ctid;


int NUM_THREADS=1;
pthread_attr_t attr;
pthread_mutex_t vout, nout, add;
pthread_t threads[37];
struct thread_data_t thread_data[37];

bool COLORED=true;   // use color in intermediate matching result
bool VERBOSE=false;  // print intermediate matching result
bool DUMPBAM=false;  // dump selected reads
bool HEADER=false;   // has HEADER been printed
int MAXDISTANCE=(int)1E6; // distance to search
int MINIMUMS=11;     // minimum S number
int MINIMUMOVERLAP=25;     // minimum length of common string
int MINIMUMMAPQ=10;  // minumum mapping quality
int MINIMUMBASEQ=0; // minumum base quality
int ERR_EDGE=1;      // ignore error in first EDGE and last EDGE bases in a read
int ERR_MATCH=3;     // allowed mismatch when matching ms sm reads
bool COMPACT=false;
bool RDFILTER=true;

int maxMismatch=2;
long int BUFFERSIZE=(int)600E6; // buffer size to hold reads, two buffers used
string SEQBUFFER="";

//! commonString : common string with a tolorance of ERR_MATCH with ERR_EDGE ignored
//! p1 : position of the beginning of match in readMS, 0-based
//! p2 : position of the beginning of match in readSM, 0-based; 0 is matched 
void readsMSandSMmatchfuzzy(const string& readMS, int Sms, 
			    const string& readSM, int Ssm, 
			    string& commonString, int& p1, int& p2) 
{
  p1=p2=-1;
  commonString="";
  
  int i,j,k;
  int maxerror=ERR_MATCH;
  
  bool match=false;
  int ndiff=0;
  
  //! look for deletion and tandem duplication
  //! for those, overlap is equal or longer than Sms + Ssm
  int minOver=Sms+Ssm;
  for(i=readMS.length()-minOver+1-1;i>=0;--i) {
    ndiff=0;
    for(k=i+ERR_EDGE,j=ERR_EDGE;
	k<(int)readMS.length()-ERR_EDGE && j<(int)readSM.length();
	++k,++j) {
      ndiff += ( readMS[k]!=readSM[j] );
      if (ndiff>maxerror) break;
    }
    if (ndiff<=maxerror) { 
      match=true;
      p1=i;
      p2=0;
      //commonString=readMS.substr(p1);
      string comMS=readMS.substr(p1);
      string comSM=readSM.substr(0,comMS.length());
      int c2=(comMS.length()-Sms-Ssm)/2;
      int L_MS=Ssm+c2;
      //int L_SM=comMS.length()-L_MS;
      if ( L_MS < 0 ) L_MS=comMS.length()/2;
      commonString=comMS.substr(0,L_MS)+comSM.substr(L_MS);
      //commonString=comSM.substr(0,L_MS)+comMS.substr(L_MS);
      break;
    }
  }
  if ( ndiff<=maxerror ) return;
  
  if ( MAXDISTANCE>2000 ) return;
  
  //! check if overlap is short random insertion
  //! for those overlap is shorter than Sms+Ssm
  //! we require the reads overlap at least 30 bases
  minOver=30;
  int maxOver=Sms+Ssm;
  for(i=readMS.length()-minOver+1-1;
      i>=0 && i>=(int)readMS.length()-maxOver-1;
      --i) 
    {
      ndiff=0;
      for(k=i+ERR_EDGE,j=ERR_EDGE;
	  k<(int)readMS.length()-ERR_EDGE && j<(int)readSM.length();
	  ++k,++j) {
	ndiff += ( readMS[k]!=readSM[j] );
	if (ndiff>maxerror) break;
      }
      if (ndiff<=maxerror) { 
	match=true;
	p1=i;
	p2=0;
	//commonString=readMS.substr(p1);
	string comMS=readMS.substr(p1);
	string comSM=readSM.substr(0,comMS.length());
	int c2=(comMS.length()-Sms-Ssm)/2;
	int L_MS=Ssm+c2;
	//int L_SM=comMS.length()-L_MS;
	if ( L_MS < 0 ) L_MS=comMS.length()/2;
	commonString=comMS.substr(0,L_MS)+comSM.substr(L_MS);
	break;
      }
    }
  
  if ( ndiff>maxerror ) {
    p1=p2=-1;
    commonString="";
  };
  return;
  
}


void matching_reads_thread( int thread_id,
			    int NUM_THREADS,
			    string& SEQBUFFER,
			    vector<RSAI_st>& MSREAD,
			    vector<RSAI_st>& SMREAD,
			    string& FASTA,
			    string& RNAME,
			    vector<BREAKSEEK_st>& BPSEEK )
{
  if ( MSREAD.size()==0 || SMREAD.size()==0 ) return;
  
  size_t i,k;
  int CL,p1,p2;
  string readMS, readSM;
  size_t istart=0;
  
  string color_reset = RESET;
  string color_red = RED ;
  string color_yellow = YELLOW ;
  string color_blue = BLUE ;
  
  if ( ! COLORED ) {
    color_reset = "";
    color_red = "" ;
    color_yellow = "" ;
    color_blue = "" ;
  }
  
  if ( VERBOSE ) {
    pthread_mutex_lock(&nout);
    cerr << "thread id : " << thread_id << " / " << NUM_THREADS << "\n"
	 << "MS reads / buffer size : " 
	 << MSREAD.size() << "\t" << SEQBUFFER.size() << "\n"
	 << "SM reads / buffer size : " 
	 << SMREAD.size() << "\t" << SEQBUFFER.size() << endl;
    if ( ERR_EDGE>=2 ) 
      cerr << "warning : allowed number of mismatches maybe too large " 
	   << ERR_EDGE << endl;
    pthread_mutex_unlock(&nout);
  }
  
  int REFPOS_com_beg;
  int REFPOS_com_end;
  string MATCHED;
  string MATCHEDREF;
  string commonString;
  
  BREAKSEEK_st ABP;
  vector<BREAKSEEK_st> INSERT(0);
  vector<BREAKSEEK_st> CNV(0);
  BREAKSEEK_st Empty;
  BREAKSEEK_st iindel=Empty;
  
  int imm=MSREAD[0].pos/1000000;
  for(i=thread_id; i<MSREAD.size(); i+=NUM_THREADS) {
    readMS=SEQBUFFER.substr(MSREAD[i].p1, MSREAD[i].len);
    iindel=Empty;
    if ( thread_id==0 && MSREAD[i].pos/1000000 > imm ) {
      imm=MSREAD[i].pos/1000000;
      cerr << "#" << RNAME << "@" << MSREAD[i].pos << '\xD'; 
    }
    for(k=istart; k<SMREAD.size(); ++k) {
      if ( SMREAD[k].pos < MSREAD[i].pos - MAXDISTANCE ) {
	istart=k+1;
	continue;
      }
      if ( SMREAD[k].pos > MSREAD[i].pos + MAXDISTANCE ) break;
      
      readSM=SEQBUFFER.substr(SMREAD[k].p1, SMREAD[k].len);
      
      readsMSandSMmatchfuzzy(readMS, MSREAD[i].S, 
			     readSM, SMREAD[k].S, 
			     MATCHED, p1, p2);
      
      // not MS-SM matched 
      if ( p1<0 || p2<0 ) continue;
      
      CL=MATCHED.length();
      p2=0;
      
      // short random insertion
      if ( CL < (MSREAD[i].S+SMREAD[k].S) && CL>30 ) {
	int insbp1=MSREAD[i].Spos-1;
	int insbp2=SMREAD[k].pos;
	int insgap=insbp2-insbp1;
	int inslen=MSREAD[i].S+SMREAD[k].S-CL ;
	string ins="",gap="";
	
	string MERGED=readMS+readSM.substr(CL);
	string MERGEDNOINS=MERGED.substr(0,readMS.length()-MSREAD[i].S)
	  +MERGED.substr(readMS.length()-MSREAD[i].S+inslen);
	string MERGEDNOINSREF=
	  FASTA.substr(MSREAD[i].pos-MSREAD[i].Mrpos-1,readMS.length()-MSREAD[i].S)+
	  FASTA.substr(SMREAD[k].pos-1, MERGED.size()-readMS.length()+MSREAD[i].S-inslen);
	
	if ( MERGEDNOINS.length() != MERGEDNOINSREF.length() ) {
	  cerr << "length error for insert\n"; 
	  exit(0);
	}
	
	int med=edit_distance( MERGEDNOINS,  MERGEDNOINSREF );
	if ( med > (int)MERGEDNOINS.length()*0.08 && med>=4 ) continue;
	if ( med >= (int)SMREAD[k].len*0.06 ) continue;
	
	if ( inslen > 2.5*abs(insgap) ) {
	  if ( MSREAD[i].S >= inslen ) 
	    ins=readMS.substr(readMS.size()-MSREAD[i].S, inslen);
	  else if ( SMREAD[k].S >= inslen )
	    ins=readSM.substr(SMREAD[k].S-inslen,inslen);
	  else {
	    ins=readMS.substr(readMS.size()-MSREAD[i].S)
	      +readSM.substr(CL,SMREAD[k].S-CL);
	  }
	  if ( insbp2>insbp1 ) gap=FASTA.substr(insbp1-1,insbp2-insbp1-1);
	  else gap=FASTA.substr(insbp2-1,insbp1-insbp2+1);
	  
	  if ( insbp1==insbp2 ) {
	    insbp2++;
	    inslen++;
	    ins+=FASTA[insbp1-1];
	    insgap=insbp2-insbp1;
	  }
	  
	  iindel.T = (insbp1 < insbp2) ? 'I' : 'J';
	  iindel.tid=MSREAD[i].tid;
	  iindel.P1=insbp1;
	  iindel.P2=insbp2;
	  iindel.UN=0;
	  iindel.L1=MSREAD[i].len;
	  iindel.S1=MSREAD[i].S;
	  iindel.L2=SMREAD[k].len;
	  iindel.S2=SMREAD[k].S;
	  iindel.q1=MSREAD[i].q1;
	  iindel.q2=SMREAD[k].q1;
	  iindel.n1=MSREAD[i].mm;
	  iindel.n2=SMREAD[k].mm;
	  iindel.CL=CL;
	  iindel.ED=med;
	  // iindel.len= (insbp1<insbp2) ? insbp2-insbp1-1:insbp1-insbp2+1;
	  iindel.len= ins.length();
	  iindel.INSEQ=ins;
	  iindel.MERGE=readMS+readSM.substr(CL);
	  iindel.count=1;
	  
	  string SV30="123456789012345678901234567890123456789012345678900";
	  string tmps="NA";
	  // using reference
	  int svlen=50;
	  if ( insbp1-svlen/2>0 && insbp2+svlen/2<(int)FASTA.size() )
	    tmps=FASTA.substr(insbp1-svlen/2+1-1,svlen/2)+FASTA.substr(insbp2-1,svlen/2);
	  SV30.replace(0,tmps.size(),tmps);
	  SV30[tmps.size()]=0;
	  
	  strcpy(iindel.SVSEQ,SV30.substr(0,svlen+1).c_str());
	  
	  bool NEW=false;
	  if ( INSERT.size()==0 ) {
	    INSERT.push_back(iindel) ;
	    NEW=true;
	  }
	  else if ( iindel.P1 != INSERT.back().P1 || 
		    iindel.P2 != INSERT.back().P2 || 
		    iindel.INSEQ != INSERT.back().INSEQ ) {
	    INSERT.push_back(iindel) ;
	    NEW=true;
	  }
	  else {
	    INSERT.back().count++;
	    if (iindel.CL > INSERT.back().CL) {
	      iindel.count=INSERT.back().count++;
	      INSERT.back()=iindel;
	    }
	    NEW=false;
	  }
	  
	  if ( NEW ) {
	    cerr << INSERT.back().S1 << " + " << INSERT.back().S2 << "\t" 
		 << INSERT.back().CL << "\tins="
		 << INSERT.back().INSEQ.size() << "\t"
		 << INSERT.back().P1 << "\t" 
		 << INSERT.back().P2 << "\tgap="
		 << insgap << "\n"
		 << INSERT.back().INSEQ << "\t" << gap << "\t" 
		 << INSERT.back().MERGE << endl;
	  }
	  
	}

	continue;

      }
      
      //! calculate break point from read position
      //! 5' end position of common string on reference
      REFPOS_com_beg=p1-MSREAD[i].Mrpos+MSREAD[i].pos;
      // 3' end position of common string on reference
      REFPOS_com_end=CL-1-SMREAD[k].Mrpos+SMREAD[k].pos;
      
      //! calculate break point from mismatched position
      //! 5' end position of common string on reference
      int REFPOS_com_beg1=MSREAD[i].Spos-CL+MSREAD[i].S;
      //! 3' end position of common string on reference
      int REFPOS_com_end1=CL-1-SMREAD[k].S+SMREAD[k].pos;
      //! the difference is due to short indel in CIGAR
      //! an alternative is to calculate based on the CIGAR strings
      //! will be implemented later
      int dxindel_beg=(int)abs(REFPOS_com_beg1-REFPOS_com_beg);
      int dxindel_end=(int)abs(REFPOS_com_end1-REFPOS_com_end);
      /*
      if ( REFPOS_com_beg1 != REFPOS_com_beg ) 
	cerr << "REFPOS_com_beg1\t" 
	     << REFPOS_com_beg1 << "\t" << REFPOS_com_beg << "\t"
	     << REFPOS_com_beg1-REFPOS_com_beg <<  endl;
      if ( REFPOS_com_end != REFPOS_com_end ) 
	cerr << "REFPOS_com_end1\t" 
	     << REFPOS_com_end1 << "\t" << REFPOS_com_end << "\t" 
	     << REFPOS_com_end1 - REFPOS_com_end << endl;
      */
      REFPOS_com_beg=REFPOS_com_beg1;
      REFPOS_com_end=REFPOS_com_end1;
      
      //! check if _beg and _end positions are accurate
      string left10=MATCHED.substr(0,MINIMUMS);
      int ndiff=0;
      for(int i1=0;i1<(int)left10.size();++i1) 
	ndiff+= ( left10[i1] != FASTA[REFPOS_com_beg-1+i1] );
      if ( ndiff > 1 ) {
	dxindel_beg=max(10, dxindel_beg);
	int minndiff=ndiff;
	int mindx=0;
	for(int dx=-dxindel_beg;dx<=dxindel_beg;++dx) {
	  int ndiffdx=0;
	  for(int i1=0;i1<(int)left10.size();++i1) {
	    ndiffdx+= ( left10[i1] != FASTA[REFPOS_com_beg+dx-1+i1] );
	    if ( ndiffdx>minndiff ) break;
	  }
	  if ( ndiffdx<minndiff ) {
	    minndiff=ndiffdx;
	    mindx=dx;
	  }
	  if ( minndiff <=1 ) break;
	}
	if ( minndiff < ndiff ) {
	  REFPOS_com_beg+=mindx;
	  //cerr << "REFPOS_com_beg += " << mindx << "\t" 
	  //     << ndiff << "\t" << minndiff << endl;
	}
      }
      
      string right10=MATCHED.substr(MATCHED.size()-MINIMUMS);
      ndiff=0;
      for(int i1=right10.size()-1, k1=REFPOS_com_end-1; i1>=0;--i1,--k1) 
	ndiff+= ( right10[i1] != FASTA[k1] );
      if ( ndiff > 1 ) {
	int minndiff=ndiff;
	int mindx=0;
	dxindel_end=max(10, dxindel_end);
	for(int dx=-dxindel_end;dx<=dxindel_end;++dx) {
	  int ndiffdx=0;
	  for(int i1=right10.size()-1, k1=REFPOS_com_end-1+dx; i1>=0;--i1,--k1) {
	    ndiffdx+= ( right10[i1] != FASTA[k1] );
	    if ( ndiffdx>minndiff ) break;
	  }
	  if ( ndiffdx<minndiff ) {
	    minndiff=ndiffdx;
	    mindx=dx;
	  }
	  if ( minndiff <=1 ) break;
	}
	if ( minndiff < ndiff ) {
	  REFPOS_com_end+=mindx;
	  //cerr << "REFPOS_com_end += " << mindx << "\t" 
	  //     << ndiff << "\t" << minndiff << endl;
	}
      }
      
      //! find break points by minimizing edit distance between MATCHED and ref
      int mdis=CL;
      int mdx1=CL-p2-1;
      int uncertainty=0;
      int premdis=mdis;
      int premdx1=mdx1;
      int dx1,dx2,bp1,bp2,edistance;
      // break point 1, break point 2
      for(dx1=SMREAD[k].S;dx1<=CL-MSREAD[i].S;++dx1) {
	dx2=CL-dx1;
	bp1=REFPOS_com_beg+dx1-1;
	bp2=REFPOS_com_end-dx2+1;
	MATCHEDREF=FASTA.substr(REFPOS_com_beg-1, dx1)+FASTA.substr(bp2-1, dx2);
	if ( MATCHEDREF.length() != MATCHED.length() ) {
	  cerr << "wrong MATCHEDREF length" << endl;
	  exit(0);
	}
	edistance=edit_distance(MATCHEDREF, MATCHED);
	if (  edistance < mdis ) {
	  premdis=edistance;
	  premdx1=dx1;
	}
	if (  edistance <= mdis ) {
	  mdis=edistance;
	  mdx1=dx1;
	}
      }
      // right position
      //dx1=mdx1;
      //dx2=CL-mdx1;
      // left position
      dx1=premdx1;
      dx2=CL-premdx1;
      bp1=REFPOS_com_beg+dx1-1;
      bp2=REFPOS_com_end-dx2+1;
      //if ( premdis==mdis ) uncertainty=mdx1-premdx1;
      //if ( CL-SMREAD[k].S-MSREAD[i].S>0 ) uncertainty=mdx1-premdx1;
      uncertainty=CL-SMREAD[k].S-MSREAD[i].S;
      if ( uncertainty > abs(bp2-bp1) ) uncertainty=mdx1-premdx1;
      edistance=mdis;
      MATCHEDREF=FASTA.substr(REFPOS_com_beg-1, dx1)+FASTA.substr(bp2-1, dx2);
      
      //! according to the equation in paper; not so accurate
      // bp1=MSREAD[i].Spos-CL+MSREAD[i].S+SMREAD[k].S+1;
      // bp2=SMREAD[k].pos;
      // uncertainty=CL-MSREAD[i].S-SMREAD[k].S;
      
      //! extend uncertainty according to repeats on reference
      //! does not seem to help much, but helps
      if ( bp1+1!=bp2 ) {
	dx1=0;
	int i5=bp1+uncertainty;
	int i3=bp2-1+uncertainty;
	if ( i5<(int)FASTA.size() && i3 <(int)FASTA.size() ) {
	  while(FASTA[i5+dx1]==FASTA[i3+dx1] && FASTA[i5+dx1]!='N' ) {
	    dx1++ ;
	    if ( (i5+dx1)>(int)FASTA.size()-10 || (i3+dx1)>(int)FASTA.size()-10 ) break; 
	    if ( dx1>CL ) break;
	  }
	}
	//if ( dx1>0 ) cerr << "3' extended by " << dx1 << endl;
	
	dx2=0;
	i5=bp1-1;
	i3=bp2-2;
	if ( i5>0 && i3>0 ) {
	  while(FASTA[i5+dx2]==FASTA[i3+dx2] && FASTA[i5+dx2]!='N' ) {
	    dx2--;
	    if ( (i5+dx2)<1 || (i3+dx2)<1 ) break; 
	    if ( abs(dx2)>CL ) break;
	  }
	}
	//if ( dx2<0 ) cerr << "5' extended by " << dx2 << endl;
	bp1+=dx1;
	bp2+=dx2;
	uncertainty+=dx1-dx2;
      }
      
      //! this is the main filter that removes false positives
      //! matched sequence and reference have too many mismatched bases
      //if ( edistance > CL*0.08 ) continue;
      //if ( edistance > CL*0.08 && edistance>=4 ) continue;
      if ( edistance >= (MSREAD[i].len+SMREAD[k].len)*0.03 ) continue;
      
      //! check if merged string belong to either regions
      //! this filter does not seem to do much
      string MERGED=readMS+readSM.substr(CL);
      string REF1=FASTA.substr(REFPOS_com_beg-readMS.size()+CL-1, MERGED.size());
      string REF2=FASTA.substr(REFPOS_com_end-readMS.size(), MERGED.size());
      int e1=edit_distance(MERGED,REF1);
      int e2=edit_distance(MERGED,REF2);
      //int mergedmismatch=int((MERGED.size())*0.08+0.5);
      //if ( e1<=mergedmismatch || e2<=mergedmismatch ) continue;
      
      string INFER;
      char TYPE;
      int LENGTH;
      // continuous; could be short insertion
      if ( bp2-bp1==1 ) {
	INFER="INS@"+RNAME+":"+to_string(bp1)+"-"+to_string(bp2);
	LENGTH=bp1-bp2+1;
	TYPE='I';
	//continue;
      }
      else if ( bp1>=bp2 ) {
	INFER="DUP@"+RNAME+":"+to_string(bp2)+"-"+to_string(bp1);
	LENGTH=bp1-bp2+1;
	TYPE='A';
      }
      else if ( bp1<bp2 ) {
	INFER="DEL@"+RNAME+":"+to_string(bp1)+"-"+to_string(bp2);
	LENGTH=bp2-bp1-1;
	TYPE='D';
      }
      
      string SV30="123456789012345678901234567890123456789012345678900";
      string tmps="NA";
      // using reference
      int svlen=50;
      if ( bp1-svlen/2>0 && bp2+svlen/2<(int)FASTA.size() )
	tmps=FASTA.substr(bp1-svlen/2+1-1,svlen/2)+FASTA.substr(bp2-1,svlen/2);
      SV30.replace(0,tmps.size(),tmps);
      SV30[tmps.size()]=0;
      
      ABP.tid=MSREAD[i].tid;
      ABP.P1=bp1;
      ABP.P2=bp2;
      ABP.UN=uncertainty;
      ABP.L1=MSREAD[i].len;
      ABP.S1=MSREAD[i].S;
      ABP.n1=MSREAD[i].mm;
      ABP.q1=MSREAD[i].q1;
      ABP.L2=SMREAD[k].len;
      ABP.S2=SMREAD[k].S;
      ABP.n2=SMREAD[k].mm;
      ABP.q2=SMREAD[k].q1;
      ABP.CL=CL;
      ABP.len=LENGTH;
      ABP.ED=edistance;
      ABP.ED5p=e1;
      ABP.ED3p=e2;
      ABP.T=TYPE;
      ABP.count=1;
      strcpy(ABP.SVSEQ,SV30.substr(0,svlen+1).c_str());
      ABP.MERGE=MERGED;
      
      bool NEW=true;
      if ( CNV.size()==0 ) NEW=true;
      else if ( ABP.P1==CNV.back().P1 && 
		ABP.P2==CNV.back().P2 &&
		ABP.T==CNV.back().T ) NEW=false;
      if ( NEW ) CNV.push_back(ABP);
      else {
	if ( ABP.CL > CNV.back().CL ) {
	  ABP.MERGE=CNV.back().MERGE;
	  CNV.back()=ABP;
	}
	CNV.back().count++;
      }
      
      if ( VERBOSE ) {
	pthread_mutex_lock(&nout);
	string pad1;
	if ( p2-p1>=0 ) pad1=string(p2-p1,' ');
	else pad1="";
	string pad2;
	if ( p1-p2>=0 ) pad2=string(p1-p2,' ');
	else pad2="";
	
	cout << "INFER\t" 
	     << INFER << "\t"
	     << LENGTH
	     << "\n"
	  
	     << RNAME << ":" << MSREAD[i].pos << "\t" 
	     << MSREAD[i].len << "\t"
	     << MSREAD[i].S << "\t"
	     << p1 << "\t" 
	     << MSREAD[i].len-p1-CL << "\t"
	     << RNAME << ":" << SMREAD[k].pos << "\t" 
	     << SMREAD[k].len << "\t"
	     << SMREAD[k].S << "\t" 
	     << p2 << "\t"
	     << CL 
	     << "\n"
	  
	     << "READ1\t"
	     << pad1 
	     << readMS.substr(0, MSREAD[i].len-MSREAD[i].S) 
	     << color_yellow
	     << readMS.substr(MSREAD[i].len-MSREAD[i].S) 
	     << color_reset
	     << "\n"
	  
	     << "REFE1\t"
	     << pad1 
	     << string(MSREAD[i].Mrpos,' ') 
	     << color_blue
	     << FASTA.substr(MSREAD[i].pos-1, 1)
	     << color_reset
	     << FASTA.substr(MSREAD[i].pos, MSREAD[i].len-MSREAD[i].Mrpos-1)
	     << "\n"
	  
	     << "READ2\t"
	     << pad2 
	     << color_yellow
	     << readSM.substr(0, SMREAD[k].S) 
	     << color_reset
	     << readSM.substr(SMREAD[k].S) 
	     << "\n"
	  
	     << "REFE2\t"
	     << pad2 
	     << FASTA.substr(SMREAD[k].pos-SMREAD[k].Mrpos-1, SMREAD[k].Mrpos)
	     << color_blue
	     << FASTA.substr(SMREAD[k].pos-1, 1)
	     << color_reset
	     << FASTA.substr(SMREAD[k].pos, SMREAD[k].len-SMREAD[k].Mrpos-1)
	     << "\n"
	  
	     << "MERGE\t"
	     << pad1 << readMS.substr(0,p1) 
	     << color_red
	     << readMS.substr(p1, CL) 
	     << color_reset
	     << readSM.substr(p2+CL)
	     << "\n"
	  
	     << "REF11\t"
	     << pad1 << string(p1,' ') 
	     << "1" 
	     << string(dx1-2,'-') 
	     << "1" 
	     << "  " << RNAME << ":" << REFPOS_com_beg << "-" << bp1 
	     << "\n"
	  
	     << "REF22\t"
	     << pad2 
	     << string(dx1,' ') 
	     << "2"
	     << string(dx2-2,'-') 
	     << "2"
	     << "  " << RNAME << ":" << bp2 << "-" << REFPOS_com_end  
	     << "\n"
	  
	     << "REFMG\t"
	     << pad1 << string(p1,' ') 
	     << MATCHEDREF << "\tMLEN:" <<  MATCHEDREF.length() << "\tED:" << mdis 
	     << "\n"
	  
	     << "BREAK\t" << bp1 << "\t" << bp2 
	     << "\n"
	  
	     << endl;
	pthread_mutex_unlock(&nout);
      }
      
    } // for(k=istart; k<SMREAD.size(); ++k) {
    
  } //   for(i=thread_id; i<MSREAD.size(); i+=NUM_THREADS) {
  
  if ( thread_id==0 ) cerr << "\n"; 

  //! lock threads and send back results to memory
  pthread_mutex_lock(&nout);
  BPSEEK.insert(BPSEEK.end(),CNV.begin(),CNV.end());
  BPSEEK.insert(BPSEEK.end(),INSERT.begin(),INSERT.end());
  for(size_t i1=0; i1<INSERT.size(); ++i1) { 
    cerr << INSERT[i1].S1 << " + " << INSERT[i1].S2 << "\t" 
	 << INSERT[i1].CL << "\tins="
	 << INSERT[i1].INSEQ.size() << "\t"
	 << INSERT[i1].P1 << "\t" 
	 << INSERT[i1].P2 << "\tgap="
	 << INSERT[i1].P2-INSERT[i1].P1 << "\n" 
	 << INSERT[i1].INSEQ << "\t"  
	 << INSERT[i1].MERGE << "\t" 
	 << INSERT[i1].count << endl;
  }
  pthread_mutex_unlock(&nout);
  
}

//! pass pointers to threads
void* matching_reads_thread_wrapper(void* threadarg)
{
  struct thread_data_t *my_data = (struct thread_data_t *) threadarg;
  
  int thread_id = my_data->thread_id;
  int NUM_THREADS = my_data->NUM_THREADS;
  string* SEQBUFFER = my_data->SEQBUFFER;
  vector<RSAI_st>* MSREAD = my_data->MSREAD;
  vector<RSAI_st>* SMREAD = my_data->SMREAD;
  string* FASTA = my_data->FASTA;
  string* RNAME = my_data->RNAME;
  vector<BREAKSEEK_st>* BPSEEK = my_data->BPSEEK;
  
  matching_reads_thread(thread_id, NUM_THREADS,
			(*SEQBUFFER), (*MSREAD), (*SMREAD), 
			(*FASTA), (*RNAME),
			(*BPSEEK) );
  
  pthread_exit((void*) 0);
}

//! print out BP
string BP_format(BREAKSEEK_st &BP, bool h)
{
  std::stringstream ss;

  if ( h ) { // header?
    ss << "#RNAME\t" 
       << "START\t"
       << "END\t"
       << "TYPE\t"
       << "UNCERTAINTY\t"
       << "LENGTH\t"
       << "INSERT\t"
       << "MATCHINFO\t"
       << "MATCHCOUNT\t"
       << "VALID\t" 
       << "SVSEQ\t"
       << "MERGE\t"
       << "BP";
  }
  else {
    string BPMARKER="BP_"+to_string(BP.T)+"_"+
      BP.RNAME()+":"+
      to_string(BP.P1)+"-"+to_string(BP.P2); 
    ss << BP.RNAME() << "\t" 
       << BP.P1 << "\t" 
       << BP.P2 << "\t" 
       << BP.T << "\t" 
       << BP.UN << "\t" 
       << BP.len << "\t"
       << BP.INSEQ << "\t"
       << "CL" << BP.CL << ","
       << "ED" << BP.ED << ","
       << "L" << BP.L1 << ","
       << "S" << BP.S1 << ","
       << "N" << BP.n1 << ","
       << "Q" << BP.q1 << ","
       << "L" << BP.L2 << ","
       << "S" << BP.S2 << ","
       << "N" << BP.n2 << ","
       << "Q" << BP.q2 << ";"
       << BP.TAG << "\t"
       << BP.count << "\t"
       << BP.isit << "\t"
       << BP.SVSEQ << "\t"
       << BPMARKER;
    //   << BP.MERGE << "\t"
    //	 << BPSEEK[i].ED5p << "\t"
    //	 << BPSEEK[i].ED3p << "\t"
  }
  
  return ss.str() ;
}

//string S1="ATTTCCCAGGTGCCGGTCCATCCTTGTTGT";
//string S2="CATTTCCCAGGTGCCGGTCCATCCTTGTTG";
//cout << edit_distance(S2,S1) << endl;
//exit(0);
//result is 2
bool is_exact_overlap_svseq(string& S1, string& S2)
{
  int i,k,dx;
  bool is_same=true;
  for(dx=-2; dx<=2; ++dx ) {
    is_same=true;
    for(i=0, k=dx; i<(int)S1.size() && k<(int)S2.size() ; ++i, ++k) {
      if ( k<0 ) continue;
      if ( S1[i]!=S2[k] ) { is_same=false; break;}
    }
    if ( is_same ) break;
  }
  
  return( is_same );
}
//! allow 2 mismatch for edit distance
bool is_same_svseq(string& S1, string& S2)
{
  if ( S1.size()<5 || S2.size()<5 ) return( false );
  if ( S1.find("NA") != string::npos ) return(false);
  if ( S2.find("NA") != string::npos ) return(false);
  
  if ( S1.size()==S2.size() ) {
    int ndiff=0;
    for(size_t i=0;i<S1.length();++i) ndiff+=( S1[i]!=S2[i] );
    if ( ndiff<=1 ) return(true);
    return( edit_distance(S1,S2)<=2 );
  }
  
  if ( S1.size()< S2.size() ) return( S2.find(S1)!=string::npos );
  return( S1.find(S2)!=string::npos );
}

void SVSEQ_filter( vector<BREAKSEEK_st>& BP )
{
  if ( BP.size()<=1 ) return;
  
  vector<int> SVID(0);   // index for same SVSEQ
  vector<bool> is_keep(BP.size(),true);
  for(size_t i=0;i<BP.size();++i) {
    
    if ( BP[i].CL<MINIMUMOVERLAP ) is_keep[i]=false;
    if ( !is_keep[i] ) continue;
    
    string SVSEQ=BP[i].SVSEQ;
    SVID.clear();
    SVID.push_back(i);
    for(size_t k=i+1;k<BP.size();++k) {
      if ( !is_keep[k] ) continue;
      string SVSEQ2=BP[k].SVSEQ;
      if ( is_same_svseq(SVSEQ, SVSEQ2 ) ) SVID.push_back(k);
    }
    if ( SVID.size()==1 ) continue;
    
    //cerr << "------------------------\n";
    //for(size_t k=0;k<SVID.size();++k) 
    //  cerr << BP_format(BP[ SVID[k] ], 0) << endl;
    
    int kept_count=0;
    
    //! select one with maximum count 
    int max_count=0;
    for(size_t k=0;k<SVID.size();++k) 
      if ( BP[ SVID[k] ].count > max_count ) max_count=BP[ SVID[k] ].count;
    for(size_t k=0;k<SVID.size();++k) 
      if ( BP[ SVID[k] ].count < max_count -3 ) is_keep[ SVID[k] ]=false;
    
    kept_count=0;
    for(size_t k=0;k<SVID.size();++k) kept_count+=is_keep[ SVID[k] ];
    if ( kept_count<=1 ) continue;
    
    if ( max_count<=3 ) {
      for(size_t k=0;k<SVID.size();++k) {
	if ( BP[ SVID[k] ].ED > BP[ SVID[k] ].CL*8/100 ) is_keep[ SVID[k] ]=false;
	if ( BP[ SVID[k] ].ED > BP[ SVID[k] ].L1*5/100 ) is_keep[ SVID[k] ]=false;
	if ( BP[ SVID[k] ].ED > BP[ SVID[k] ].L2*5/100 ) is_keep[ SVID[k] ]=false;
      }
    }
    kept_count=0;
    for(size_t k=0;k<SVID.size();++k) kept_count+=is_keep[ SVID[k] ];
    if ( kept_count<=1 ) continue;
    
    //! select one with minimum edit distance/mlen for common string 
    float min_ER= 1.0;
    for(size_t k=0;k<SVID.size();++k) {
      float ER=(float)(BP[ SVID[k] ].ED+0.001) / (float)BP[ SVID[k] ].CL;
      if ( ER < min_ER ) min_ER=ER;
    }
    for(size_t k=0;k<SVID.size();++k) {
      if ( BP[ SVID[k] ].ED == 0 ) continue;
      float ER=(float)(BP[ SVID[k] ].ED+0.001) / (float)BP[ SVID[k] ].CL;
      if ( ER > min_ER+0.000001 ) is_keep[ SVID[k] ]=false;
    }
    kept_count=0;
    for(size_t k=0;k<SVID.size();++k) kept_count+=is_keep[ SVID[k] ];
    if ( kept_count<=1 ) continue;
    
  }
  
  vector<BREAKSEEK_st> BP0(0);
  for(size_t i=0;i<BP.size();++i) if ( is_keep[i] ) BP0.push_back( BP[i] );
  
  BP=BP0;
  return;
}

void compact_breakpoint( vector<BREAKSEEK_st>& BPSEEK,
			 vector<BREAKSEEK_st>& BP)
{
  size_t i,k;
  
  BP.clear();
  if ( BPSEEK.size()<=1 ) {
    BP=BPSEEK;
    return;
  }
  
  int count0=0;
  for(i=0;i<BPSEEK.size();++i) count0+=BPSEEK[i].count;
  
  vector<BREAKSEEK_st> BPGROUP;
  vector<BREAKSEEK_st> BPTMP;
  BPGROUP.clear();
  
  for(i=0;i<BPSEEK.size();++i) 
    if ( BPSEEK[i].P1 > BPSEEK[i].P2 ) swap(BPSEEK[i].P1, BPSEEK[i].P2);
  
  /* sort according to the p1 positions */
  vector<int> p1(BPSEEK.size());
  vector<int> p1_idx(BPSEEK.size());
  for(i=0;i<BPSEEK.size();++i) p1[i]=BPSEEK[i].P1;
  arrayindex(p1,p1_idx,1);
  BPGROUP.resize(BPSEEK.size());
  for(i=0;i<BPSEEK.size();++i) BPGROUP[i]=BPSEEK[p1_idx[i]];
  BPSEEK=BPGROUP;
  for(i=1;i<BPSEEK.size();++i) 
    if ( BPSEEK[i].P1 < BPSEEK[i-1].P1 ) { 
      cerr << "Error sorting according to P1\n"; 
      exit(0); 
    }
  
  int maxUN=2;
  vector<bool> iskept(BPSEEK.size(),true);
  vector<int> gid(0);
  // merge according to both positions
  // allowed distance is 2.
  BP.clear();
  BPGROUP.clear();  
  BPGROUP.push_back(BPSEEK[0]);
  gid.push_back(0);
  for(i=1;i<BPSEEK.size();++i) {
    if ( (int)abs( BPSEEK[i].P1 - BPSEEK[ gid[0] ].P1 ) <= maxUN &&
	 (int)abs( BPSEEK[i].P2 - BPSEEK[ gid[0] ].P2 ) <= maxUN &&
	 BPSEEK[i].T == BPSEEK[ gid[0] ].T &&
	 BPSEEK[i].tid == BPSEEK[ gid[0] ].tid ) {
      gid.push_back(i);
      continue;
    }                      // group closeby 
    else  {                // select one
      int gcount=0;
      int ic_max=0;
      int im_max=0;
      for(k=0;k<gid.size();++k) {
	gcount+=BPSEEK[ gid[k] ].count;
	if ( BPSEEK[ gid[k] ].count > BPSEEK[ gid[ic_max] ].count ) ic_max=k;
	if ( BPSEEK[ gid[k] ].CL > BPSEEK[ gid[ic_max] ].CL ) im_max=k;
      }
      BP.push_back( BPSEEK[ gid[ic_max] ] );
      BP.back().count=gcount;
      
      gid.clear(); gid.push_back(i);
    }
  }
  if ( gid.size()>0 ) {
    int gcount=0;
    int ic_max=0;
    int im_max=0;
    for(k=0;k<gid.size();++k) {
      gcount+=BPSEEK[ gid[k] ].count;
      if ( BPSEEK[ gid[k] ].count > BPSEEK[ gid[ic_max] ].count ) ic_max=k;
      if ( BPSEEK[ gid[k] ].CL > BPSEEK[ gid[ic_max] ].CL ) im_max=k;
    }
    BP.push_back( BPSEEK[ gid[ic_max] ] );
    BP.back().count=gcount;
    
    gid.clear(); gid.push_back(i);
  }
  BPSEEK=BP;
  
  //! merge according to the p2 positions
  //! allowed distance is 3.
  BP.clear();  
  BPGROUP.clear();  
  BPGROUP.push_back(BPSEEK[0]);
  for(i=1;i<BPSEEK.size();++i) {
    if ( BPSEEK[i].P1 == BPGROUP.back().P1 &&
	 BPSEEK[i].T == BPGROUP.back().T &&
	 BPSEEK[i].tid == BPGROUP.back().tid ) {
      BPGROUP.push_back(BPSEEK[i]);
      continue;
    }
    else  { 
      /* sort according to the p2 positions */
      p1.resize(BPGROUP.size());
      p1_idx.resize(BPGROUP.size());
      for(k=0;k<BPGROUP.size();++k) p1[k]=BPGROUP[k].P2;
      arrayindex(p1,p1_idx,1);
      BPTMP.resize(BPGROUP.size());
      for(k=0;k<BPGROUP.size();++k) BPTMP[k]=BPGROUP[p1_idx[k]];
      BPGROUP=BPTMP;
      BPTMP.clear();
      
      BP.push_back(BPGROUP[0]);
      for(k=1;k<BPGROUP.size();++k) {
	if ( BPGROUP[k].P2 - BP.back().P2 < 3 ) {
	  BP.back().count+=BPGROUP[k].count;
	  if ( BPGROUP[k].CL > BP.back().CL ) {
	    BPGROUP[k].count=BP.back().count;
	    BP.back()=BPGROUP[k];
	  }
	  continue;
	}
	BP.push_back(BPGROUP[k]);
      }
      
      BPGROUP.clear();
      BPGROUP.push_back(BPSEEK[i]);
    }
  }
  if ( BPGROUP.size()>0 ) {
    p1.resize(BPGROUP.size());
    p1_idx.resize(BPGROUP.size());
    for(k=0;k<BPGROUP.size();++k) p1[k]=BPGROUP[k].P2;
    arrayindex(p1,p1_idx,1);
    BPTMP.resize(BPGROUP.size());
    for(k=0;k<BPGROUP.size();++k) BPTMP[k]=BPGROUP[p1_idx[k]];
    BPGROUP=BPTMP;
    
    BP.push_back(BPGROUP[0]);
    for(k=1;k<BPGROUP.size();++k) {
      if ( abs(BPGROUP[k].P2 - BP.back().P2) < 3 ) {
	BP.back().count+=BPGROUP[k].count;
	if ( BPGROUP[k].CL > BP.back().CL ) {
	  BPGROUP[k].count=BP.back().count;
	  BP.back()=BPGROUP[k];
	}
	continue;
      }
      BP.push_back(BPGROUP[k]);
    }
  }
  
  //! merge according to the p1 positions
  //! allowed distance is 3. also check SVSEG
  if ( BP.size()>1 ) {
    BPTMP.clear();
    BPTMP.push_back(BP[0]);
    for(i=1;i<BP.size();++i) {
      if ( abs(BP[i].P1-BPTMP.back().P1) < 3 && 
	   abs(BP[i].P2-BPTMP.back().P2) < 3 ) {
	string SVSEQ1=BPTMP.back().SVSEQ;
	string SVSEQ2=BP[i].SVSEQ;
	//if ( edit_distance(SVSEQ2,SVSEQ1) <= 2*(BP[i].P1-BPTMP.back().P1) ) {
	if ( is_same_svseq(SVSEQ2,SVSEQ1)  ) {
	  int tmpcount=BP[i].count+BPTMP.back().count;
	  if ( BP[i].ED < BPTMP.back().ED ) BPTMP.back()=BP[i];
	  BPTMP.back().count=tmpcount;
	  continue;
	}
	else BPTMP.push_back(BP[i]);
      }
      else BPTMP.push_back(BP[i]);
    }
    BP=BPTMP;
  }
  BPTMP.clear();
  
  int count=0;
  for(i=0;i<BP.size();++i) {
    count+=BP[i].count;
    if ( i>0 ) if ( BP[i].P1<BP[i-1].P1 ) {
	cerr << "order error BP\n"; 
	exit(0); 
      }
  }
  if ( count != count0 ) {
    cerr << "[compact_breakpoint] Error: Not all variations are counted\n"
	 << "expecting " << count0 << " found " << count << endl;
    for(i=0;i<BPSEEK.size();++i) 
      cerr << BPSEEK[i].RNAME()  << "\t" 
	   << BPSEEK[i].P1 << "\t" << BPSEEK[i].P2 << "\t"
	   << BPSEEK[i].SVSEQ << "\t" 
	   << BPSEEK[i].count << "\n";
    cerr << "----------" << endl;
    for(i=0;i<BP.size();++i) 
      cerr << BP[i].RNAME() << "\t" 
	   << BP[i].P1 << "\t" << BP[i].P2 << "\t" 
	   << BP[i].SVSEQ << "\t"
	   << BP[i].count << "\n";
    cerr << "----------" << endl;
  }
  
  SVSEQ_filter( BP );  
  return;
}

//! recycle buffer
void shift_buffer(string& SEQBUFFER, vector<RSAI_st>& MSREAD, vector<RSAI_st>& SMREAD)
  
{
  size_t i;
  
  int last_pos=max(SMREAD.back().pos, MSREAD.back().pos);
  int last_cutoff=last_pos-MAXDISTANCE;
  size_t last_char_ms=0;
  size_t last_char_sm=0;
  
  if ( MSREAD.size() > 0 ) {
    i=MSREAD.size()-1;
    while(MSREAD[i].pos>last_cutoff) {
      if (i==0) break;
      i--;
    }
    
    if ( i==MSREAD.size()-1 ) {
      last_char_ms=MSREAD[i].p1+MSREAD[i].len+1;
      MSREAD.clear();
    }
    else if (MSREAD[i].pos<=last_cutoff) {
      last_char_ms=MSREAD[i+1].p1;
      MSREAD.erase(MSREAD.begin(),MSREAD.begin()+i+1);
    }
  }
  
  if ( SMREAD.size() > 0 ) {
    i=SMREAD.size()-1;
    while(SMREAD[i].pos>last_cutoff) {
      if (i==0) break;
      i--;
    }
    
    if ( i==SMREAD.size()-1 ) {
      last_char_sm=SMREAD[i].p1+SMREAD[i].len+1;
      SMREAD.clear();
    }
    else if(SMREAD[i].pos<=last_cutoff) {
      last_char_sm=SMREAD[i+1].p1;
      SMREAD.erase(SMREAD.begin(),SMREAD.begin()+i+1);
    }
  }
  
  size_t last_char=max(last_char_ms,last_char_sm);
  if ( last_char>SEQBUFFER.size() ) last_char=SEQBUFFER.size();
  SEQBUFFER.erase(0, last_char);
  for(i=0;i<MSREAD.size();++i) MSREAD[i].p1-=last_char;
  for(i=0;i<SMREAD.size();++i) SMREAD[i].p1-=last_char;
  
  if ( VERBOSE ) {
    cerr << "buffers shifted down to " 
	 << MSREAD.size() << "\t" 
	 << SMREAD.size() << "\t" 
	 << SEQBUFFER.size() << endl;
  }
  
  return;
}


void reduce_readsidx(string& BUFFER, vector<RSAI_st>& READIDX)
{
  int maxMismatch=1;
  
  if ( READIDX.size()<20000 ) return;
  size_t i,j,k;
  
  vector<bool> iskept(READIDX.size(), false);
  
  for(i=0;i<READIDX.size();++i) {
    int nOverlap=0;
    for(k=i+1;k<READIDX.size();++k) {
      int com_beg=max(READIDX[i].sbeg, READIDX[k].sbeg);
      int com_end=min(READIDX[i].send, READIDX[k].send);
      if ( com_beg<com_end ) nOverlap+=1;
      if ( READIDX[k].sbeg >  READIDX[i].send+READIDX[i].len ) break;
    }
    if ( nOverlap==0 ) iskept[i]=true;
  }
  
  size_t ipre=0;
  for(i=0;i<READIDX.size();++i) {
    string SEQi=BUFFER.substr(READIDX[i].p1,READIDX[i].len);
    string CLIPi= READIDX[i].sbeg<READIDX[i].pos ? 
      SEQi.substr(0,READIDX[i].S) :
      SEQi.substr(READIDX[i].len-READIDX[i].S) ;
    ipre=i;
    //int nOverlap=0;
    for(k=i+1;k<READIDX.size();++k) {
      int com_beg=max(READIDX[i].sbeg, READIDX[k].sbeg);
      int com_end=min(READIDX[i].send, READIDX[k].send);
      if ( com_beg<com_end ) {
	string SEQk=BUFFER.substr(READIDX[k].p1,READIDX[k].len);
	string CLIPk= READIDX[k].sbeg<READIDX[k].pos ? 
	  SEQk.substr(0,READIDX[k].S) :
	  SEQk.substr(READIDX[k].len-READIDX[k].S) ;
	
	string comi=CLIPi.substr(com_beg-READIDX[i].sbeg, com_end-com_beg+1);
	string comk=CLIPk.substr(com_beg-READIDX[k].sbeg, com_end-com_beg+1);
	int ndiff=0;
	for(j=0;j<comi.size();++j) {
	  if ( comi[j]!=comk[j] ) ndiff++;
	  if ( ndiff>maxMismatch ) break;
	}
	if ( ndiff<=maxMismatch ) {
	  iskept[i]=true;
	  iskept[k]=true;
	  ipre=k;
	}
      }
      if ( READIDX[k].sbeg >  READIDX[i].send+READIDX[i].len ) break;
    }
    i=ipre;
  }
  
  vector<RSAI_st> READIDX0(0);
  for(i=0;i<READIDX.size();++i) 
    if ( iskept[i] ) READIDX0.push_back( READIDX[i] );
  if ( READIDX.size() != READIDX0.size() ) 
    cerr << "#read compact:\t" << READIDX.size() << "\t" << READIDX0.size() << endl;
  READIDX=READIDX0;
  return;
}


//! main calculation block
//! start multi thread 
void process_data_block( vector<BREAKSEEK_st>& BPSEEK,
			 string& SEQBUFFER,
			 vector<RSAI_st>& MSREAD, vector<RSAI_st>& SMREAD,
			 string& FASTA, string& RNAME)
{
  if ( COMPACT ) {
    // reduce_readsidx(MSBUFFER, MSREAD);
    // reduce_readsidx(SMBUFFER, SMREAD);
  }
  
  size_t i;
  for(i=0;i<(size_t) NUM_THREADS;++i) {
    thread_data[i].thread_id=i;
    thread_data[i].NUM_THREADS=NUM_THREADS;
    thread_data[i].SEQBUFFER = &SEQBUFFER;
    thread_data[i].MSREAD = &MSREAD;
    thread_data[i].SMREAD = &SMREAD;
    thread_data[i].FASTA = &FASTA;
    thread_data[i].RNAME = &RNAME;
    thread_data[i].BPSEEK = &BPSEEK;
    int rc = pthread_create(&threads[i], 
			    &attr, 
			    matching_reads_thread_wrapper, 
			    &thread_data[i] );
    if (rc) {
      cerr << "ERROR; return code from pthread_create() is " << rc << endl; 
      exit(-1);
    }
  }
  for (i=0; i<(size_t)NUM_THREADS; i++) pthread_join(threads[i], NULL);
  //! wait for threads to stop
  
  vector<BREAKSEEK_st> BPSEEK0;
  compact_breakpoint(BPSEEK, BPSEEK0);
  BPSEEK.clear();
  BPSEEK=BPSEEK0;
  return;
}

void write_cnv_to_file(vector<BREAKSEEK_st>& BPSEEK)
{
  if ( ! HEADER ) {   // has header been printed?
    cout << BP_format(BPSEEK[0],1) << "\n";
    HEADER=true;
  }
  for(size_t i=0;i<BPSEEK.size();++i) cout << BP_format(BPSEEK[i],0) << "\n";
  
  return;
}

void load_cnv_from_file(vector<BREAKSEEK_st>& BPSEEK, string cnvFile)
{
  size_t i;
  
  BPSEEK.clear();
  BREAKSEEK_st ABP;
  ifstream FIN(cnvFile.c_str());
  if ( !FIN ) { cerr << "Can't open file " << cnvFile << endl; exit(0); }
  
  while ( !FIN.eof() ) {
    string chr,tmps,comment;
    getline(FIN,tmps);
    if ( tmps.size()<3 ) continue;
    if ( tmps[0]=='#' ) continue;
    
    istringstream iss(tmps);
    
    iss >> chr
	>> ABP.P1
	>> ABP.P2
	>> ABP.T
	>> ABP.UN
	>> ABP.len
	>> ABP.INSEQ
	>> comment
	>> ABP.count
	>> ABP.isit
	>> ABP.SVSEQ
	>> ABP.MERGE;
    
    if ( ABP.RNAMEINT.count(chr)==0 ) 
      cerr << "#" << chr << " not found in bam file" << endl;
    ABP.tid=ABP.RNAMEINT[chr];
    
    for(i=0;i<comment.size();++i) {
      if ( comment[i]==',' ) comment[i]='\t';
      if ( comment[i]=='C' ) comment[i]=' ';
      if ( comment[i]=='L' ) comment[i]=' ';
      if ( comment[i]=='E' ) comment[i]=' ';
      if ( comment[i]=='D' ) comment[i]=' ';
      if ( comment[i]=='S' ) comment[i]=' ';
      if ( comment[i]=='N' ) comment[i]=' ';
      if ( comment[i]=='Q' ) comment[i]=' ';
    }
    iss.clear();
    iss.str(comment);
    
    iss >> ABP.CL
	>> ABP.ED
	>> ABP.L1
	>> ABP.S1
	>> ABP.n1
	>> ABP.q1
	>> ABP.L2
	>> ABP.S2
	>> ABP.n2
	>> ABP.q2;
    
    BPSEEK.push_back(ABP);
  }
  
  // cerr << "laoded " << BPSEEK.size() << " cnv" << endl;
  FIN.close();
  return;
}


void filters(vector<BREAKSEEK_st>& BP, string bamFile)
{
  
  bam_pileup_api_wrapper_t tmp;  
  tmp.in = NULL;
  tmp.bamidx = NULL;
  
  if ( is_binary(bamFile) ) {
    tmp.in = samopen(bamFile.c_str(), "rb", 0);
    tmp.bamidx = bam_index_load(bamFile.c_str()); // load BAM index
    if (tmp.bamidx == 0) {  
      cerr << "BAM indexing file for " << bamFile << " is not available." << endl ;  
      return ;  
    }  
  }
  else return;
  
  // cerr << "Filtering " << BP.size() << " cnvs" << endl;
  
  size_t i, k;
  string bamRegion="";
  
  for(i=0;i<BP.size();++i) {
    size_t dp1=0,dp2=0,din=0,dout=0;
    int dx=max(BP[i].L1,BP[i].L2)*2;
    if ( dx<(int)abs(BP[i].P2-BP[i].P1) ) dx=(int)abs(BP[i].P2-BP[i].P1);
    if ( dx>400) dx=400;
    int p1=min(BP[i].P1,BP[i].P2) -dx;
    int p2=max(BP[i].P1,BP[i].P2) +dx;
    if ( p1<1 ) p1=1;
    int bp1=BP[i].P1+BP[i].UN;
    int bp2=BP[i].P2;
    if ( BP[i].UN > BP[i].len/2 ) {
      bp1=BP[i].P1;
      bp2=BP[i].P2;
    }
    int Lout=bp1-p1-1 + p2-bp2-1;
    int Lin=bp2-bp1-1;
    if ( Lin<1 ) Lin=1;
    
    int klen=kmerlength(BP[i].MERGE);
    
    //if ( BP[i].count>10 && klen<40 ) goto FILTERDONE;
    //if ( BP[i].count==1 && BP[i].CL<30 && BP[i].len>100000 ) BP[i].isit=false;
    //if ( BP[i].isit==false ) goto FILTERDONE;
    
    tmp.T=BP[i].T;
    tmp.bp1=BP[i].P1;
    tmp.bp2=BP[i].P2;
    if ( tmp.bp1 > tmp.bp2 ) swap(tmp.bp1,tmp.bp2);
    bampileup(tmp, BP[i].RNAME(), p1, p2);
    
    for(k=p1,dout=0; (int)k<bp1; ++k) dout+=tmp.n[k-p1];
    for(k=bp2+1; (int)k<=p2; ++k) dout+=tmp.n[k-p1];
    for(k=bp1+1, din=0; (int)k<BP[i].P2; ++k) din+=tmp.n[k-p1];
    dp1=tmp.n[bp1-p1];
    dp2=tmp.n[bp2-p1];
    
    dout /= Lout;
    din  /=  Lin;
    
    float q0=(float)tmp.c0/(float)(tmp.c0+tmp.c1);
    
    //if ( BP[i].count==1 && BP[i].CL<30 ) BP[i].isit=false;
    //if ( BP[i].count==2 && BP[i].CL<30 ) BP[i].isit=false;
    
    if ( BP[i].count==1 || klen>40 ) {
      if ( BP[i].T=='D' && din>=dout*3/4 ) BP[i].isit=false;
      if ( BP[i].T=='A' && din<=dout*5/4 ) BP[i].isit=false;
    }
    
    if ( BP[i].count==2 ) {
      if ( BP[i].T=='D' && din>=dout ) BP[i].isit=false;
      if ( BP[i].T=='A' && din<=dout ) BP[i].isit=false;
    }
    
    if (dout>30 || din>30)  {
      if ( BP[i].T=='D' && din>=dout ) BP[i].isit=false;
      if ( BP[i].T=='A' && din<=dout ) BP[i].isit=false;
    }
    if (dout>50 || din>50)  {
      if ( BP[i].T=='D' && din>=dout*3/4 ) BP[i].isit=false;
      if ( BP[i].T=='A' && din<=dout*5/4 ) BP[i].isit=false;
    }

    if ( klen > 33 && BP[i].count<5 ) BP[i].isit=false;
    
    //FILTERDONE:    
    BP[i].TAG="DO="+to_string(dout)+","
      +"DI="+to_string(din)+","
      +"D1="+to_string(dp1)+","
      +"D2="+to_string(dp2)+";"
      +"PR="+to_string(tmp.pair)+";"
      +"KM="+to_string(klen)+";"
      +"Q0="+to_string((int)(q0*100+0.499))+"%";
    cerr << BP[i].RNAME() << "\t" 
	 << BP[i].P1 << "\t"
	 << BP[i].P2 << "\t"
 	 << BP[i].T << "\t" 
 	 << BP[i].UN << "\t" 
 	 << BP[i].len << "\t" 
	 << (float)tmp.c0/(float)(tmp.c0+tmp.c1) << "\t"
      //<< dout << "," 
      //<< din << ","
      //<< dp1 << ","  
      //<< dp2 << "\t"
	 << BP[i].TAG << "\t"
	 << BP[i].count << "\t"
	 << klen << "\t" 
	 << BP[i].isit << endl; 
    
  }
  
  samclose(tmp.in);
  bam_index_destroy(tmp.bamidx);  
  
  return;
}


int usage_match_MS_SM_reads(int argc, char* argv[]) {
  cerr << "Usage:\n" 
       << "  " << argv[0] << " <options> -f REFFILE -b BAMFILE [REGION]\n"
       << "\nOptions:\n"
       << "  -t  INT  number of threads [1,36], INT=1 \n"
       << "  -ee INT  maximum allowed mismatches at 3' edges, INT=1 \n"
       << "  -me INT  max allowed mismatches when matching strings, INT=1 \n"
       << "  -s  INT  minimum number of soft clipped bases, INT=10 \n"
       << "  -q  INT  minimum mapping score, INT=10 \n"
       << "  -Q  INT  minimum base read quality, INT=0 \n"
       << "  -2       use secondary alignment \n"
       << "  -L  INT  check reads before and after INT bases, INT=1E6 \n"
       << "  -nc      no colored display for intermediate info, colored=true \n"
       << "  -o  STR  outputfile, STR=STDOUT \n"
       << "  -dump    only dump reads with S part larger than threshold(-s) \n"
       << "  -v       print intermediate matching information \n"
       << "   REGION  if given should be in samtools's region format \n"
       << "\nExamples:\n"
       << "  " << argv[0] << " -f human_g1k_v37.fasta -o bwap40X.bam.bp -b bwap40X.bam \n"
       << "  " << argv[0] << " -f human_g1k_v37.fasta -o bwap40X.bam.bp -q 30 -b bwap40X.bam chr19 -t 4 \n"
       << "\nNote:\n"
       << "  REGION must be understood by samtools. But, internally, chr1 and 1 are same.\n"
       << "  The default setting should work for most current platforms as of 2013\n"
       << "  If reads are longer than 200, you may want to use -me 2\n"
       << "  If reads qualities are poor at 3' end, you may try -ee 2 or more.\n"
       << "  Or, truncate reads before mapping.\n"
       << "  If you got nothing, run with -v to check mapq and S number stats.\n"
       << endl;
  
  return(0);
}

void match_MS_SM_reads(int argc, char* argv[])
{
  if ( argc<3 )  exit( usage_match_MS_SM_reads(argc, argv) );
  
  string mycommand="";
  size_t i, k, count=0;
  string QNAME,RNAME,POS,MAPQ,CIGAR,MRNM,MPOS,ISIZE,SEQ,QUAL,OPT;
  string QSEQ;
  string CIGAR1,SEQ1,QUAL1,OPT1;
  string REFSEQ=string(200,'C');
  string REFCLIP=string(200,'C');
  string CLIPPEDSEQ=string(200,'C');
  string CLIPPEDQUAL=string(200,'C');
  string CLIPPEDCIGAR=string(200,'C');
  string samtoolsOpt1="", samtoolsOpt2="", samtoolsCommand="";
  vector<size_t>mismatchCount(10000,0);
  vector<size_t>mapqCount(10000,0);
  string oper;
  //FILE *fp;
  bool NOSECONDARY=true;
  string bamFile="",bamRegion="",fastaFile="",outputFile="STDOUT";
  string cnvFile="";
  vector<string> inputArgv;
  vector<BREAKSEEK_st> BPSEEK(0);
  vector<BREAKSEEK_st> BPALL(0);
  vector<BREAKSEEK_st> BPSEEK0(0);
  
  //  string S1="ATTTCCCAGGTGCCGGTCCATCCTTGTTGT";
  //  string S2="ATTTCCCAGGTGCCGGTCCATCCTTGTTGT";
  //  cout << edit_distance(S1,S2) << endl;
  //  cout << is_same_svseq(S1,S2) << endl;
  //  exit(0);
  
  for(i=0;i<(size_t) argc;++i) inputArgv.push_back(string(argv[i]));
  mycommand=inputArgv[0];
  for(i=1;i<inputArgv.size();++i) mycommand+=" "+inputArgv[i];
  for(i=1;i<inputArgv.size();++i) {
    if ( inputArgv[i]=="-b" ) {  // input bam file
      bamFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
      k=i+2;
      if ( k<inputArgv.size() ) {
	if ( inputArgv[k][0] != '-' ) {
	  bamRegion=inputArgv[k];
	  inputArgv[k]="";
	}
      }
    }
    if ( inputArgv[i]=="-f" ) {  // input reference file
      fastaFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-o" ) {  // output file
      outputFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    }
    if ( inputArgv[i]=="-t" ) {   // number of threads
      NUM_THREADS=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-ee" ) {   // number of bases at edge that can be ignored
      ERR_EDGE=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-me" ) {   // err_rate when fuzzy matching
      ERR_MATCH=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-s" ) {   // number of softclipped bases
      MINIMUMS=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-q" ) {   // min mapq score
      MINIMUMMAPQ=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-Q" ) {   // min mapq score
      MINIMUMBASEQ=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-2" ) {   // min mapq score
      NOSECONDARY=false;
      inputArgv[i]="";
    } 
    if ( inputArgv[i]=="-L" ) {   // max distance between two reads
      MAXDISTANCE=atoi(inputArgv[i+1].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-nc" ) {   // no colored output
      COLORED=false;
      inputArgv[i]="";
    } 
    if ( inputArgv[i]=="-v" ) {   // print intenediate results
      VERBOSE=true;
      inputArgv[i]="";
    } 
    if ( inputArgv[i]=="-dump" ) {   // only dump reads
      DUMPBAM=true;
      inputArgv[i]="";
    } 
    if ( inputArgv[i]=="-cnv" ) {   // load cnv file
      cnvFile=inputArgv[i+1];
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
    if ( inputArgv[i]=="-noc" ) {   // don't compact the reads
      COMPACT=false;
      inputArgv[i]="";
    } 
    if ( inputArgv[i]=="-nof" ) {   // don't compact the reads
      RDFILTER=false;
      inputArgv[i]="";
    } 
    if ( inputArgv[i]=="-ce" ) {   // when mismatch when compacting reads
      maxMismatch=atoi(inputArgv[i].c_str());
      inputArgv[i]="";
      inputArgv[i+1]="";
    } 
  }
  if ( bamFile=="" ) {
    cerr << "Need bam file\n";
    exit( usage_match_MS_SM_reads(argc, argv) );
  }
  if ( fastaFile=="" ) {
    cerr << "Need reference file\n";
    exit( usage_match_MS_SM_reads(argc, argv) );
  }
  for(i=1;i<inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) cerr << "unknown argument:\t" << inputArgv[i] << endl;
  for(i=1;i<inputArgv.size();++i) 
    if ( inputArgv[i]!="" ) exit( usage_match_MS_SM_reads(argc, argv) );
  
  
  streambuf* sbuf = cout.rdbuf();
  ofstream FOUT;
  if ( outputFile != "STDOUT" && !DUMPBAM ) {
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
  }
  
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  string comment="#Command Line     : " + mycommand + "\n" +
    "#Input bamfile    : " + bamFile + "\n" +
    "#Bamfile region   : " + bamRegion + "\n" +
    "#FASTAfile        : " + fastaFile + "\n" +
    "#Num of S bases   : " + to_string(MINIMUMS) + "\n" +
    "#Minumum overlap  : " + to_string(MINIMUMOVERLAP) + "\n" +
    "#Minimum mapq     : " + to_string(MINIMUMMAPQ) + "\n" +
    "#3' bases ignored : " + to_string(ERR_EDGE) + "\n" +
    "#Allowed mismatch : " + to_string(ERR_MATCH) + "\n" +
    "#Maximum distance : " + to_string(MAXDISTANCE) + "\n" +
    "#Output           : " + outputFile + "\n";
  
  cerr << comment << endl;
  if ( !DUMPBAM && VERBOSE ) cout <<  comment << endl;
  
  vector<RSAI_st> MSREAD(0);         // info used to find the read
  vector<RSAI_st> SMREAD(0);
  RSAI_st iread;
  string FASTA;
  string fastaname="";
  
  SEQBUFFER.reserve(BUFFERSIZE+1); // buffer for reads with MS type of CIGAR
  BREAKSEEK_st::TARGET.clear();
  
  //is_binary(bamFile); exit(0);
  
  samfile_t *fp_in = NULL;
  samfile_t *fp_out = NULL;
  bam1_t *b=NULL;   
  bam_index_t *bamidx=NULL;
  bam_iter_t iter=0;
  b = bam_init1();
  int ref=0, beg=0, end=0x7fffffff;
  if ( is_binary(bamFile) ) {
    fp_in = samopen(bamFile.c_str(), "rb", 0);
    bamidx = bam_index_load(bamFile.c_str()); // load BAM index
    for(int i=0;i<fp_in->header->n_targets;++i) {
      BREAKSEEK_st::TARGET.push_back( string(fp_in->header->target_name[i]) );
      BREAKSEEK_st::RNAMEINT[ string(fp_in->header->target_name[i]) ]=i;
    }
    if ( BREAKSEEK_st::TARGET.size() == 0 ) 
      cerr << "#no header in bam file? " << bamFile << endl;
  }
  
  vector<string> REGIONS(0);
  if ( bamRegion!="" ) {
    REGIONS.clear();
    REGIONS.push_back(bamRegion);
  }
  else REGIONS=BREAKSEEK_st::TARGET;
  
  if ( cnvFile!="" ) goto TESTCNV;
  
  //  for(size_t iRegion=0;iRegion<REGIONS.size();++iRegion) {
  //    bamRegion=REGIONS[iRegion];
  //    preprocess_bam_file(bamFile, bamRegion, 
  //			MINIMUMS, MINIMUMMAPQ, MINIMUMBASEQ, fastaFile);
  //  }
  //  exit(0);
  
  
  if ( DUMPBAM ) {
    fp_out = outputFile=="STDOUT" ? 
      samopen("-", "w", fp_in->header) :
      samopen(outputFile.c_str(), "wb", fp_in->header) ;
  }
  
  bam_parse_region(fp_in->header, bamRegion.c_str(), &ref, &beg, &end); 
  iter = bam_iter_query(bamidx, ref, beg, end);
  
  while( bamread(fp_in, iter, b, bamRegion) > 0 ) {
    count++;
    if ( count%1000000==0 ) cerr << "#processed " << count << " reads" << endl;
    
    if ( (int)b->core.n_cigar <=1 ) continue;
    if ( (int)b->core.qual < MINIMUMMAPQ ) continue;
    if ( (int)b->core.tid < 0 ) continue;
    if ( b->core.mtid != b->core.tid && b->core.mtid>0 ) continue;
    //if ( NOSECONDARY && (b->core.flag & 0x100)>0 ) continue;
    
    POSCIGAR_st bm;
    resolve_cigar_pos(b, bm);  

    // S part before the start of reference
    if ( bm.cop[0]+bm.nop[0] < bm.cop[0] ) continue;
    
    if ( bm.pos==0 ) continue;
    if ( bm.iclip<0 ) continue;
    if ( (int)bm.nop[bm.iclip] < MINIMUMS ) continue;
    
    uint8_t *t = bam1_qual(b);
    // minumum base quality
    if ( MINIMUMBASEQ>1 && t[0] != 0xff) {
      int bQ=0, bQCount=0;
      for(i=0;i<bm.nop[bm.iclip];++i) {
	bQ=t[ bm.qop[bm.iclip]+i ];
	if ( bQ < MINIMUMBASEQ ) ++bQCount;
      }
      if ( bQCount > (int)bm.nop[bm.iclip]/3 ) continue;
    }
    
    get_rname(fp_in->header, b, RNAME);
    if ( RNAME != BREAKSEEK_st::TARGET[b->core.tid] ) 
      cerr << "#TARGET read error" << endl;
    get_qseq(b, SEQ);  
    get_cigar(b, CIGAR);  
    
    /* load reference if necessary and process data for previous chromosome*/
    if ( RNAME != fastaname ) {
      if ( MSREAD.size()>0 && SMREAD.size()>0 && FASTA.size()>1 ) {
	cerr << "#" << fastaname << "\t" << MSREAD.size() << "\t" << SMREAD.size() << endl;
	process_data_block(BPSEEK, SEQBUFFER, MSREAD, SMREAD, FASTA, fastaname);
	shift_buffer(SEQBUFFER, MSREAD, SMREAD);
	write_cnv_to_file(BPSEEK);
	if ( BPSEEK.size()>0 ) BPALL.insert(BPALL.end(),BPSEEK.begin(),BPSEEK.end());
	cerr << "#" << fastaname << "\t" << BPSEEK.size() << endl;
      }
      
      string RNAMEtmp=RNAME;
      if ( ci_find(RNAME,"chr") != string::npos ) RNAMEtmp=RNAME.substr(3); 
      read_fasta(fastaFile,RNAMEtmp,FASTA);
      fastaname=RNAME;
      cerr << "#" << fastaFile << "\t" 
	   << fastaname << "\t" 
	   << FASTA.size() << " loaded" << endl;
      
      if ( BREAKSEEK_st::RNAMEINT.count(RNAME) == 0 ) {
	BREAKSEEK_st::RNAMEINT[RNAME] = BREAKSEEK_st::TARGET.size();
	BREAKSEEK_st::TARGET.push_back( RNAME );
      }
      
      // reset for the new chromosome
      BREAKSEEK_st::ctid = BREAKSEEK_st::RNAMEINT[RNAME];      
      BPSEEK.clear();
      MSREAD.clear();
      SMREAD.clear();
      SEQBUFFER="";
    }
    if ( FASTA.size()<1 ) continue;
    
    //int nAdjust=calibrate_resolved_cigar_pos(FASTA, SEQ, bm);  
    //if ( VERBOSE && nAdjust>50 ) cerr << "nAdjust=" << nAdjust << endl;
    // S part beyond reference
    // cerr << bm.pos << "\t" << CIGAR << endl;
    if ( bm.cop.back()+bm.nop.back() >= FASTA.size() ) continue;
    if ( bm.cop[0]<1 ) continue;
    calibrate_resolved_cigar_pos(FASTA, SEQ, bm);  
    
    if ( bm.pos==0 ) continue;
    if ( bm.iclip<0 ) continue;
    if ( (int)bm.nop[bm.iclip] < MINIMUMS ) continue;
    if ( bm.nop[bm.iclip]*1.25 > bm.l_qseq ) continue;
    //if ( bm.nop[bm.iclip]*2 > bm.l_qseq ) continue;
    
    int ndiff_m=0;
    for (k = 0; k < bm.op.size(); ++k) {
      if ( bm.op[k]!=BAM_CMATCH && bm.op[k]!=BAM_CEQUAL ) continue; 
      string CLIPPEDSEQ=SEQ.substr(bm.qop[k], bm.nop[k]);
      string REFCLIP=FASTA.substr(bm.cop[k]-1, bm.nop[k]);
      for(i=0; i<REFCLIP.size(); ++i) if ( REFCLIP[i]!=CLIPPEDSEQ[i] ) ++ndiff_m;
    }
    int ndiff_s=0, nN=0;
    if ( bm.iclip>=0 ) {
      k = bm.iclip;
      string CLIPPEDSEQ=SEQ.substr(bm.qop[k], bm.nop[k]);
      string REFCLIP=FASTA.substr(bm.cop[k]-1, bm.nop[k]);
      for(i=0; i<REFCLIP.size(); ++i) {
	if ( REFCLIP[i]!=CLIPPEDSEQ[i] ) ++ndiff_s;
	if ( CLIPPEDSEQ[i]=='N' ) ++nN;
      }
    }
    if ( ndiff_s <= 2 ) continue;    
    if ( ndiff_s <= (int)bm.nop[bm.iclip]/4 ) continue;
    if ( ndiff_m >= (int)bm.l_qseq*8/100  ) continue;
    if ( nN >= MINIMUMS/2 ) continue;
    int nS=0,nIndel=0;
    for(i=0;i<bm.op.size();++i) {
      if (bm.op[i]==BAM_CSOFT_CLIP ) ++nS;
      if (bm.op[i]==BAM_CINS || bm.op[i]==BAM_CDEL || bm.op[i]==BAM_CREF_SKIP || bm.op[i]==BAM_CPAD ) ++nIndel;
    }
    if ( nS>3 || nIndel>3 ) continue;    // CIGAR too complicated
    int Snum_adjust=bm.nop[bm.iclip];
    int Sotherend=0;
    if ( bm.iclip==0 && bm.op.back()==BAM_CSOFT_CLIP ) Sotherend=bm.nop.back();
    if ( bm.iclip>0  && bm.op[0]==BAM_CSOFT_CLIP ) Sotherend=bm.nop[0];
    //if ( Sotherend > SEQ.size()*0.08 ) continue;
    
    mismatchCount[ndiff_s]++;
    mapqCount[b->core.qual]++;
    
    if ( DUMPBAM ) {
      bam_aux_append(b, "ns", 'i', 4, (uint8_t*)&Snum_adjust);
      bam_aux_append(b, "nm", 'i', 4, (uint8_t*)&ndiff_m);
      samwrite(fp_out, b);
      continue;
    }
    
    iread.tid=BREAKSEEK_st::ctid;
    iread.pos=bm.pos;
    iread.q1=b->core.qual;
    iread.len=bm.l_qseq;
    iread.M=bm.nop[bm.anchor];
    iread.Mrpos=bm.qop[bm.anchor];
    iread.S=bm.nop[bm.iclip];;
    iread.Spos=bm.cop[bm.iclip];
    iread.sbeg=bm.cop[bm.iclip];
    iread.send=bm.cop[bm.iclip]+bm.nop[bm.iclip]-1;
    iread.mm=ndiff_s;
    iread.p1=SEQBUFFER.size();
    SEQBUFFER+=SEQ+"\t";
    if ( bm.iclip > bm.anchor ) MSREAD.push_back(iread); // type M...S
    else  SMREAD.push_back(iread);  // type S...M
    
    if ( SEQBUFFER.size() > BUFFERSIZE-5*SEQ.length() ) {
      cerr << "#" << fastaname << "\t" << MSREAD.size() << "\t" << SMREAD.size() << endl;
      process_data_block(BPSEEK, SEQBUFFER, MSREAD, SMREAD, FASTA, RNAME);
      shift_buffer(SEQBUFFER, MSREAD, SMREAD);
    }
    
  }
  samclose(fp_in);
  bam_destroy1(b);
  bam_iter_destroy(iter);
  bam_index_destroy(bamidx);
  if ( fp_out ) samclose(fp_out);
  if ( DUMPBAM ) {
    if ( outputFile!="" ) bam_index_build(outputFile.c_str());
    cerr << "#[matchclips] exit\n";
    cerr << comment << endl;
    return;
  }
  
  if ( MSREAD.size()>0 && SMREAD.size()>0 && FASTA.size()>5 ) {
    cerr << "#" << fastaname << "\t" << MSREAD.size() << "\t" << SMREAD.size() << endl;
    process_data_block(BPSEEK, SEQBUFFER, MSREAD, SMREAD, FASTA, fastaname);
    shift_buffer(SEQBUFFER, MSREAD, SMREAD);
    write_cnv_to_file(BPSEEK);
    if ( BPSEEK.size()>0 )
      BPALL.insert(BPALL.end(),BPSEEK.begin(),BPSEEK.end());
    cerr << "#" << fastaname << "\t" << BPSEEK.size() << endl;
  }
  
 TESTCNV:
  if ( cnvFile!="" ) load_cnv_from_file(BPALL, cnvFile);
  //if ( RDFILTER ) filters(BPALL, bamFile);
  if ( RDFILTER ) RD_PAIR_filters(BPALL, bamFile);
  
  if ( outputFile != "STDOUT" ) {
    FOUT.close();
    remove(outputFile.c_str());
    FOUT.open(outputFile.c_str());
    cout.rdbuf(FOUT.rdbuf());
    HEADER=false;
    write_cnv_to_file(BPALL);
  }
  FOUT.close();
  cout.rdbuf(sbuf);
  
  if ( VERBOSE ) {
    cerr << "\n#Soft clipped Summary:" << endl;
    for(i=0;i<mismatchCount.size();++i) 
      if ( mismatchCount[i]>0 ) cerr << "#\t" << i << "\t" << mismatchCount[i] << "\n";
    cerr << endl;
    cerr << "\n#MAPQ Summary:" << endl;
    for(i=0;i<mapqCount.size();++i) 
      if ( mapqCount[i]>0 ) cerr << "#\t" << i << "\t" << mapqCount[i] << "\n";
    cerr << endl;
  }
  
  cerr << "#[matchclips] exit\n";
  cerr << comment << endl;
  
  return;
}

