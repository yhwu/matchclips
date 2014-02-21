#ifndef _MATCH_READS_H
#define _MATCH_READS_H

#include <map>

struct RSAI_st {        // read suffic array index
  int tid;              // target id as in TARGET
  int pos;              // 1-based position on reference
  size_t p1;          // 0-based position in memory buffer
  int q1;
  int len;              // length of read
  int len_cigar;        // length of CIGAR
  int M;                // number of M bases
  int Mrpos;            // position of first matched base in the read
  int S;                // number of S bases
  int Spos;             // 1-based pos
  int sbeg;             // 1-based pos begin of S, sbeg=Spos
  int send;             // 1-based pos end of S
  int pos_beg;          // 1-based pos of first base of read on reference
  int pos_end;          // 1-based pos of last base of read on reference
  int mm;               // mismatched wrt REF in S part of read
  RSAI_st():tid(-1),
	    pos(0),
	    p1(0),
	    q1(0),
	    len(0), 
	    len_cigar(0), 
	    M(0),
	    Mrpos(0),
	    S(0),
	    Spos(0),
	    pos_beg(0),
	    pos_end(0),
	    mm(0){};
};

struct BREAKSEEK_st {        // describes break point discovery
  static vector<string> TARGET;
  static map<string, int> RNAMEINT;
  static int ctid;      // current tid
  int tid;              // target id as in TARGET
  int P1;               // first break point position in chromosome
  int P2;               // second break point position in chromosome
  int L1;               // length of readMS
  int S1;               // clipped bases of readMS
  int n1;               // mismatched bases in clipped port of readMS
  int L2;               // length of readSM
  int S2;               // clipped bases of readSM
  int n2;               // mismatched bases in clipped port of readMS
  int q1;               // map quality of read1
  int q2;               // map quality of read2
  int CL;               // common string length of reads
  int len;              // length of variation
  int ED;               // edit_distance between matched and reference
  int ED5p;             // edit_distance between merged and reference at 5'
  int ED3p;             // edit_distance between merged and reference at 3'
  int UN;               // uncertainty of break point
  char T;               // D for DEL, A for ADD
  char SVSEQ[51];       // 30 chars befroe and after the SV, includes BP 
  int count;            // number of supportive pairs
  bool isit;            // number of supportive pairs
  string INSEQ;         // inserted sequence for Insertion only
  string MERGE;         // sequence of merged reads
  string TAG;           // comment on cnv
  string RNAME() { return ( TARGET[tid] ) ;}
  BREAKSEEK_st():tid(-1),
		 P1(0),
		 P2(0),
		 L1(0),
		 n1(0),
		 L2(0),
		 n2(0),
		 CL(0),
		 len(0),
		 ED(0),
		 UN(0),
		 count(0),
		 isit(true),
		 INSEQ("."),
		 MERGE(200,'C'),
		 TAG(30,'d'){};
};
//vector<string> BREAKSEEK_st::TARGET;
//map<string, int> BREAKSEEK_st::RNAMEINT;
//int BREAKSEEK_st::ctid;

struct thread_data_t {
  int thread_id;
  int NUM_THREADS;
  string* SEQBUFFER;
  vector<RSAI_st>* MSREAD;
  vector<RSAI_st>* SMREAD;
  string* FASTA;
  string* RNAME;
  vector<BREAKSEEK_st>* BPSEEK;
};

void match_MS_SM_reads(int argc, char* argv[]);

#endif
