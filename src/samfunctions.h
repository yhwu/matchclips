#ifndef _SAMFUNCTIONS_H
#define _SAMFUNCTIONS_H

/**** samtools headers ****/
#include <string>
#include <vector>
#include <bam.h>
#include <sam.h>

struct POSCIGAR_st {     
  int tid;
  int pos;              // 1-based position in chromosome
  int qual;
  int anchor;              // 1-based position in chromosome
  int iclip;
  unsigned int l_qseq;              // length of qseq
  vector<unsigned int> op;
  vector<unsigned int> nop;
  vector<unsigned int> cop;
  vector<unsigned int> qop;
  string cigar_s;
  string seq;
  POSCIGAR_st(): tid(-1),
		 pos(0),
		 qual(0),
		 anchor(0),
		 iclip(-1),
		 l_qseq(0),
		 op(0),
		 nop(0),
		 cop(0),
		 qop(0),
		 cigar_s(),
		 seq(){};
  
};

/*! 
  @abstract a wrapper for bam api
  
  @field  ref   reference id in target 0-based
  @field  beg   position on reference 0-based
  @field  end   position on reference 1-based
  @field  bp1   position of CNV, 1-based
  @field  bp2   position of CNV, 1-based, bp1<bp2
  @field  T     type of CNV
  @field  c0    number of reads with mapping quality <=1
  @field  c1    number of reads with mapping quality >1
  @field  pair  number of reads whose pairs cover bp1 and bp2
                for deletion, pairs cover [bp1,bp2]
		for duplication, pairs with negative insert size
  @field  in    bam input
  @field  bamidx bam index
  @field  buf   buffer
  @field  n     vector of depths
*/
typedef struct {  
  int ref, beg, end;  // for bam api, ref and bed are 0-based, end is 1-based
  int bp1, bp2;       // 1-based break point
  char T;
  size_t c0, c1;  
  size_t pair;  
  samfile_t *in;  
  bam_index_t *bamidx;
  bam_plbuf_t *buf;  
  std::vector<int> n;
} bam_pileup_api_wrapper_t;  

int bamread(samfile_t *fp_in, bam_iter_t& iter, bam1_t *b, std::string &bamRegion);
void resolve_cigar_pos(const bam1_t *b,  POSCIGAR_st& m);
void resolve_cigar_string(int POS, std::string& CIGAR, POSCIGAR_st& m);
int calibrate_resolved_cigar_pos(string& FASTA, string& SEQ, POSCIGAR_st& m);
void get_cigar(const bam1_t *b,  string& cigar);  
void get_qseq(const bam1_t *b,  std::string& seq);
void get_rname(const bam_header_t *header, const bam1_t *b,  std::string& rname);  
void bampileup(bam_pileup_api_wrapper_t& tmp, string RNAME, int p1, int p2);

#endif
