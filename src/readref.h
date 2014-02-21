#ifndef _READ_REF_H
#define _READ_REF_H

void read_fasta(string fastaFile, string chr, string& ref);
void get_N_regions(string& fasta, vector<int>& beg, vector<int>& end);

#endif
