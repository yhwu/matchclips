/*************************************************************************
BEGIN OF LICENSE
Copyright (c) Yinghua Wu and Hongzhe Li (rsicnv project).

This program is free software; you can redistribute and/or modify
the codes written by the authors under the terms of the GNU General 
Public License as published by the Free Software Foundation 
(www.fsf.org); either version 2 of the License, or (at your option) 
any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; See the GNU General Public License for 
more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
END OF LICENSE

CONTACT: wu_yinghua@hotmail.com; hongzhe@upenn.edu
*************************************************************************/
/**** system headers ****/
#include <iostream>
#include <time.h>
#include <cstdlib>
using namespace std;

#include "functions.h"
#include "matchreads.h"

int main(int argc, char* argv[])
{
  time_t begin_T,end_T;
  int pid = getpid();

  time(&begin_T);
  match_MS_SM_reads(argc, argv);  
  
  cerr << procpidstatus(pid,"VmPeak") ;
  time(&end_T);
  cerr << "#Time elapsed: " << difftime(end_T,begin_T) << " seconds\n";
  exit(0);
} 

