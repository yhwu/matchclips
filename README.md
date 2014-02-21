![The MATCHCLIP algorithm](https://dl.dropboxusercontent.com/u/39236968/match.svg)

## FREEZED
Please check the updated version [matchclips2](https://github.com/yhwu/matchclips2). It is much faster and reliable for large files.

This version of matchclips has been abandoned. The codes are left here for result reproduction purpose. The reason is that it is rather slow and too sensitive for real whole genome sequencing data. Real data is usually much noisier and has weird characteristics. What is really useful biologically to our collaborators are results reproducible in wet lab. So, matchclips is going fast and conservative by incorporating read pair and read depth information. The strategy is to find long variations using read pair information, and find shorter variations with reads matching. Our goal is to process high coverage(~70X) whole genome sequencing data such as on TCGA within hours. Please check the [matchclips2](https://github.com/yhwu/matchclips2) project.

## Synopsis
MATCHCLIPS detects the precise break points of CNVs through a fuzzy string matching algorithm using both CIGAR and POS information. In case the two break points of a CNV are in repeated regions and the break points are not unique, MATCHCLIPS reports the range where the break points can slide. A full description will be available after our paper is accepted. Nonetheless, the algorithm is shown in the illustration.

TAGCNV tags some properties to CNVs, such as number of heterozygous sites within a CNV region, read depths inside and outside of a CNV regions, and portion of ambiguous reads in a CNV region. These qualities may indicate the qualities of detected CNVs. For example, in a deletion region, there should not be many heterozygous sites and the read depth inside should be lower than outside. 

Please cite MATCHCLIPS by ```doi: 10.3389/fgene.2013.00157```.

Related work, [RSICNV](https://github.com/yhwu/rsicnv), read depth method based on robust statistics and machine learning.

## Code Example
```
matchclips                                                           #help
matchclips -f human_g1k_v37.fasta -o bwap40X.bam.bp -b bwap40X.bam   #default setting
matchclips -f human_g1k_v37.fasta -o bwap40X.bam.bp -q 30 -b bwap40X.bam chr1:1-1000000  #using reads with mapq>=30 in chr1:1-1000000 

tagcnv                                                               #help
tagcnv pair -v bwap40X.bam.bp -b bwap40X.bam                         #check evident pairs
tagcnv het  -v bwap40X.bam.bp -b bwap40X.bam -f human_g1k_v37.fasta  #check zygosity
tagcnv rd   -v bwap40X.bam.bp -b bwap40X.bam -q 10                   #check read depth
tagcnv q0   -v bwap40X.bam.bp -b bwap40X.bam                         #check map quality
tagcnv sub  -v bwap40X.bam.bp -r chr1:1-10000000                     #subset CNVs in chr1:1-1000000
tagcnv tab  -v bwap40X.bam.bp bt2p40X.bam.bp                         #tabulate CNV overlap
```

## Installation
If you have GIT installed, you can download and compile the codes with
```
git clone git@github.com:yhwu/matchclips.git
cd matchclips 
make
```
or you can download the files directly with
```
mkdir matchclips
cd matchclips
wget http://github.com/yhwu/matchclips/zipball/master -O matchclips.zip 
unzip matchclips.zip 
mv yhwu-matchclips-*/* .
rm -R yhwu-matchclips-*
make
```

Note: you will need g++ and gcc compilers and -pthread -lm -lz libs to compile the programs. These are installed by default in most cluster systems. Samtools is included in the download as matchclips needs samtools' API to read bam files. It is not needed after compilation. Tagcnv needs working samtools and bcftools binaries to read bam files and make zygosity calls. 

## Tests
This program requires a BAM file and the corresponding reference fasta file. The BAM file must be sorted and both must be indexed by samtools. Our program has been tested to work well with BAM files produced by bwa, bowtie2, and novoalign, either single or paired-end mode. It should work for other mapping tools too, as long as the reads are aligned with local mapping. The test BAM files are aligned against hg19 using paired-end bwa and bowtie2.  
```
#download hg19 reference if needed
wget ftp://ftp.sanger.ac.uk/pub/1000genomes/tk2/main_project_reference/human_g1k_v37.fasta.fai
wget ftp://ftp.sanger.ac.uk/pub/1000genomes/tk2/main_project_reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz

#download test data
wget http://dl.dropbox.com/u/39236968/matchclips_data.tar
tar xf matchclips_data.tar

#detect from sample data; should finish in less than 2 minutes
matchclips -f human_g1k_v37.fasta -b bwap40X.bam -o bwap40X.bam.bp
matchclips -f human_g1k_v37.fasta -b bt2p40X.bam -o bt2p40X.bam.bp 
matchclips -f human_g1k_v37.fasta -q 30 -b bwap40X.bam chr1

#check some properties of cnv
tagcnv pair -v bwap40X.bam.bp -b bwap40X.bam
tagcnv het  -v bwap40X.bam.bp -b bwap40X.bam -f human_g1k_v37.fasta
tagcnv rd   -v bwap40X.bam.bp -b bwap40X.bam -q 10
tagcnv q0   -v bwap40X.bam.bp -b bwap40X.bam
tagcnv sub  -v bwap40X.bam.bp -r chr1:1-10000000
tagcnv tab  -v bwap40X.bam.bp bt2p40X.bam.bp

#check accuracy with intended CNV set
checkbp.pl -p 3 NA12878hg19.txt bwap40X.bam.bp
checkbp.pl -p 3 NA12878hg19.txt bt2p40X.bam.bp
```

Note test BAM files can be generated with, e.g., 10X coverage,
```
simbp -f human_g1k_v37.fasta -cnv NA12878hg19.txt -X 10 -s $RANDOM -o tmp
bwa aln hs37m.fa tmp_1.fq > r1.sai && bwa aln hs37m.fa tmp_2.fq > r2.sai \  
    && bwa sampe hs37m r1.sai r2.sai tmp_1.fq tmp_2.fq > bwap10X.sam  
bowtie2 -x hs37m --local -X 650 -q -1 tmp_1.fq -2 tmp_2.fq -S bt2p10X.sam  
samtools view -uS bwap10X.sam | samtools sort - bwap10X
samtools index bwap10X.bam
samtools view -uS bt2p10X.sam | samtools sort - bt2p10X
samtools index bt2p10X.bam
```
or using other mapping software as described on the page http://lh3lh3.users.sourceforge.net/alnROC.shtml .

## Output fields
### matchclips output:
<table>
  <tr><td> 1 </td><td> RNAME </td><td> chromosome name </td></tr>
  <tr><td> 2 </td><td> START </td><td> start(lower) position, 1 based </td></tr>
  <tr><td> 3 </td><td> END </td><td>  end(higher) position, 1 based </td></tr>
<tr></tr>
  <tr><td> 4 </td><td> TYPE </td><td>  A for tandem duplication(bases from START to END are repeated) , D for deletion(bases from START+1 to END-1 are missing), I for insertion(bases from START+1 to END-1 are replaced by INSERT), J for insertion(INSERT is inserted between tandem duplication ) </td></tr>
  <tr><td> 5 </td><td> UNCERTAINTY </td><td>  break point pairs, (START+i,END+i) with 0&lt;=i&lt;=UNCERTAINTY, are equivalent </td></tr>
  <tr><td> 6 </td><td> LENGTH </td><td> length of variation, for deletion END-START-1, for duplication END-START+1, for insertion the length of INSERT </td></tr>

  <tr><td> 7 </td><td> INSERT </td><td> inserted sequence, dot if not available </td></tr>

  <tr><td> 8 </td><td> MATCHINFO </td><td> information for the matched reads with the longest overlap in the format "CL#,ED#,L#,S#,N#,Q#,L#,S#,N#,Q#", where CL is length of common sequence, ED is edit distance between the common sequence and the corresponding reference, L is length of read, S is the number of soft clipped bases in the read, N is the number of mismatched bases for the soft clipped bases, Q is mapping quality, the first set of LSNQ is for the left read(5') and the second set for the right(3') side; read depth information for CNV: DO=read depth outside, DI=read depth inside, D1=read depth on 1st break point, D2=on the second; paired-end information: PR=number of pairs that enclose the CNV; KM=length of kmers for the matched reads; Q0=percentage of ambiguous q=0 reads.</td></tr>
  <tr><td> 9 </td><td> MATCHCOUNT </td><td>  number of pairs of reads that support this CNV </td></tr>
  <tr><td> 10 </td><td> SVSEQ </td><td>  25 bases before and 25 bases after CNV on the reference </td></tr>
  <tr><td> 11 </td><td> BP </td><td>  CNV identifier </td></tr>
</table>
	
### tagcnv output:
<table>
  <tr><td> . </td><td> PAIRS </td><td> $n1,$n2,$n3, where n1=# of evident reads, n2=# of paired reads, n3=# of total reads. This field indicates whether there are pairs that support the variation. For deletion, evident pairs should cover both break points; for duplication, we look for pairs with negative insert within the break points. Use awk '$NF~/^0,/' to grep CNVs with no pair support. </td></tr>
  <tr><td> . </td><td> ZYGOSITY </td><td> HOM or HET=$T($MSRQ):$POS($Q)[,$POS($Q),...] T=total num of heterozyguous sites, MSRQ=mean square root of quality scores of zygosity as called by samtools/bcftools, POS=position, Q=quality score for each site </td></tr>
  <tr><td> . </td><td> RD[DIS=+-1000] </td><td> RDO=$n1;RDI=$n2;RD1=$n3;RD2=$n4  n1=average RD of outter region, n2=average RD inside CNV region, n3=RD of the lower end, n4=RD of the upper end </td></tr>
  <tr><td> . </td><td> MAP </td><td> [019][019]. This field indicates whether the merged read indeed spans the break points of a variation. The first digit is for the 5' break point and the second 3'. The purpose of this field is to test whether the merged read belongs to a continuous region somewhere else. The whole pierce of MERGE, left 70%, and right 70% are mapped using bwa bwasw. 0 is given if the break point is not covered, and 1 if it is. 9 is given if the whole MERGE is mapped with less than 8% of soft clipped bases. </td></tr>
  <tr><td> . </td><td> NOSEQ </td><td> PASS or FAIL </td></tr>
  <tr><td> . </td><td> MAPQ0 </td><td> proportion of low quality reads </td></tr>
</table>

## Parameters
For usage, type ```matchclips``` and return. The default parameters should be good for most current platforms as of 2013. If your reads are longer than 200 bases or the read quality is low, change maximum allowed mismatched bases to 2 or 3 with option ```-me```. 

Working on large real data, e.g., whole genome 20X coverage data, one should use as many threads as there are in the node with t. It might be slow and you may want to use ```-v``` option to check the intermediate result. For targeted sequencing, you can process the reads by REGIONs or use ```-L``` option to limit the search range. For exon sequencing, use ```-L 20000```; for short indel use ```-L 200```.

```
[yhwu@debian matchclip]$ ./matchclips
Usage:
  ./matchclips <options> -f REFFILE -b BAMFILE [REGION]

Options:
  -t  INT  number of threads [1,36], INT=1
  -h  INT  if set, only process first INT reads
  -ee INT  maximum allowed mismatches at 3' edges, INT=1
  -me INT  max allowed mismatches when matching strings, INT=1
  -s  INT  minimum number of soft clipped bases, INT=10
  -q  INT  minimum mapping score, INT=10
  -L  INT  check reads before and after INT bases, INT=3E6
  -nc      no colored display for intermediate info, colored=true
  -o  STR  outputfile, STR=STDOUT
  -dump    only dump reads with S part larger than threshold(-s)
  -v       print intermediate matching information
   REGION  if given should be in samtools's region format

Examples:
  ./matchclips -f human_g1k_v37.fasta -o bwap40X.bam.bp -b bwap40X.bam
  ./matchclips -f human_g1k_v37.fasta -o bwap40X.bam.bp -q 30 -b bwap40X.bam chr19 -t 4

Note:
  REGION must be understood by samtools. But, internally, chr1 and 1 are same.
  The default setting should work for most current platforms as of 2013
  If reads are longer than 200, you may want to use -me 2
  If reads qualities are poor at 3' end, you may try -ee 2 or more.
  Or, truncate reads before mapping.
  If you got nothing, run with -v to check mapq and S number stats.
```

## AUTHOR
Yinghua Wu, 
Department of Biostistics and Epidemiology, University of Pennsylvania, Philadelphia, PA 19104 

## License
GNU General Public License, version 3.0 (GPLv3)
