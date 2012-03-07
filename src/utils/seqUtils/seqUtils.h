/*****************************************************************************
  seqUtils.h
  Last-modified: 07 Mar 2012 11:58:04 AM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#ifndef SEQUTILS_H
#define SEQUTILS_H

#include <string>

#include "common.h"
#include "fileUtils.h"

using namespace std;

class Fasta
{
public:
    string id;
    string seq;

public:
    // Construction
	Fasta(){};
    Fasta(const string &tid,const string &tseq);

    // Destruction
    ~Fasta(void){};

    // Length of the seq
    CHRPOS length(void)
    {
        return seq.size();
    }

	// Format output
	void print(CHRPOS lineSize=80);
};

class Fastq: public Fasta
{
public:
    string phred;
    
public:
    // Construction
    Fastq(const string &tid, const string  &tseq, const string &tphred);

    // Destruction
    ~Fastq(void){};

    // Get the mean Q(Phred) score
    float meanQ(int offset=33)
    {
        int psum,lseq;
        for(int i=0;i<lseq;i++)
        {
            psum+=phred.c_str()[i]-offset;
        }
        return (float)psum/lseq;
    }
        
}; 

class SeqReader: public Reader
{

public:
    // Construction
    SeqReader(const string &fname="stdin"); // Since we have the filename, it's better to open it here.
	
	// Destruction
    ~SeqReader(void); // Remember to close the file handle.

    // Get next Fasta record
    bool getNext(Fasta &fa);
    
    // Get next Fastq record
    bool getNext(Fastq &fq);

};


namespace SeqUtils
{
	string toUpper (const string &seq); // DNA/RNA to uppercase
	string toLower (const string &seq); // DNA/RNA to lowercase
	string toRNA   (const string &seq); // DNA to RNA
	string toDNA   (const string &seq); // RNA to DNA
	string rc      (const string &seq); // reverse complementary sequence
	float  GC      (const string &seq); // GC content
	float  TM      (const string &seq); // TM value
	float  MW      (const string &seq); // Molecular weigth
}

#endif //SEQUTILS_H
