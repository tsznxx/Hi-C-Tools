#include <iostream>
#include <list>
#include <algorithm>
#include <cctype>

#include "seqUtils.h"

using namespace std;

/*************************************
 Fasta
*************************************/


Fasta::Fasta(const string &tid, const string &tseq):id(tid),seq(tseq)
{
}

// Format output (lineSize=0: print the whole sequence in one line
void Fasta::print(CHRPOS lineSize)
{
	int i,curlen;

	curlen=seq.size()-lineSize;
	cout << ">" << id << endl;
	for(i=0; i< curlen;i+=lineSize)
	{
		cout << seq.substr(i,lineSize) << endl;
	}
	cout << seq.substr(i) << endl;
}

/*************************************
 Fastq
*************************************/

// Constructor
Fastq::Fastq(const string &tid, const string &tseq, const string &tphred):Fasta(tid,tseq)
{
	phred=tphred;
}

/*************************************
Genome
*************************************/


/*************************************
 SeqReader
*************************************/

// Constructor
SeqReader::SeqReader(const string &fname):Reader(fname)
{
}

// Destructor
SeqReader::~SeqReader(void)
{
	close();
}


// Get next fasta record
bool SeqReader::getNext(Fasta &fa)
{
	// Initiation
	string curline,tseq;
	static string tid;
	list<string> seqlst;
	bool flag=true;
	
	// Get next Fasta
	while (flag)
	{
		if(_inStream->good())
			getline(*_inStream,curline);
		else // The last run
		{
			curline=">";
			flag=false;
		}

		if (*(curline.begin())=='>')
		{
			if (tid!="")
			{
				fa.id=tid;
				tseq.reserve(seqlst.size()* seqlst.begin()->size());
				for(list<string>::iterator it=seqlst.begin(); it!=seqlst.end();it++)
					tseq.append(*it);
				fa.seq=tseq;
				tid=curline.substr(1);
				return true;
			}
			else // The first run
				tid=curline.substr(1);
		}
		else
		{
			transform(curline.begin(), curline.end(), curline.begin(), (int (*)(int))toupper);
			seqlst.push_back(curline);
		}
	}
	if (tid!="" && tseq!="") // The last run
	{
		fa.id=tid;
		tseq.reserve(seqlst.size()* seqlst.begin()->size());
		for(list<string>::iterator it=seqlst.begin(); it!=seqlst.end();it++)
			tseq.append(*it);
		fa.seq=tseq;
		return true;
	}
	return false;
}

// Get next fastq record
bool SeqReader::getNext(Fastq &fq)
{
	// Initiations
	string curline;

	// Read new record
	if (_inStream->good())
	{
		getline(*_inStream,curline);
		if (curline.size())
		{
			fq.id=curline.substr(1);
		}
		else
			return false;
		getline(*_inStream,fq.seq);
		getline(*_inStream,curline);
		getline(*_inStream,fq.phred);
		return true;
	}
	else
		return false;
}

/*************************************
 Namespace SeqUtils
*************************************/

// DNA/RNA to uppercase
string SeqUtils::toUpper(const string &seq)
{
	string upper;
	transform(seq.begin(),seq.end(),upper.begin(),(int (*)(int)) toupper );
	return upper;
}

// DNA/RNA to lowercase
string SeqUtils::toLower(const string &seq)
{
	string lower;
	transform(seq.begin(),seq.end(),lower.begin(),(int (*)(int)) tolower );
	return lower;
}

// DNA to RNA
string SeqUtils::toRNA  (const string &seq)
{
	string rna=seq;
	size_t found;

	found=rna.find_first_of("tT");
	while(found!=string::npos)
	{
		rna[found++]++; // 't' to 'u' or 'T' to 'U'
		found=rna.find_first_of("tT",found);
	}
	return rna;
}


// RNA to DNA
string SeqUtils::toDNA  (const string &seq)
{
	string dna=seq;
	size_t found;

	found=dna.find_first_of("uU");
	while(found!=string::npos)
	{
		dna[found++]--; // 'u' to 't' or 'U' to 'T'
		found=dna.find_first_of("uU",found);
	}
	return dna;
}

// reverse complermentary table
char transtable(char i)
{
	static string trans="ACGNTacgntTGCNAtgcna";
	static size_t found;

	found=trans.find_first_of(i);
	if(found!=string::npos)
		return trans[found+10];
	else
		return 'N';
}

// reverse complementary sequence
string SeqUtils::rc (const string &seq)
{
	string rcseq=seq;
	reverse(rcseq.begin(),rcseq.end());
	transform(rcseq.begin(),rcseq.end(),rcseq.begin(),transtable);
	return rcseq;
}

// GC content
float SeqUtils::GC (const string &seq)
{
	size_t gc=0;
	size_t found;
	found=seq.find_first_of("CGcg");
	while(found!=string::npos)
	{
		gc++;
		found=seq.find_first_of("CGcg",found+1);
	}
	return 100.0*gc/seq.size();
}

// TM value
float SeqUtils::TM(const string &seq)
{
	// Add code here.
	return 60;
}

// Molecular weight
float SeqUtils::MW(const string &seq)
{
	// Add code here.
	return 3000;
}






