/*****************************************************************************
  genomeToREFrags.cpp
  Last-modified: 07 Mar 2012 12:02:34 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>
#include <cstdlib>

#include "common.h"
#include "stringUtils.h"
#include "fileUtils.h"
#include "bedUtils.h"

using namespace std;

//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: genomeToREFrags (v" << VERSION << ")" << endl; 
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Split genome file into Resitriction Enzyme fragments." << endl;
    cerr << "Usage:   genomeToREFrags [OPTIONS] -g <genome.fa> -c <AAGCTT> -b <REname.bed> -f <REends.fa>" << endl << endl;
	
	cerr << "Options:" << endl;
	
	cerr << "Input:"   << endl;
	cerr << "\t-g:--genome           Genome fasta format file." << endl << endl;
	cerr << "\t-c:--cut_seq          Restriction enzyme cut sequence. (default HindIII <AAGCTT>)." << endl << endl;
	
	cerr << "Output:" << endl;
	cerr << "\t-b:--bed_output       Bed file of the RE fragments. (default <stdout>)" << endl << endl;
	cerr << "\t-f:--fa_output        Fasta file of the RE end fragments. (default <ends.fa>)" << endl << endl;
	
	cerr << "Parameters:" << endl;
	cerr << "\t-s:--split            Split RE fragments to two ends. (0 or 1, default <0>)" << endl << endl;
	cerr << "\t-r:--read_len         Read length (RE end length) from the cutting position. (default <20>)" << endl << endl;
	cerr << "\t-h:--help             Show help information." << endl << endl;
    
}

int main( int argc, char *argv[])
{
	// Options
	bool     showHelp          = false;
	string   cutSeq   = "AAGCTT";
	string   genomeFile;
	string   bedFile           = "stdout";
	string   faFile            = "ends.fa";
	bool     split             = false;
	CHRPOS   readLen           = 20;

	// Show help when has no options
	if(argc <= 1)
	{
		Help();
		return 0;
	}

	// Parsing options
	for(int i = 1; i < argc; i++)
	{
		int parameterLength = (int)strlen(argv[i]);
		if((PARAMETER_CHECK("-h", 2, parameterLength)) || (PARAMETER_CHECK("--help", 5, parameterLength))) 
			showHelp=true;
		else if((PARAMETER_CHECK("-g", 2, parameterLength)) || (PARAMETER_CHECK("--genome", 8, parameterLength)))
		{
			if ((++i) < argc) 
				genomeFile = argv[i];
		}
		else if((PARAMETER_CHECK("-c", 2, parameterLength)) || (PARAMETER_CHECK("--cut_seq", 9, parameterLength)))
		{
			if ((++i) < argc)
				cutSeq = argv[i];
		}
		else if ((PARAMETER_CHECK("-b", 2, parameterLength)) || (PARAMETER_CHECK("--bed_output", 12, parameterLength)))
		{
			if ((++i) < argc)
				bedFile=argv[i];
		}
		else if ((PARAMETER_CHECK("-f", 2, parameterLength)) || (PARAMETER_CHECK("--fa_output", 11, parameterLength)))
		{
			if ((++i) < argc)
				faFile=argv[i];
		}
		else if ((PARAMETER_CHECK("-s", 2, parameterLength)) || (PARAMETER_CHECK("--split", 7, parameterLength)))
		{
			if ((++i) < argc)
			{
				if (argv[i][0]=='1')
					split=true;
			}

		}
		else if ((PARAMETER_CHECK("-r", 2, parameterLength)) || (PARAMETER_CHECK("--read_len", 10, parameterLength)))
		{
			if ((++i) < argc)
				readLen = ToValue<CHRPOS>(argv[i]);
		}
		else
		{
			cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}
	}
	
	// Show help if no proper auguments.
	if (showHelp)
	{
		Help();
		return 0;
	}

	// Statistical variables
	map <string, int, less<string> > bedCount;
	map <string, CHRPOS, less<string> > bedSum;
	map <string, CHRPOS, less<string> > faSize;

    // Variables
	CHRPOS lindex,rindex;
	int siteLen=cutSeq.size();
	bool flag;
	Fasta curFa;
	SeqReader fhfa(genomeFile);
	Writer bedOutput(bedFile);
	Writer faOutput(faFile);

	
	// Open files
	fhfa.Open();
	bedOutput.Open();
	faOutput.Open();

	// Read the genome file.
	while (fhfa.getNext(curFa))	
	{
		// Statistics
		bedCount[curFa.id]=0;
		bedSum[curFa.id]=0;
		faSize[curFa.id]=curFa.length();
		
		// Find next recognition site.
		lindex=rindex=0;
		flag=true;
		while (flag)
		{
			rindex=curFa.seq.find(cutSeq,lindex);
			if(rindex==CHRPOS(string::npos))
			{
				rindex=curFa.seq.size();
				flag=false;
			}

			if(split)
			{
				(*(bedOutput.Printer())) << curFa.id << "\t" << lindex << "\t" << (lindex+rindex)/2 << "\t" << (3*lindex+rindex)/4 << endl;
				(*(bedOutput.Printer())) << curFa.id << "\t" << (lindex+rindex)/2 << "\t" << rindex << "\t" << (lindex+3*rindex)/4 << endl;
				if (rindex - lindex >= readLen)
				{
					bedCount[curFa.id]++;
					bedSum[curFa.id]+=rindex-lindex;
					(*(faOutput.Printer())) << ">" << curFa.id << "_" << (3*lindex+rindex)/4 << endl;
					(*(faOutput.Printer())) << curFa.seq.substr(lindex,readLen) << endl;
					(*(faOutput.Printer())) << ">" << curFa.id << "_" << (lindex+3*rindex)/4 << endl;
					(*(faOutput.Printer())) << curFa.seq.substr(rindex-readLen,readLen) << endl;
				}
			}
			else
			{
				(*(bedOutput.Printer())) << curFa.id << "\t" << lindex << "\t" << rindex << "\t" << (lindex+rindex)/2 << endl;
				if (rindex - lindex >= 2*readLen)
				{
					bedCount[curFa.id]++;
					bedSum[curFa.id]+=rindex-lindex;
					(*(faOutput.Printer())) << ">" << curFa.id << "_" << (lindex+rindex)/2 << "_L" << endl;
					(*(faOutput.Printer())) << curFa.seq.substr(lindex,readLen) << endl;
					(*(faOutput.Printer())) << ">" << curFa.id << "_" << (lindex+rindex)/2 << "_R" << endl;
					(*(faOutput.Printer())) << curFa.seq.substr(rindex-readLen,readLen) << endl;
				}
			}
			lindex=rindex+siteLen;
		}
	}
	
	// Close files
	fhfa.Close();
	bedOutput.Close();
	faOutput.Close();

	// print statistics into log file
	Writer log(cutSeq+".log");
	log.Open();
	log.Close();
	
	return 0;
}

