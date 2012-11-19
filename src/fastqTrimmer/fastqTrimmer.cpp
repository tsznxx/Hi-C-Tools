/*****************************************************************************
  fastqTrimmer.cpp
  Last-modified: 17 Sep 2012 03:17:23 PM

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
#include "seqUtils.h"

using namespace std;

//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: fastqTrimmer (v" << VERSION << ")" << endl; 
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Trim fastq format files (both ends) according to average Q scores in given window size. Paired-end support is in progress." << endl;

	cerr << "         trimmed reads : prefix+_trimmed.fq" << endl;
	cerr << "         low quality: prefix+_untrimmed.fq" << endl;
	
    cerr << "Usage:   fastqTrimmer [OPTIONS] -i reads.fq -p prefix" << endl << endl;
	
	cerr << "Options:" << endl;
	
	cerr << "Input:"   << endl;
	cerr << "\t-i:--input            Fastq format file." << endl << endl;
//	cerr << "\t-1:--first            First reads in Fastq format." << endl << endl;
//	cerr << "\t-2:--second           Second reads in Fastq format." << endl << endl;
	
	cerr << "Output:" << endl;
	cerr << "\t-p:--prefix           Prefix of output files." << endl << endl;
	
	cerr << "Parameters:" << endl;
	cerr << "\t-t:--threshold        Threshold for average Q scores. (default 30)" << endl << endl;
	cerr << "\t-w:--winsize          Window size. (default 6)" << endl << endl;
	cerr << "\t-o:--offset           Offset for Q scores. (default 33)" << endl << endl;
	cerr << "\t-m:--minlen           Minimum reads length. (default 30)" << endl << endl;
	cerr << "\t-h:--help             Show help information." << endl << endl;
    
}

int main( int argc, char *argv[])
{
	// Options
	string   readfile;
	string   prefix            = "result";
	int      winsize           = 6;
	int      minlen            = 30;
	int      threshold         = 30;
	int      offset            = PHREDOFFSET;
	bool     showHelp          = false;

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
		else if((PARAMETER_CHECK("-i", 2, parameterLength)) || (PARAMETER_CHECK("--input", 7, parameterLength)))
		{
			if ((++i) < argc) 
				readfile = argv[i];
		}
		else if((PARAMETER_CHECK("-w", 2, parameterLength)) || (PARAMETER_CHECK("--winsize", 9, parameterLength)))
		{
			if ((++i) < argc)
				winsize = StringUtils::toValue<int>(argv[i]);
		}
		else if ((PARAMETER_CHECK("-m", 2, parameterLength)) || (PARAMETER_CHECK("--minlen", 8, parameterLength)))
		{
			if ((++i) < argc)
				minlen = StringUtils::toValue<int>(argv[i]);
		}
		else if ((PARAMETER_CHECK("-p", 2, parameterLength)) || (PARAMETER_CHECK("--prefix", 8, parameterLength)))
		{
			if ((++i) < argc)
				prefix = argv[i];
		}
		else if ((PARAMETER_CHECK("-o", 2, parameterLength)) || (PARAMETER_CHECK("--offset", 8, parameterLength)))
		{
			if ((++i) < argc)
				offset = StringUtils::toValue<int>(argv[i]);
		}
		else if ((PARAMETER_CHECK("-t", 2, parameterLength)) || (PARAMETER_CHECK("--threshold", 11, parameterLength)))
		{
			if ((++i) < argc)
			{
				threshold = StringUtils::toValue<int>(argv[i]);
				if (threshold <0 || threshold >50)
				{
					cerr << endl << "*****ERROR: Threshold should be within 0~50, reset to default: 30. *****" << endl << endl;
					threshold = 30;
				}
			}
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

    // Variables
	int counter=0,lowcounter=0,i;
	int qsum=0,sumThreshold;
	char   buf[MAXSEQREADLEN]; // Store the sequence for trimming
	int start,seqlen;
	Fastq  seq;
	SeqReader fq(readfile);
	Writer trimmed(prefix+"_trimmed.fq");
	Writer lowqual(prefix+"_untrimmed.fq");
	
	// open files
	fq.open();
	trimmed.open();
	lowqual.open();
	
	// Read the fastq file.
	sumThreshold=(threshold+offset)*winsize;
	while (fq.getNext(seq))	
	{
		if (++counter%100000==0)
			cerr << "Finished  " << counter << " reads.\r";
		strcpy(buf,seq.phred.c_str());
		seqlen=seq.length();
		qsum=0;
		start=0;

		// Trim 5' ends
		for(i=0;i<winsize;i++)
			qsum+=buf[i];
		while(qsum<sumThreshold)
		{
			start++;
			if (seqlen-start< minlen)
				break;
			qsum+=buf[start+winsize-1]-buf[start-1];
		}
		// Trim 3' ends
		if (seqlen-start >= minlen)
		{
			qsum=0;
			for( i=seqlen-1;i>=seqlen-winsize;i--)
				qsum+=buf[i];
			while (qsum<sumThreshold)
			{
				seqlen--; // trim the last nucleotide.
				if (seqlen-start < minlen)
					break;
				qsum+=buf[seqlen-winsize-1]-buf[seqlen-1];
			}
		}
		if (seqlen-start>=minlen)
		{
			seq.seq=seq.seq.substr(start,seqlen-start);
			seq.phred=seq.phred.substr(start,seqlen-start);
			(*(trimmed.Printer())) << seq;
		}
		else
		{
			(*(lowqual.Printer())) << seq;
			lowcounter++;
		}
	}

	// close files
	fq.close();
	trimmed.close();
	lowqual.close();

	// Print summary
	cerr << endl;
	cerr << "Parameters: \n\tthreshold = " << threshold << "\n\twindow size = " << winsize << "\n\tminimum length = " << minlen << endl; 
	cerr << "Summary:" << endl;
	cerr << "\tTotal reads  : " << counter << endl;
	cerr << "\tTrimmed reads: " << counter-lowcounter << endl;
	cerr << "\tShort reads  : " << lowcounter <<endl;

	return 0;
}

