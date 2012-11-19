/*****************************************************************************
  fastqFilter.cpp
  Last-modified: 18 Jul 2012 05:43:56 PM

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

    cerr << "Program: fastqFilter (v" << VERSION << ")" << endl; 
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Filter fastq format files according to average Q scores." << endl;
	cerr << "         single end : prefix+_se.fq" << endl;
	cerr << "         paired end : prefix+_pe_read1/2.fq" << endl;
	cerr << "         low quality: prefix+_low_qual.fq" << endl;
	
    cerr << "Usage:   fastqFilter [OPTIONS] -1 read1.fq [-2 read2.fq] -p prefix" << endl << endl;
	
	cerr << "Options:" << endl;
	
	cerr << "Input:"   << endl;
	cerr << "\t-1:--read1            First read file." << endl << endl;
	cerr << "\t-2:--read2            Second read file if paired (optional)" << endl << endl;
	
	cerr << "Output:" << endl;
	cerr << "\t-p:--prefix           Prefix of output files." << endl << endl;
	
	cerr << "Parameters:" << endl;
	cerr << "\t-t:--threshold        Threshold for Q scores. (default 30)" << endl << endl;
	cerr << "\t-o:--offset           Offset for Q scores. (default 33)" << endl << endl;
	cerr << "\t-h:--help             Show help information." << endl << endl;
    
}

int main( int argc, char *argv[])
{
	// Options
	bool     isPaired          = false;
	string   read1;
	string   read2;
	string   prefix            = "result";
	int      threshold         = 30;
	int      offset            = 33;
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
		else if((PARAMETER_CHECK("-1", 2, parameterLength)) || (PARAMETER_CHECK("--read1", 7, parameterLength)))
		{
			if ((++i) < argc) 
				read1 = argv[i];
		}
		else if((PARAMETER_CHECK("-2", 2, parameterLength)) || (PARAMETER_CHECK("--read2", 7, parameterLength)))
		{
			if ((++i) < argc)
			{
				read2 = argv[i];
				isPaired = true;
			}
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
	if (read1=="")
	{
		if (read2!="")
		{
			read1    = read2;
			read2    = "";
			isPaired = false;
		}
		else
		{
			cerr << endl << "*****ERROR: No fastq file provided! *****" << endl << endl;
			Help();
			return 0;
		}
	}


    // Variables
	int counter=0,lowcounter=0;
	float  qual1,qual2;
	Fastq  seq1,seq2;
	SeqReader fq1(read1);
	SeqReader fq2(read2);
	Writer single(prefix+"_filtered_se.fq");
	Writer paired1(prefix+"_filtered_R1.fq");
	Writer paired2(prefix+"_filtered_R2.fq");
	Writer lowqual(prefix+"_filetered_lowqual.fq");
	
	// open files
	fq1.open();
	single.open();
	lowqual.open();
	if (isPaired)
	{
		fq2.open();
		paired1.open();
		paired2.open();
	}
	
	// Read fastq files.
	while (fq1.getNext(seq1))	
	{
		if (++counter%10000==0)
			cerr << "Finished  " << counter << " reads." << endl;
		qual1=seq1.meanQ(offset);
		if (isPaired) // paired end
		{
			fq2.getNext(seq2);
			qual2=seq2.meanQ(offset);
			if (qual1 > threshold && qual2 > threshold) // both good
			{
				(*(paired1.Printer())) << seq1;
				(*(paired2.Printer())) << seq2;
			}
			else if (qual1 > threshold)  // read1 good
			{
				seq1.id+="/1";
				(*(single.Printer()))  << seq1;
				seq2.id+="/2";
				(*(lowqual.Printer())) << seq2;
				lowcounter++;
			}
			else if (qual2 > threshold) // read2 good
			{
				seq1.id+="/1";
				(*(lowqual.Printer())) << seq1;
				lowcounter++;
				seq2.id+="/2";
				(*(single.Printer()))  << seq2;
			}
			else // both bad
			{
				seq1.id+="/1";
				(*(lowqual.Printer())) << seq1;
				seq2.id+="/2";
				(*(lowqual.Printer())) << seq2;
				lowcounter+=2;
			}
		}
		else // single end
		{
			if (qual1 > threshold)
				(*(single.Printer()))  << seq1;
			else
			{
				(*(lowqual.Printer())) << seq1;
				lowcounter++;
			}
		}
	}

	// close files
	fq1.close();
	single.close();
	lowqual.close();
	if (isPaired)
	{
		fq2.close();
		paired1.close();
		paired2.close();
	}

	// Print summary
	cerr << endl;
	cerr << "Summary:" << endl;
	cerr << "Total reads  : " << counter << endl;
	cerr << "Passed filter: " << counter-lowcounter << "\treads in " << prefix << (isPaired ? "_pe_read1/2.fq": "_se.fq") << endl;
	cerr << "low quality  : " << lowcounter << "\treads in " << prefix << "_low_qual.fq" << endl;

	return 0;
}


