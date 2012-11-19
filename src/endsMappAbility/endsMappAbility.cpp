/*****************************************************************************
  endsMappAbility.cpp
  Last-modified: 13 Mar 2012 04:49:18 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>
#include <iomanip>

#include "common.h"
#include "bedUtils.h"
#include "fileUtils.h"

using namespace std;

//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: endsMapAbility (v" << VERSION << ")" << endl; 
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Ends map ability of the Resitriction Enzyme fragments. The mapability is 1 over # of best hits." << endl;
    cerr << "Usage:   endsMapAbility [options] -b <REname.bed> -m <RE_ends.map>" << endl << endl;
        
    cerr << "Options:" << endl;
        
    cerr << "Input:"   << endl;
    cerr << "    -b:--bed              Bed file of the RE fragments." << endl << endl;
    cerr << "    -m:--map              Bowtie mapped ends of the RE end fragments." << endl;
	cerr << "                           - Default: stdin" << endl << endl;
        
    cerr << "Output:"  << endl;
	cerr << "    -o:--output           Bed format with revised scores." << endl << endl;
	
	cerr << "Parameters:" << endl;
    cerr << "    -h:--help             Show help information." << endl << endl;
        
}
typedef map<CHRPOS, Bed, std::less<CHRPOS> > BedChr;
typedef map<string, BedChr, std::less<string> > BedGenome;

int main(int argc, char* argv[])
{
	bool          showHelp        = false;
	string        hitsfile        = "stdin";
	string        bedfile;
	string        outfile         = "stdout";
	ColumnReader  bedreader;
	ColumnReader  hitsreader;
	Writer        outwriter;

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
		else if ((PARAMETER_CHECK("-b", 2, parameterLength)) || (PARAMETER_CHECK("--bed", 5, parameterLength)))
		{
			if ((++i) < argc)
				bedfile=argv[i];
		}
		else if ((PARAMETER_CHECK("-m", 2, parameterLength)) || (PARAMETER_CHECK("--map", 5, parameterLength)))
		{
			if ((++i) < argc)
				hitsfile=argv[i];
		}
		else if ((PARAMETER_CHECK("-o", 2, parameterLength)) || (PARAMETER_CHECK("--output", 8, parameterLength)))
		{
			if ((++i) < argc)
				outfile=argv[i];
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
	Bed tbed;
	BedGenome bedgnms;
	LINE_ELEMS elems;

	// Read the bed file
	bedreader.setFileName(bedfile);
	bedreader.open();
	while(bedreader.getNext(elems)!=LINE_INVALID)
	{
		Bed tbed(elems);
		bedgnms[tbed.chrom][(tbed.start+tbed.end)/2]=tbed;
	}
	bedreader.close();

	// Read the mapping file
	hitsreader.setFileName(hitsfile);
	hitsreader.open();
	vector<string> keys;

	while(hitsreader.getNext(elems)!=LINE_INVALID)
	{
		StringUtils::tokenize(elems[0],keys,"_");
		bedgnms[keys[0]][StringUtils::toValue<CHRPOS>(keys[1])].score++;
		keys.clear();
	}

	// close file
	hitsreader.close();

	// Write into file
	outwriter.setFileName(outfile);
	outwriter.open();

	for(BedGenome::iterator it=bedgnms.begin();it!=bedgnms.end();it++)
		for(BedChr::iterator it2=it->second.begin();it2!=it->second.end();it2++)
		{
			it2->second.score = it2->second.score ? 2/it2->second.score : 0; 
			(*(outwriter.Printer())) << setprecision(3) << it2->second << endl;
		}
	
	// close writer handle
	outwriter.close();
	
    return 0;
}

