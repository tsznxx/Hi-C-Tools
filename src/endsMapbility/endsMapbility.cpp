/*****************************************************************************
  endsMapbility.cpp
  Last-modified: 07 Mar 2012 12:04:05 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>

#include "common.h"
#include "bedUtils.h"
#include "fileUtils.h"

using namespace std;

//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: endsMapbility (v" << VERSION << ")" << endl; 
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Split genome file into Resitriction Enzyme fragments." << endl;
    cerr << "Usage:   endsMapbility [options] -b <REname.bed> -f <REends.fa>" << endl << endl;
        
    cerr << "Options:" << endl;
        
    cerr << "Input:"   << endl;
    cerr << "\t-b:--bed              Bed file of the RE fragments. (default <stdout>)" << endl << endl;
    cerr << "\t-f:--fa               Fasta file of the RE end fragments. (default <ends.fa>)" << endl << endl;
        
    cerr << "Parameters:" << endl;
    cerr << "\t-s:--split            Split RE fragments to two ends. (0 or 1, default <0>)" << endl << endl;
    cerr << "\t-r:--read_len         Read length (RE end length) from the cutting position. (default <20>)" << endl << endl;
    cerr << "\t-h:--help             Show help information." << endl << endl;
        
}
typedef map<CHRPOS, Bed, std::less<CHRPOS> > BedChr;
typedef map<string, BedChr, std::less<string> > BedGenome;

int main(int argc, char* argv[])
{
	bool          showHelp        = false;
	string        hitsfile      = "stdin";
	string        bedfile         = "stdin";
	ColumnReader     bedreader;
	ColumnReader  hitsreader;

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
		else if ((PARAMETER_CHECK("-m", 2, parameterLength)) || (PARAMETER_CHECK("--maped", 7, parameterLength)))
		{
			if ((++i) < argc)
				hitsfile=argv[i];
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
	BedGenome bedvecs;
	LINE_ELEMS elems;

	// Read the bed file
	bedreader.setFileName(bedfile);
	bedreader.Open();
	while(bedreader.getNext(elems)!=LINE_INVALID)
	{
		Bed tbed(elems);
		bedvecs[tbed.chrom][ToValue<CHRPOS>(tbed.name)]=tbed;
	}
	bedreader.Close();

	// Read the mapping file
	hitsreader.setFileName(hitsfile);
	hitsreader.Open();
	vector<string> keys;

	while(hitsreader.getNext(elems)!=LINE_INVALID)
	{
		Tokenize(elems[0],keys,"_");
		bedvecs[keys[0]][ToValue<CHRPOS>(keys[1])].score++;
		keys.clear();
	}

	hitsreader.Close();

	for(BedGenome::iterator it=bedvecs.begin();it!=bedvecs.end();it++)
		for(BedChr::iterator it2=it->second.begin();it2!=it->second.end();it2++)
			cout << it2->second << endl;

    return 0;
}

