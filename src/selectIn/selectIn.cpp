/*****************************************************************************
  fastqTrimmer.cpp
  Last-modified: 20 May 2012 09:23:52 PM

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
#include <algorithm>

#include "common.h"
#include "stringUtils.h"
#include "fileUtils.h"

using namespace std;

//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: selectIn (v" << VERSION << ")" << endl; 
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Select rows whose nth column has items in the given list." << endl;
	
    cerr << "Usage:   selectIn [OPTIONS] -i input.txt -l items.lst -n 1 -p prefix " << endl << endl;
	
	cerr << "Options:" << endl;
	
	cerr << "Input:"   << endl;
	cerr << "\t-i:--input            Column file." << endl << endl;
	cerr << "\t-l:--list             List of items." << endl << endl;
	
	cerr << "Output:" << endl;
	cerr << "\t-p:--prefix           Prefix of output files." << endl << endl;
	
	cerr << "Parameters:" << endl;
	cerr << "\t-n:--nthcol           The nth column for selection (default 1)." << endl << endl;
	cerr << "\t-v:--reverse          Reverse selection." << endl << endl;
	cerr << "\t-h:--help             Show help information." << endl << endl;
    
}

int main( int argc, char *argv[])
{
	// Options
	string   infile;
	string   listfile;
	string   outfile           = "stdout";
	bool     revsel            = false;
	int      nthcol            = 1;
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
				infile = argv[i];
		}
		else if((PARAMETER_CHECK("-l", 2, parameterLength)) || (PARAMETER_CHECK("--list", 6, parameterLength)))
		{
			if ((++i) < argc)
				listfile = argv[i];
		}
		else if((PARAMETER_CHECK("-n", 2, parameterLength)) || (PARAMETER_CHECK("--nthcol",8, parameterLength)))
		{
			if ((++i) < argc)
				nthcol = StringUtils::toValue<int>(argv[i]);
		}
		else if ((PARAMETER_CHECK("-v", 2, parameterLength)) || (PARAMETER_CHECK("--revsel", 8, parameterLength)))
		{
			revsel = true ;
		}
		else if ((PARAMETER_CHECK("-o", 2, parameterLength)) || (PARAMETER_CHECK("--output", 8, parameterLength)))
		{
			if ((++i) < argc)
				outfile = argv[i];
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
	vector<string>::iterator it;
	LINE_ELEMS elems;
	vector<string> items;
	ColumnReader fhi(infile);
	ColumnReader fhl(listfile);
	Writer       fho(outfile);
	
	// Read items from list file
	fhl.open();
	while (fhl.getNext(elems)!=LINE_INVALID)
		items.push_back(elems[0]);
	sort(items.begin(),items.end());
	fhl.close();
	
	// Read the column file for selection
	fhi.open();
	fho.open();
	nthcol--;
	if (nthcol<0)
	{
		cerr << "*****ERROR: Column number should be greater than 0. Set to " << endl << endl;
		exit(0);
	}
	while (fhi.getNext(elems)!=LINE_INVALID)
	{
		if ( revsel != binary_search(items.begin(),items.end(),elems[nthcol]) )
		{
			for( it=elems.begin();it<elems.end()-1;it++)
				(*(fho.Printer())) << *it << "\t";
			(*(fho.Printer())) << *(elems.end()-1) << endl;
		}
	}
	
	// Close files
	fhi.close();
	fho.close();
	
	return 0;
}

