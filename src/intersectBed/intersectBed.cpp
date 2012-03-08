/*****************************************************************************
  intersectBed.cpp
  Last-modified: 08 Mar 2012 10:04:53 AM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>

#include "fileUtils.h"
#include "bedUtils.h"

using namespace std;

//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: endsMapbility (v" << VERSION << ")" << endl; 
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Split genome file into Resitriction Enzyme fragments." << endl;
    cerr << "Usage:   intersect [options] -a <a.bed> -b <b.bed>" << endl << endl;
            
    cerr << "Options:" << endl;
            
    cerr << "Input:"   << endl;
    cerr << "\t-a:--beda             Bed file of the RE fragments. (default <stdout>)" << endl << endl;
    cerr << "\t-b:--bedb             Fasta file of the RE end fragments. (default <ends.fa>)" << endl << endl;
            
    cerr << "Parameters:" << endl;
    cerr << "\t-s:--forcestrand      Split RE fragments to two ends. (0 or 1, default <0>)" << endl << endl;
    cerr << "\t-f:--fraction         Read length (RE end length) from the cutting position. (default <20>)" << endl << endl;
    cerr << "\t-h:--help             Show help information." << endl << endl;
            
}


int main(int argc, char* argv[])
{
    bool          showHelp        = false;
    string        bedafile;
    string        bedbfile;
	bool          forcestrand     = false;
	float         fraction        = 0.8;

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
        else if ((PARAMETER_CHECK("-a", 2, parameterLength)) || (PARAMETER_CHECK("--beda", 6, parameterLength)))
        {
            if ((++i) < argc)
                bedafile=argv[i];
        }
        else if ((PARAMETER_CHECK("-b", 2, parameterLength)) || (PARAMETER_CHECK("--bedb", 6, parameterLength)))
        {       
            if ((++i) < argc)
                bedbfile=argv[i];
        }
		else if ((PARAMETER_CHECK("-s", 2, parameterLength)) || (PARAMETER_CHECK("--forcestrand", 12, parameterLength)))
		{
			if ((++i) < argc)
				forcestrand=ToValue<int>(argv[i]);
		}
		else if ((PARAMETER_CHECK("-f", 2, parameterLength)) || (PARAMETER_CHECK("--fraction", 10, parameterLength)))
		{
			if ((++i) < argc)
				fraction = ToValue<float>(argv[i]);
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
	
	LINE_ELEMS elems;
	Hits hits;
	BedMap  bedmap;
	ColumnReader beds(bedafile);
	beds.Open();
	
	BedUtils::loadBedFileToMap(bedmap,bedbfile);
	while(beds.getNext(elems)!=LINE_INVALID)
	{
		Bed tbed(elems);
		BedUtils::intersectBed(tbed,bedmap,hits,forcestrand);
		cout << tbed ;
		if(hits.size())
			cout << "\t" << hits.begin()->second[0] << "\t" << hits.begin()->first << endl;
		else
			cout << "\t.\t.\t.\t.\t.\t.\t0" << endl;
		hits.clear();
	}
	
	beds.Close();


    return 0;
}

