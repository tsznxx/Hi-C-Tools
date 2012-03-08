/*****************************************************************************
  readsAnnotation.cpp
  Last-modified: 08 Mar 2012 11:10:27 AM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>

#include "bedUtils.h"
#include "fileUtils.h"

using namespace std;

//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: endsMapbility (v" << VERSION << ")" << endl;
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Annotate Bed with multiple annotation files." << endl;
    cerr << "Usage:   readsAnnotation [options] -i <a.bed> -a <anno1.bed>[,anno2.bed,....]" << endl << endl;
                  
    cerr << "Options:" << endl;   
                  
    cerr << "Input:"   << endl;
    cerr << "\t-i:--input             input file. (Bed format)" << endl << endl;
    cerr << "\t-a:--annotations       Annotation files seperated with \",\". (Bed, Tab or bowtie format)" << endl << endl;
            
    cerr << "Parameters:" << endl;
    cerr << "\t-s:--forcestrand      Force strandedness for overlap check. (default false <0>)" << endl << endl;
    cerr << "\t-h:--help             Show help information." << endl << endl;

}     


int main(int argc, char* argv[])
{
	//options
	// -i inputfile -a annotationfiles -f forcestrand

	bool           showHelp        = false;
    string         infile;
    vector<string> annofiles;
    bool           forcestrand     = false;
    
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
        else if ((PARAMETER_CHECK("-i", 2, parameterLength)) || (PARAMETER_CHECK("--input", 7, parameterLength)))
        {
            if ((++i) < argc)
                infile=argv[i];
        }
		else if ((PARAMETER_CHECK("-a", 2, parameterLength)) || (PARAMETER_CHECK("--annotations", 13, parameterLength)))
        {       
            if ((++i) < argc)
			{
				string fs=argv[i];
				cout << fs << endl;
                Tokenize(fs,annofiles,",");
			}
        }
        else if ((PARAMETER_CHECK("-s", 2, parameterLength)) || (PARAMETER_CHECK("--forcestrand", 13, parameterLength)))
        {
            if ((++i) < argc)
                forcestrand= !strcmp("0",argv[i]);
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
	
	// definition
	int annocount=annofiles.size();    // File count
	vector<BedMap> bedmaps(annocount); // BedMap 
	vector<double> annos(annocount+1,0);   // Count reads for each kind of annotaions. The last one counts the remained reads.
	LINE_ELEMS elems;                  // Column reader buffer
	Hits hits;
	ColumnReader beds(infile);

	// Load annotations into BedMaps
	for(int i=0;i<annocount;i++)
		BedUtils::loadBedFileToMap(bedmaps[i],annofiles[i]);

	// Annotation
	beds.Open();
	while(beds.getNext(elems)!=LINE_INVALID)
	{
		Bed tbed(elems);
		tbed.score++;
		cout << tbed;
		bool intersectflag=true;
		for(int i=0;i<annocount;i++)
		{
			if(intersectflag)
			{
				BedUtils::intersectBed(tbed,bedmaps[i],hits,forcestrand);
				if(hits.size())
				{
					annos[i]+=tbed.score;
					cout << "\t" << tbed.score;
					intersectflag=false; // if found hits, suppress intersection for the following annotations.
				}
				else
					cout << "\t0";
				hits.clear();
			}
			else
				cout << "\t0";
		}
		if(intersectflag) // No overlas with annotations.
			annos[annocount]+=tbed.score;
		cout << endl;
	}


	for(int i=0;i<annocount;i++)
		cerr << annofiles[i] << "\t" << annos[i] << endl;
	cerr << "Remained\t" << annos[annocount] << endl;

    return 0;
}

