/*****************************************************************************
  seqToContact.cpp
  Last-modified: 16 Mar 2012 10:42:52 AM

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
#include "stringUtils.h"
#include "fileUtils.h"
#include "bedUtils.h"

using namespace std;

vector<string> ids;
map<string, vector< BedVec> > hitsmap;

//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: seqToContact (v" << VERSION << ")" << endl; 
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Bowtie mapped files to RE fragment contacts." << endl;
    cerr << "Usage:   seqToContact [options] -b <REname.bed> -1 <seq_1.map> -2 <seq_2.map>" << endl << endl;
            
    cerr << "Options:" << endl;
            
    cerr << "Input:"   << endl;
    cerr << "    -b:--bed             Bed file of the RE fragments." << endl << endl;
    cerr << "    -1:--first            First Bowtie mapped file." << endl;
    cerr << "    -2:--second           Second Bowtie mapped file." << endl;
            
    cerr << "Output:"  << endl;
    cerr << "    -o:--output           Bed contacts with frequency." << endl << endl;
        
    cerr << "Parameters:" << endl;
	cerr << "    -r:--readlen          Read length of RNA-Seq." << endl << endl;
    cerr << "    -h:--help             Show help information." << endl << endl;
            
}

bool findPairs(ColumnReader &cr1, ColumnReader &cr2)
{
	int counts;
	Bed tbed;
	LINE_ELEMS elems;

	// Read hits from the first file
	if (cr1.getNext(elems)!=LINE_INVALID)
	{
		counts=StringUtils::toValue<int>(elems[6]);
		if(hitsmap[elems[0]].size()==0)  // if this id (elems[0]) not in the hitsmap
			hitsmap[elems[0]]=vector<BedVec>(2);
		// Put the line into hitsmap.
		BedUtils::bowtieToBed(tbed,elems);
		hitsmap[elems[0]][0].push_back(tbed);
		// Reportable if second vector has elements.
		if(hitsmap[elems[0]][1].size()) 
			ids.push_back(elems[0]);
		// read additional lines with the same id
		for( int i=0;i< counts; i++)
		{
			cr1.getNext(elems);
			BedUtils::bowtieToBed(tbed,elems);
			hitsmap[elems[0]][0].push_back(tbed);
		}
	}
	// Read hits from the second file
	else if(cr2.getNext(elems)!=LINE_INVALID)
	{
		counts=StringUtils::toValue<int>(elems[6]);
		// if this id (elems[0]) not in the hitsmap.
		if(hitsmap[elems[0]].size()==0)
			hitsmap[elems[0]]=vector<BedVec>(2);
		// Put the line into hitsmap.
		BedUtils::bowtieToBed(tbed,elems);
		hitsmap[elems[0]][1].push_back(tbed);
		// Reportable if second vector has elements.
		if(hitsmap[elems[0]][0].size())
			ids.push_back(elems[0]);
		// read additional lines with the same id
		for( int i=0;i< counts; i++)
		{
			cr2.getNext(elems);
			BedUtils::bowtieToBed(tbed,elems);
			hitsmap[elems[0]][1].push_back(tbed);
		}
	}
	else
		return false; // if both file are empty.
	return true;
}

int main(int argc, char* argv[])
{
	bool   showHelp    = false;
	int    readLen     = 20;
	string REFile;
	string firstFile;
	string secondFile;

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
				REFile=argv[i];
		}
		else if ((PARAMETER_CHECK("-1", 2, parameterLength)) || (PARAMETER_CHECK("--first", 7, parameterLength)))
		{
			if ((++i) < argc)
				firstFile=argv[i];
		}
		else if ((PARAMETER_CHECK("-2", 2, parameterLength)) || (PARAMETER_CHECK("--second", 8, parameterLength)))
		{
			if ((++i) < argc)
				secondFile=argv[i];
		}
		else if ((PARAMETER_CHECK("-r", 2, parameterLength)) || (PARAMETER_CHECK("--readlen", 9, parameterLength)))
		{
			if ((++i) < argc)
				readLen=StringUtils::toValue<int>(argv[i]);
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

	// RE bed file
	BedMap REMap;
	LINE_ELEMS elems;
	ColumnReader REReader(REFile);

	// Read RE fragment ends into BedMap.
	REReader.open();
	while(REReader.getNext(elems)!=LINE_INVALID)
	{
		Bed tbed(elems);
		Bed tbed2(tbed.chrom, tbed.start,       tbed.start+readLen, tbed.name, tbed.score);
		Bed tbed3(tbed.chrom, tbed.end-readLen, tbed.end,           tbed.name, tbed.score);
		REMap[tbed.chrom][tbed2.getBIN()].push_back(tbed2);
		REMap[tbed.chrom][tbed3.getBIN()].push_back(tbed3);
	}
	REReader.close();

	// File handles
	ColumnReader firstHitsReader(firstFile);
	ColumnReader secondHitsReader(secondFile);

	// Open files
	firstHitsReader.open();
	secondHitsReader.open();

	// Find Pair from pair-end hits
	BedVec *vec1,*vec2;
	pBedVec REVecs1,REVecs2;
	OverlapMap overlaps; // <CHRPOS, pBedVec>
	map<string,float> contactMap;
	
	while(ids.size() || findPairs(firstHitsReader, secondHitsReader))
	{
		if (ids.size())
		{
			vec1=&(hitsmap[ids[0]][0]);
			vec2=&(hitsmap[ids[0]][1]);
			for(BedVec::iterator it=vec1->begin();it!=vec1->end();it++)
			{
				BedUtils::intersectBed(*it,REMap,overlaps);
				if(overlaps.size())
					REVecs1.push_back(&(*(overlaps.begin()->second[0])));
			}
			
			for(BedVec::iterator it=vec2->begin();it!=vec2->end();it++)
			{
				BedUtils::intersectBed(*it,REMap,overlaps);
				if(overlaps.size())
					REVecs2.push_back(&(*(overlaps.begin()->second[0])));
			}
			
			for(pBedVec::iterator it1=REVecs1.begin();it1!=REVecs1.end();it1++)
				for(pBedVec::iterator it2=REVecs2.begin();it2!=REVecs2.end();it2++)
				{
					contactMap[(*it1)->chrom+"\t"+(*it1)->name+"\t"+(*it2)->chrom+"\t"+(*it2)->name]+=((*it1)->score+(*it2)->score)/2;
				}
			REVecs1.clear();
			REVecs2.clear();
			ids.erase(ids.begin());
		}
	}
	
	// Close files
	firstHitsReader.close();
	secondHitsReader.close();

	// Print contact map
	for(map<string,float>::iterator it=contactMap.begin();it!=contactMap.end();it++)
			cout << it->first << "\t" << it->second << endl;

    return 0;
}
