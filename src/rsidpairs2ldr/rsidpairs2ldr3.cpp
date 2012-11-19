/*****************************************************************************
  rsidpairs2ldr.cpp
  Last-modified: 12 Nov 2012 03:06:30 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include <cmath>

#include "common.h"
#include "stringUtils.h"
#include "fileUtils.h"

using namespace std;


float P11(float A, float B, float C, float D, float x1, float x2, float C_OK, int cnt)
{
	float p11,R;
	p11 = (x1+x2)/2;
	R = A*pow(p11,3)+B*pow(p11,2)+C*p11+D;
	
	if(fabs(x1-x2)<C_OK)
	{
		return p11;
	}
	else
	{
		if(R==0)
		{
			return p11;
		}
		else if(R>0)
		{
			return P11(A, B, C, D, x1, p11, C_OK, cnt+1);
		}
		else
		{
			return P11(A, B, C, D, p11, x2, C_OK, cnt+1);
		}
	} 
	return 0;
}

void LDR2( const vector<string> &array1, const vector<string> &array2, float &ldr, float &errcode)
{
	static float cutoff = 0.1;

	// different array size
	if (array1.size() != array2.size() )
	{
		ldr=-1;
		errcode=-1;
		return;
	}

	// Count AA pairs
	map< string, map<string,int> > hapmapA;
	map< string, map<string,int> > hapmapB;
	vector<string>::const_iterator it1=++array1.begin(); // shift the first element.
	vector<string>::const_iterator it2=++array2.begin();

	map<string,int> keysA,keysB;

	int  cnt=0;
	while(it1!=array1.end())
	{
		if( it1->find("N") == string::npos && it2->find("N") == string::npos ) // No 'N's
		{
			if (hapmapA.find(*it1)==hapmapA.end() || hapmapA[*it1].find(*it2)==hapmapA[*it1].end()) // key not exist
				hapmapA[*it1][*it2]=0;
			if (hapmapB.find(*it2)==hapmapB.end() || hapmapB[*it2].find(*it1)==hapmapB[*it2].end())
				hapmapB[*it2][*it1]=0;
			hapmapA[*it1][*it2]++;
			hapmapB[*it2][*it1]++;
			if ( keysA.find(*it1)==keysA.end())
				keysA[*it1]=0;
			if ( keysB.find(*it2)==keysB.end())
				keysB[*it2]=0;
			keysA[*it1]++;
			keysB[*it2]++;
			cnt++;
		}
		it1++;
		it2++;
	}
	// cnt/length < cutoff
	if (float(cnt)/array1.size() <= cutoff)
	{
		ldr=-1;
		errcode=-2;
		return;
	}

	// dominant AA counts
	map<string,int>::iterator it;
    string A1A1="";
	string B1B1="";
	string A1A2="";
	string B1B2="";
	int minc=0;

	for(it=keysA.begin();it!=keysA.end();it++)
	{
		if ( it->first[0] != it->first[1])
		{
			A1A2=it->first;
		}
		else if (it->second > minc)
		{
			minc=it->second;
			A1A1=it->first;
		}
	}
	minc=0;

	for(it=keysB.begin();it!=keysB.end();it++)
	{
		if ( it->first[0] != it->first[1])
		{
			B1B2=it->first;
		}
		else if (it->second > minc)
		{
			minc=it->second;
			B1B1=it->first;
		}
	}
	

	float p1,p2;
	p1 = (keysA[A1A1]*2 + keysA[A1A2]) /2.0/cnt;
	p2 = (keysB[B1B1]*2 + keysB[B1B2]) /2.0/cnt;

	if (p1 == 1 || p2 ==1)
	{
		ldr=-1;
		errcode=-3;
		return;
	}

	// 
	float n11=0,n12=0,n21=0,n22=0;
	if (hapmapA.find(A1A1)!=hapmapA.end() && hapmapA[A1A1].find(B1B1)!=hapmapA[A1A1].end())
		n11 = hapmapA[A1A1][B1B1];
	if (hapmapA.find(A1A1)!=hapmapA.end() && hapmapA[A1A1].find(B1B2)!=hapmapA[A1A1].end())
		n12 = hapmapA[A1A1][B1B2];
	if (hapmapA.find(A1A2)!=hapmapA.end() && hapmapA[A1A2].find(B1B1)!=hapmapA[A1A2].end())
		n21 = hapmapA[A1A2][B1B1];
	if (hapmapA.find(A1A2)!=hapmapA.end() && hapmapA[A1A2].find(B1B2)!=hapmapA[A1A2].end())
		n22 = hapmapA[A1A2][B1B2];

	float a_p3, b_p2, c_p1, d_p0, p11, D;
	a_p3 = 1;
	b_p2 = (2*cnt*(1-2*p1-2*p2)-n22-2*(2*n11+n12+n21))/(4*cnt);
	c_p1 = (2*cnt*p1*p2-n22*(1-p1-p2)-(2*n11+n12+n21)*(1-2*p1-2*p2))/(4*cnt);
	d_p0 = p1*p2*(2*n11+n12+n21)/((0-4)*cnt);

	p11 = P11(a_p3, b_p2, c_p1, d_p0, (2*n11+n12+n21)/(2*cnt), (2*n11+n12+n21+n22)/(2*cnt), 0.001, 1);
	D = p11 - p1*p2;
	ldr = D*D/(p1*(1-p1)*p2*(1-p2));

	errcode=float(cnt)/array1.size();
	return;
}


//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: fastqFilter (v" << VERSION << ")" << endl;
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com)" << endl;
    cerr << "Summary: Filter fastq format files according to average Q scores." << endl;
    cerr << "         single end : prefix+_filtered.fq" << endl;
    cerr << "         paired end : prefix+_filtered_R1/2.fq" << endl;
    cerr << "         low quality: prefix+_filetered_lowqual.fq" << endl;
        
    cerr << "Usage:   fastqFilter [OPTIONS] -1 read1.fq [-2 read2.fq] -p prefix" << endl << endl;
        
    cerr << "Options:" << endl;
        
    cerr << "Input:"   << endl;
    cerr << "\t-i:--rsidfile         rsid file." << endl << endl;
    cerr << "\t-h:--hapmap           Hapmap file." << endl << endl;
        
    cerr << "Output:" << endl;
    cerr << "\t-o:--output           output file." << endl << endl;
}


int main(int argc, char* argv[])
{
	// files
	string rsidfile;
	string hapmapfile;
	string outputfile;
	bool showHelp = false;


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
		else if((PARAMETER_CHECK("-i", 2, parameterLength)) || (PARAMETER_CHECK("--rsid", 6, parameterLength)))
		{
			if ((++i) < argc) 
				rsidfile = argv[i];
		}
		else if((PARAMETER_CHECK("-H", 2, parameterLength)) || (PARAMETER_CHECK("--Hapmap", 8, parameterLength)))
		{
			if ((++i) < argc)
				hapmapfile = argv[i];
		}
		else if ((PARAMETER_CHECK("-o", 2, parameterLength)) || (PARAMETER_CHECK("--output", 8, parameterLength)))
		{
			if ((++i) < argc)
				outputfile = argv[i];
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
	
	// Initiation
	LINE_ELEMS elems;
	set<string> rsids; // unique rsids
	vector< pair<string,string> > rsidpairs; //rsid pairs
	ColumnReader rsidfh(rsidfile);

	// read pairs
	cerr << "Reading " << rsidfile << endl;
	rsidfh.open();
	while(rsidfh.getNext(elems)!=LINE_INVALID)
	{
		rsids.insert(elems[0]);
		rsids.insert(elems[1]);
		rsidpairs.push_back(make_pair(elems[0],elems[1]));
	}
	rsidfh.close();
	cerr << "Finish reading " << rsidfile << endl;

	// Hapmap initiation
	map<string,vector<string> > hapmap;
	ColumnReader hapmapfh(hapmapfile);
	
	// Read hapmap file
	cerr << "Reading " << hapmapfile << endl;
	hapmapfh.open();
	while(hapmapfh.getNext(elems)!=LINE_INVALID)
	{
		if (rsids.count(elems[0]))
			hapmap[elems[0]]=elems;
	}
	hapmapfh.close();
	cerr << "Finish reading " << hapmapfile << endl;

	// Output
	cerr << "Reading " << outputfile << endl;
	Writer outputfh(outputfile);

	// Calculate LDR
	float ldr=0,errcode=0;
	outputfh.open();
	for (vector< pair<string,string> >::iterator it=rsidpairs.begin(); it!=rsidpairs.end();it++)
	{
		if (hapmap.find(it->first)!=hapmap.end() && hapmap.find(it->second)!=hapmap.end() )
		{
			LDR2(hapmap[it->first],hapmap[it->second],ldr,errcode);
			(*(outputfh.Printer())) << it->first << "\t" << it->second << "\t" << ldr << "\t" << errcode << endl;
		}
	}
	outputfh.close();
	cerr << "Finish reading " << outputfile << endl;

    return 0;
}
