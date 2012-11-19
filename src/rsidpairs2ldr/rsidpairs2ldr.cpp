/*****************************************************************************
  rsidpairs2ldr.cpp
  Last-modified: 19 Nov 2012 02:15:48 PM

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
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <bitset>

#include "common.h"
#include "stringUtils.h"
#include "gzstream.h"

using namespace std;

#define ARRYLEN 409
#define LINELEN ARRYLEN*2 // twobit
#define BITLEN  824 // ceiling(409*2/8)*8 = 824 
#define CUTOFF  0.1


double p11(double A, double B, double C, double D, double x1, double x2, double C_OK, int cnt)
{
	double p11r = (x1+x2)/2;
	double R = A*pow(p11r,3)+B*pow(p11r,2)+C*p11r+D;
	
	if(fabs(x1-x2)<C_OK){
		return p11r;
	}else{
		if(R==0){
			return p11r;
		}else if(R>0){
			return p11(A, B, C, D, x1, p11r, C_OK, cnt+1);
		}else{
			return p11(A, B, C, D, p11r, x2, C_OK, cnt+1);
		}
	}
}

void LDR2(const bitset<BITLEN> &arry1, const bitset<BITLEN> &arry2, char &ldr, char &errcode)
{
	static int ARRY[3][3];
	static int ARRY1[3];
	static int ARRY2[3];
	static float cnt;
	static int i,j, n11, n12, n21, n22, n1a, n2a, na1, na2;
	static double p1, p2, b_p2, c_p1, d_p0, p11r, D, r2;
	static char tmp1, tmp2;

	cnt = 0;
	for(i=0;i<3;i++)
	{
		ARRY1[i]=0;
		ARRY2[i]=0;
		for(j=0;j<3;j++)
			ARRY[i][j]=0;
	}

	cout << arry1.to_string() << endl;
	for(i=0; i<ARRYLEN; i++)
	{
		cout << arry1[i*2+1] << arry1[i*2];
		tmp1 = arry1[i*2]+arry1[i*2+1]*2;
		tmp2 = arry2[i*2]+arry2[i*2+1]*2;

		if(tmp1==3 || tmp2==3){
			continue;
		}else{
			cnt++;
		}
		ARRY[tmp1][tmp2]++;
		ARRY1[tmp1]++;
		ARRY2[tmp2]++;
	}
	cout << endl;

	if(cnt/ARRYLEN < CUTOFF)
	{
		ldr=-1;
		errcode=cnt/ARRYLEN;
		return;
	}
	
	n11 = ARRY[0][0];
	n12 = ARRY[0][1];
	n21 = ARRY[1][0];
	n22 = ARRY[1][1];
	n1a = ARRY1[0];
	n2a = ARRY1[1];
	na1 = ARRY2[0];
	na2 = ARRY2[1];
	p1 = (n1a*2.0 + n2a)/(cnt*2.0);
	p2 = (na1*2.0 + na2)/(cnt*2.0);
	
	if(p1 == 1 || p2 == 1 || p1 == 0 || p2 == 0)
	{
		ldr=-1;
		errcode=-2;
		return;
	}

	b_p2 = (2*cnt*(1-2*p1-2*p2)-n22-2*(2*n11+n12+n21))/(4*cnt);
	c_p1 = (2*cnt*p1*p2-n22*(1-p1-p2)-(2*n11+n12+n21)*(1-2*p1-2*p2))/(4*cnt);
	d_p0 = p1*p2*(2*n11+n12+n21)/(-4*cnt);
	p11r = p11(1, b_p2, c_p1, d_p0, ((2.0*n11+n12+n21)/(2.0*cnt)), ((2.0*n11+n12+n21+n22)/(2.0*cnt)), 0.001, 1);
	D = p11r - p1*p2;
	r2 = D*D/(p1*(1-p1)*p2*(1-p2));
	
	if(r2<0.1)
	{
		ldr=0;
		errcode=0;
		return;
	}
	ldr=char(int(r2*100+0.5));
	errcode = char(int(cnt/ARRYLEN*100+0.5));
	return;
}


//Help function
void Help()
{
    cerr << endl;

    cerr << "Program: rsidpairs2ldr (v" << VERSION << ")" << endl;
    cerr << "Author:  Yunfei Wang (tszn1984@gmail.com), Changning Liu" << endl;
        
    cerr << "Usage:   rsidpairs2ldr  -i i1 i2 -j j1 j2 -H hapmapfile" << endl << endl;
	cerr << "         Output file format: i1_i2_j1_j2_count_data.gz"  << endl << endl;
        
    cerr << "Options:" << endl;
        
    cerr << "Input:"   << endl;
    cerr << "\t-i:--irange           \"i1 i2\". Index range of i. Make sure mini >= minj." << endl << endl;
	cerr << "\t-j:--jrange           \"j1 j2\". Index range of j." << endl << endl;
    cerr << "\t-H:--Hapmap           Hapmap file." << endl << endl;
        
}


int main(int argc, char* argv[])
{
	// parameters
	int mini,minj,maxi,maxj;
	string hapmapfile;
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
		else if((PARAMETER_CHECK("-i", 2, parameterLength)) || (PARAMETER_CHECK("--irange", 8, parameterLength)))
		{
			if ((i+2) < argc) 
			{
				mini = StringUtils::toValue<int>(argv[++i]);
				maxi = StringUtils::toValue<int>(argv[++i]);
			}
		}
		else if((PARAMETER_CHECK("-j", 2, parameterLength)) || (PARAMETER_CHECK("--jrange", 8, parameterLength)))
		{
			if ((i+2) < argc)
			{
				minj = StringUtils::toValue<int>(argv[++i]);
				maxj = StringUtils::toValue<int>(argv[++i]);
			}
		}
		else if((PARAMETER_CHECK("-H", 2, parameterLength)) || (PARAMETER_CHECK("--Hapmap", 8, parameterLength)))
		{
			if ((++i) < argc)
				hapmapfile = argv[i];
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
	
	// Hapmap initiation
	map<int, bitset<BITLEN> > hapmap;
	ifstream hapmapfh;
	string line;
	
	// Read hapmap file
	cerr << "Reading " << hapmapfile << endl;
	hapmapfh.open(hapmapfile.c_str());
	int lnum;

	// Read minj -> maxj ( minj <= mini )
	lnum=minj;  
	hapmapfh.seekg((LINELEN+1)*(lnum-1)); // Go to the first line. (minj<mini)
	for(;lnum<maxi;lnum++)
	{
		getline(hapmapfh,line);
		bitset<BITLEN> bits(line);
		hapmap[lnum]=bits;
	}
	// Read max(maxj,mini) -> maxi
	lnum = mini > maxj ? mini:maxj; // max(maxj,mini)
	hapmapfh.seekg((LINELEN+1)*(lnum-1)); // Go to first of the rest lines
	for(; lnum<maxi; lnum++)
	{
		getline(hapmapfh,line);
		bitset<BITLEN> bits(line);
		hapmap[lnum]=bits;
	}
	hapmapfh.close();
	cerr << "Finish reading " << hapmapfile << endl;

	// Output
	ogzstream outputfh;
	string fname,newfname;
	fname = StringUtils::toString(mini) + "_" + StringUtils::toString(maxi) + "_" + StringUtils::toString(minj) + "_" + StringUtils::toString(maxj);
	cerr << "Starting write " << fname << endl;
	outputfh.open(fname.c_str());
	char ldr,errcode;
	int i,j;
	long int count=0;
	// Calculate LDR
	if (mini>minj)
	{
		for(i=mini;i<maxi;i++)
		{
			for(j=minj;j<maxj;j++)
			{
				LDR2(hapmap[i],hapmap[j],ldr,errcode);
				if(errcode>0)
				{
					count++;
					outputfh.write( (char*)&i, sizeof(i));
					outputfh.write( (char*)&j, sizeof(j));
					outputfh.write( (char*)&ldr,sizeof(ldr));
					outputfh.write( (char*)&errcode,sizeof(errcode));
				}
			}
		}
	}
	else if(mini==minj)
	{
		for(i=mini;i<maxi;i++)
		{
			for(j=minj;j<=i;j++)
			{
				LDR2(hapmap[i],hapmap[j],ldr,errcode);
				if(errcode>0)
				{
					count++;
					outputfh.write( (char*)&i, sizeof(i));
					outputfh.write( (char*)&j, sizeof(j));
					outputfh.write( (char*)&ldr,sizeof(ldr));
					outputfh.write( (char*)&errcode,sizeof(errcode));
				}
			}
		}
	}
	outputfh.close();
	cerr << "Finish writing" << fname << endl;
	newfname = fname+"_"+StringUtils::toString(count)+"_data.gz";
	cerr << "Rename file from " << fname << " to " << newfname << endl;
	rename(fname.c_str(),newfname.c_str());

    return 0;
}
