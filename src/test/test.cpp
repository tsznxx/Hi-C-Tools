#include "test.h"
#include "seqUtils.h"
#include "bedUtils.h"

#include <vector>
#include <fstream>
#include <iostream>
#include "BamReader.h"
#include <iostream>
#include "gzstream.h"
#include <bitset>

using namespace std;
/*
bool comp(const Bed &A, const Bed &B)
{
	return (A.chrom <B.chrom || A.start<B.start || A.end<B.end || A.score <B.score);
}*/

void help(void)
{
	cout << "test *.bed" <<endl;
}


int main( int argc, char *argv[] )
{
	// test file fseek.
	ifstream in;
	ofstream out;
	string line;
	char buf[103]; // (412+3)/4
	int i;
	in.open("tmp");
	out.open("hapmap_binary.file");
	map<int,bitset<818> > hapmap;
	i=0;
	while (getline(in,line))
	{
		bitset<818> bits((string)line);
		hapmap[i]=bits;
		i++;
/*		// covert to binary
		for(i=0;i<102;i++)
			buf[i] = (line[i*4]-'0')*64 + (line[i*4+1]-'0')*16 + (line[i*4+2]-'0')*4 + line[i*4+3] - '0';
		// last one
		buf[102]= (line[408]-'0')*64 + 3*16 + 3*4 +3; // X333
		out.write((char*)&buf,103);*/
	}
	for(map<int,bitset<818> >::iterator it=hapmap.begin();it!=hapmap.end();it++)
		cout << it->second.to_string() << endl;

	in.close();
	out.close();
/*	// read binary file
	in.open("hapmap_binary.file");
	while(in.good())
	{
		in.read(buf,103);
		for (i=0;i<103;i++)
		{
			cout << (buf[i]%4+4)%4;
			buf[i] >> 2;
			cout << (buf[i]%4+4)%4;
			buf[i] >> 2;
			cout << (buf[i]%4+4)%4;
			buf[i] >> 2;
			cout << (buf[i]%4+4)%4;
		}
		cout << endl;
	}
	in.close();*/
/*  // Test reading and writing binary and gzip file.
    if (argc<=1)
	{
		cout << argv[0] << " rsid.txt.gz start linenum" << endl;
		exit(0);
	}
	igzstream in;
	ogzstream out;
	int i,j;
	char ldr,errcode;
	ldr='a';
	for (int k=0;k<10;k++)
	{
		bout.Write((char*)&k,sizeof(k));
		bout.Write((char*)&k,sizeof(k));
		bout.Write((char*)&ldr,sizeof(ldr));
		bout.Write((char*)&ldr,sizeof(ldr));
	}
	bout.Close();
	int cnt=0;
	int start=0,end=0;
	if (argc >2)
		start=atoi(argv[2]);
	if (argc >3)
		end=atoi(argv[2])+atoi(argv[3]);
	bin.Open(argv[1],"rb");
	bin.Read((char*)&header1,sizeof(header1));
	cout << header1.mini <<"," << header1.count << endl;
	for(cnt=0;cnt<header1.count;cnt++)
	{
		in.read((char*)&i,sizeof(i));
		in.read((char*)&j,sizeof(j));
		in.read((char*)&ldr,sizeof(ldr));
		in.read((char*)&errcode,sizeof(errcode));

		bin.Read((char*)&i,sizeof(i));
		bin.Read((char*)&j,sizeof(j));
		bin.Read((char*)&ldr,sizeof(ldr));
		bin.Read((char*)&errcode,sizeof(errcode));
		cnt++;
		if (cnt>=start && cnt<end)
		{
			cout << i << "\t" << j << "\t" << int(ldr) << "\t" << int(errcode) << endl;
		}
		else if (cnt>end)
		{
			break;
		}
	}
	bin.Close();*/

/*	if(argc==1)
	{
		help();
		return 0;
	}
	string infile=argv[1];
	cout << infile << endl;
	ColumnReader bedfile(infile);
	bedfile.Open();
	LINE_ELEMS elems;
	LineStatus ls;
	BedVec beds;
	cout << "Starting reading ..." <<endl;
	while( (ls = bedfile.getNext(elems)) != LINE_INVALID)
	{
		Bed tbed(elems);
		//cout << *tbed << endl;
		beds.push_back(tbed);
		//cout << beds.size() << endl;
		elems.clear();
	}
	cout << "Total: " << beds.size() <<endl;

	cout << "Sorting ...." <<endl;

	sort(beds.begin(),beds.end());
	BedVec::iterator it;
	for (it=beds.begin();it!=beds.end()-1;it++)
	{
		cout << *it << endl;
	}

	bedfile.Close();
*/
/* Fasta test
	Fasta fa("chr01","AAACCACACCGTGTTGCAcacacgagccgatcgatcgtagcgcaact");
	fa.print();
	cout << SeqUtils::toRNA(fa.seq) << endl;
	cout << SeqUtils::rc(fa.seq) << endl;
	cout << SeqUtils::GC(fa.seq) << endl;
	cout << SeqUtils::TM(fa.seq) << endl;
	*/

// BedMap test
	// Load File to BedMap
	//BedMap bedmap;
	//BedUtils::loadBedFileToMap(bedmap,argv[1]);
	/*for(BedMap::iterator it=bedmap.begin();it!=bedmap.end();it++)
		for(binBeds::iterator it2=it->second.begin();it2!=it->second.end();it2++)
			for(BedVec::iterator it3=it2->second.begin();it3!=it2->second.end();it3++)
				cout << *it3 << endl;*/
	/*
	string bedfile=argv[2];
	BedVec bedvec; // bed file
	OverlapMap overlaps;

	BedUtils::loadBedFileToVec(bedvec,bedfile);

	for(BedVec::iterator it=bedvec.begin();it!=bedvec.end();it++)
	{
		BedUtils::intersectBed(*it,bedmap,overlaps);

		cout << (*it) << "\t";
		if(overlaps.size())
		{
			cout << "\t" << *(overlaps.begin()->second[0]) << "\t" << overlaps.begin()->first << endl;
		}
		else
			cout << "\t.\t.\t.\t.\t.\t.\t0" << endl;
		overlaps.clear();
	}*/

/*// BedVec test
	BedVec bv;
	Bed tbed;
	for(int i=0;i<10;i++)
	{
		tbed.score++;
		bv.push_back(tbed);
	}
	for(BedVec::iterator it=bv.begin();it!=bv.end();it++)
		cout << *it << endl;
*/

/*// Bamtools test
	BamTools::BamReader br;
	br.Open("test.bam");
	cout << br.IsOpen() <<endl;
	br.Close();
	cout << br.IsOpen();*/
return 0;
}


