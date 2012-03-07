#include "test.h"
#include "seqUtils.h"

#include <vector>

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
	BedVector beds;
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
	BedVector::iterator it;
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
	BedMap bedmap;
	BedUtils::loadBedFileToMap(bedmap,argv[1]);
	for(BedMap::iterator it=bedmap.begin();it!=bedmap.end();it++)
		for(binBeds::iterator it2=it->second.begin();it2!=it->second.end();it2++)
			for(BedVector::iterator it3=it2->second.begin();it3!=it2->second.end();it3++)
				cout << *it3 << endl;
	
	return 0;
}


