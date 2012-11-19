#include <cstdlib>
#include <vector>


#include "common.h"
#include <iostream>
#include "gzstream.h"
#include "stringUtils.h"

using namespace std;

int main( int argc, char *argv[] )
{
	if (argc<=1)
	{
		cout << argv[0] << " i1_i2_j1_j2_data.gz startline linenum" << endl;
		exit(0);
	}
	igzstream in;
	ogzstream out;
	int i,j;
	char ldr,errcode;
	ldr='a';
	int cnt=0;
	int start=0,end=0;

	if (argc >2)
		start=atoi(argv[2]);
	if (argc >3)
		end=atoi(argv[2])+atoi(argv[3]);
	in.open(argv[1]);
	vector<string> elems;
	StringUtils::tokenize(argv[1],elems,"_");
	cnt=StringUtils::toValue<int>(elems[4])+1;
	end = end <= cnt ? end : cnt;
	for(int k=1;k<end;k++)
	{
		in.read((char*)&i,sizeof(i));
		in.read((char*)&j,sizeof(j));
		in.read((char*)&ldr,sizeof(ldr));
		in.read((char*)&errcode,sizeof(errcode));
		if (k>=start)
		{
			cout << i << "\t" << j << "\t" << int(ldr) << "\t" << int(errcode) << endl;
		}
	}
	in.close();
	return 0;
}
