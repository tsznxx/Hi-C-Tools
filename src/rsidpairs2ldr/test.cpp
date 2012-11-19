/*****************************************************************************
  test.cpp
  Last-modified: 09 Nov 2012 03:22:50 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

int main(int argc, char* argv[])
{
	map< int, vector<int> > myvec;
	vector<int> ints;

	for (int i=0;i<3;i++)
	{
		ints.push_back(i);
		ints.push_back(i*i);
		myvec[i]=ints;
		ints.clear();
	}
	for( map< int, vector<int> >::iterator it=myvec.begin();it!=myvec.end();it++)
	{
		cout<< it->first << "\t" << it->second[0] << "\t" << it->second[1] << endl;
	}
    return 0;
}

