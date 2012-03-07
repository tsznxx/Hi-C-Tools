/*****************************************************************************
  test.cpp
  Last-modified: 05 Mar 2012 04:15:06 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>
#include "seqUtils.h"

using namespace std;

int main(int argc, char* argv[])
{
	Fasta fa("chr01","AAACCACACCGTGTTGCAcacacgagccgatcgatcgtagcgcaact");
	fa.print();
	cout << SeqUtils::toRNA(fa.seq) << endl;
	cout << SeqUtils::rc(fa.seq) << endl;
	
    return 0;
}

