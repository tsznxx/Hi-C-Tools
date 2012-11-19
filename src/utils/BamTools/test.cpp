/*****************************************************************************
  test.cpp
  Last-modified: 29 May 2012 04:39:43 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>

#include "BamReader.h"

using namespace std;

int main(int argc, char* argv[])
{
	BamReader br;
	br.Open("test.bam");
	cout << br.IsOpen() << endl;
	br.Close();
    return 0;
}

