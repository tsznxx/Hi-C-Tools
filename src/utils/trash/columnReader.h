/*****************************************************************************
  columnReader.h
  Last-modified: 06 Mar 2012 03:47:55 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#ifndef COLUMNREADER_H
#define COLUMNREADER_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include "gzstream.h"
#include "reader.h"
#include "stringUtils.h"

using namespace std;

// Line status defination for ColumnReader
enum LineStatus
{
    LINE_INVALID=-1,
    LINE_HEADER=0,
    LINE_BLANK=1,
    LINE_VALID=2
};

typedef vector<string> LINE_ELEMS;

///////////////////////////////////////
// ColumnReader Class
///////////////////////////////////////

class ColumnReader: public Reader
{
public:
    // Constructor
    ColumnReader(const string &fname="stdin", const string delimiter = "\t");
    
    // Destructor
    ~ColumnReader();

	// Get next column line
	LineStatus getNext(LINE_ELEMS &elems, bool withheader = false);

protected:
	// Data
	int _colNum;
	string _delimiter;
};

#endif //COLUMNREADER_H

