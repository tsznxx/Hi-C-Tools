/*****************************************************************************
  fileUtils.h
  Last-modified: 15 Nov 2012 04:28:05 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <cstdlib>

#include "gzstream.h"
#include "stringUtils.h"

using namespace std;

//************************************************
// File reader class
//************************************************

class Reader
{
public:
    // Constructor
    Reader(const string &fname="stdin");

    // Destructor
    ~Reader();

    // Set file name
    void setFileName(const string &filename="stdin");
        
    // Open
    void open(void);

    // Close
    void close(void);

protected:
    istream *_inStream;
    string _inFile;
        
};

//************************************************
// File writer class
//************************************************

class Writer
{
public:
    // Constructor
    Writer(void);
    Writer(const string &fname);

    // Destructor
    ~Writer(void);

    // Set file name;
    void setFileName(const string &fname);

    // Open
    void open(void); // may extended to binary open

    // Close
    void close(void);

    // Printer
    ostream * Printer(void);

	// overload <<
/*	template<class T>
	friend ostream& operator<<(ostream& os, const T &item)
	{
		os = _outStream;
		os << item;
		return os;
	}*/

protected:
    string   _outFile;
    ostream* _outStream;

};


//************************************************
// Column file reader
//************************************************

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


#endif //FILEUTILS_H

