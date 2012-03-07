/*****************************************************************************
  reader.h
  Last-modified: 06 Mar 2012 11:26:30 AM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#ifndef READER_H
#define READER_H

#include <string>
#include <sstream>
#include <cstdlib>
using namespace std;


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
    void Open(void);

    // Close
    void Close(void);

protected:
    istream *_inStream;
    string _inFile;
    
};

#endif // READER_H

