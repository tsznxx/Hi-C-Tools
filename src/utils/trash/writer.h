/*****************************************************************************
  writer.h
  Last-modified: 05 Mar 2012 01:43:59 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#ifndef WRITER_H
#define WRITER_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

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
	void Open(void); // may extended to binary open

	// Close
	void Close(void);

	// Printer
	ostream * Printer(void);

protected:
	string   _outFile;
	ostream* _outStream;

};

#endif //WRITER_H

