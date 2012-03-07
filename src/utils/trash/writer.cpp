/*****************************************************************************
  writer.cpp
  Last-modified: 05 Mar 2012 01:44:54 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>

#include "writer.h"

using namespace std;

// Constructor
Writer::Writer(void): _outStream(NULL)
{
}

Writer::Writer(const string &outfile):_outFile(outfile),_outStream(NULL)
{
}

// Destructor
Writer::~Writer(void)
{
	Close();
}

// Set file name
void Writer::setFileName(const string &fname)
{
	_outFile=fname;
	Close();
}

//// Open
void Writer::Open(void)
{
	if(_outFile == "stdout")
		_outStream = &cout;
	else
	{
		ofstream tryopen(_outFile.c_str(), ios::out);
		if(!tryopen)
		{
			cerr << "Error: The requested file (" << _outFile << ") could not be opened. Exiting!" << endl;
			cerr << "Take stdout as default output." << endl;
			_outStream = &cout;
		}
		else
		{
			tryopen.close();
			_outStream = new ofstream(_outFile.c_str(), ios::out);
		}
	}
}

// Close
void Writer::Close(void)
{
	if( _outStream != &cout && _outStream!=NULL )
	{
		((ofstream *) _outStream) -> close(); // Close the file with ofstream function
		delete _outStream;
		_outStream=NULL;
	}
}

// Printer
ostream * Writer::Printer(void)
{
	return _outStream;
}

/****************************************************
//Test program:
int main(int argc, char* argv[])
{
	Writer test("test.txt");
	if(test.Open())
	{
		(*test._outStream) << "fsff" <<endl;
	}
    return 0;
}
****************************************************/

