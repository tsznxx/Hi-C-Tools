/*****************************************************************************
  fileUtils.cpp
  Last-modified: 07 Mar 2012 11:46:42 AM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>

#include "fileUtils.h"
using namespace std;


//************************************************
// File reader class
//************************************************

// Constructor
Reader::Reader(const string &fname):_inFile(fname)
{}

// Destructor
Reader::~Reader(void)
{
	Close();
}

// Set file name
void Reader::setFileName(const string &fname)
{
	if(_inStream) // close existing file handle
		Close();
	_inFile=fname;
}
	
// Open
void Reader::Open(void)
{
	// stdin
	if (_inFile == "stdin")
	{
		_inStream = &cin;
	}
	else
	{
		size_t foundPos;
		foundPos = _inFile.find_last_of(".gz");
		// GZ file
		if (foundPos == _inFile.size() - 1)
		{
			igzstream tryopen(_inFile.c_str(), ios::in);
			if ( !tryopen )
			{
				cerr << "Error: The requested file (" << _inFile << ") could not be opened. Exiting!" << endl;
				exit(1);
			}
			else
			{
				tryopen.close();
				_inStream = new igzstream(_inFile.c_str(), ios::in);
			}
		}
		// Not GZ file
		else
		{
			ifstream tryopen(_inFile.c_str(), ios::in);
			if ( !tryopen )
			{
				cerr << "Error: The requested file (" << _inFile << ") could not be opened. Exiting!" << endl;
				exit(1);
			}
			else
			{
				tryopen.close();
				_inStream = new ifstream(_inFile.c_str(), ios::in);
			}
		} // GZ file 
	} // stdin
}

// Close
void Reader::Close(void)
{
	if (_inFile != "stdin" )
	{
		if (_inStream!=NULL)
		{
			((ifstream *) _inStream)-> close();
			delete _inStream;
			_inStream=NULL;
		}
	}
}


//************************************************
// File writer class
//************************************************

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

//************************************************
// Column file reader
// ************************************************

// Constructor
ColumnReader::ColumnReader(const string &fname, const string delimiter): Reader(fname),_colNum(20),_delimiter(delimiter)
{}

// Destructor : call the destructor of Reader automatically
ColumnReader::~ColumnReader(void)
{}

LineStatus ColumnReader::getNext(LINE_ELEMS &elems,bool withheader)
{
	string     curline;

	elems.clear(); // clear the elems before use.
	elems.reserve(_colNum);

	while(_inStream->good())
	{
		getline(*_inStream,curline);
		if (curline.find_first_of("#")==0 && withheader) // suppress header information by set withheader = false
		{
			elems.push_back(curline);
			return LINE_HEADER;
		}
		else
		{
			Tokenize(curline,elems,_delimiter);
			if (elems.size())
			{
				_colNum=elems.size();
				return LINE_VALID;
			}
		}
	}
	return LINE_INVALID;
}
