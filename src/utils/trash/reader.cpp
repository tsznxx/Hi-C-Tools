#include <iostream>
#include "gzstream.h"
#include "reader.h"

using namespace std;

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


