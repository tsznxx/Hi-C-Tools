#include "columnReader.h"

using namespace std;

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

