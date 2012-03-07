/*****************************************************************************
  stringUtils.h
  Last-modified: 16 Feb 2012 08:56:19 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>

using namespace std;

// ToString
template <typename T>
inline
string ToString(const T &value)
{
    stringstream ss;
    ss<< value;
    return ss.str();
}

// ToValue
template <typename T>
inline
T ToValue( const string &str)
{
    stringstream stream; //for temporal storage
    T value; // for output value
    stream<<str; // read into t
    stream>>value;
    return value;
}

// Convert for basic types
template<class outT,class inT>
outT convert(const inT & t)
{
    stringstream stream;
    outT result;//这里存储转换结果
    stream<<t;//向流中传值
    stream>>result;//向result中写入值
    return result;

}

// Tokenize
inline 
void Tokenize(const string &str, vector<string> &elems, const string delimiter="\t")
{
    char *tok;
    char cchars[str.size()+1]; // chars to store str
    char *cstr= &cchars[0];    // pointer to cchars
    strcpy(cstr,str.c_str());

    tok=strtok(cstr,delimiter.c_str());
    while( tok!=NULL )
    {
        elems.push_back(tok);
        tok=strtok(NULL,delimiter.c_str());
    }
}

template<typename T>
inline 
void Tokenize(const string &str, vector<T> &elems, const string delimiter="\t")
{
    char *tok;
    char cchars[str.size()+1];
    char *cstr=&cchars[0];
    strcpy(cstr,str.c_str());

    tok=strtok(cstr,delimiter.c_str());
    while (tok!=NULL)
    {
        elems.push_back(convert<T>(tok));
        tok=strtok(NULL,delimiter.c_str());
    }
}


#endif //STRINGUTILS_H
