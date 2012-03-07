#include <iostream>
#include <string>

#include "reader.h"
#include "gzstream.h"

using namespace std;

int main()
{
	Reader f("reader.h");
	f.Open();
	f.Close();
    return 0;
}

