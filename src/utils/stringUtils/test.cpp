#include <iostream>
#include <string>
#include "stringUtils.h"

using namespace std;

int main()
{
	string a="10.2";
	string b="100";
	string c="10\t20\t30";
	vector<int> lst;
	int i=8;
	float j=0.24;
	cout <<ToString(i) << ToString(j) <<endl;
	i=ToValue<int>(b);
	j=ToValue<float>(a);
	cout << i << j << ToString(i) << ToString(j) << endl;
    Tokenize(c,lst);
	cout <<lst[0] << "," <<lst[1] <<endl;
	return 0;
}

