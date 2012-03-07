#include <iostream>
#include <string>
#include <sstream>
using namespace std;

template<typename out_type,typename in_value>

out_type convert(const in_value & t)

{

	stringstream stream;

	stream<<t;//向流中传值

	out_type result;//这里存储转换结果

	stream>>result;//向result中写入值

	return result;

}

int main()
{
	double d;

	string salary;

	string s="12.56";

	d=convert<double>(s);//d等于12.56

	salary=convert<string>(9000.0);//salary等于”9000”
	cout << d << "," <<salary <<endl;
    return 0;
}

