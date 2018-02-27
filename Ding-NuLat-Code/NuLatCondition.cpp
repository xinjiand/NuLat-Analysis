#include <iostream>
#include <fstream>
using namespace std;
int main()
{
	ofstream fout;
	fout.open("condition.txt",fstream::trunc);
	fout << "lowthreshold[10][10]=" << endl;
	for (int i=0; i<10; i++)
	{
		fout << "0\t0\t0\t0\t0\t0\t0\t0\t0\t0" << endl;
	}
	fout << "highthreshold[10][10]=	"<< endl;
	for (int i=0; i<10; i++)
	{
		fout << "500000\t500000\t500000\t500000\t500000\t500000\t500000\t500000\t500000\t500000" << endl;
	}
	fout << "calibration[10][10]=" << endl;
	for (int i=0; i<10; i++)
	{
		fout << "0\t0\t0\t0\t0\t0\t0\t0\t0\t0" << endl;
	}
	return 0;

}
