#include <iostream>
#include <fstream>
#include <cstdlib> // for exit()
#include <string.h>
#include <vector>
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"
#include "anafunction.h"
using namespace std;
typedef vector<int> Row;
typedef vector<Row> Matrix;
typedef vector<double> RowPsd;
typedef vector<RowPsd> MatrixPsd;


int main(int argc, char* argv[])
{
	/*initalizae the variable using to reading different files*/
	int lineMaximum=5000;
	char str[lineMaximum];
	const char *d=" '\t'";
	char* p;
	/*initialize the option this analysis needs to do
		
		T: Just for test, will creat 100 event 2D energy map
		C: Threshold cut and then output 2D historgram of energy map and 1D spectrum for each channel
		RC: Recalibration the energy just be reasonable number (needs to modified incase the number out of range)	
	*/
	/*string option;
	cout << " Choose which command you need : AB pick (AB) ; B event analysis(B); Energy Mapping (m) ; Timing check for AB (T); AB pick out with timing check (ABT) " << endl;
	cin >> option;*/
	/* open different file ready to be input
		
		energy	psd	peak	timing	cubeID	event
		
	*/
	string option_test;
	char analysis_new;

	string analysisname=argv[1];
	ifstream fincondition;
	fincondition.open ("condition.txt");
	ifstream finenergyA;
	finenergyA.open (analysisname+" energyA.txt");
	ifstream finpsdA;
	finpsdA.open (analysisname+" psdA.txt");
	ifstream finpeakA;
	finpeakA.open (analysisname+" peakA.txt");
	ifstream fintimingA;
	fintimingA.open (analysisname+" timingA.txt");
	ifstream finenergyB;
	finenergyB.open (analysisname+" energyB.txt");
	ifstream finpsdB;
	finpsdB.open (analysisname+" psdB.txt");
	ifstream finpeakB;
	finpeakB.open (analysisname+" peakB.txt");
	ifstream fintimingB;
	fintimingB.open (analysisname+" timingB.txt");	
	ifstream fincube;
	fincube.open (analysisname+" cubeID.txt");
	ifstream finevent;
	finevent.open (analysisname+" event.txt");
	ifstream fintimingAB;
	fintimingAB.open (analysisname+" timeAB.txt");
	/* reading event file and create a event vector to be used and event cateory in the future, also deine the number of line needs to be read from other file by counting event number*/	
	vector<int> eventID;
	int n=0; 
	while (finevent.getline(str,lineMaximum))
	{
		p=strtok (str,d);		
		eventID.push_back (atoi(p));
		//cout << p << endl;
		n++;
	}
	/*Fill the map matrix*/
	Matrix map(10,Row(10));
	ifstream finmap;
	finmap.open ("map.txt");
	int a=0;
	int b=0;

	while (finmap.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			map[a][b]=atoi(p);
			b++;
			p=strtok(NULL,d);
		}
		b=0;
		a++;
	}

	/* output a matrix mapping being reading here*/
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			cout << map[i][j]<< "\t" ;
		}
		cout << endl;
	}
	/*create different matrix to store data in the future and also 2D matrix as for calibration and threshold cut*/
	Matrix energyA(n*10,Row(10));	
	MatrixPsd psdA(n*10,RowPsd(10));
	Matrix peakA(n*10,Row(10));
	Matrix timingA(n*10,Row(10));
	Matrix energyB(n*10,Row(10));	
	MatrixPsd psdB(n*10,RowPsd(10));
	Matrix peakB(n*10,Row(10));
	Matrix timingB(n*10,Row(10));	
	Matrix cube(n*10,Row(10));
	Matrix timingAB(n*10,Row(10));
	
	int lowthreshold[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};
	
	int highthreshold[10][10]={	{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000},
								{2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000,	2000}};	

	float calibration[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};	
	
	/*Start reading data from different file*/
	int i=0;	
	int j=0;
	while (fincondition.getline(str,lineMaximum))
	{
		
		p=strtok (str,d);
		if (i<11 && i>0)			
		{	while (p)
			{	
				lowthreshold[i-1][j]=atoi(p);
				j++;			
				p=strtok(NULL,d);	
			}
		}
		else if (i>11 && i < 22)			
		{	while (p)
			{	
				highthreshold[i-12][j]=atoi(p);
				j++;			
				p=strtok(NULL,d);	
			}
		}
		else if (i >22 && i<33)			
		{	while (p)
			{	
				calibration[i-23][j]=atoi(p);
				j++;			
				p=strtok(NULL,d);	
			}
		}
		j=0;
		i++;
		//cout << i << endl;
	}
	cout << "finish condition" << endl;

	i=0;
	while (finenergyA.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{	
			energyA[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpeakA.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{			
			peakA[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpsdA.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			psdA[i][j]=atof(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (fintimingA.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			timingA[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
		i=0;
	while (finenergyB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{	
			energyB[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpeakB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{			
			peakB[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpsdB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			psdB[i][j]=atof(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (fintimingB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			timingB[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (fincube.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			cube[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
       	i=0;
	while (fintimingAB.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			timingAB[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
        //printMatrix (timingAB);
	
	/*choice of using function*/
	do
	{
		cout<<"input the function you want to call : AB(AB-Event), B (B-Event), ABT (ABT-Event), one (One-Event), Radon(Radon calibration)"<<endl;
		cin >> option_test;
		if (option_test=="AB")
			ABEvent (argv[1], eventID, n, energyA, energyB, peakA, peakB, timingA, timingB, psdA, psdB, timingAB);
		else if (option_test=="B")
			BEvent (argv[1], eventID, n, energyA, energyB, peakA, peakB, timingA, timingB, psdA, psdB, timingAB);
		else if (option_test=="ABT")
			ABtEvent (argv[1], eventID, n, energyA, energyB, peakA, peakB, timingA, timingB, psdA, psdB, timingAB);
                else if (option_test=="one")
                        OneEvent (argv[1], eventID, n, energyA, energyB, peakA, peakB, timingA, timingB, psdA, psdB, timingAB);
		else if (option_test=="Radon")
                        RadonCalibration (argv[1], eventID, n, energyA, energyB, peakA, peakB, timingA, timingB, psdA, psdB, timingAB);


		cout << "Do you need to do another analysis (y/n)" << endl;
		cin >> analysis_new ;
	} while (analysis_new=='y');

        /*Choice AB: doing event search requires both A and B are geometrically well defined.*/
	
        
       /*Choice ABT: doing event search requires both A and B are geometrically well defined and also their timing information also being used as veto requirement.*/



    
	/*create a root file that recreate a root file hold data into different histogram*/
	cout << "good analysis " << n << " event being"<<endl;
	fincondition.clear();
	fincondition.close();
	finenergyA.clear();
	finenergyA.close();
	finpsdA.clear();
	finpsdA.close();
	finpeakA.clear();
	finpeakA.close();
	fintimingA.clear();
	fintimingA.close();
	finenergyB.clear();
	finenergyB.close();
	finpsdB.clear();
	finpsdB.close();
	finpeakB.clear();
	finpeakB.close();
	fintimingB.clear();
	fintimingB.close();
	fintimingAB.clear();
	fintimingAB.close();
	fincube.clear();
	fincube.close();
	finevent.clear();
	finevent.close();
	return 0;
}
