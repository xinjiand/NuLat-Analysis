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
using namespace std;

vector<int> adjust (vector<int> &l);
int maxfind (vector<int>&l);
int minfind(vector<int> &l);
int sum (vector<int>&l, int number1, int number2);
double psd (vector<int>&l,int number1, int number2, int number3);
bool pick(int a[10][10],int b[10][10],int c[10][10]);
vector<int> CutShock (vector<int> &l);


int main(int argc, char* argv[])
{
	if (argc==1)
	{
		cerr << "Usage:" << argv[0] << "filename[s]\n";
		exit(EXIT_FAILURE);
	}
	ifstream fin;
	TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",argv[argc-1]);
	string analysisname=argv[argc-1];
	TFile* f=new TFile (analysisfile+rootap,"recreate");
	cout << "root file being created here" << endl;	
	cout << "Input event number needs to be drawed" << endl;
	int eventtestnumber=0;
	cin >> eventtestnumber;
	TString erowchar, echchar,ecolchar;
	/* variable to hold content for each line from data*/
	int lineMaximum=1000;
	char str[lineMaximum];
	const char *d=" '\t'";
	char* p;
	/*holding event information for the current line being analyzed*/
	int countline=0;
	int tempdatacount=0;
	int eventnumber=0;
	int eventrow=0;
	int eventchanl=0;
	int eventcol=0;
	int tempcondition[4]={0}; 	// {eventnumber rolnum colnum channelnum}
	/*vector to hold valid analysis information, being processed event by event and clear after reasonable histogram being created and data storage each event*/	
	vector<int> pulse;
	vector<int> adjustedpulse;
	vector<int> cuttedpulse;
	bool ChBool=false;
	for (int i=1; i<argc-1; i++)
	{
		fin.open (argv[i]);
		countline=0;
		if (!fin.is_open())
		{
			cerr << "Could not open "<< argv[i]<< endl;
			fin.clear();
		}
		cout << argv[i] << "\t file being analysised here" << endl;
		while (fin.getline(str,lineMaximum))
		{	
			/* line by line processing */
			if (countline<3)
			cout << "note output"<<endl;
			else
			{
				p=strtok (str,d);
				while (p)
				{
					/* word by word frome one line processing*/
					if (tempdatacount==1)
					{
						tempcondition[0]=atoi(p);
						if (countline==3)
						{
							eventnumber=tempcondition[0];							
							cout << "event number being read here with value \t" << eventnumber << endl;
						}
						
					}
					else if (tempdatacount==3)
					{
						tempcondition[1]=atoi(p);
						if (countline==3)
						{
							eventrow=tempcondition[1];
							cout << "event row being initialized here with value \t" << eventrow << endl;
						}
					}
					else if (tempdatacount==4)
                                        {
						tempcondition[2]=atoi(p);
						if (countline==3)
						{
							eventcol=tempcondition[2];
							cout << "event col being initialized here with value \t" << eventcol << endl;
						}
                                        }
					else if (tempdatacount==5)
                                        {
						tempcondition[3]=atoi(p);
						if (countline==3)
						{
							eventchanl=tempcondition[3];
							cout << "event channel being initialized here with value \t" << eventchanl << endl;
						}
                                        }
					else if (tempdatacount==8)
					{
						/*After all the pre-process of the data, start deal with the issue from last event*/	
						ChBool = (tempcondition[1]!=eventrow || tempcondition[2]!=eventcol || tempcondition[3]!=eventchanl	);				
						if (eventnumber==eventtestnumber && ChBool)
						{				
							cuttedpulse=CutShock(pulse);					
							adjustedpulse=adjust(cuttedpulse);
							int *pulsey=new int[adjustedpulse.size()];
							int *xaxis=new int[adjustedpulse.size()];						
							for (int j=0; j<adjustedpulse.size(); j++)
							{
								pulsey[j]=adjustedpulse[j];
								xaxis[j]=j;
							}
							erowchar.Form("%d",eventrow) ;
							echchar.Form("%d",eventchanl) ;
							ecolchar.Form("%d",eventcol) ; 
							TGraph *g=new TGraph(pulse.size(),xaxis,pulsey);
							g->SetTitle (erowchar+" "+ecolchar+" "+echchar);								
							g->Write();
							delete g; 								
							pulse.clear();
							adjustedpulse.clear();
							delete []pulsey;
							delete []xaxis;	
							/*new event condition being initialized here*/
							eventnumber=tempcondition[0];
							eventrow=tempcondition[1];
							eventcol=tempcondition[2];
							eventchanl=tempcondition[3];	
						}
						else if(tempcondition[1]!=eventrow || tempcondition[2]!=eventcol || tempcondition[3]!=eventchanl || tempcondition[0]!=eventnumber)
						{
							pulse.clear();
							eventnumber=tempcondition[0];
							eventrow=tempcondition[1];
							eventcol=tempcondition[2];
							eventchanl=tempcondition[3];	
						}
						/*After all the pre-process of the data, start deal with the issue from different channel but same event*/
					}
					else if (tempdatacount>8)
					pulse.push_back(atoi(p));
					p=strtok(NULL,d);
					tempdatacount++;
				}
				tempdatacount=0;
			}
			countline++;
		}
		fin.clear();
		fin.close();
	}
	f->Write();
	f->Close();
	return 0;
}
/* function to be called to flip the pulse and ajust the pulse to go to the base line*/
vector<int> adjust (vector<int> &l)
{
	/*Current algorithum only take the left side of the pulse into account, find the average of them as the reference of basline shift.
	  define the limit of leftzeropose to be 50, so if pulse is actually just noise, the peakpos might be very small, then base line will compute at least 50 point.
	  minpos is actually the position of the pulse peak	*/	
	vector<int> adjust;
	int *array = new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];	
	}
	int minpos=minfind(l);
	//cout << "minimum find" << array[minpos] << endl;
	int maxpos=maxfind(l);
	int maxamp=array[maxpos];
	//cout << "maximum find" << maxamp << endl;
	int threshold=20;
	int leftzeropos=0;
	/*if (minpos<51)
	leftzeropos=50;
	else
	leftzeropos=minpos-50;*/
	int adjustoff=0;
	for (int i=0; i<minpos;i++)
	{
		if (array[minpos-i]>(maxamp-threshold))
		{
			leftzeropos=minpos-i;
			break;
		}		
	}
	if (leftzeropos<50)
	leftzeropos=50;	
	int totalcount=leftzeropos;
	adjustoff=sum(l,0,leftzeropos);
	int adjustoffset=int (adjustoff/totalcount);
	for (int i=0; i<l.size(); i++)
	{
		adjust.push_back(-(l[i]-adjustoffset)); //adjust the baseline and flip the pulse
	}
	delete []array;
	return adjust;
}

vector<int> CutShock (vector<int> &l)
{
	/*Function being called to cut all the shock signal which will mimic the pulse or confused peak finding mechanism*/	
	vector<int> cutpulse;
	int *array = new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];	
	}
	int threshold=50; // creat a threshold , so all the data cross +50 and -50 will take into account.
	int arraypos=0;
	int lefttest=0;
	int righttest=0;
	vector<int> spik;
	for (int i=0; i<l.size(); i++)
	{
		if (array[i]>threshold || array[i] < -threshold)
		spik.push_back(i);
	} // store all the position that pulse higher the threhold
	for (int i=0; i<spik.size(); i++)
	{
		arraypos=spik[i];
		lefttest=abs (array[arraypos]-array[arraypos-1]);
		righttest=abs (array[arraypos]-array[arraypos+1]);	// search all the potential spike, identify them by compare to the neighbour, if both neighbour have 50 different, then identify them as spike
		if (lefttest>50&&righttest>50)
		array[arraypos]=array[arraypos-1]; //rewrite the spike by make them equal to their neighbour.
	}
	for (int i=0; i<l.size(); i++)
	{
		cutpulse.push_back(array[i]); //rewrite the cutted spike pulse back to the pulse
	}
	delete []array;
	return cutpulse;
}


int maxfind (vector<int>&l)
{
	int* array=new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];
	}
	int tempvalue=0;
	int maxpos=0;
	for (int i=1; i<l.size(); i++)
	{
		if (tempvalue<array[i])
		{
			tempvalue=array[i];
			maxpos=i;
		}
		
	}
	delete []array;
	return maxpos;
}

int minfind(vector<int> &l)
{
	int* array=new int[l.size()];
	int minpos=0;
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];
	}
	int tempdata=array[0];
	for (int i=0; i<l.size(); i++)
	{
		if (tempdata>array[i])
		{
			tempdata=array[i];
			minpos=i;
		}
	}
	delete []array;
	return minpos;
}

int sum (vector<int>&l, int number1, int number2)
{
	int* array=new int[l.size()];
        for (int i=0; i<l.size(); i++)
        {
                array[i]=l[i];
        }
	int tempsum=0;
	if (number1==number2)
	{
		delete []array;		
		return 0;
	}
	else 
	{
		for (int i=0; i< number2-number1 ; i++)
		{
			tempsum=tempsum+array[number1+i];
		}
		delete []array;
		return tempsum;
	}
}

double psd (vector<int>&l,int number1,int number2, int number3)
{
	int* array=new int[l.size()];
        for (int i=0; i<l.size(); i++)
        {
                array[i]=l[i];
        }
	int peakleft=number1;
	int peak=number2;
	int peakright=number3;
	int totalenergy=sum(l,peakleft,peakright);
	int tailenergy=sum(l,peak,peakright);
	double psdana=0.;
	psdana=(double (tailenergy)/double (totalenergy));
	return psdana;	
}

bool pick(int a[10][10],int b[10][10], int c[10][10])
{
	bool eventveto=false;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (c[i][j]<a[i][j] || c[i][j]>b[i][j])
				eventveto=true;
		}
	} 
	return eventveto;
}




/*
if (j==0)	
					{
						cout << event[j] << "\t" <<xID[j] <<"\t" << yID[j] << "\t" << cubeID[j] << "\t" << energyspec[j] << "\t" << energypeak[j] << endl;
						fout << event[j] << "\t" <<xID[j] <<"\t" << yID[j] << "\t" << cubeID[j] << "\t" << energyspec[j] << "\t" << energypeak[j] << endl;
					}
								else if (j==cubeID.size()-1)	
					{
						cout << xID[j] <<"\t" << yID[j] << "\t"  << event[j] << "\t"<< cubeID[j] << "\t" << energyspec[j] << "\t" << energypeak[j] << endl;
						fout << xID[j] <<"\t" << yID[j] << "\t"  << event[j] << "\t"<< cubeID[j] << "\t" << energyspec[j] << "\t" << energypeak[j] << endl;
					}
								else
					{
						cout << xID[j] <<"\t" << yID[j] << "\t"  << event[j] << "\t"<< cubeID[j] << "\t" << energyspec[j] << "\t" << energypeak[j] << endl;
						fout << xID[j] <<"\t" << yID[j] << "\t"  << event[j] << "\t"<< cubeID[j] << "\t" << energyspec[j] << "\t" << energypeak[j] << endl;
					}								
								}


*/
