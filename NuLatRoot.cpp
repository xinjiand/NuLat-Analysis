#include <iostream>
#include <fstream>
#include <cstdlib> // for exit()
#include <string.h>
#include <vector>
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TH2D.h"
using namespace std;

vector<int> adjust (vector<int> &l);
int maxfind (vector<int>&l);
int minfind(vector<int> &l);
int sum (vector<int>&l, int number1, int number2);
double psd (vector<int>&l);

int main(int argc, char* argv[])
{
	if (argc==1)
	{
		cerr << "Usage:" << argv[0] << "filename[s]\n";
		exit(EXIT_FAILURE);
	}
	ifstream fin;
	TFile* f=new TFile ("analysis.root","recreate");
	int lineMaximum=500;
	char str[lineMaximum];
	const char *d=" '\t'";
	char* p;
	int countline=0;
	int tempdatacount=0;
	int eventnumber=0;
	int eventrow=0;
	int eventchanl=0;
	int eventcol=0;
	bool flag=true;
	int tempcondition[4]={0};
	int totalenergy;
	double psdratio;
	TNtuple ntuple("ntuple","energy spectrum","totalenergy:psd");
	TH2D* psdana=new TH2D("psdana","PSD Analysis",1000,0,2,100,0,1000);
	vector<int> pulse;
	for (int i=1; i<argc; i++)
	{
		fin.open (argv[i]);
		countline=0;
		if (!fin.is_open())
		{
			cerr << "Could not open "<< argv[i]<< endl;
			fin.clear();
		}
		while (fin.getline(str,lineMaximum))
		{	
			if (countline<3)
			cout << "note output"<<endl;
			else
			{
				p=strtok (str,d);
				while (p)
				{
					if (tempdatacount==1)
					{
						tempcondition[0]=atoi(p);
						
					}
					else if (tempdatacount==3)
					{
						tempcondition[1]=atoi(p);
					}
					else if (tempdatacount==4)
                                        {
						tempcondition[2]=atoi(p);
                                        }
					else if (tempdatacount==5)
                                        {
						tempcondition[3]=atoi(p);
                                        }
					else if (tempdatacount==8)
					{
						if (tempcondition[0]!=eventnumber)
						{
							vector<int> adjustedpulse=adjust(pulse);
							int *pulsey=new int[adjustedpulse.size()];
							int *xaxis=new int[adjustedpulse.size()];
							totalenergy=sum(adjustedpulse,0,adjustedpulse.size());
							psdratio = psd (adjustedpulse);
							ntuple.Fill(totalenergy,psdratio);
							psdana->Fill(psdratio,totalenergy);
							for (int j=0; j<adjustedpulse.size(); j++)
							{
								pulsey[j]=adjustedpulse[j];
								xaxis[j]=j;
							}
							TGraph *g=new TGraph(pulse.size(),xaxis,pulsey);
							g->Write();
							delete g;
							pulse.clear();
							eventnumber=tempcondition[0];
							eventrow=tempcondition[1];
							eventcol=tempcondition[2];
							eventchanl=tempcondition[3];
						}
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
		psdana->Write();
		f->Write();
		f->Close();
		fin.clear();
		fin.close();      
	}
	return 0;
}

vector<int> adjust (vector<int> &l)
{
	vector<int> adjust;
	int *array = new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];	
	}
	int minpos=minfind(l);
	int threshold=30;
	int leftzeropos;
	int rightzeropos;
	int adjustoff=0;
	for (int i=0; i<minpos;i++)
	{
		if (-array[minpos-i]<threshold)
		{
			leftzeropos=minpos-i;
			break;
		}		
	}
	for (int i=0; i< l.size()-minpos; i++)
	{
		if (-array[minpos+i]<threshold)
		{
			rightzeropos=minpos+i;
			break;
		}
	}
	int totalcount= l.size()-rightzeropos+leftzeropos;
	adjustoff=(sum(l,0,leftzeropos)+sum(l,rightzeropos,l.size()));
	cout << totalcount <<"\t" << adjustoff <<endl;
	if (totalcount<1)
	totalcount=1;
	int adjustoffset=int (adjustoff/totalcount);
//	cout << adjustoffset << endl;
	for (int i=0; i<l.size(); i++)
	{
		adjust.push_back(-(l[i]-adjustoffset));
	}
	delete []array;
	//cout << "adjust being called" << endl;
	return adjust;
}

int maxfind (vector<int>&l)
{
	int* array=new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];
	}
	int tempvalue=0;
	int maxpos;
	for (int i=1; i<l.size(); i++)
	{
		if (tempvalue<array[i])
		{
			tempvalue=array[i];
			maxpos=i;
		}
		
	}
	delete []array;
//	cout << "max being called" << endl;
	return maxpos;
}

int minfind(vector<int> &l)
{
	int* array=new int[l.size()];
	int minpos;
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
//	cout << "min being called" << endl;
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
	for (int i=0; i< number2-number1+1 ; i++)
	{
		tempsum=tempsum+array[number1+i];
	}
	delete []array;
//	cout << "sum being called" << endl;
	return tempsum;
}

double psd (vector<int>&l)
{
	int maxpos=maxfind(l);
	int totalenergy=sum(l,0,l.size());
	int tailenergy=sum(l,maxpos,l.size());
	cout <<"temp output"<< totalenergy << "\t" << tailenergy << endl;
	double psdana=0.;
	if( totalenergy==0)
	cout << "totol energy become 0" << endl;
	else
	{
//		cout << "psdana put number" << endl;
		psdana=(double (tailenergy)/double (totalenergy));
	}
	return psdana;	
}
