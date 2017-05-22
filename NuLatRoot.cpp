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
#include "TH1D.h"
#include "TString.h"
using namespace std;

vector<int> adjust (vector<int> &l);
int maxfind (vector<int>&l);
int minfind(vector<int> &l);
int sum (vector<int>&l, int number1, int number2);
double psd (vector<int>&l,int number1, int number2, int number3);

int main(int argc, char* argv[])
{
	if (argc==1)
	{
		cerr << "Usage:" << argv[0] << "filename[s]\n";
		exit(EXIT_FAILURE);
	}
	ifstream fin;
	TFile* f=new TFile ("analysis.root","recreate");
	ofstream fout;
	fout.open("analysis.txt",fstream::trunc);

	int lineMaximum=1000;
	{
		int i=10;
		TString test;
		int x=13/8;
		cout << x << endl;
		test.Form ("%d",i);
		TString test2="babababa";
		test2=test2+" "+test;
		cout << test2 <<endl;
	}
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
	int tempcondition[4]={0}; // {eventnumber rolnum colnum channelnum}
	int totalenergy=0;
	int peakamp=0;
	int valleyamp=0;
	double psdratio=0.;
	int graphnumber=0;
	int channelnum[64]={0};
	TString rowstr="row ";	
	TString colstr="col ";
	TString chanlstr="channel ";
	TString energystr="energy ";
	TString peakstr="peak ";
	TString psdstr="psd ";
	TH1D peakhist[64];
	TH1D integralhist[64];
	TH1D psdhist[64];	
	for (int i=0; i<64; i++)
	{
		TString rowname,colname,chanlname;
		int rowcount=i/64;
		int colcount=i/8;
		int channelcount=i%8;
	//	cout << rowcount <<" "<<colcount < " " << channelcount <<endl;
		cout << rowcount << colcount << channelcount << endl;	
		rowname.Form ("%d",rowcount);
		colname.Form ("%d",colcount);
		chanlname.Form ("%d",channelcount);
		TString titlename = rowstr+" "+rowname+" "+colstr+" "+colname+" "+chanlstr+" "+chanlname;
		TString histname = rowname+" "+colname + " "+chanlname;
		peakhist[i]=TH1D(peakstr+histname,titlename,2500,0,2500);	
		integralhist[i]=TH1D(energystr+histname,titlename,100000,0,100000);
		psdhist[i]=TH1D(psdstr+histname,titlename,1000,0,1);
	}
	TString test="energy spectrum";
	TNtuple ntuple("ntuple","energy spectrum","totalenergy:psd:peak");
	TH2D* psdana=new TH2D("psdana","PSD Analysis",1000,0,2,100,0,1000);
	TH1D* peakspec=new TH1D("peakspec","Peak amp",2500,0,2500);
	vector<int> pulse;
	vector<int> event;
	vector<int> energyspec;
	vector<int> energypeak;
	vector<int> psdanalysis;
	vector<int> row;
	vector<int> col;
	vector<int> channel;
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
						if (countline==3)
						{
							eventnumber=tempcondition[0];
							cout << "event number being initialized here with valuse \t" << eventnumber << endl;
						}
						
					}
					else if (tempdatacount==3)
					{
						tempcondition[1]=atoi(p);
						if (countline==3)
						{
							eventrow=tempcondition[1];
							cout << "event row being initialized here with valuse \t" << eventrow << endl;
						}
					}
					else if (tempdatacount==4)
                                        {
						tempcondition[2]=atoi(p);
						if (countline==3)
						{
							eventcol=tempcondition[2];
							cout << "event col being initialized here with valuse \t" << eventcol << endl;
						}
                                        }
					else if (tempdatacount==5)
                                        {
						tempcondition[3]=atoi(p);
						if (countline==3)
						{
							eventchanl=tempcondition[3];
							cout << "event channel being initialized here with valuse \t" << eventchanl << endl;
						}
                                        }
					else if (tempdatacount==8)
					{
						if (tempcondition[0]!=eventnumber || tempcondition[1]!=eventrow || tempcondition[2]!=eventcol || tempcondition[3]!=eventchanl)
						{
						//	cout << "tempcondition now is" << tempcondition[0] << "\t while event number is "<< eventnumber << endl;
							vector<int> adjustedpulse=adjust(pulse);
							int *pulsey=new int[adjustedpulse.size()];
							int *xaxis=new int[adjustedpulse.size()];
							int peakpos=maxfind(adjustedpulse);
							int valleypos=minfind(adjustedpulse);
							for (int j=0; j<adjustedpulse.size(); j++)
							{
								pulsey[j]=adjustedpulse[j];
								xaxis[j]=j;
							}
							int threshold=15;
							int leftzeropos=0;
							int rightzeropos=0;
							for (int j=0; j<peakpos;j++)
							{
								if (pulsey[peakpos-j]<threshold)
								{
									leftzeropos=peakpos-j;
									break;
								}		
							}
							for (int j=0; j< adjustedpulse.size()-peakpos; j++)
							{
								if (pulsey[peakpos+j]<threshold)
								{
									rightzeropos=peakpos+j;
									break;
								}
							}
							fout << "peakleft=" << leftzeropos <<"\t peak=" << peakpos << "\t peakright=" << rightzeropos << endl;
							peakamp=pulsey[peakpos];
							valleyamp=pulsey[valleypos];
							int peakfullwidth=rightzeropos-leftzeropos;
							int deltapeak=peakamp-valleyamp;
							if (deltapeak < 40 || peakfullwidth <10)
							{
								pulse.clear();
								cout << "Event number=" << eventnumber << "\t just noise" << endl;
								eventnumber=tempcondition[0];
								eventrow=tempcondition[1];
								eventcol=tempcondition[2];
								eventchanl=tempcondition[3];								
								break;
							}
							else
							{
								int histcount=eventchanl+eventcol*8+eventrow*64;								
								channelnum[histcount]++;								
								totalenergy=sum(adjustedpulse,leftzeropos,rightzeropos);
								psdratio = psd (adjustedpulse,leftzeropos,peakpos,rightzeropos);
								psdana->Fill(psdratio,totalenergy);								
								ntuple.Fill(totalenergy,psdratio,peakamp);	
								peakhist[histcount].Fill(peakamp);
								integralhist[histcount].Fill(totalenergy);
								psdhist[histcount].Fill(psdratio);					
								TGraph *g=new TGraph(pulse.size(),xaxis,pulsey);
								g->Write();
								delete g;
								pulse.clear();
								event.push_back(eventnumber);
								energyspec.push_back(totalenergy);
								energypeak.push_back(peakamp);
								psdanalysis.push_back(psdratio);
								row.push_back(eventrow);
								col.push_back(eventcol);
								channel.push_back(eventchanl);								
								eventnumber=tempcondition[0];
								eventrow=tempcondition[1];
								eventcol=tempcondition[2];
								eventchanl=tempcondition[3];
								fout << "Event=" << eventnumber<<"\t Graph number="<<graphnumber <<" \t Peak amp=" << peakamp << "\t pulse size=" << adjustedpulse.size() << "\t total energy=" << totalenergy <<endl;	
								graphnumber++;
							}
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
		
		for (int z=0; z<64; z++)
		{	
			if (channelnum[z]!=0)			
			{
				peakhist[z].Write();
				integralhist[z].Write();
				psdhist[z].Write();
			}
		}		
		psdana->Write();
		f->Write();
		f->Close();
		fin.clear();
		fin.close();      
		fout.close();
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
	int threshold=5;
	int leftzeropos=0;
	int rightzeropos=0;
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
//	cout << totalcount <<"\t" << adjustoff <<endl;
	if (totalcount<1)
	{
		totalcount=1;
		cout << "have to adjust total count to be 1 here" << endl;	
	}
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
//	cout << "max being called" << endl;
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
	for (int i=0; i< number2-number1 ; i++)
	{
		tempsum=tempsum+array[number1+i];
	}
	delete []array;
//	cout << "sum being called" << endl;
	return tempsum;
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
	// cout <<"temp output"<< totalenergy << "\t" << tailenergy << endl;
	double psdana=0.;
//	cout << "psdana put number" << endl;
	psdana=(double (tailenergy)/double (totalenergy));
	return psdana;	
}
