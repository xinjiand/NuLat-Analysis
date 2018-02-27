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
int main(int argc, char* argv[])
{
	/*program running structure with n term
	0	 	1---n-2		n-1
	exe 	data file  	output-file-name
	*/
	if (argc==1)
	{
		cerr << "Usage:" << argv[0] << "filename[s]\n";
		exit(EXIT_FAILURE);
	}
	ifstream fin;
	TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",argv[argc-1]);
	TFile* f=new TFile (analysisfile+rootap,"recreate"); //summary root file being created, every time code being runned, file being recreated.
	cout << "root file being created here" << endl;
	
	bool eventcheck=false;
	bool channelcheck=false;
	bool thresholdcheck=false;
	string graphoption;
	//cout << "input event-channel type or channel option---eg----c" << endl;
	//cin >> graphoption;
	graphoption="f";
	vector<int> eventplot;
	int eventplotinput;
	char plotnew;
	vector<int> channelplot;
	int channelplotinput;
	char channelnew;
	bool drawoption=false;	
	int threshold;
	char thresholdflag;
	if (graphoption=="eg")
	{
		do
		{
			cout<<"input the event you need to plot"<<endl;
			cin >> eventplotinput;
			eventplot.push_back(eventplotinput);
			cout << "Do you need input another event to be plot (y/n)" << endl;
			cin >> plotnew ;
		} while (plotnew=='y');

		do
		{
			cout<<"input the channel you need to plot"<<endl;
			cin >> channelplotinput;
			channelplot.push_back(channelplotinput);
			cout << "Do you need input another channel to be plot (y/n)" << endl;
			cin >> channelnew;
		} while (channelnew=='y');	

		cout << "input threshold flag (y/n)";
		cin >> thresholdflag;
		if (thresholdflag=='y')
		{
			cout << "input threhold you need to set " << endl;
			cin >> threshold;
		}
	}	
	else if (graphoption=="c")
	{
		do
		{
			cout<<"input the channel you need to plot"<<endl;
			cin >> channelplotinput;
			channelplot.push_back(channelplotinput);
			cout << "Do you need input another channel to be plot (y/n)" << endl;
			cin >> channelnew;
		} while (channelnew=='y');	
		
		cout << "input threshold flag (y/n)";
		cin >> thresholdflag;
		if (thresholdflag=='y')
		{
			cout << "input threhold you need to set " << endl;
			cin >> threshold;
		}
	}
	else if (graphoption=="f")
	{
		int chanplot=atoi(argv[argc-3]);
		int thres=atoi(argv[argc-2]);
		cout << chanplot << "\t" << thres << endl;
		channelplot.push_back(chanplot);
		threshold=int(thres);
		thresholdflag='y';
	}
	/* variable to hold content for each line from data*/
	int lineMaximum=1000;
	char str[lineMaximum];
	const char *d=" '\t'";
	char* p;
	/*holding event information for the current line being analyzed*/
	int countline=0;
	int chan=0;
	int tempdatacount=0;
	int eventnumber=0;
	int eventrow=0;
	int eventchanl=0;
	int eventcol=0;
	int eventfirst=0; // variable for initialize the analysis for the first event	 
	int tempcondition[4]={0}; 	// {eventnumber rolnum colnum channelnum}
	
	vector<int> pulse;	
	vector<int> cuttedpulse;					
	vector<int> adjustedpulse;
	vector<double> pulseinfo;
	
	TString eeventchar, erowchar, echchar,ecolchar;
	for (int i=1; i<argc-1; i++)
	{
		/* reading data file, the variable input when run the program has the struture : file1 , file2, file3,.... , filelast, output-file-name*/
		fin.open (argv[i]);
		countline=0; // initialize countline=0; data file is read line by line, use this number to give reading process proper tag.
		if (!fin.is_open())
		{
			cerr << "Could not open "<< argv[1]<< endl;
			fin.clear();
		}
		cout << argv[i] << "\t file being analysised here" << endl; // Output the file being analysised at current run
		while (fin.getline(str,lineMaximum))
		{	
			/* line by line processing */
			if (countline<3)
			cout << "note output"<<endl;	// the first three line just file information and structure
			else
			{
				p=strtok (str,d);
				while (p)
				{
					/* word by word frome one line processing, the first 8 number is the inforamtion belong to basic hardware information*/
					if (tempdatacount==1)
					{
						tempcondition[0]=atoi(p);
						if (countline==3)
						{
							eventnumber=tempcondition[0];
							eventfirst=tempcondition[0];
							cout << "event number being initialized here with value \t" << eventfirst << endl;
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
						/* if event number is not the same, means new event starts here, at this point, analysis all the information analysised before for the previous event.*/
											
						if (tempcondition[0]!=eventnumber || tempcondition[1]!=eventrow || tempcondition[2]!=eventcol || tempcondition[3]!=eventchanl) 
						{
							channelcheck=false;
							eventcheck=false;						
							if (graphoption=="eg")
							{
								for (int j=0; j<eventplot.size(); j++)
								{
									if (eventnumber==eventplot[j])
									eventcheck=true;
								}
								chan=eventrow*100+eventcol*10+eventchanl;
								for (int j=0; j<channelplot.size(); j++)
								{
									if (chan==channelplot[j])
									channelcheck=true;	
								}
								drawoption=eventcheck&&channelcheck;	
								//cout << "drawoption being decided"<<endl;	
							}
							else if (graphoption=="c")
							{
								//cout << "drawoption being prepared"<<endl;
								chan=eventrow*100+eventcol*10+eventchanl;
								for (int j=0; j<channelplot.size(); j++)
								{
									if (chan==channelplot[j])
									channelcheck=true;
								}
								drawoption=channelcheck;	
								//cout << "drawoption being decided"<<endl;	
							}
							else if (graphoption=="f")
							{
								//cout << "drawoption being prepared"<<endl;
								chan=eventrow*100+eventcol*10+eventchanl;
								for (int j=0; j<channelplot.size(); j++)
								{
									if (chan==channelplot[j])
									channelcheck=true;
								}
								drawoption=channelcheck;	
								//cout << "drawoption being decided"<<endl;	
							}
							if (drawoption)	
							{						
								cuttedpulse=CutShock(pulse);					
								adjustedpulse=adjust(cuttedpulse);
								int *pulsey=new int[pulse.size()];
								int *xaxis=new int[pulse.size()];						
								for (int j=0; j<adjustedpulse.size(); j++)
								{
									pulsey[j]=adjustedpulse[j];
									xaxis[j]=j;
								}
								if (thresholdflag=='y')
								{
									pulseinfo=pulseProcess(pulse);
									if (pulseinfo[1]>threshold)
										thresholdcheck=true;
									else
										thresholdcheck=false;
								}
								else 
								thresholdcheck=true;
								echchar.Form("%d",eventchanl) ;
								ecolchar.Form("%d",eventcol) ;
								eeventchar.Form("%d",eventnumber);
								erowchar.Form("%d",eventrow) ;
								if (thresholdcheck)
								{
									TGraph *g=new TGraph(pulse.size(),xaxis,pulsey);
									g->SetTitle (eeventchar+""+erowchar+" "+ecolchar+" "+echchar+" A");								
									g->Write();
									//cout << "Graph input here" << endl;
									delete g; 								
								} 
								delete []pulsey;
								delete []xaxis;
							}
							pulse.clear();
							/*new event condition being initialized here*/
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
		fin.clear();
		fin.close();
	}
	
	f->Write();
	f->Close();
	return 0;
}

