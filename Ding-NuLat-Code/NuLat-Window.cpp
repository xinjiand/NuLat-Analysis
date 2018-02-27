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
	//cout << "input event-channel type or channel option---eg----c----f" << endl;
	//cin >> graphoption;
	int eventplotinput;
	int channelplotinput;
	cout<<"input the event you need to plot"<<endl;
	cin >> eventplotinput;
	cout<<"input the channel you need to plot"<<endl;
        cin >> channelplotinput;
					
	
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
	int tempcondition[5]={0}; 	// {eventnumber rolnum colnum channelnum}
	
	vector<int> pulse;	
	vector<int> pulseA;
	vector<int> pulseB;
	vector<int> cuttedpulseA;					
	vector<int> adjustedpulseA;
	vector<int> cuttedpulseB;					
	vector<int> adjustedpulseB;
	vector<double> pulseinfoA;
	vector<double> pulseinfoB;
	
	TString eeventchar, erowchar, echchar,ecolchar, ewindow;

	for (int i=1; i<argc-1; i++)
	{
		/* reading data file, the variable input when run the program has the struture : file1 , file2, file3,.... , filelast, output-file-name*/
		fin.open (argv[i]);
		countline=0; // initialize countline=0; data file is read line by line, use this number to give reading process proper tag.
		if (!fin.is_open())
		{
			cerr << "Could not open "<< argv[i]<< endl;
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
                                        else if (tempdatacount==7)
				        {
				   		tempcondition[4]=atoi(p);
						timewindow.push_back(tempcondition[4]);
                                                eventwindow=tempcondition[4];
						if (countline==3)
						{
							cout << "event window being initialized here with value \t" << eventwindow << endl;
						}
	        		        }
					else if (tempdatacount==8)
					{
						/*After all the pre-process of the data, start deal with the issue from last event*/						
						/* if event number is not the same, means new event starts here, at this point, analysis all the information analysised before for the previous event.*/
                                                chan=eventrow*100+eventcol*10+eventchanl;						
                                                if (eventnumber==eventplotinput && chan==channelplotinput)					
						{
                                                        while (p)
                                                        {
                                                        pulse.push_back(atoi(p));
                					p=strtok(NULL,d);
                					}

						      
						        int *pulsey=new int[pulse.size()];
						        int *xaxis=new int[pulse.size()];	
						       				
						        for (int j=0; j<pulse.size(); j++)
						        {
							        pulsey[j]=pulse[j];
							        xaxis[j]=j;
						        }
						        
                                                        eeventchar.Form("%d",eventnumber);                                                        
                                                        echchar.Form("%d",eventchanl) ;
							ecolchar.Form("%d",eventcol) ;
							erowchar.Form("%d",eventrow) ;
                                                        ewindow.Form("%d", eventwindow) ;
						        TGraph *gr=new TGraph(pulse.size(),xaxis,pulsey);
						        gr->SetTitle (eeventchar+""+erowchar+" "+ecolchar+" "+echchar+""+ewindow);		
						        gr->Write();
						        delete gr; 								    
						        delete []pulsey;
						        delete []xaxis;
						        pulse.clear();
		        			}
					}
					
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

