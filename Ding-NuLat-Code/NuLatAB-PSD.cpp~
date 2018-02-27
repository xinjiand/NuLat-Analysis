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
#include <ctime>
using namespace std;
typedef vector<int> Row;
typedef vector<Row> Matrix;

bool eventtest (int eventid, vector<int> &event);


int main(int argc, char* argv[])
{
        clock_t start;
        double duration;
        start=clock();
	if (argc==1)
	{
		cerr << "Usage:" << argv[0] << "filename[s]\n";
		exit(EXIT_FAILURE);
	}
	ifstream fin;
	TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",argv[argc-1]);
	TFile* f=new TFile (analysisfile+"-PSD"+rootap,"recreate"); //summary root file being created, every time code being runned, file being recreated.
	cout << "root file being created here" << endl;
        ofstream fout;
        ofstream foutdata;
        string filename=argv[argc-1];
	foutdata.open (filename+"-PSD.csv");
	foutdata << "Event Num, push, energyA, psdA, energyB, psdB, \n";
	
	/* variable to hold content for each line from data*/
	int lineMaximum=1000;
	char str[lineMaximum];
	const char *d=" '\t'";
	char* p;
	ifstream finevent;
	finevent.open (argv[argc-2]);
	vector<int> eventIDpick;
	while (finevent.getline(str,lineMaximum))
	{
		p=strtok (str,d);		
		eventIDpick.push_back (atoi(p));
	}
        /*input channel needs to do PSD search test*/
	int chantest=atoi(argv[argc-3]);
        cout << chantest << " being tested for PSD" << endl;
	/*holding event information for the current line being analyzed*/
	int countline=0;
	int chan=0;
	int tempdatacount=0;
        int eventscrodnumber=0;
        int eventscrod=0;
	int eventnumber=0;
	int eventrow=0;
	int eventchanl=0;
	int eventcol=0;
	int eventfirst=0; // variable for initialize the analysis for the first event	 
	int tempcondition[5]={0}; 	// {eventnumber rolnum colnum channelnum}
	/* Tstring name and 64 histogram being created, right now the histogram is PSD,Energy by both peak and Integral, might get rid off peak method*/
	TString eventIDA="eventA ";
	TString eventIDB="eventB ";
	TString eventID="eventA+eventB ";
	TH2D psdhistA[200];
	TH2D psdhistB[200];
	TH2D psdhist[200];
	for (int i=0; i<100; i++)
	{
		TString psdname ;
		psdname.Form ("%d",i);
		psdhistA[i]=TH2D(eventIDA+psdname,eventIDA+psdname,260,0,70000,50,-0.3,1.2);
		psdhistB[i]=TH2D(eventIDB+psdname,eventIDB+psdname,260,0,70000,50,-0.3,1.2);
		psdhist[i]=TH2D(eventID+psdname,eventID+psdname,260,0,70000,50,-0.3,1.2);
	}
	/*vector to hold valid analysis information, being processed event by event and clear after reasonable histogram being created and data storage each event*/	
	vector<int> pulse;
	vector<int> pulseA;
	vector<int> pulseB;
	vector<double> pulseinfoA;
	vector<double> pulseinfoB;
	vector<int> cuttedpulseA;					
	vector<int> adjustedpulseA;
	vector<int> cuttedpulseB;					
	vector<int> adjustedpulseB;
	bool eventpick=true;
        double psdratioA, psdratioB;
        double tailenergyA, tailenergyB;
        double totalenergyA, totalenergyB;
        int leftzeroposA, rightzeroposA, leftzeroposB, rightzeroposB, peakposA, peakposB;   
        int threshold=10;
        int psdposA, psdposB;

								
	/* reading data file, the variable input when run the program has the struture : file1 , file2, file3,.... , filelast, output-file-name*/
	for (int i=1; i<argc-3; i++)
	{
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
                                        else if (tempdatacount==2)
					{
						tempcondition[4]=atoi(p);
						if (countline==3)
						{
							eventscrodnumber=tempcondition[4];
							cout << "event row being initialized here with value \t" << eventrow << endl;
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
                                                chan=tempcondition[4]*1000+tempcondition[1]*100+tempcondition[2]*10+tempcondition[3]; 
                                                              					
						if (tempcondition[0]!=eventnumber || tempcondition[1]!=eventrow || tempcondition[2]!=eventcol || tempcondition[3]!=eventchanl || tempcondition[4]!=eventscrod) 
						{							
                                                        chan=eventscrod*1000+eventrow*100+eventcol*10+eventchanl;							
                                                        if (chan==chantest)
							{							
								//if (eventtest(eventnumber,eventIDpick))
								//{						
									pulseA=pulseAB(pulse,0);
									pulseB=pulseAB(pulse,1);
									cuttedpulseA=CutShock(pulseA);					
									adjustedpulseA=adjust(cuttedpulseA);
									cuttedpulseB=CutShock(pulseB);					
									adjustedpulseB=adjust(cuttedpulseB);									
									int *pulseyA=new int[adjustedpulseA.size()];
									int *xaxisA=new int[adjustedpulseA.size()];
									for (int j=0; j<adjustedpulseA.size(); j++)
									{
										pulseyA[j]=adjustedpulseA[j];
										xaxisA[j]=j;
									}
									int *pulseyB=new int[adjustedpulseB.size()];
									int *xaxisB=new int[adjustedpulseB.size()];
									for (int j=0; j<adjustedpulseB.size(); j++)
									{
										pulseyB[j]=adjustedpulseB[j];
										xaxisB[j]=j;
									}
									peakposA=maxfind(adjustedpulseA);
                                                                        peakposB=maxfind(adjustedpulseB);
                                                                        for (int j=0; j<peakposA;j++)
	                                                                {
		                                                                if (pulseyA[peakposA-j]<threshold)
		                                                                {
			                                                                leftzeroposA=peakposA-j;
			                                                                break;
		                                                                }		
	                                                                }
                                                                        for (int j=peakposA; j<adjustedpulseA.size();j++)
	                                                                {
		                                                                if (pulseyA[j]<threshold)
		                                                                {
			                                                                rightzeroposA=j;
			                                                                break;
		                                                                }		
	                                                                }
                                                                        for (int j=0; j<peakposB;j++)
	                                                                {
		                                                                if (pulseyB[peakposB-j]<threshold)
		                                                                {
			                                                                leftzeroposB=peakposB-j;
			                                                                break;
		                                                                }		
	                                                                }
                                                                        for (int j=peakposB; j<adjustedpulseB.size();j++)
	                                                                {
		                                                                if (pulseyB[j]<threshold)
		                                                                {
			                                                                rightzeroposB=j;
			                                                                break;
		                                                                }		
	                                                                }
                                                                        totalenergyA=sum(adjustedpulseA,leftzeroposA,adjustedpulseA.size());
                                                                        totalenergyB=sum(adjustedpulseB,leftzeroposB,adjustedpulseB.size());

									for (int j=0; j<100; j++)
									{
										psdposA=peakposA+10+j*5;
                                                                                psdposB=peakposB+10+j*5;
                                                                                if (psdposA>adjustedpulseA.size())
                                                                                psdposA=adjustedpulseA.size()-2;  
                                                                                if (psdposB>adjustedpulseB.size())
                                                                                psdposB=adjustedpulseB.size()-2;                           
                                                                                tailenergyA=sum(adjustedpulseA,psdposA,adjustedpulseA.size());
										tailenergyB=sum(adjustedpulseB,psdposB,adjustedpulseB.size());
                                                                                psdratioA=tailenergyA/totalenergyA;
                                                                                psdratioB=tailenergyB/totalenergyB;
                                                                                //foutdata << "Event Num, push, energyA, psdA, energyB, psdB, \n";
                                                                                foutdata << eventnumber << " ,"<< j << " ,"<< totalenergyA << " ,"<<psdratioA << " ,"<<totalenergyB << " ," << psdratioB << endl;  
                                                                                psdhistA[j].Fill(totalenergyA,psdratioA);
										psdhistB[j].Fill(totalenergyB,psdratioB);
										psdhist[j].Fill(totalenergyA,psdratioA);
										psdhist[j].Fill(totalenergyB,psdratioB);
									}
                                                                        //cout << eventnumber << eventrow << eventcol << eventchanl << eventscrod << " being processed" << endl;
                                                                        delete []pulseyA;
									delete []xaxisA;
									delete []pulseyB;
                                                                        delete []xaxisB;                                      
                                                                        pulseA.clear();
									pulseB.clear();
								//}
							}
							pulse.clear();
							/*new event condition being initialized here*/
							eventnumber=tempcondition[0];
							eventrow=tempcondition[1];
							eventcol=tempcondition[2];
							eventchanl=tempcondition[3];
                                                        eventscrod=tempcondition[4];	
							
						}
					}
					else if (tempdatacount>8)
                                        {                                                                      
                                                if (chan==chantest)					        
                                                pulse.push_back(atoi(p));
                                        }
					p=strtok(NULL,d);
					tempdatacount++;
				}
				tempdatacount=0;
			}
			countline++;
		}
                duration=(clock()-start)/(double) CLOCKS_PER_SEC;
                cout << "program running time is " << duration << endl;
		fin.clear();
		fin.close();
	}
        duration=(clock()-start)/(double) CLOCKS_PER_SEC;
        cout << "program running time is " << duration << endl;
        foutdata.clear();
        foutdata.close();           
	f->Write();
	f->Close();
	return 0;
}


bool eventtest (int eventid, vector<int> &event)
{
	bool eventpick=false;	
        int low=0, high=event.size()-1, middle;
        while (low <= high)
        {
                middle=(high+low)/2;
                if (event[middle]==eventid)
                {
                        eventpick=true;
                        break;
                }
                else if (event[middle]>eventid)
                {
                        high=middle-1;
                }
                else
                {
                        low=middle+1;                
                }
        } 
       
        /*
	for(int i=0; i<event.size(); i++)
	{
		if (eventid==event[i])
		eventpick=true;
	}
        */

	return eventpick;
}



