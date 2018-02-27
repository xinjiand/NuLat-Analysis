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
#include <unordered_map>
#include <ctime>
using namespace std;
typedef vector<int> Row;
typedef vector<Row> Matrix;

int main(int argc, char* argv[])
{
	/*program running structure with n term
	0	 	1---n-2		n-1
	exe 	data file  	output-file-name
	*/
        clock_t start;
        double duration;
        start=clock();
	if (argc==1)
	{
		cerr << "Usage:" << argv[0] << "filename[s]\n";
		exit(EXIT_FAILURE);
	}
	ifstream fin;
	/*TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",argv[argc-1]);
	TFile* f=new TFile (analysisfile+rootap,"recreate"); //summary root file being created, every time code being runned, file being recreated.
	cout << "root file being created here" << endl;*/
	/* The output txt file will be listed as follows
	 * A event-> peakA.txt, energyA.txt, psdA.txt, timingA.txt
	 * B event-> peakB.txt, energyB.txt, psdB.txt, timingB.txt
	 * summary-> event.txt, summary.txt, cubeID.txt
	*/
	string analysisname=argv[argc-1];	
	ofstream peakfileMuon;
	peakfileMuon.open (analysisname+" peakMuon.txt",fstream::app);
	ofstream energyfileMuon;
	energyfileMuon.open (analysisname+" energyMuon.txt",fstream::app);
	ofstream psdfileMuon;
	psdfileMuon.open (analysisname+" psdMuon.txt",fstream::app);
	//ofstream timingfileMuon;
	//timingfileA.open (analysisname+" timingMuon.txt",fstream::app);
	ofstream cubeIDfile;
	cubeIDfile.open (analysisname+" cubeID.txt",fstream::app);	
	ofstream eventfile;
	eventfile.open (analysisname+" event.txt",fstream::app);	
	ofstream anasummary;
	anasummary.open (analysisname+" summary.txt",fstream::app);
        ofstream fcsv;
	fcsv.open (analysisname+".csv",fstream::app);
        fcsv << "Event Num, ChanID, Trig, energyA, energyB, peakA, peakB, timingA, timingB, timingAB, psdA, psdB, \n";
	
	/* variable to hold content for each line from data*/
	int lineMaximum=1000;
	char str[lineMaximum];
	const char *d=" '\t'";
	char* p;
	/*holding event information for the current line being analyzed*/
	int countline=0;
	int chan=0;
	int tempdatacount=0;
        int eventscrodnumber=0;
	int eventnumber=0;
        int eventscrod=0;
	int eventrow=0;
	int eventchanl=0;
	int eventcol=0;
	int eventwindow=0;
	int eventfirst=0;               // variable for initialize the analysis for the first event	 
	int tempcondition[6]={0}; 	// {eventnumber rolnum colnum channelnum timewindow scrod}
	int totalenergysumA=0;
	int totalenergysumB=0;
	int timewindowend=0;
	int timewindowstart=0;
	int timeinterval=0;
	int timewindowsize=0;
	int histcount=0;
        bool trigstatus=false;

	/*vector to hold valid analysis information, being processed event by event and clear after reasonable histogram being created and data storage each event*/	
	vector<int> pulse;
	vector<int> timewindow;	
	//vector<int> timeAB;
	vector<double> pulseinfoMuon;
	vector<int> event;
	vector<double> energyspecMuon;
	vector<double> energypeakMuon;
	vector<double> psdanalysisMuon;
	//vector<double> timingA;
	//vector<double> timingB;
        vector<int> scrod;
	vector<int> row;
	vector<int> col;
	vector<int> channel;
	vector<int> cubeID;
        vector<int> TriggeredChan;
        bool TrigChanCheck=false;

	unordered_map<int,int> channelmapXY;
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
                        channelmapXY[atoi(p)]=a*10+b;
                        b++;			
			p=strtok(NULL,d);	
		}
		//cout<< a << b << "\t" << map[a][b] << endl;
                b=0;
		a++;
                
	}
        printMatrix(map);
	
	bool pulsecheck=false; // bool variable used to check whether signal is from valid channel on the map
        std::unordered_map<int,int>::const_iterator got;
	int xID=0;
	int yID=0;
	
	cout << "Matrix initialize here" << endl;
	/*Matrix initialization 
		EnergyPeakMatrix	EnergyInteMatrix	CubeMatrix	PsdMatrix       TimingMatrix    TimingABMatrix
	*/
	int TimingMatrixMuon[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};
	
	int EnergyPeakMatrixMuon[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};
	
	int EnergyInteMatrixMuon[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};

	double PsdMatrixMuon[10][10] = {	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};
	


	int CubeMatrix[10][10] = {	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
								{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};

	int timeABMatrix[10][10] = {	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
									{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};

	
	/* reading data file, the variable input when run the program has the struture : file1 , file2, file3,.... , filelast, output-file-name*/
	for (int i=1; i<argc-1; i++)
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
						tempcondition[5]=atoi(p);
						if (countline==3)
						{
							eventscrodnumber=tempcondition[5];
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

				        else if (tempdatacount==7)
				        {
				   		tempcondition[4]=atoi(p);
                                                timewindow.push_back(tempcondition[4]);	
						if (countline==3)
						{
							eventwindow=tempcondition[4];
							cout << "event window being initialized here with value \t" << eventwindow << endl;
						}
	        		        }
                                        
                                  
					else if (tempdatacount==8)
					{
                                                /*deal with trig information first*/
                                                if (atoi(p)==1)
                                                {
                                                        TrigChanCheck=true;                                                        
                                                        chan=tempcondition[5]*1000+tempcondition[1]*100+tempcondition[2]*10+tempcondition[3];
                                                        if (TriggeredChan.size()!=0)
                                                        {
                                                                for (int j=0; j<TriggeredChan.size(); j++)
                                                                {
                                                                     if (TriggeredChan[j]==chan)
                                                                     {   
                                                                        TrigChanCheck=false;
                                                                        break;
                                                                     }
                                                                }
                                                        }
                                                        if (TrigChanCheck)
                                                             TriggeredChan.push_back(chan);   
                                                }
                                        
                                                						
                                                /*After all the pre-process of the data, start deal with the issue from last event*/						
						/* if event number is not the same, means new event starts here, at this point, analysis all the information analysised before for the previous event.*/				
						if (tempcondition[0]!=eventnumber || tempcondition[1]!=eventrow || tempcondition[2]!=eventcol || tempcondition[3]!=eventchanl || tempcondition[5]!=eventscrod) 
						{
							chan=eventscrod*1000+eventrow*100+eventcol*10+eventchanl;
							pulsecheck=false;
							/*loop being calledl to check whether the channel is match with map, pick up all the signal which is matched the map and only doing analysis for them*/
                                                        got=channelmapXY.find (chan);						
                                                        /*for (int j=0; j<10; j++)
							{
								for (int k=0; k<10; k++)
								{
									if (chan==map[j][k])
									{
										pulsecheck=true;
                                                                                break;
									}
								}
							}*/							
							if (got != channelmapXY.end())
                                                        //if (pulsecheck)
							{	
								//cout << "pulse check function running" << endl;						
								//pulseA=pulseAB(pulse,0);
								//pulseB=pulseAB(pulse,1);
								pulseinfoMuon=pulseProcess(pulse);
								histcount=eventchanl+eventcol*8+eventrow*32+1;
								cubeID.push_back(histcount);
								/*pulseinfo has the structure that
								 * 	pulseinfo.[0](totalenergy);
									pulseinfo.[1](peakamp);
									pulseinfo.[2](psdratio);
									pulseinfo.[3](peakpos);	 // for the peakpos as timing, we need to redo this use CFD */
                                                                /**/
								//timeAB.push_back(timeinterval*WINDOW-pulseinfoA[3]*SAMPLE-(pulse.size()-pulseinfoB[3])*SAMPLE);
                                                              //  cout << eventscrod << eventrow << eventcol << eventchanl << endl;
                                                              //  cout << histcount << "\t" << pulseinfoB[0] << endl;
                                                                //cout << timeinterval*WINDOW-pulseinfoA[3]*SAMPLE-(pulse.size()-pulseinfoB[3])*SAMPLE << "\t" << pulseinfoA[3]*SAMPLE << "\t" << (pulse.size()-pulseinfoB[3])*SAMPLE << endl; 
                                                                
								//cout << "Function after timeAB push" << endl;
								energyspecMuon.push_back(pulseinfoMuon[0]);
								energypeakMuon.push_back(pulseinfoMuon[1]);
								psdanalysisMuon.push_back(pulseinfoMuon[2]);
								//timingA.push_back(pulseinfoA[3]);
								event.push_back(eventnumber);
                                                                scrod.push_back(eventscrod);
								row.push_back(eventrow);
								col.push_back(eventcol);
								channel.push_back(eventchanl);
                                                                eventwindow=timewindow[timewindowsize-1];
                                                                timewindow.clear();
                                                                timewindow.push_back(eventwindow);
                                                                //cout <<eventrow << eventcol << eventchanl << "pulse being processed here" << endl;	
								//cout << "pulse check finished" << endl;
							}
							if (tempcondition[0]!=eventnumber)
							{				
								//cout << "event=" << eventnumber << " being processed" << endl;
								for (int j=0; j<cubeID.size(); j++) // create a proper mapping for different channel to a 2D Matrix
								{								
									/*call the function to find the position certain channel should be located*/
									//xID= mapXID (map , scrod[j], row[j], col[j], channel[j]);
									//yID= mapYID (map , scrod[j], row[j], col[j], channel[j]);
                                                                        chan=scrod[j]*1000+row[j]*100+col[j]*10+channel[j];
                                                                        xID=channelmapXY[chan]%10;
                                                                        yID=(channelmapXY[chan]/10)%10;
                                                                        //cout << eventnumber << "\t" << xID << "\t" << yID<< "\t" << scrod[j] << row[j]<<col[j]<<channel[j]<< "\t"<< energyspecA[j] << "\t" <<cubeID[j] << endl;	
									// write certain information into proper matrix
									EnergyPeakMatrixMuon[yID][xID]=energypeakMuon[j];		
									EnergyInteMatrixMuon[yID][xID]=energyspecMuon[j];
									PsdMatrixMuon[yID][xID] = psdanalysisMuon[j];
									//TimingMatrixA[yID][xID]=timingA[j];	
									CubeMatrix[yID][xID]=cubeID[j]; // might just get rid of this matrix and the txt file and just use the map
									//timeABMatrix[yID][xID]=timeAB[j];
                                                                        /*fcsv << "Event Num, ChanID, trigstatus, energyA, energyB, peakA, peakB, timingA, timingB, timingAB, psdA, psdB, \n";
	*/                                                              
									trigstatus=false;                                                                        
									for (int k=0; k<TriggeredChan.size(); k++)
									{
										 if (TriggeredChan[k]==chan)
										 {   
										    trigstatus=true;
										    break;
										 }
									}
									//fcsv << event[0] <<" ," << chan <<" ,"<< trigstatus << " ," << energyspecA[j] <<" ," << energyspecB[j]<<" ," << energypeakA[j]<<" ," << energypeakB[j]<<" ," << timingA[j]<<" ," << timingB[j]<<" ," << timeAB[j]<<" ," << psdanalysisA[j]<<" ," << psdanalysisB[j] <<" ," << "\n";

            								if (xID<5 && yID<5)									
									{
										totalenergysumA+=energyspecMuon[j]; // totalenergy currently add whole blue face which is from 04-0							
									}
								
								}
								//cout << "information being stored into Matrix" << endl;
								/*put sum into matrix [8][9], leave [9][9] for event number*/
								//EnergyInteMatrixA[8][9]=totalenergysumA;
								//EnergyInteMatrixB[8][9]=totalenergysumB;
								//CubeMatrix[8][9]=128; 
								/*write result into a preanalysis root file so we can check briefly how the data looks like*/
								/*for(int j=0;j<10;j++)
								{
									for(int k=0; k<10; k++)								
									{
										if (CubeMatrix[j][k]>0 && CubeMatrix[j][k]<129)
										{
                                                                                        cout << CubeMatrix[j][k] << endl;			
											peakhistA[CubeMatrix[j][k]-1].Fill(EnergyPeakMatrixA[j][k]);
											integralhistA[CubeMatrix[j][k]-1].Fill(EnergyInteMatrixA[j][k]);
											psdhistA[CubeMatrix[j][k]-1].Fill(EnergyInteMatrixA[j][k],PsdMatrixA[j][k]);
											peakhistB[CubeMatrix[j][k]-1].Fill(EnergyPeakMatrixB[j][k]);
											integralhistB[CubeMatrix[j][k]-1].Fill(EnergyInteMatrixB[j][k]);
											psdhistB[CubeMatrix[j][k]-1].Fill(EnergyInteMatrixB[j][k],PsdMatrixB[j][k]);
										}
																			
									}
								}*/
								/*writing the pre-analysis result to txt file for future analysis*/
								//cout << "write info into txt file" << endl;
								for(int j=0;j<10;j++)
								{
								
									for(int k=0; k<10; k++)
									{			
										//cout << k ;					
										if ( j==9 && k==9)
										{
											/*event number being record every matrix at [9][9], currently it is easier for people to check*/
											peakfileMuon << event[0]  << "\t";
											energyfileMuon << event[0] << "\t";
											psdfileMuon << event[0] << "\t";	
											//timingfileA << event[0] << "\t";
											cubeIDfile << event[0] << "\t" ;
									//		timeABfile << event[0] << "\t" ;
											eventfile << event[0] << endl;	
										}
										else
										{
											//cout << j << k << "\t" << CubeMatrix[j][k] <<endl;	
											peakfileMuon << EnergyPeakMatrixMuon[j][k] << "\t" ;
											energyfileMuon << 	EnergyInteMatrixMuon[j][k] << "\t";
											psdfileMuon << PsdMatrixMuon[j][k] << "\t";
											//timingfileA << TimingMatrixA[j][k] << "\t";
											cubeIDfile << CubeMatrix[j][k] << "\t";
											//timeABfile << timeABMatrix[j][k] << "\t";
										}	
																			
									}
									energyfileMuon << endl;
									peakfileMuon << endl;
									psdfileMuon << endl;	
									//timingfileA<<endl;
									cubeIDfile << endl;		
									//timeABfile << endl;								
								}
                                                                anasummary << event[0] << "\t" ;
                                                                for (int j=0; j<TriggeredChan.size(); j++)
                                                                {
                                                                        anasummary << TriggeredChan[j] << "\t";
                                                                }
                                                                anasummary << endl;
								/*clear all the information from last event*/
								for (int j=0; j<10; j++)
								{
									for (int k=0; k<10; k++)
									{									
										EnergyPeakMatrixMuon[j][k]=0;
										EnergyInteMatrixMuon[j][k]=0;
										CubeMatrix[j][k]=0;
										PsdMatrixMuon[j][k]=0; 										
										//TimingMatrixA[j][k]=0;		
									}
								}
								totalenergysumA=0;
								totalenergysumB=0;
								//timeAB.clear();		
								cubeID.clear();
								energyspecMuon.clear();
								energypeakMuon.clear();
								psdanalysisMuon.clear();
								//timingA.clear();
								event.clear();
                                                                scrod.clear();
								row.clear();
								col.clear();
								channel.clear();
                                                                TriggeredChan.clear();
								//cout << "event process finished" << endl;
							} // end of last event processing											
							pulse.clear();
							/*new event condition being initialized here*/
							//cout << "pulse information being cleared"<< endl;
							eventnumber=tempcondition[0];
							eventrow=tempcondition[1];
							eventcol=tempcondition[2];
							eventchanl=tempcondition[3];
                                                        eventscrod=tempcondition[5];	
							//cout << "end of information processing" << endl;
						}
					}
					else if (tempdatacount>8)
					{
						pulse.push_back(atoi(p));
						//cout << "test" << "\t" ;
					}
					p=strtok(NULL,d);
					tempdatacount++;
				}
				tempdatacount=0;
			}
			countline++;
		}
		anasummary << "analysis file is=" << argv[i] << "\t first event is=" << eventfirst << "\t last event is=" << eventnumber << endl;
		fin.clear();
		fin.close();
	}
	duration=(clock()-start)/(double) CLOCKS_PER_SEC;
	cout << "program running time is " << duration << endl;                
	anasummary.clear();	
	anasummary.close();
	cubeIDfile.clear();
	cubeIDfile.close();
	eventfile.clear();
	eventfile.close();
	//timingfileA.clear();
	//timingfileA.close();
	peakfileMuon.clear();
	peakfileMuon.close();
	energyfileMuon.clear();
	energyfileMuon.close();
	psdfileMuon.clear();
	psdfileMuon.close();
        fcsv.clear();
        fcsv.close();
	/*f->Write();
	f->Close();*/
	return 0;
}
