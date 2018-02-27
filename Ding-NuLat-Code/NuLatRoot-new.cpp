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
typedef vector<int> Row;
typedef vector<Row> Matrix;

vector<int> adjust (vector<int> &l);
int maxfind (vector<int>&l);
int minfind(vector<int> &l);
int sum (vector<int>&l, int number1, int number2);
double psd (vector<int>&l,int number1, int number2, int number3);
vector<int> CutShock (vector<int> &l);
int mapXID (Matrix &l , int row, int col,int chanl);
int mapYID (Matrix &l , int row, int col,int chanl);


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
	TFile* f=new TFile (analysisfile+rootap,"recreate"); //summary root file being created, every time code being runned, file being recreated.
	cout << "root file being created here" << endl;
	/*Current output file have 7 different txt file
	5 of them store different information by Matrix format, Peak, Energy, Psd, Timing, CubeID
	event.txt file is create here so it can be used to initial event tag during analysis
	summary .txt file needs to be thinking how much summary information will be useful
	*/
	ofstream peakfile;
	peakfile.open (analysisname+" peak.txt",fstream::app);
	ofstream energyfile;
	energyfile.open (analysisname+" energy.txt",fstream::app);
	ofstream psdfile;
	psdfile.open (analysisname+" psd.txt",fstream::app);
	ofstream timingfile;
	timingfile.open (analysisname+" timing.txt",fstream::app);	
	ofstream cubeIDfile;
	cubeIDfile.open (analysisname+" cubeID.txt",fstream::app);	
	ofstream eventfile;
	eventfile.open (analysisname+" event.txt",fstream::app);	
	ofstream anasummary;
	anasummary.open (analysisname+" summary.txt",fstream::app);
	
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
	int eventfirst=0;
	/* this flag is used to deal with the condtion when the last pulse of one event is just noise, then it will be throw away at the start of the new event, however the information of the event will not be organized, so the flag goes up, then the next good puse from other event will deal with this event although at that time, the eventnumber is matching at that time.*/
	int firsteventflag=0;		 
	int tempcondition[4]={0}; 	// {eventnumber rolnum colnum channelnum}
	/*pulse information being processed for each analysis*/	
	int totalenergy=0;
	int peakamp=0;
	int valleyamp=0;
	double psdratio=0.;
	int totalenergysum=0;
	/* noise graph number counnting , just to verify reasonable being throwed away, might get rid off*/
	int graphnumber=1;
	/* Tstring name and 64 histogram being created, right now the histogram is PSD,Energy by both peak and Integral, might get rid off peak method*/
	TString rowstr="row ";	
	TString colstr="col ";
	TString chanlstr="channel ";
	TString energystr="energy ";
	TString peakstr="peak ";
	TString psdstr="psd ";
	TString eventstr="event";
	TH1D peakhist[128];
	TH1D integralhist[128];
	TH2D psdhist[128];
	int histcount=0;
	for (int i=1; i<129; i++)
	{
				
		TString rowname,colname,chanlname;
		int rowcount=(i-1)/32;
		int colcount=0;		
		if (i<33)
		colcount=(i-1)/8;
		else if (i<65)
		colcount=(i-33)/8;
		else if (i<97)
		colcount=(i-65)/8;
		else if (i<129)
		colcount=(i-97)/8;
		int channelcount=(i-1)%8;
		rowname.Form ("%d",rowcount);
		colname.Form ("%d",colcount);
		chanlname.Form ("%d",channelcount);
		TString titlename = rowstr+" "+rowname+" "+colstr+" "+colname+" "+chanlstr+" "+chanlname;
		TString histname = rowname+" "+colname + " "+chanlname;
		peakhist[i-1]=TH1D(peakstr+histname,titlename,2000,0,2000);	
		integralhist[i-1]=TH1D(energystr+histname,titlename,1500,0,500000); // needs to be modified the size and maximum
		psdhist[i-1]=TH2D(psdstr+histname,titlename,1500,0,500000,100,0,1);
	}
	TH1D* totalenergyhist=new TH1D("totalenergyhist","total energy spec",1500,0,500000);
	/*vector to hold valid analysis information, being processed event by event and clear after reasonable histogram being created and data storage each event*/	
	vector<int> pulse;
	vector<int> event;
	vector<int> energyspec;
	vector<int> energypeak;
	vector<int> timing;
	vector<double> psdanalysis; //007-trig
	/*Geometry information being prcessed current version
				x=0	x=1	x=2	x=3	x=4	x=5	x=6	x=7	x=8	x=9			
		y=0		130	131	132	133	134	030	031	032	033	034
		y=1		135	137	124	125	126	035	037/106	020	021	022
		y=2		127	110	111	112	113	023	024	025	026	027			
		y=3		114	115	116	117	100	010	011	012	013	014
		y=4		101	102	103	104	105	015	016	017	106?	107?
		y=5		
		y=6		
		y=7										107	037	
		y=8		000	001	002	003	004	005	006	007	136	036	
		y=9		120	121	122	123	sum+				extra	eventID
		031 033 035 037 021 023 are the channel 1/10 attenuation 
		007 is the channel being trigged		
	*/
	
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
	
	vector<int> row;
	vector<int> col;
	vector<int> channel;
	vector<int> cubeID;
	int xID[64]={0};
	int yID[64]={0};
	int energyinteg[64]={0};
	cout << "Matrix initialize here" << endl;
	/*Matrix initialization 
		EnergyPeakMatrix	EnergyInteMatrix	CubeMatrix	PsdMatrix
	*/
	int TimingMatrix[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};
	
	int EnergyPeakMatrix[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};
	
	int EnergyInteMatrix[10][10]={	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
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

	double PsdMatrix[10][10] = {	{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
					{0,	0,	0,	0,	0,	0,	0,	0,	0,	0}};

	vector<int> adjustedpulse;
	vector<int> cuttedpulse;
	int peakpos=0;
	int valleypos=0;	
	int threshold=0;
	int leftzeropos=0;
	int rightzeropos=0;
	int peakfullwidth=0;
	int deltapeak=0;
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
											
						if (tempcondition[0]!=eventnumber) 
						{
							
							if (firsteventflag==1)
							firsteventflag=0;	 // firstevent flag initialize to be zero so the default treat it as good pulse	
							cuttedpulse=CutShock(pulse);					
							adjustedpulse=adjust(cuttedpulse);
							int *pulsey=new int[adjustedpulse.size()];
							int *xaxis=new int[adjustedpulse.size()];
							peakpos=maxfind(adjustedpulse);
							valleypos=minfind(adjustedpulse);							
							for (int j=0; j<adjustedpulse.size(); j++)
							{
								pulsey[j]=adjustedpulse[j];
								xaxis[j]=j;
							}
							peakamp=pulsey[peakpos];
							valleyamp=pulsey[valleypos];
							threshold=0;
							leftzeropos=0;
							rightzeropos=adjustedpulse.size(); // if pulse is really long, the end of the graph still not long enough for the pulse to reach zero, then the default end of the pulse will be the end of graph
							/*find the left size zero point of the pulse*/
							for (int j=0; j<peakpos;j++)
							{
								if (pulsey[peakpos-j]<threshold)
								{
									leftzeropos=peakpos-j;
									break;
								}		
							}
							/*find the right side zero point of the pulse*/
							for (int j=0; j< adjustedpulse.size()-peakpos; j++)
							{
								if (pulsey[peakpos+j]<threshold)
								{
									rightzeropos=peakpos+j;
									break;
								}
							}
							/*compute the peakfullwidth and deltapeak, use both to veto the noise*/
							peakfullwidth=rightzeropos-leftzeropos; 
							if (peakamp <20 || peakamp >1500)
							{
								firsteventflag=1;	
							}
							/*Good pulse being selectd and deal with here*/
							if (firsteventflag==0)
							{	
								/*analyze information for the new pulse*/
								histcount=eventchanl+eventcol*8+eventrow*32+1;
								if (histcount>64)
								{
									cout << "Unexpected signal come throught" << endl;
									//fout << "Unexpected signal come throught" << endl;
									eventnumber=tempcondition[0];
									eventrow=tempcondition[1];
									eventcol=tempcondition[2];
									eventchanl=tempcondition[3];
									break;	
								}																
								totalenergy=sum(adjustedpulse,leftzeropos,adjustedpulse.size());
								psdratio = psd (adjustedpulse,leftzeropos,peakpos+30,adjustedpulse.size());
								/*event  information storage*/
								cubeID.push_back(histcount);
								timing.push_back(peakpos);
								event.push_back(eventnumber);
								energyspec.push_back(totalenergy);
								energypeak.push_back(peakamp);
								psdanalysis.push_back(psdratio);
								row.push_back(eventrow);
								col.push_back(eventcol);
								channel.push_back(eventchanl);
							}	
							/* cubeID.size >0 means last event has useful data, store the useful data here*/					
							if (cubeID.size()>0)
							{																
								for (int j=0; j<cubeID.size(); j++) // create a proper mapping for different channel to a 2D Matrix
								{								
									xID[j]=	mapXID (map , row[j], col[j], channel[j]);
									yID[j]= mapYID (map , row[j], col[j], channel[j]);								
									energyinteg[j]=energyspec[j];
									EnergyPeakMatrix[yID[j]][xID[j]]=energypeak[j];
									EnergyInteMatrix[yID[j]][xID[j]]=energyspec[j];
									CubeMatrix[yID[j]][xID[j]]=cubeID[j];
									PsdMatrix[yID[j]][xID[j]] = psdanalysis[j];	
									TimingMatrix[yID[j]][xID[j]]=timing[j];		
									if (xID[j]<5 && yID[j]<5)									
									{
										totalenergysum+=energyinteg[j]; // totalenergy currently add whole blue face which is from 04-04
									}	
														
								}
								totalenergyhist->Fill(totalenergysum);				
								for(int j=0;j<10;j++)
								{
									for(int k=0; k<10; k++)
									{
																							
										if (j==9 && k==4)
										{
											/*store program sum at row 9, col 4 assign cubeID 1 for this for current purpose*/	
											EnergyInteMatrix[j][k]=totalenergysum;
											CubeMatrix[j][k]=1; // currently overlap with 000, assign another number when fully functional board
										}											
										else if (CubeMatrix[j][k]>0)
										{											
											peakhist[CubeMatrix[j][k]-1].Fill(EnergyPeakMatrix[j][k]);
											integralhist[CubeMatrix[j][k]-1].Fill(EnergyInteMatrix[j][k]);
											psdhist[CubeMatrix[j][k]-1].Fill(EnergyInteMatrix[j][k],PsdMatrix[j][k]);				
										}
																				
									}
								}
								/*writing the pre-analysis result to txt file for future analysis*/
								for(int j=0;j<10;j++)
								{
									for(int k=0; k<10; k++)
									{
											
										if ( j==9 && k==9)
										{
											peakfile << event[0]  << "\t";
											energyfile << event[0] << "\t";
											psdfile << event[0] << "\t";	
											timingfile << event[0] << "\t";
											cubeIDfile << event[0] << "\t" ;
											eventfile << event[0] << endl;									
										}
										else
										{
											peakfile << EnergyPeakMatrix[j][k] << "\t" ;
											energyfile << 	EnergyInteMatrix[j][k] << "\t";
											psdfile << PsdMatrix[j][k] << "\t";
											timingfile << TimingMatrix[j][k] << "\t";
											cubeIDfile << CubeMatrix[j][k] << "\t";
										}											
									}
									peakfile << endl;
									energyfile << endl;
									psdfile << endl;	
									timingfile<<endl;					
									cubeIDfile << endl;				
								}	
								/*clear all the information after event information being stored so they can be ready for next event*/	
								for (int j=0; j<64; j++)
								{
									xID[j]=0;
									yID[j]=0;
									energyinteg[j]=0;
								}		
								for (int j=0; j<10;j++)
								{
									for (int k=0; k<10; k++)
									{									
										EnergyPeakMatrix[j][k]=0;
										EnergyInteMatrix[j][k]=0;
										CubeMatrix[j][k]=0;
										PsdMatrix[j][k] = 0;	
										TimingMatrix[j][k]=0;		
									}
								}
								totalenergysum=0;		
								cubeID.clear();
								timing.clear();
								event.clear();
								energyspec.clear();
								energypeak.clear();
								psdanalysis.clear();
								row.clear();
								col.clear();
								channel.clear();										
							}
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
						/*After all the pre-process of the data, start deal with the issue from different channel but same event*/
						/* event number is the same, it means it's new channel of the same event, analysis the information and stored into vector*/
						else if (tempcondition[1]!=eventrow || tempcondition[2]!=eventcol || tempcondition[3]!=eventchanl)
						{								
							cuttedpulse=CutShock(pulse);					
							adjustedpulse=adjust(cuttedpulse);
							int *pulsey=new int[adjustedpulse.size()];
							int *xaxis=new int[adjustedpulse.size()];
							peakpos=maxfind(adjustedpulse);
							valleypos=minfind(adjustedpulse);							
							for (int j=0; j<adjustedpulse.size(); j++)
							{
								pulsey[j]=adjustedpulse[j];
								xaxis[j]=j;
							}
							peakamp=pulsey[peakpos];
							valleyamp=pulsey[valleypos];
							threshold=0;
							leftzeropos=0;
							rightzeropos=adjustedpulse.size(); // if pulse is really long, the end of the graph still not long enough for the pulse to reach zero, then the default end of the pulse will be the end of graph
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
							//fout << "peakleft=" << leftzeropos <<"\t peak=" << peakpos << "\t peakright=" << rightzeropos << endl;
							peakfullwidth=rightzeropos-leftzeropos;
							deltapeak=peakamp-valleyamp;
							if (peakamp < 20 || peakamp >1500)
							{		
								graphnumber++;	
							}
							else
							{
								histcount=eventchanl+eventcol*8+eventrow*32+1;
								if (histcount>64)
								{
									cout << "Unexpected signal come throught" << endl;
								//	fout << "Unexpected signal come throught" << endl;
									eventnumber=tempcondition[0];
									eventrow=tempcondition[1];
									eventcol=tempcondition[2];
									eventchanl=tempcondition[3];
									break;	
								}																
								totalenergy=sum(adjustedpulse,leftzeropos,adjustedpulse.size());
								psdratio = psd (adjustedpulse,leftzeropos,peakpos+30,adjustedpulse.size());
								/*event  information storage*/							
								cubeID.push_back(histcount);
								timing.push_back(peakpos);
								energyspec.push_back(totalenergy);
								energypeak.push_back(peakamp);
								psdanalysis.push_back(psdratio);
								event.push_back(eventnumber);
								row.push_back(eventrow);
								col.push_back(eventcol);
								channel.push_back(eventchanl);	
							}	
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
		anasummary << "analysis file is=" << argv[i] << "\t first event is=" << eventfirst << "\t last event is=" << eventnumber << endl;
		fin.clear();
		fin.close();
	}
	anasummary.clear();
	cubeIDfile.clear();
	timingfile.clear();
	peakfile.clear();
	energyfile.clear();
	psdfile.clear();
	eventfile.clear();
	eventfile.close();
	cubeIDfile.close();
	timingfile.close();
	peakfile.close();
	energyfile.close();
	psdfile.close();
	anasummary.close();
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

int mapXID (Matrix &l , int row, int col,int chanl)
{
	int xid=0;
	int chan=row*100+col*10+chanl;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (chan==l[i][j])
			xid=j;
		}
	}
	return xid;
}

int mapYID (Matrix &l , int row, int col,int chanl)
{
	int yid=0;
	int chan=row*100+col*10+chanl;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (chan==l[i][j])
			yid=i;
		}
	}
	return yid;
}



