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

void recalibration (float a[10][10], int rl,int rh,int cl, int ch);

bool pick(int a[10][10],int b[10][10], Matrix &l , int n);

void print2Darray (int a[10][10]);

void printMatrix (Matrix &l);

void printfloat2Darray ( float a[10][10]);

bool compare(int energycompare, Matrix &l, int row, int col,int n);

bool compareTwo(int energycompare, Matrix &l, int row1, int col1, int row2, int col2, int n);

bool compete (Matrix &l , int r, int c, int rowl, int coll,int n);




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
	string option;
	cout << " Choose which command you need : Test (t) ; Cube (c) ; Energy Mapping (m)" << endl;
	cin >> option;
	/* open different file ready to be input
		
		energy	psd	peak	timing	cubeID	event
		
	*/
	string analysisname=argv[1];
	ifstream fincondition;
	fincondition.open ("condition.txt");
	ifstream finenergy;
	finenergy.open (analysisname+" energyMuon.txt");
	ifstream finpsd;
	finpsd.open (analysisname+" psdMuon.txt");
	ifstream finpeak;
	finpeak.open (analysisname+" peakMuon.txt");
	ifstream fintiming;
	fintiming.open (analysisname+" timingMuon.txt");
	ifstream fincube;
	fincube.open (analysisname+" cubeID.txt");
	ifstream finevent;
	finevent.open (analysisname+" event.txt");
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
	/* output a matrix mapping being reading here*/
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
	/*create different matrix to store data in the future and also 2D matrix as for calibration and threshold cut*/
	Matrix energyMuon(n*10,Row(10));	
	MatrixPsd psdMuon(n*10,RowPsd(10));
	Matrix peakMuon(n*10,Row(10));
	Matrix timingMuon(n*10,Row(10));
	Matrix cube(n*10,Row(10));
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
	i=0;
	while (finenergy.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{	
			energyMuon[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpeak.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{			
			peakMuon[i][j]=atoi(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (finpsd.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			psdMuon[i][j]=atof(p);
			j++;			
			p=strtok(NULL,d);	
		}
		j=0;
		i++;
	}
	i=0;
	while (fintiming.getline(str,lineMaximum))
	{
		p=strtok (str,d);
		while (p)
		{
			timingMuon[i][j]=atoi(p);
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

	/*test, just veto suspective muon event*/
	if (option=="test, just veto muon event" || option=="t")
	{
		cout << "function test pick" << endl;		
		/*create root file to hold spectrum*/
		TString rootap=".root";
		TString analysisfile;
		analysisfile.Form ("%s",argv[1]);
		TFile* f=new TFile (analysisfile+"-vetomuon-pickout"+rootap,"recreate");
		ofstream foutevent;
		foutevent.open ("cube-Event.txt");
		int eventcount=0;
		/*create category of histogram to hold spectrum*/
		TString rowstr="row ";
		TString colstr="col ";
		TString chanlstr="channel ";
		TString energystr="energy ";
		TString peakstr="peak ";
		TString psdstr="psd ";
		TString eventstr="event";
		TH1D peakhistMuon[128];
		TH1D integralhistMuon[128];
		TH2D psdhistMuon[128];
		TString titleMuon=" Muon";
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
			peakhistMuon[i-1]=TH1D(peakstr+histname+titleMuon,titlename,2000,0,2000);
			integralhistMuon[i-1]=TH1D(energystr+histname+titleMuon,titlename,1500,0,350000); // needs to be modified the size and maximum
			psdhistMuon[i-1]=TH2D(psdstr+histname+titleMuon,titlename,1500,0,350000,100,0,1);
		}
		bool vetoMuonpick=true;
		bool noisevo=true;
		int trigch;		
		cout << "input which channel work as trigger" << endl;
		cin >> trigch;
		int xtrig=0;
		int ytrig=0;
		xtrig=mapXIDtrig(map,trigch);
		ytrig=mapYIDtrig(map,trigch);
		string MapOption;
		cout << "Do you need 2D events mapping? (y/n)" << endl;
		cin >> MapOption;

		for (int i=0; i<n; i++)
		{
			noisevo=noiseveto (peakMuon , xtrig, ytrig, i);
			if (!noisevo) 
			{
				foutevent << eventID[i] << endl;
				eventcount++;					
				for(int j=0;j<10;j++)
				{
					for(int k=0; k<10; k++)
					{
						if (j<5||k<5)
						{
							if (cube[j][k]>0)
							{								
								peakhistMuon[cube[j][k]-1].Fill(peakMuon[j+10*i][k]);
								integralhistMuon[cube[j][k]-1].Fill(energyMuon[j+10*i][k]);
								psdhistMuon[cube[j][k]-1].Fill(energyMuon[j+10*i][k],psdMuon[j+10*i][k]);
								//cout << "fill one channel" << endl;
							}
						}
					}
				}
				if (MapOption=="y")
				{			
					TString eventchar;						
					eventchar.Form ("%d",eventID[i]);	
					TH2D* temp2dhisMuon=new TH2D("eventMuon "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5 || k<5)
							{
								for (int l=0;l<energyMuon[j+i*10][k];l++)
								{
									temp2dhisMuon->Fill(k,9-j);						
								}
							}
						}
					}
					cout << "wrtie the 2D map" << endl;
					temp2dhisMuon->Write();
					delete temp2dhisMuon;
				}
				
			}
		}
		cout << "total "<< eventcount << "\t events being found" << endl;
		foutevent.clear();
		foutevent.close();
		f->Write();
		f->Close();
	}
	



	/*Muon one cube pick*/
	if (option=="Muon one cube event pick" || option=="c")
	{
		cout << "function cube pick" << endl;		
		/*create root file to hold spectrum*/
		TString rootap=".root";
		TString analysisfile;
		analysisfile.Form ("%s",argv[1]);
		TFile* f=new TFile (analysisfile+"-cube-pickout"+rootap,"recreate");
		ofstream foutevent;
		foutevent.open ("cube-Event.txt");
		int eventcount=0;
		/*create category of histogram to hold spectrum*/
		TString rowstr="row ";
		TString colstr="col ";
		TString chanlstr="channel ";
		TString energystr="energy ";
		TString peakstr="peak ";
		TString psdstr="psd ";
		TString eventstr="event";
		TH1D peakhistMuon[128];
		TH1D integralhistMuon[128];
		TH2D psdhistMuon[128];
		TString titleMuon=" Muon";
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
			peakhistMuon[i-1]=TH1D(peakstr+histname+titleMuon,titlename,2000,0,2000);
			integralhistMuon[i-1]=TH1D(energystr+histname+titleMuon,titlename,1500,0,350000); // needs to be modified the size and maximum
			psdhistMuon[i-1]=TH2D(psdstr+histname+titleMuon,titlename,1500,0,350000,100,0,1);
		}
		/*
			X=0					Z=0
		Y=0	137	101	107	102	103	030	031	032	033	034
			104	105	106	120	121	035	036	037	020	021
			122	123	124	125	126	022	023	024	025	026
			127	110	111	112	113	027	010	011	012	013
			014	015	016	017	000	114	115	116	117	100
			201	231	202	233	234	999	999	999	999	999
			235	236	237	220	221	999	999	999	999	999
			222	223	224	225	226	999	999	999	999	999
			227	210	211	212	213	999	999	999	999	999
			214	215	216	217	200	999	230	231	004	999
		*/
		int xidMuon, yidMuon, zidMuon;
                double fracVeto=0.;
                cout << "input factor for veto" << endl;
		cin >> fracVeto;	
		string MapOption;
		int faceveto;
		cout << "input xid for Muon event" << endl;
		cin >> xidMuon;
		cout << "input yid for Muon event" << endl;
		cin >> yidMuon;
		cout << "input zid for Muon event" << endl;
		cin >> zidMuon;
		cout << "Do you need 2D events mapping? (y/n)" << endl;
		cin >> MapOption;
		cout << "input how much face needs to be veto 2/3" << endl;
		cin >> faceveto;
		bool vetoMuon=true;
		bool vetoMuonpick=true;
		bool noisevo=true;
		int trigch;		
		cout << "input which channel work as trigger" << endl;
		cin >> trigch;
		int xtrig=0;
		int ytrig=0;
		xtrig=mapXIDtrig(map,trigch);
		ytrig=mapYIDtrig(map,trigch);
		for (int i=0; i<n; i++)
		{
			vetoMuonpick=pick(lowthreshold,highthreshold,energyMuon,i);
			noisevo=noiseveto (peakMuon , xtrig, ytrig, i);
			if (true) //!vetoMuonpick && !noisevo
			{
				//cout << eventID[i] << "Muon 3 data"<< endl;
				if (faceveto==2)
	 			vetoMuon=twosideveto(energyMuon, xidMuon,yidMuon,zidMuon,fracVeto,i);
				else
				vetoMuon=cubeveto(energyMuon, xidMuon,yidMuon,zidMuon,fracVeto,i);
				if (!vetoMuon)
				{
					foutevent << eventID[i] << endl;
					eventcount++;					
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5||k<5)
							{
								if (cube[j][k]>0)
								{								
									peakhistMuon[cube[j][k]-1].Fill(peakMuon[j+10*i][k]);
									integralhistMuon[cube[j][k]-1].Fill(energyMuon[j+10*i][k]);
									psdhistMuon[cube[j][k]-1].Fill(energyMuon[j+10*i][k],psdMuon[j+10*i][k]);
									//cout << "fill one channel" << endl;
								}
							}
						}
					}
					if (MapOption=="y")
					{			
						TString eventchar;						
						eventchar.Form ("%d",eventID[i]);	
						TH2D* temp2dhisMuon=new TH2D("eventMuon "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
						for(int j=0;j<10;j++)
						{
							for(int k=0; k<10; k++)
							{
								if (j<5 || k<5)
								{
									for (int l=0;l<energyMuon[j+i*10][k];l++)
									{
										temp2dhisMuon->Fill(k,9-j);						
									}
								}
							}
						}
						cout << "wrtie the 2D map" << endl;
						temp2dhisMuon->Write();
						delete temp2dhisMuon;
					}
				}
			}
		}
		cout << "total "<< eventcount << "\t events being found" << endl;
		foutevent.clear();
		foutevent.close();
		f->Write();
		f->Close();
	}

	if (option=="Event Energy mapping " || option=="m")
	{
		cout << "function mapping pick" << endl;		
		/*create root file to hold spectrum*/
		TString rootap=".root";
		TString analysisfile;
		analysisfile.Form ("%s",argv[1]);
		TFile* f=new TFile (analysisfile+"-energy-map"+rootap,"recreate");
		/*ask for event to be plot*/
		vector<int> eventIDplot;
		int eventid=0;
		string eventNew;
		do
		{
			cout<<"input the event you need to map"<<endl;
			cin >> eventid;
			eventIDplot.push_back(eventid);
			cout << "Do you need input another channel to be plot (y/n)" << endl;
			cin >> eventNew;
		} while (eventNew=='y');

		bool eventcheck=false;
		for (int i=0; i<n; i++)
		{	
			eventcheck=false;			
			for (int j=0; j<eventIDplot.size(); j++)
			{
				if (eventID[i]==eventIDplot[j])
					eventcheck=true;
			}
			if (eventcheck)
			{
				
				TString eventchar;						
				eventchar.Form ("%d",eventID[i]);	
				TH2D* temp2dhisMuon1=new TH2D("eventA1 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
				TH2D* temp2dhisMuon2=new TH2D("eventA2 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
				TH2D* temp2dhisMuon3=new TH2D("eventA3 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
				for(int j=0;j<5;j++)
				{
					for(int k=0; k<5; k++)
					{
						{
							for (int l=0;l<energyMuon[j+i*10][k];l++)
							{
								temp2dhisMuon1->Fill(k,4-j);						
							}
						}
					}
				}
				for(int j=0;j<5;j++)
				{
					for(int k=5; k<10; k++)
					{
						{
							for (int l=0;l<energyMuon[j+i*10][k];l++)
							{
								temp2dhisMuon2->Fill(k-5,4-j);						
							}
						}
					}
				}
				for(int j=5;j<10;j++)
				{
					for(int k=0; k<5; k++)
					{
						{
							for (int l=0;l<energyMuon[j+i*10][k];l++)
							{
								temp2dhisMuon3->Fill(k,9-j);						
							}
						}
					}
				}				

				cout << "wrtie the 2D map" << endl;
				temp2dhisMuon1->Write();
				temp2dhisMuon2->Write();
				temp2dhisMuon3->Write();
				delete temp2dhisMuon1;
				delete temp2dhisMuon2;
				delete temp2dhisMuon3;
			}
		}
		f->Write();
		f->Close();
	}						

	cout << "good analysis " << n << " event being"<<endl;
	finenergy.clear();
	finpsd.clear();
	finpeak.clear();
	fintiming.clear();
	fincube.clear();
	finenergy.close();
	finpsd.close();
	finpeak.close();
	fintiming.close();
	fincube.close();
	return 0;
}



