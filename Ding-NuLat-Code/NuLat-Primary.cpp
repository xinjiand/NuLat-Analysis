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
	int totalenergysumMuon=0;
	/* Tstring name and 64 histogram being created, right now the histogram is PSD,Energy by both peak and Integral, might get rid off peak method*/
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
	TString title=" muon";
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
		peakhistMuon[i-1]=TH1D(peakstr+histname+title,titlename,2000,0,2000);	
		integralhistMuon[i-1]=TH1D(energystr+histname+title,titlename,1500,0,500000); // needs to be modified the size and maximum
		psdhistMuon[i-1]=TH2D(psdstr+histname+title,titlename,1500,0,500000,100,0,1);
		}
	/*vector to hold valid analysis information, being processed event by event and clear after reasonable histogram being created and data storage each event*/	
	vector<int> pulse;	
	vector<double> pulseinfoMuon;
	
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
											
						if (tempcondition[0]!=eventnumber || tempcondition[1]!=eventrow || tempcondition[2]!=eventcol || tempcondition[3]!=eventchanl) 
						{
							pulseinfoMuon=pulseProcess(pulse);
							histcount=eventchanl+eventcol*8+eventrow*32+1;
							if (histcount < 129)
							peakhistMuon[histcount-1].Fill(pulseinfoMuon[1]);
							integralhistMuon[histcount-1].Fill(pulseinfoMuon[0]);
							psdhistMuon[histcount-1].Fill(pulseinfoMuon[0],pulseinfoMuon[2]);
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
		fin.clear();
		fin.close();
	}
	f->Write();
	f->Close();
	return 0;
}
