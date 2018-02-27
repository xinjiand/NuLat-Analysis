#ifndef anafunction_H_
#define anafunction_H_
#include <stdlib.h>
#include <math.h>
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"

#define WINDOW 23.5
#define SAMPLE 0.37
using namespace std;
typedef vector<int> Row;
typedef vector<Row> Matrix;
typedef vector<double> RowPsd;
typedef vector<RowPsd> MatrixPsd;
/*function generall used when deal with the raw data */
vector<int> adjust (vector<int> &l);
vector<int> CutShock (vector<int> &l);
int maxfind (vector<int>&l);
int minfind(vector<int> &l);
int sum (vector<int>&l, int number1, int number2);
double psd (vector<int>&l,int number1, int number2, int number3, int number4);
int mapXID (Matrix &l , int row, int col,int chanl);
int mapYID (Matrix &l , int row, int col,int chanl);
vector<double> pulseProcess (vector<int>&l);
/*fuction being used directly for the AB trig event*/
vector<int> pulseAB (vector<int> &l, int a);
/*fuction being used for post anslysis*/
void recalibration (float a[10][10], int rl,int rh,int cl, int ch);
bool pick(int a[10][10],int b[10][10], Matrix &l , int n);
void print2Darray (int a[10][10]);
void printMatrix (Matrix &l);
void printfloat2Darray ( float a[10][10]);
int mapXID (Matrix &l , int channum);
int mapYID (Matrix &l , int channum);
bool pmtsideveto (Matrix &l , int a, int b, double c,int side, int n);
bool cubeveto (Matrix &l , int x, int y, int z, double c,int n);
vector<double> pulsePSDProcess (vector<int>&l,int push);
int timedif (int timebegin, int timeend);
bool twosideveto (Matrix &l , int x, int y, int z,double c, int n);
bool noiseveto (Matrix &l , int x, int y, int n);
int mapXIDtrig (Matrix &l , int ch);
int mapYIDtrig (Matrix &l , int ch);
vector<int> CFDpulse (vector<int> &l);
int BEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, Matrix &psdA, Matrix &psdB, Matrix &timingAB);
int ABEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, Matrix &psdA, Matrix &psdB, Matrix &timingAB);
int ABtEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, MatrixPsd &psdA, MatrixPsd &psdB, Matrix &timingAB);
int MapEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB );
int PulseCFD (vector<int>&l);
int OneEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, MatrixPsd &psdA, MatrixPsd &psdB, Matrix &timingAB);
int RadonCalibration (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, MatrixPsd &psdA, MatrixPsd &psdB, Matrix &timingAB);
int lowerCubeBound (int i);
int upperCubeBound (int i);
int lowerBinBound (int i, int j);
int upperBinBound (int i, int j);
/* function generall used for both NuLat-AB and NuLat-Muon and NuLat-PSD 
	adjust -> flip the negative pulse into positive so when do the integral, it directly becomes positive energy valuse
				It also readjust the base line by average the point before the pulse, treat that as baseline
	CutShock-> Function being called to cut all the shock signal which will mimic the pulse or confused peak finding mechanism
	maxfine ; minfind -> return the position where peak and valley of pulse located
	sum -> add point of pulse from certain range together, it used for different purpose like adjust baseline, find energy, and deal with psd	
	psd -> return the psd ratio value
	mapXID and mapYID -> return the xid and yid basid on the map of the experiment 
	pulseProcess-> deal with a pulse being cutshocked and adjusted, return a vector that contain meaningful informaiton about energy, peak-height, timing, psd
*/
vector<int> flip (vector<int> &l)
{
	vector<int> flip;
	for (int i=0; i<l.size(); i++)
	{
		flip.push_back(-l[i]); //adjust the baseline and flip the pulse
	}
	return flip;
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
	int averageOfbegin=sum(l,0,30)/30; // treat the average of the first 30 point as the brief base point, another algorithum used before is find the maximum point, but that might come from accident shock and very sharp 
	int threshold=5;
	int leftzeropos=50; // set the inital value of leftzeropos as 50
	int adjustoff=0;
        int rightzeropos=l.size()-1;
	for (int i=0; i<minpos;i++)
	{
		if (array[minpos-i]>(averageOfbegin-threshold))
		{
			leftzeropos=minpos-i;
			break;
		}		
	}
        for (int i=minpos; i<l.size();i++)
	{
		if (array[minpos+i]>(averageOfbegin-threshold))
		{
			rightzeropos=minpos-i;
			break;
		}		
	}
	int totalcount=leftzeropos+l.size()-rightzeropos;
	if (leftzeropos<5)
		leftzeropos=5;// to avoid condition that a purly noise pulse make the first point as the pulse peak
	adjustoff=sum(l,0,leftzeropos)+sum(l,rightzeropos, l.size());
	int adjustoffset=int (adjustoff/totalcount);
	for (int i=0; i<l.size(); i++)
	{
		adjust.push_back(-(l[i]-adjustoffset)); //adjust the baseline and flip the pulse
	}
	delete []array;
	return adjust;
}

/*Funciton used to cut the shock of the signal, the shock will win the competition to the real pulse peak and cause the result confusing, the algorithum is search and tag all the sample higher then 50, then compare them to the neighbours, if the difference larger then 50, just rewrite them to their neighbours*/
vector<int> CutShock (vector<int> &l)
{
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
		if (abs(array[i])>1500) // deal with the condition that the begining of the pulse is really low , then a hugh step to real pulse
		array[i]=0;
		if (abs(array[i])>threshold)
		spik.push_back(i);
		
	} // store all the position that pulse higher the threhold
	for (int i=0; i<spik.size(); i++)
	{
		arraypos=spik[i];
		lefttest=abs (array[arraypos]-array[arraypos-1]);
		righttest=abs (array[arraypos]-array[arraypos+1]);	// search all the potential spike, identify them by compare to the neighbour, if both neighbour have 50 different, then identify them as spike
		if (lefttest>50&&righttest>50)     // currently make it cut at 50, so it will not cut the shock less then 50
		array[arraypos]=array[arraypos-1]; //rewrite the spike by make them equal to their neighbour.
	}
	for (int i=0; i<l.size(); i++)
	{
		cutpulse.push_back(array[i]); //rewrite the cutted spike pulse back to the pulse
	}
	delete []array;
	return cutpulse;
}
/*Function to implemented CFD method to find timing information*/
int timeCFD (vector<int>&l)
{
	vector<int> modifiedpulse;
	modifiedpulse=CFDpulse(l);

    int minpos=minfind(modifiedpulse);
    int maxpos=maxfind(modifiedpulse);
    int FirstNegative=minpos;
    int FirstPositive=maxpos;
    double zeropos=0;
    //cout << minpos << "\t" << maxpos << "\t" << modifiedpulse[minpos] << "\t" <<  modifiedpulse[maxpos] <<endl;
    if (minpos < maxpos && modifiedpulse[minpos]<0 && modifiedpulse[maxpos]>0)
    {        
            for (int i=0; i<abs(maxpos-minpos); i++)
            {
                    if (modifiedpulse[maxpos-i]<0)
                    {
                            FirstNegative=maxpos-i;
                            break;
                    }                
            }
            for (int i=0; i<abs(maxpos-minpos); i++)
            {
                    if (modifiedpulse[minpos+i]>0)
                    {
                            FirstPositive=minpos+i;
                            break;
                    }                
            }
            //zeropos=((modifiedpulse[FirstPositive]-modifiedpulse[FirstNegative])*FirstPositive-(FirstPositive-FirstNegative)*modifiedpulse[FirstPositive])/(modifiedpulse[FirstPositive]-modifiedpulse[FirstNegative]);
            zeropos=(FirstPositive+FirstNegative)/2;
    }

    return zeropos;
                     
}

/*Function use fraction directly from pulse shape*/

int PulseCFD (vector<int>&l)
{
	int *array = new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];	
	}
	double fraction = 0.6;
    int peakpos=maxfind(l);
    int leftzeropos=0;
    for (int i=0; i<peakpos;i++)
	{
		if (array[peakpos-i]<10)
		{
			leftzeropos=peakpos-i;
			break;
		}		
	}
    int FirstNegative=leftzeropos;
    int FirstPositive=peakpos;
    double zeropos=0;
	double PulseFrac = fraction * l[peakpos];
    //cout << minpos << "\t" << maxpos << "\t" << modifiedpulse[minpos] << "\t" <<  modifiedpulse[maxpos] <<endl;
    for (int i=0; i<abs(peakpos-leftzeropos); i++)
    {
            if (l[peakpos-i]<PulseFrac)
            {
                    FirstNegative=peakpos-i;
                    break;
            }                
    }
    for (int i=0; i<abs(peakpos-leftzeropos); i++)
    {
            if (l[leftzeropos+i]>PulseFrac)
            {
                    FirstPositive=leftzeropos+i;
                    break;
            }                
    }
	zeropos=(FirstPositive+FirstNegative)/2;
    return zeropos;
                     
}


/*function to creat a CFD vector*/

vector<int> CFDpulse (vector<int> &l)
{
	vector<int> cfd;
	int* arrayA=new int[l.size()];
	int* arrayB=new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		if (i<30)
		arrayA[i]=0;
		else	        
		arrayA[i]=l[i-30];
	}
	for (int i=0; i<l.size(); i++)
	{
		arrayB[i]=-l[i]*0.4;
	}
	
	for (int i=0; i<l.size(); i++)
	{
		cfd.push_back(arrayA[i]+arrayB[i]); //adjust the baseline and flip the pulse
	}

	delete []arrayA;
	delete []arrayB;
	return cfd;
}


/*Function used ot find the position where the sample is maximum*/
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

/*Function used ot find the position where the sample is minmum*/
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

/*Function used to add sample together, majorly used as integral of the pulse to find the energy, variable will be the vector hold the pulse sample information, and the window for integration*/
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

/*Fuction used to find psd, majorly call the function sum to get result with two different integration window and find the ratio of them*/
double psd (vector<int>&l,int number1,int number2, int number3, int number4)
{
	/*find the ratio between peakpush(number2), peakleftzero(number1), peakrightzero(number3 currently define as the end of pulse)*/	
	int* array=new int[l.size()];
        for (int i=0; i<l.size(); i++)
        {
                array[i]=l[i];
        }
	int peakleft=number1;
	int peakpush=number2;
	int peakright=number3;
        int peakend=number4;
	int totalenergy=sum(l,peakleft,peakright);
	int tailenergy=sum(l,peakpush,peakend);
	double psdana=0.;
	psdana=(double (tailenergy)/double (totalenergy));
	return psdana;	
}

/*Function to map the channel to get related X position*/
int mapXID (Matrix &l , int scrod, int row, int col,int chanl)
{
	int xid=0;
	int chan=scrod*1000+row*100+col*10+chanl;
        //cout << chan << endl;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (chan==l[i][j])
			{
				xid=j;
				break;
			}
		}
	}
	return xid;
}

/*Function to map the channel to get related Y position*/
int mapYID (Matrix &l ,int scrod, int row, int col,int chanl)
{
	int yid=0;
	int chan=scrod*1000+row*100+col*10+chanl;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (chan==l[i][j])
			{
				yid=i;
				break;
			}
		}
	}
	return yid;
}

/*Function to map the trig channel to get related X position, the difference between this and former function is this one read the channel number directly*/
int mapXIDtrig (Matrix &l , int ch)
{
	int xid=0;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (ch==l[i][j])
			{
				xid=j;
			}
		}
	}
	return xid;
}

/*Function to map the trig channel to get related Y position, the difference between this and former function is this one read the channel number directly*/
int mapYIDtrig (Matrix &l , int ch)
{
	int yid=0;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (ch==l[i][j])
			{
				yid=i;
			}
		}
	}
	return yid;
}

int lowerCubeBound (int i)
{
	if (i>0)
	return i-1;
	
	else 
	return 0;
}

int upperCubeBound (int i)
{
	if (i<4)
	return i+1;
	else 
	return 4;
}

int lowerBinBound (int i, int j)
{
	if (i-j>0)
	return i-j;
	
	else 
	return 0;
}

int upperBinBound (int i, int j)
{
	if (i+j<6000)
	return i+j;
	else 
	return 6000;
}


/*Function to be used to deal with the pulse, get useful informaiton for futher anlaysis
	First: Call the CutShock and adjust function, this two funciton basically modified the pulse shape, cut the shock, readjust the baseline and flip the pulse
	Second: Call the math function, majorly sum() to get integral, and also find peak information
	Third: store all the information into certain vector named pulseinfo;
		pulseinfo [0] -> Energy
		pulseinfo [1] -> Peak
		pulseinfo [2] -> PSD
		pulseinfo [3] -> peakpos ---realted to timing issue, needs more work to add CFD ability
*/
vector<double> pulseProcess (vector<int>&l)
{
        //cout<< "start " << endl;
	vector<double> pulseinfo;	
	vector<int> cuttedpulse=CutShock(l);					
	vector<int> adjustedpulse=flip(cuttedpulse);
	int *pulsey=new int[adjustedpulse.size()];
	int *xaxis=new int[adjustedpulse.size()];
	int totalenergy=0;		
	int peakamp=0;
	double psdratio=0.;	
	int peakpos=0;	
	int threshold=5;
	int leftzeropos=0;
	int rightzeropos=adjustedpulse.size();
        int pulsetiming=0;	
	for (int j=0; j<adjustedpulse.size(); j++)
	{
		pulsey[j]=adjustedpulse[j];
		xaxis[j]=j;
	}
	peakpos=maxfind(adjustedpulse);
	peakamp=pulsey[peakpos];
	for (int j=0; j<peakpos;j++)
	{
		if (pulsey[peakpos-j]<threshold)
		{
			leftzeropos=peakpos-j;
			break;
		}		
	}
	for (int j=peakpos; j<adjustedpulse.size();j++)
	{
		if (pulsey[j]<threshold)
		{
			rightzeropos=j;
			break;
		}		
	}
        int psdpos=peakpos+95;
	if (psdpos>adjustedpulse.size())
	psdpos=peakpos; 
	totalenergy=sum(adjustedpulse,leftzeropos,rightzeropos);
	psdratio = psd (adjustedpulse,leftzeropos,psdpos,rightzeropos,adjustedpulse.size());
	pulsetiming= PulseCFD (adjustedpulse);
	/*event  information storage*/						
	pulseinfo.push_back(totalenergy);
	pulseinfo.push_back(peakamp);
	pulseinfo.push_back(psdratio);
	pulseinfo.push_back(pulsetiming);
	adjustedpulse.clear();
        //cout << "end" << endl;
	delete []pulsey;
	delete []xaxis;	
	return pulseinfo;						
}

/*Function to separate AB events, when variable a=0, A event being filled, otherwise fill B event*/
vector<int> pulseAB (vector<int> &l, int a)
{
	int *array = new int[l.size()];
	for (int i=0; i<l.size(); i++)
	{
		array[i]=l[i];	
	}
	
	int halfcount=l.size()/2;
	/*when a=0, fill the A event, when a=1, fill the B event, cut half and fill into both A and B*/
	vector<int> pulse;	
	if (a==0)
	{
		for (int i=0; i<halfcount; i++)
		{
			pulse.push_back(array[i]);
			
		}
	}
	else 
	{
		for (int i=halfcount; i<l.size(); i++)
		{
			pulse.push_back(array[i]);
		}	
	}
	delete []array;
	return pulse;
}

/*Function define the time difference between AB event time window*/
int timedif (int timebegin, int timeend)
{
	int timeinterval=0;
	if (timeend>timebegin)
		timeinterval=timeend-timebegin+1;
	else
		timeinterval=512-abs(timeend-timebegin)+1;
	return timeinterval;
}

/*function to be used as post analysis*/

void print2Darray (int a[10][10])
{
	for (int i=0; i< 10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (j==9)
			{
				cout << a[i][j] << endl;
			}
			else
			{
				cout << a[i][j] << "\t" ;
			}
		}
	}
}

void printfloat2Darray ( float a[10][10])
{
	for (int i=0; i< 10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (j==9)
			{
				cout << a[i][j] << endl;
			}
			else
			{
				cout << a[i][j] << "\t" ;
			}
		}
	}
}

void printMatrix (Matrix &l)
{
	for (int i=0; i< l.size(); i++)
	{
		for (int j=0; j<10; j++)
		{
			if (j==9)
			{
				cout << l[i][j] << endl;
			}
			else
			{
				cout << l[i][j] << "\t" ;
			}
		}
	}
}

void recalibration (float a[10][10], int rl,int rh,int cl, int ch)
{
	float multiply=100;
	for (int i=rl ; i< rh+1 ; i++)
	{
		for (int j=cl; j<ch+1; j++)
		{
			multiply=multiply*a[i][j]/10;
			cout << "multiply=" << multiply << endl;
		}
		
	}
	for (int i=rl ; i< rh+1 ; i++)
	{
		for (int j=cl; j<ch+1; j++)
		{
			a[i][j]=multiply/a[i][j]/1000000000;
		}
	}
}


bool pick(int a[10][10],int b[10][10], Matrix &l , int n)
{
	/*deal with the threshold*/
	bool eventveto=false;
	for (int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)
		{
			if (i<5|| j<5)
			{			
		
				if (l[i+10*n][j]>b[i][j] || l[i+10*n][j] < a[i][j] )
				{			
					//cout << "veto from pick \t" << l[i+10*n][j] << endl;			
					eventveto=true;
				}
			}
		}
	} 
	return eventveto;
}

bool noiseveto (Matrix &l , int x, int y, int n)
{
	bool eventveto=true;
	if (l[y+10*n][x]>1000)
	eventveto=false;
	return eventveto;
	
}

bool cubeveto (Matrix &l , int x, int y, int z, double c, int n)
{
	bool eventveto=false;
	bool sideone=pmtsideveto(l,x,y,c,0,n);
	bool sidetwo=pmtsideveto(l,x,z,c,1,n);
	bool sidethree=pmtsideveto(l,z,y,c,2,n);
	if (sideone || sidethree || sidetwo)	
		eventveto=true;
	return eventveto;
}

bool twosideveto (Matrix &l , int x, int y, int z,double c, int n)
{
	bool eventveto=false;
	bool sideone=pmtsideveto(l,x,y,c,0,n);
	bool sidethree=pmtsideveto(l,z,y,c,2,n);
	if (sideone || sidethree)	
		eventveto=true;
	return eventveto;
}


bool pmtsideveto (Matrix &l , int a, int b, double c, int side, int n)
{
	bool eventveto=false;

	int matrixpickx=0;
	int matrixpicky=0;
	int xid=0;
	int yid=0;
	/*
	side=0 -> x=0-4, y=0-4;
	side=1 -> x=0-4, y=5-9;
	side=2 -> x=5-9, y=0-4; 
	*/
	if (side==0)
	{
		matrixpickx=0;
		matrixpicky=0;
		xid=a;
		yid=b;
	}
	else if (side==1)
	{
		matrixpickx=5;
		matrixpicky=0;
		xid=a;
		yid=b+5;
	}
	else
	{
		matrixpickx=0;
		matrixpicky=5;
		xid=a+5;
		yid=b;
	}
	int trigchan, testchan;	
	//cout << l[yid+n*10][xid] << endl;
	for (int i=matrixpickx; i<matrixpickx+5; i++)
	{
		for (int j=matrixpicky; j<matrixpicky+5; j++)
		{
			
			if (i!=yid || j!=xid)
			{		
				if (l[i+10*n][j]>l[yid+n*10][xid]*c)
				{
					eventveto=true;
					//cout << "veto "<<testchan << "\t" << trigchan<<"\t" << l[yid+n*10][xid] <<endl;
				}
			}
				
		}
	} 
	return eventveto;
}

/*PSD*/

vector<double> pulsePSDProcess (vector<int>&l,int push)
{
	vector<double> pulseinfo;	
	vector<int> cuttedpulse=CutShock(l);					
	vector<int> adjustedpulse=flip(cuttedpulse);
	int *pulsey=new int[adjustedpulse.size()];
	int *xaxis=new int[adjustedpulse.size()];
	int totalenergy=0;		
	int peakamp=0;
	double psdratio=0.;	
	int peakpos=0;
	int threshold=0;
	int leftzeropos=0;
        int rightzeropos=adjustedpulse.size();
	//cout << "pulse being called" << endl;	
	for (int j=0; j<adjustedpulse.size(); j++)
	{
		pulsey[j]=adjustedpulse[j];
		//cout << pulsey[j] ;
		xaxis[j]=j;
	}
	peakpos=maxfind(adjustedpulse);
	int psdpos=peakpos+push;
	if (psdpos>adjustedpulse.size())
	psdpos=adjustedpulse.size()-2; 
	for (int j=0; j<peakpos;j++)
	{
		if (pulsey[peakpos-j]<threshold)
		{
			leftzeropos=peakpos-j;
			break;
		}		
	}
        for (int j=peakpos; j<adjustedpulse.size();j++)
	{
		if (pulsey[j]<threshold)
		{
			rightzeropos=j;
			break;
		}		
	}
	totalenergy=sum(adjustedpulse,leftzeropos,adjustedpulse.size());
	//cout << peakamp << endl;
	psdratio = psd (adjustedpulse,leftzeropos,psdpos, rightzeropos, adjustedpulse.size());
	/*event  information storage*/						
	pulseinfo.push_back(totalenergy);
	pulseinfo.push_back(psdratio);
	//cout << "totalenergy=" << totalenergy << "\t" << "psdratio=" << psdratio << endl;
	adjustedpulse.clear();
	delete []pulsey;
	delete []xaxis;	
	return pulseinfo;						
}

/*	Matrix energyA(n*10,Row(10));	
	MatrixPsd psdA(n*10,RowPsd(10));
	Matrix peakA(n*10,Row(10));
	Matrix timingA(n*10,Row(10));
	Matrix energyB(n*10,Row(10));	
	MatrixPsd psdB(n*10,RowPsd(10));
	Matrix peakB(n*10,Row(10));
	Matrix timingB(n*10,Row(10));	
	Matrix cube(n*10,Row(10));
	Matrix timingAB(n*10,Row(10));
	
*/

int ABEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, MatrixPsd &psdA, MatrixPsd &psdB, Matrix &timingAB)
{
	cout << "function AB pick" << endl;
	int xidA, yidA, zidA, xidB, yidB, zidB,faceid;
	string MapOption;
	TString xAid, yAid, zAid, xBid, yBid, zBid, Aname, Bname, fracname;
	Aname="Event-A-";
	Bname="Event-B-";
	cout << "input xid for A event" << endl;
	cin >> xidA;
	cout << "input yid for A event" << endl;
	cin >> yidA;
	cout << "input zid for A event" << endl;
	cin >> zidA;
	cout << "input xid for B event" << endl;
	cin >> xidB;
	cout << "input yid for B event" << endl;
	cin >> yidB;
	cout << "input zid for B event" << endl;
	cin >> zidB;
	cout << "Input the face pmt being trigged" << endl;
	cin >> faceid;
	cout << "Do you need 2D events mapping? (y/n)" << endl;
	cin >> MapOption;
	int facepick;
	double fracVeto=0.;
	cout << "input factor for veto" << endl;
	cin >> fracVeto;	
	fracname.Form("%f",fracVeto);	
	xAid.Form("%d",xidA);
	xBid.Form("%d",xidB);
	yAid.Form("%d",yidA);
	yBid.Form("%d",yidB);
	zAid.Form("%d",zidA);
	zBid.Form("%d",zidB);
	TString rootid=Aname+xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+"-"+fracname;
	
	/*create root file to hold spectrum*/
	TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",filename);
	TFile* f=new TFile (analysisfile+"-"+rootid+rootap,"recreate");
	
	/*record picked event number for future analysis purpose*/
	ofstream foutevent;
	foutevent.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".txt");

	/*Output a csv file for graph purpose*/
	ofstream foutdata;
	foutdata.open (xAid+"-"+yAid+"-"+zAid+"-"+Aname+xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+"-"+fracname+".csv");
	foutdata << "Event Num, Event Count, Orange A, Blue A, Green A, Orange B, Blue B, Green B, Orange A PSD, Blue A PSD, Green A PSD, Orange B PSD, Blue B PSD, Green B PSD,timingAB-Orange, timingAB-Blue, timingAB-Green, timingA-Orange. timingA-Blue, timingA-Green, timingB-Orange, timingB-Blue, timingB-Green, \n";
	
	int eventcount=0;
	
	/*create category of histogram to hold spectrum*/
	TString energystr="energy ";
	TString peakstr="peak ";
	TString psdstr="psd ";
	TH1D peakhistA[100];
	TH1D integralhistA[100];
	TH1D peakhistB[100];
	TH1D integralhistB[100];
	TH2D psdhistA[100];
	TH2D psdhistB[100];
    TH2D psdhistAB[100];
	TString titleA=" A";
	TString titleB=" B";
	int xid=0;
	int yid=0;
	for (int i=0; i<100; i++)
	{
		TString xidname,yidname;
		xid=i%10;
		yid=i/10;
		xidname.Form ("%d",xid);
		yidname.Form ("%d",yid);
		TString titlename = yidname+" "+xidname;
		TString histname = yidname+" "+xidname;
		peakhistA[i]=TH1D(peakstr+histname+titleA,titlename,2000,0,2000);
		integralhistA[i]=TH1D(energystr+histname+titleA,titlename,6000,0,350000); // needs to be modified the size and maximum
		psdhistA[i]=TH2D(psdstr+histname+titleA,titlename,1500,0,350000,100,0,1);
		peakhistB[i]=TH1D(peakstr+histname+titleB,titlename,2000,0,2000);
		integralhistB[i]=TH1D(energystr+histname+titleB,titlename,6000,0,350000); // needs to be modified the size and maximum
		psdhistB[i]=TH2D(psdstr+histname+titleB,titlename,1500,0,350000,100,0,1);
		psdhistAB[i]=TH2D(psdstr+histname+titleA+titleB,titlename,1500,0,350000,100,0,1);
	}

	TH1D *timeA12 = new TH1D("timeA12", "timeA-1-2-dif", 100, -50, 50);
	TH1D *timeA13 = new TH1D("timeA13", "timeA-1-3-dif", 100, -50, 50);
	TH1D *timeA23 = new TH1D("timeA23", "timeA-2-3-dif", 100, -50, 50);
	TH1D *timeB12 = new TH1D("timeB12", "timeB-1-2-dif", 100, -50, 50);
	TH1D *timeB13 = new TH1D("timeB13", "timeB-1-3-dif", 100, -50, 50);
	TH1D *timeB23 = new TH1D("timeB23", "timeB-2-3-dif", 100, -50, 50);
	TH1D *timeABtrig= new TH1D("timeAB", "timeAB-dif", 3000, -10, 30000);

	int TimeAFace1, TimeAFace2, TimeAFace3, TimeBFace1, TimeBFace2,TimeBFace3, TimeABFace;
				

	int histTag=0; // based on the x and y value to define the hist number

	bool vetoA=true;
	bool vetoB=true;
	bool vetoApick=true;
	bool vetoBpick=true;
	
	for (int i=0; i<n; i++)
	{
		/*currently disable the threshold ability, should be added when needed*/
		//vetoApick=pick(lowthreshold,highthreshold,energyA,i);
		//vetoBpick=pick(lowthreshold,highthreshold,energyB,i);

		if (true)
		{
			//cout << eventID[i] << "A 3 data"<< endl;
			vetoA=cubeveto(energyA, xidA,yidA,zidA, 1, i);
			//cout << eventID[i] << "B 3 data"<< endl;
			vetoB=cubeveto(energyB, xidB,yidB,zidB, fracVeto, i);
			//cout << vetoA << "\t" << vetoB << endl;
			if (!vetoA && !vetoB)
			{
				foutevent << eventID[i] << endl;
				eventcount++;		
				foutdata << eventID[i] << "," << i << "," << energyA[yidA+i*10][xidA] << "," <<  energyA[yidA+i*10][zidA+5] << "," << energyA[5+i*10+zidA][xidA] << "," << energyB[yidB+i*10][xidB] << "," <<  energyB[yidB+i*10][zidB+5] << "," << energyB[5+i*10+zidB][xidB] << "," << psdA[yidA+i*10][xidA] << "," <<  psdA[yidA+i*10][zidA+5] << "," << psdA[5+i*10+zidA][xidA] << "," << psdB[yidB+i*10][xidB] << "," <<  psdB[yidB+i*10][zidB+5] << "," << psdB[5+i*10+zidB][xidB]<<timingAB[yidB+i*10][xidB] << "," << timingAB[yidB+i*10][5+zidB] << "," << timingAB[5+zidB+i*10][xidB] << "," << timingA[yidA+i*10][xidA] << "," << timingA[yidA+i*10][5+zidA] << "," << timingA[5+zidA+i*10][xidA] << "," << timingB[yidB+i*10][xidB] << "," << timingB[yidB+i*10][5+zidB] << "," << timingB[5+zidB+i*10][xidB] << "," << "\n";
				for(int j=0;j<10;j++)
				{
					for(int k=0; k<10; k++)
					{
						if (j<5||k<5)
						{
							histTag=j*10+k;								
							peakhistA[histTag].Fill(peakA[j+10*i][k]);
							integralhistA[histTag].Fill(energyA[j+10*i][k]);
							psdhistA[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
							peakhistB[histTag].Fill(peakB[j+10*i][k]);
							integralhistB[histTag].Fill(energyB[j+10*i][k]);
							psdhistB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
							psdhistAB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
							psdhistAB[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
						
						}
					}
				}

				if (MapOption=="y")
				{			
					TString eventchar;						
					eventchar.Form ("%d",eventID[i]);	
					TH2D* temp2dhisA=new TH2D("eventA "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5 || k<5)
							{
								for (int l=0;l<energyA[j+i*10][k];l++)
								{
									temp2dhisA->Fill(k,9-j);						
								}
							}
						}
					}
					TH2D* temp2dhisB=new TH2D("eventB "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5 || k<5)
							{
								for (int l=0;l<energyB[j+i*10][k];l++)
								{
									temp2dhisB->Fill(k,9-j);						
								}
							}
						}
					}
					//cout << "wrtie the 2D map" << endl;
					temp2dhisA->Write();
					temp2dhisB->Write();
					delete temp2dhisA;
					delete temp2dhisB;
				}

				TimeAFace1=timingA[yidA+i*10][xidA]; 
				TimeAFace2=timingA[yidA+i*10][5+zidA];
				TimeAFace3=timingA[5+zidA+i*10][xidA];
				TimeBFace1=timingB[yidB+i*10][xidB];
				TimeBFace2=timingB[yidB+i*10][5+zidB];
				TimeBFace3=timingB[5+zidB+i*10][xidB];
				if(faceid==1)
					TimeABFace=timingAB[yidA+10*i][xidA];
				else if (faceid==2)
					TimeABFace=timingAB[yidA+i*10][5+zidA];
				else
					TimeABFace=timingAB[5+zidA+i*10][xidA];
				timeA12->Fill(TimeAFace1-TimeAFace2);
				timeA13->Fill(TimeAFace1-TimeAFace3);
				timeA23->Fill(TimeAFace2-TimeAFace3);
				timeB12->Fill(TimeBFace1-TimeBFace2);
				timeB13->Fill(TimeBFace1-TimeBFace3);
				timeB23->Fill(TimeBFace2-TimeBFace3);
				timeABtrig->Fill(TimeABFace);
			}
		}
	}
	cout << "total "<< eventcount << "\t events being found" << endl;
	foutdata.clear();
	foutdata.close();
	foutevent.clear();
	foutevent.close();
	f->Write();
	f->Close();
    return 0;
}

int BEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, MatrixPsd &psdA, MatrixPsd &psdB, Matrix &timingAB)
{
	cout << "function B pick" << endl;
	int xidA, yidA, zidA, xidB, yidB, zidB;
	string MapOption;
	int facepick;
	double fracVeto=0.;
	TString xAid, yAid, zAid, xBid, yBid, zBid,Aname, Bname, fracname;
	Aname="Event-A-";
	Bname="Event-B-";
    cout << "input factor for veto" << endl;
	cin >> fracVeto;
	cout << "input xid for B event" << endl;
	cin >> xidB;
	cout << "input yid for B event" << endl;
	cin >> yidB;
	cout << "input zid for B event" << endl;
	cin >> zidB;
	cout << "Do you need 2D events mapping? (y/n)" << endl;
	cin >> MapOption;
	fracname.Form("%f",fracVeto);
	xAid.Form("%d",xidA);
	xBid.Form("%d",xidB);
	yAid.Form("%d",yidA);
	yBid.Form("%d",yidB);
	zAid.Form("%d",zidA);
	zBid.Form("%d",zidB);
	TString rootid=Bname+xBid+"-"+yBid+"-"+zBid+"-"+fracname;
	/*create root file to hold spectrum*/
	TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",filename);
	TFile* f=new TFile (analysisfile+"-"+rootid+rootap,"recreate");
	ofstream foutevent;
	foutevent.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".txt");
	ofstream foutdata;
	foutdata.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".csv");
	foutdata << "Event Num, Event Count, Orange A, Blue A, Green A, Orange B, Blue B, Green B, Orange A PSD, Blue A PSD, Green A PSD, Orange B PSD, Blue B PSD, Green B PSD,timingAB-Orange, timingAB-Blue, timingAB-Green, timingA-Orange. timingA-Blue, timingA-Green, timingB-Orange, timingB-Blue, timingB-Green, \n";
	int eventcount=0;
	/*create category of histogram to hold spectrum*/
	TString energystr="energy ";
	TString peakstr="peak ";
	TString psdstr="psd ";
	TH1D peakhistA[100];
	TH1D integralhistA[100];
	TH1D peakhistB[100];
	TH1D integralhistB[100];
	TH2D psdhistA[100];
	TH2D psdhistB[100];
	TH2D psdhistAB[100];
	TString titleA=" A";
	TString titleB=" B";
	int xid=0;
	int yid=0;
	for (int i=0; i<100; i++)
	{
		TString xidname,yidname;
		xid=i%10;
		yid=i/10;
		xidname.Form ("%d",xid);
		yidname.Form ("%d",yid);
		TString titlename = yidname+" "+xidname;
		TString histname = yidname+" "+xidname;
		peakhistA[i]=TH1D(peakstr+histname+titleA,titlename,2000,0,2000);
		integralhistA[i]=TH1D(energystr+histname+titleA,titlename,6000,0,350000); // needs to be modified the size and maximum
		psdhistA[i]=TH2D(psdstr+histname+titleA,titlename,1500,0,350000,100,0,1);
		peakhistB[i]=TH1D(peakstr+histname+titleB,titlename,2000,0,2000);
		integralhistB[i]=TH1D(energystr+histname+titleB,titlename,6000,0,350000); // needs to be modified the size and maximum
		psdhistB[i]=TH2D(psdstr+histname+titleB,titlename,1500,0,350000,100,0,1);
    	psdhistAB[i]=TH2D(psdstr+histname+titleA+titleB,titlename,1500,0,350000,100,0,1);
	}
	int histTag=0;
	bool vetoA=true;
	bool vetoB=true;
	bool vetoApick=true;
	bool vetoBpick=true;
	for (int i=0; i<n; i++)
	{
		//vetoApick=pick(lowthreshold,highthreshold,energyA,i);
		//vetoBpick=pick(lowthreshold,highthreshold,energyB,i);
		if (true)
		{
			//cout << eventID[i] << "A 3 data"<< endl;
			//vetoA=cubeveto(energyA, xidA,yidA,zidA,i);
			//cout << eventID[i] << "B 3 data"<< endl;
			vetoB=cubeveto(energyB, xidB,yidB,zidB,fracVeto,i);
			//cout << vetoA << "\t" << vetoB << endl;
			if (!vetoB)
			{
				foutevent << eventID[i] << endl;
				eventcount++;		
				foutdata << eventID[i] << "," << i << "," << energyA[yidA+i*10][xidA] << "," <<  energyA[yidA+i*10][zidA+5] << "," << energyA[5+i*10+zidA][xidA] << "," << energyB[yidB+i*10][xidB] << "," <<  energyB[yidB+i*10][zidB+5] << "," << energyB[5+i*10+zidB][xidB] << "," << psdA[yidA+i*10][xidA] << "," <<  psdA[yidA+i*10][zidA+5] << "," << psdA[5+i*10+zidA][xidA] << "," << psdB[yidB+i*10][xidB] << "," <<  psdB[yidB+i*10][zidB+5] << "," << psdB[5+i*10+zidB][xidB]<<timingAB[yidB+i*10][xidB] << "," << timingAB[yidB+i*10][5+zidB] << "," << timingAB[5+zidB+i*10][xidB] << "," << timingA[yidA+i*10][xidA] << "," << timingA[yidA+i*10][5+zidA] << "," << timingA[5+zidA+i*10][xidA] << "," << timingB[yidB+i*10][xidB] << "," << timingB[yidB+i*10][5+zidB] << "," << timingB[5+zidB+i*10][xidB] << "," << "\n";
				for(int j=0;j<10;j++)
				{
					for(int k=0; k<10; k++)
					{
						if (j<5||k<5)
						{
							histTag=j*10+k;								
							peakhistA[histTag].Fill(peakA[j+10*i][k]);
							integralhistA[histTag].Fill(energyA[j+10*i][k]);
							psdhistA[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
							peakhistB[histTag].Fill(peakB[j+10*i][k]);
							integralhistB[histTag].Fill(energyB[j+10*i][k]);
							psdhistB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
							psdhistAB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
							psdhistAB[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
						
						}
					}
				}
				if (MapOption=="y")
				{			
					TString eventchar;						
					eventchar.Form ("%d",eventID[i]);	
					TH2D* temp2dhisA=new TH2D("eventA "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5 || k<5)
							{
								for (int l=0;l<energyA[j+i*10][k];l++)
								{
									temp2dhisA->Fill(k,9-j);						
								}
							}
						}
					}
					TH2D* temp2dhisB=new TH2D("eventB "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5 || k<5)
							{
								for (int l=0;l<energyB[j+i*10][k];l++)
								{
									temp2dhisB->Fill(k,9-j);						
								}
							}
						}
					}
					//cout << "wrtie the 2D map" << endl;
					temp2dhisA->Write();
					temp2dhisB->Write();
					delete temp2dhisA;
					delete temp2dhisB;
				}
			}
		}
	}
	cout << "total "<< eventcount << "\t events being found" << endl;
	foutdata.clear();
	foutdata.close();
	foutevent.clear();
	foutevent.close();
	f->Write();
	f->Close();
	return 0;
}

int ABtEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, MatrixPsd &psdA, MatrixPsd &psdB, Matrix &timingAB)
{
	cout << "function AB pick" << endl;
	int xidA, yidA, zidA, xidB, yidB, zidB,faceid;
	string MapOption;
	TString xAid, yAid, zAid, xBid, yBid, zBid, Aname, Bname, fracname;
	Aname="Event-AT-";
	Bname="Event-BT-";
	cout << "input xid for A event" << endl;
	cin >> xidA;
	cout << "input yid for A event" << endl;
	cin >> yidA;
	cout << "input zid for A event" << endl;
	cin >> zidA;
	cout << "input xid for B event" << endl;
	cin >> xidB;
	cout << "input yid for B event" << endl;
	cin >> yidB;
	cout << "input zid for B event" << endl;
	cin >> zidB;
	cout << "Input the face pmt being trigged" << endl;
	cin >> faceid;
	cout << "Do you need 2D events mapping? (y/n)" << endl;
	cin >> MapOption;
	int facepick;
	double fracVeto=0.;
	cout << "input factor for veto" << endl;
	cin >> fracVeto;	
	fracname.Form("%f",fracVeto);	
	xAid.Form("%d",xidA);
	xBid.Form("%d",xidB);
	yAid.Form("%d",yidA);
	yBid.Form("%d",yidB);
	zAid.Form("%d",zidA);
	zBid.Form("%d",zidB);
	TString rootid=Aname+xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+"-"+fracname;
	
	/*create root file to hold spectrum*/
	TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",filename);
	TFile* f=new TFile (analysisfile+"-"+rootid+rootap,"recreate");
	
	/*record picked event number for future analysis purpose*/
	ofstream foutevent;
	foutevent.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".txt");

	/*Output a csv file for graph purpose*/
	ofstream foutdata;
	foutdata.open (xAid+"-"+yAid+"-"+zAid+"-"+Aname+xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+"-"+fracname+".csv");
	foutdata << "Event Num, Event Count, Orange A, Blue A, Green A, Orange B, Blue B, Green B, Orange A PSD, Blue A PSD, Green A PSD, Orange B PSD, Blue B PSD, Green B PSD,timingAB-Orange, timingAB-Blue, timingAB-Green, timingA-Orange. timingA-Blue, timingA-Green, timingB-Orange, timingB-Blue, timingB-Green, \n";
	
	int eventcount=0;
	
	/*create category of histogram to hold spectrum*/
	TString energystr="energy ";
	TString peakstr="peak ";
	TString psdstr="psd ";
	TH1D peakhistA[100];
	TH1D integralhistA[100];
	TH1D peakhistB[100];
	TH1D integralhistB[100];
	TH2D psdhistA[100];
	TH2D psdhistB[100];
    TH2D psdhistAB[100];
	TString titleA=" A";
	TString titleB=" B";
	int xid=0;
	int yid=0;
	for (int i=0; i<100; i++)
	{
		TString xidname,yidname;
		xid=i%10;
		yid=i/10;
		xidname.Form ("%d",xid);
		yidname.Form ("%d",yid);
		TString titlename = yidname+" "+xidname;
		TString histname = yidname+" "+xidname;
		peakhistA[i]=TH1D(peakstr+histname+titleA,titlename,2000,0,2000);
		integralhistA[i]=TH1D(energystr+histname+titleA,titlename,6000,0,350000); // needs to be modified the size and maximum
		psdhistA[i]=TH2D(psdstr+histname+titleA,titlename,1500,0,350000,100,0,1);
		peakhistB[i]=TH1D(peakstr+histname+titleB,titlename,2000,0,2000);
		integralhistB[i]=TH1D(energystr+histname+titleB,titlename,6000,0,350000); // needs to be modified the size and maximum
		psdhistB[i]=TH2D(psdstr+histname+titleB,titlename,1500,0,350000,100,0,1);
		psdhistAB[i]=TH2D(psdstr+histname+titleA+titleB,titlename,1500,0,350000,100,0,1);
	}

	TH1D *timeA12 = new TH1D("timeA12", "timeA-1-2-dif", 100, -50, 50);
	TH1D *timeA13 = new TH1D("timeA13", "timeA-1-3-dif", 100, -50, 50);
	TH1D *timeA23 = new TH1D("timeA23", "timeA-2-3-dif", 100, -50, 50);
	TH1D *timeB12 = new TH1D("timeB12", "timeB-1-2-dif", 100, -50, 50);
	TH1D *timeB13 = new TH1D("timeB13", "timeB-1-3-dif", 100, -50, 50);
	TH1D *timeB23 = new TH1D("timeB23", "timeB-2-3-dif", 100, -50, 50);
	TH1D *timeABtrig= new TH1D("timeAB", "timeAB-dif", 3000, -10, 30000);

	int TimeAFace1, TimeAFace2, TimeAFace3, TimeBFace1, TimeBFace2,TimeBFace3, TimeABFace;
				

	int histTag=0; // based on the x and y value to define the hist number

	Matrix energyAtm(n*10,Row(10));	
	Matrix energyBtm(n*10,Row(10));
	int xtrigid, ytrigid;
	cout << "input xid for trigger event" << endl;
	cin >> xtrigid;
	cout << "input yid for trigger event" << endl;
	cin >> ytrigid;
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<10; j++)
		{
			for (int k=0; k<10; k++)
			{
				if (j!=ytrigid || k!=xtrigid)
				{		
					if ( abs(timingA[j+i*10][k]-timingA[ytrigid+i*10][xtrigid])>20)
						energyAtm[j+i*10][k]=0;
					else
						energyAtm[j+i*10][k]=energyA[j+i*10][k];
					if ( abs(timingB[j+i*10][k]-timingB[ytrigid+i*10][xtrigid])>20)
						energyBtm[j+i*10][k]=0;
					else
						energyBtm[j+i*10][k]=energyB[j+i*10][k];
				}
				else
				{
    				energyAtm[j+i*10][k]=energyA[j+i*10][k];
   					energyBtm[j+i*10][k]=energyB[j+i*10][k];
				}
			}
		}
	} 
        

	bool vetoA=true;
	bool vetoB=true;
	bool vetoApick=true;
	bool vetoBpick=true;
	
	for (int i=0; i<n; i++)
	{
		/*currently disable the threshold ability, should be added when needed*/
		//vetoApick=pick(lowthreshold,highthreshold,energyA,i);
		//vetoBpick=pick(lowthreshold,highthreshold,energyB,i);

		if (true)
		{
			//cout << eventID[i] << "A 3 data"<< endl;
			vetoA=cubeveto(energyA, xidA,yidA,zidA, 1, i);
			//cout << eventID[i] << "B 3 data"<< endl;
			vetoB=cubeveto(energyB, xidB,yidB,zidB, fracVeto, i);
			//cout << vetoA << "\t" << vetoB << endl;
			if (!vetoA && !vetoB)
			{
				foutevent << eventID[i] << endl;
				eventcount++;		
				foutdata << eventID[i] << "," << i << "," << energyA[yidA+i*10][xidA] << "," <<  energyA[yidA+i*10][zidA+5] << "," << energyA[5+i*10+zidA][xidA] << "," << energyB[yidB+i*10][xidB] << "," <<  energyB[yidB+i*10][zidB+5] << "," << energyB[5+i*10+zidB][xidB] << "," << psdA[yidA+i*10][xidA] << "," <<  psdA[yidA+i*10][zidA+5] << "," << psdA[5+i*10+zidA][xidA] << "," << psdB[yidB+i*10][xidB] << "," <<  psdB[yidB+i*10][zidB+5] << "," << psdB[5+i*10+zidB][xidB]<<timingAB[yidB+i*10][xidB] << "," << timingAB[yidB+i*10][5+zidB] << "," << timingAB[5+zidB+i*10][xidB] << "," << timingA[yidA+i*10][xidA] << "," << timingA[yidA+i*10][5+zidA] << "," << timingA[5+zidA+i*10][xidA] << "," << timingB[yidB+i*10][xidB] << "," << timingB[yidB+i*10][5+zidB] << "," << timingB[5+zidB+i*10][xidB] << "," << "\n";
				for(int j=0;j<10;j++)
				{
					for(int k=0; k<10; k++)
					{
						if (j<5||k<5)
						{
							histTag=j*10+k;								
							peakhistA[histTag].Fill(peakA[j+10*i][k]);
							integralhistA[histTag].Fill(energyA[j+10*i][k]);
							psdhistA[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
							peakhistB[histTag].Fill(peakB[j+10*i][k]);
							integralhistB[histTag].Fill(energyB[j+10*i][k]);
							psdhistB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
							psdhistAB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
							psdhistAB[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
						
						}
					}
				}

				if (MapOption=="y")
				{			
					TString eventchar;						
					eventchar.Form ("%d",eventID[i]);	
					TH2D* temp2dhisA=new TH2D("eventA "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5 || k<5)
							{
								for (int l=0;l<energyA[j+i*10][k];l++)
								{
									temp2dhisA->Fill(k,9-j);						
								}
							}
						}
					}
					TH2D* temp2dhisB=new TH2D("eventB "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5 || k<5)
							{
								for (int l=0;l<energyB[j+i*10][k];l++)
								{
									temp2dhisB->Fill(k,9-j);						
								}
							}
						}
					}
					//cout << "wrtie the 2D map" << endl;
					temp2dhisA->Write();
					temp2dhisB->Write();
					delete temp2dhisA;
					delete temp2dhisB;
				}

				TimeAFace1=timingA[yidA+i*10][xidA]; 
				TimeAFace2=timingA[yidA+i*10][5+zidA];
				TimeAFace3=timingA[5+zidA+i*10][xidA];
				TimeBFace1=timingB[yidB+i*10][xidB];
				TimeBFace2=timingB[yidB+i*10][5+zidB];
				TimeBFace3=timingB[5+zidB+i*10][xidB];
				if(faceid==1)
					TimeABFace=timingAB[yidA+10*i][xidA];
				else if (faceid==2)
					TimeABFace=timingAB[yidA+i*10][5+zidA];
				else
					TimeABFace=timingAB[5+zidA+i*10][xidA];
				timeA12->Fill(TimeAFace1-TimeAFace2);
				timeA13->Fill(TimeAFace1-TimeAFace3);
				timeA23->Fill(TimeAFace2-TimeAFace3);
				timeB12->Fill(TimeBFace1-TimeBFace2);
				timeB13->Fill(TimeBFace1-TimeBFace3);
				timeB23->Fill(TimeBFace2-TimeBFace3);
				timeABtrig->Fill(TimeABFace);
			}
		}
	}
	cout << "total "<< eventcount << "\t events being found" << endl;
	foutdata.clear();
	foutdata.close();
	foutevent.clear();
	foutevent.close();
	f->Write();
	f->Close();
    return 0;
}


int MapEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB )
{
	cout << "function mapping pick" << endl;		
	
	TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",filename);
	TFile* f=new TFile (analysisfile+"-energy-map"+rootap,"recreate");
	
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
			TH2D* temp2dhisA1=new TH2D("eventA1 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
			TH2D* temp2dhisB1=new TH2D("eventB1 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
			TH2D* temp2dhisA2=new TH2D("eventA2 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
			TH2D* temp2dhisB2=new TH2D("eventB2 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
			TH2D* temp2dhisA3=new TH2D("eventA3 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
			TH2D* temp2dhisB3=new TH2D("eventB3 "+eventchar, eventchar+"Energy mapping",10,0,5,10,0,5);
			for(int j=0;j<5;j++)
			{
				for(int k=0; k<5; k++)
				{
					{
						for (int l=0;l<energyA[j+i*10][k];l++)
						{
							temp2dhisA1->Fill(k,4-j);						
						}
						for (int l=0;l<energyB[j+i*10][k];l++)
						{
							temp2dhisB1->Fill(k,4-j);						
						}
					}
				}
			}
			for(int j=0;j<5;j++)
			{
				for(int k=5; k<10; k++)
				{
					{
						for (int l=0;l<energyA[j+i*10][k];l++)
						{
							temp2dhisA2->Fill(k-5,4-j);						
						}
						for (int l=0;l<energyB[j+i*10][k];l++)
						{
							temp2dhisB2->Fill(k-5,4-j);						
						}
					}
				}
			}
			for(int j=5;j<10;j++)
			{
				for(int k=0; k<5; k++)
				{
					{
						for (int l=0;l<energyA[j+i*10][k];l++)
						{
							temp2dhisA3->Fill(k,9-j);						
						}
						for (int l=0;l<energyB[j+i*10][k];l++)
						{
							temp2dhisB3->Fill(k,9-j);						
						}
					}
				}
			}				
		

			cout << "wrtie the 2D map" << endl;
			temp2dhisA1->Write();
			temp2dhisB1->Write();
			temp2dhisA2->Write();
			temp2dhisB2->Write();
			temp2dhisA3->Write();
			temp2dhisB3->Write();
			delete temp2dhisA1;
			delete temp2dhisB1;
			delete temp2dhisA2;
			delete temp2dhisB2;
			delete temp2dhisA3;
			delete temp2dhisB3;
		
		}
	}
	f->Write();
	f->Close();
	return 0;
}

int OneEvent (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, MatrixPsd &psdA, MatrixPsd &psdB, Matrix &timingAB)
{
	cout << "function B pick" << endl;
	int xidA, yidA, zidA, xidB, yidB, zidB;
	string MapOption;
	int facepick;
	double fracVeto=0.;
	TString xAid, yAid, zAid, xBid, yBid, zBid,Aname, Bname, fracname;
	Aname="Event-A-";
	Bname="Event-B-";
    cout << "input factor for veto" << endl;
	cin >> fracVeto;
	cout << "input xid for B event" << endl;
	cin >> xidB;
	cout << "input yid for B event" << endl;
	cin >> yidB;
	cout << "input zid for B event" << endl;
	cin >> zidB;
	cout << "Do you need 2D events mapping? (y/n)" << endl;
	cin >> MapOption;
	fracname.Form("%f",fracVeto);
	xAid.Form("%d",xidA);
	xBid.Form("%d",xidB);
	yAid.Form("%d",yidA);
	yBid.Form("%d",yidB);
	zAid.Form("%d",zidA);
	zBid.Form("%d",zidB);
	TString rootid=Bname+xBid+"-"+yBid+"-"+zBid+"-"+fracname;
	/*create root file to hold spectrum*/
	TString rootap=".root";
	TString analysisfile;
	analysisfile.Form ("%s",filename);
	TFile* f=new TFile (analysisfile+"-"+rootid+rootap,"recreate");
	ofstream foutevent;
	foutevent.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".txt");
	ofstream foutdata;
	foutdata.open (xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+".csv");
	foutdata << "Event Num, Event Count, Orange A, Blue A, Green A, Orange B, Blue B, Green B, Orange A PSD, Blue A PSD, Green A PSD, Orange B PSD, Blue B PSD, Green B PSD,timingAB-Orange, timingAB-Blue, timingAB-Green, timingA-Orange. timingA-Blue, timingA-Green, timingB-Orange, timingB-Blue, timingB-Green, \n";
	int eventcount=0;
	/*create category of histogram to hold spectrum*/
	TString energystr="energy ";
	TString peakstr="peak ";
	TString psdstr="psd ";
	TH1D peakhistA[100];
	TH1D integralhistA[100];
	TH1D peakhistB[100];
	TH1D integralhistB[100];
	TH2D psdhistA[100];
	TH2D psdhistB[100];
	TH2D psdhistAB[100];
	TString titleA=" A";
	TString titleB=" B";
	int xid=0;
	int yid=0;
	for (int i=0; i<100; i++)
	{
		TString xidname,yidname;
		xid=i%10;
		yid=i/10;
		xidname.Form ("%d",xid);
		yidname.Form ("%d",yid);
		TString titlename = yidname+" "+xidname;
		TString histname = yidname+" "+xidname;
		peakhistA[i]=TH1D(peakstr+histname+titleA,titlename,2000,0,2000);
		integralhistA[i]=TH1D(energystr+histname+titleA,titlename,6000,0,350000); // needs to be modified the size and maximum
		psdhistA[i]=TH2D(psdstr+histname+titleA,titlename,1500,0,350000,100,0,1);
		peakhistB[i]=TH1D(peakstr+histname+titleB,titlename,2000,0,2000);
		integralhistB[i]=TH1D(energystr+histname+titleB,titlename,6000,0,350000); // needs to be modified the size and maximum
		psdhistB[i]=TH2D(psdstr+histname+titleB,titlename,1500,0,350000,100,0,1);
    	psdhistAB[i]=TH2D(psdstr+histname+titleA+titleB,titlename,1500,0,350000,100,0,1);
	}
	int histTag=0;
	bool vetoA=true;
	bool vetoB=true;
	bool vetoApick=true;
	bool vetoBpick=true;
	for (int i=0; i<n; i++)
	{
		//vetoApick=pick(lowthreshold,highthreshold,energyA,i);
		//vetoBpick=pick(lowthreshold,highthreshold,energyB,i);
		if (true)
		{
			//cout << eventID[i] << "A 3 data"<< endl;
			//vetoA=cubeveto(energyA, xidA,yidA,zidA,i);
			//cout << eventID[i] << "B 3 data"<< endl;
			vetoB= pmtsideveto(energyB,xidB,yidB,fracVeto,0,i);
			//cout << vetoA << "\t" << vetoB << endl;
			if (!vetoB)
			{
				foutevent << eventID[i] << endl;
				eventcount++;		
				foutdata << eventID[i] << "," << i << "," << energyA[yidA+i*10][xidA] << "," <<  energyA[yidA+i*10][zidA+5] << "," << energyA[5+i*10+zidA][xidA] << "," << energyB[yidB+i*10][xidB] << "," <<  energyB[yidB+i*10][zidB+5] << "," << energyB[5+i*10+zidB][xidB] << "," << psdA[yidA+i*10][xidA] << "," <<  psdA[yidA+i*10][zidA+5] << "," << psdA[5+i*10+zidA][xidA] << "," << psdB[yidB+i*10][xidB] << "," <<  psdB[yidB+i*10][zidB+5] << "," << psdB[5+i*10+zidB][xidB]<<timingAB[yidB+i*10][xidB] << "," << timingAB[yidB+i*10][5+zidB] << "," << timingAB[5+zidB+i*10][xidB] << "," << timingA[yidA+i*10][xidA] << "," << timingA[yidA+i*10][5+zidA] << "," << timingA[5+zidA+i*10][xidA] << "," << timingB[yidB+i*10][xidB] << "," << timingB[yidB+i*10][5+zidB] << "," << timingB[5+zidB+i*10][xidB] << "," << "\n";
				for(int j=0;j<10;j++)
				{
					for(int k=0; k<10; k++)
					{
						if (j<5||k<5)
						{
							histTag=j*10+k;								
							peakhistA[histTag].Fill(peakA[j+10*i][k]);
							integralhistA[histTag].Fill(energyA[j+10*i][k]);
							psdhistA[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
							peakhistB[histTag].Fill(peakB[j+10*i][k]);
							integralhistB[histTag].Fill(energyB[j+10*i][k]);
							psdhistB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
							psdhistAB[histTag].Fill(energyB[j+10*i][k],psdB[j+10*i][k]);
							psdhistAB[histTag].Fill(energyA[j+10*i][k],psdA[j+10*i][k]);
						
						}
					}
				}
				if (MapOption=="y")
				{			
					TString eventchar;						
					eventchar.Form ("%d",eventID[i]);	
					TH2D* temp2dhisA=new TH2D("eventA "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5 || k<5)
							{
								for (int l=0;l<energyA[j+i*10][k];l++)
								{
									temp2dhisA->Fill(k,9-j);						
								}
							}
						}
					}
					TH2D* temp2dhisB=new TH2D("eventB "+eventchar, eventchar+"Energy mapping",20,0,10,20,0,10);
					for(int j=0;j<10;j++)
					{
						for(int k=0; k<10; k++)
						{
							if (j<5 || k<5)
							{
								for (int l=0;l<energyB[j+i*10][k];l++)
								{
									temp2dhisB->Fill(k,9-j);						
								}
							}
						}
					}
					//cout << "wrtie the 2D map" << endl;
					temp2dhisA->Write();
					temp2dhisB->Write();
					delete temp2dhisA;
					delete temp2dhisB;
				}
			}
		}
	}
	cout << "total "<< eventcount << "\t events being found" << endl;
	foutdata.clear();
	foutdata.close();
	foutevent.clear();
	foutevent.close();
	f->Write();
	f->Close();
	return 0;
}

int RadonCalibration (char* filename, vector<int> eventID, int n, Matrix &energyA, Matrix &energyB, Matrix &peakA, Matrix &peakB, Matrix &timingA, Matrix &timingB, MatrixPsd &psdA, MatrixPsd &psdB, Matrix &timingAB)
{
	/*create root file to hold spectrum*/
	TString rootap="-calibration.root";
	TString analysisfile;
	analysisfile.Form ("%s",filename);
	TFile* f=new TFile (analysisfile+rootap,"recreate");
	
	/*Output a csv file for graph purpose*/

	ofstream foutdata;
	foutdata.open (analysisfile+"-radon-calibration.csv");
	foutdata << "CubeID, Orange, Blue, Green, EventCount, \n";
	
	
	//ofstream foutdata;
	//foutdata.open (xAid+"-"+yAid+"-"+zAid+"-"+Aname+xAid+"-"+yAid+"-"+zAid+"-"+Bname+xBid+"-"+yBid+"-"+zBid+"-"+fracname+".csv");
	//foutdata << "Event Num, Event Count, Orange A, Blue A, Green A, Orange B, Blue B, Green B, Orange A PSD, Blue A PSD, Green A PSD, Orange B PSD, Blue B PSD, Green B PSD,timingAB-Orange, timingAB-Blue, timingAB-Green, timingA-Orange. timingA-Blue, timingA-Green, timingB-Orange, timingB-Blue, timingB-Green, \n";
	
	
	bool vetoA=true;
	bool vetoB=true;
	bool vetoApick=true;
	bool vetoBpick=true;

	TString faceidOne="Orange";
	TString faceidTwo="Blue";
	TString faceidThree="Green";
	TString TypeEnergy="energy";
	TString TypePeak="peak";

	int energylowbound=0;
	int energyhighbound=350000;
	int bincount=6000;
	
	int binOrange=0;
	int binBlue=0;
	int binGreen=0;

	int binOrangelow, binOrangehigh, binBluelow, binBluehigh,binGreenlow, binGreenhigh;
	
	int eventnumber=0;
	int amin, amax, bmin, bmax, cmin, cmax;
	int cubeID;
	
	
	for (int a=0; a<5; a++)
	{
		for (int b=0; b<5; b++)
		{
			for (int c=0; c<5; c++)
			{
				/*a, b, c represents the cube number increment along x, y and z axis*/
				TString CubeIdName;
				cubeID=a*100+b*10+c;
				CubeIdName.Form ("%d",cubeID);
				TH1D *energyOrange = new TH1D(CubeIdName+faceidOne+TypeEnergy,CubeIdName+faceidOne+TypeEnergy,6000,0,360000);
				TH1D *energyBlue = new TH1D(CubeIdName+faceidTwo+TypeEnergy,CubeIdName+faceidTwo+TypeEnergy,6000,0,360000);
				TH1D *energyGreen = new TH1D(CubeIdName+faceidThree+TypeEnergy,CubeIdName+faceidThree+TypeEnergy,6000,0,360000);
				TH1D *peakOrange = new TH1D(CubeIdName+faceidOne+TypePeak,CubeIdName+faceidOne+TypePeak,2000,0,2000);
				TH1D *peakBlue = new TH1D(CubeIdName+faceidTwo+TypePeak,CubeIdName+faceidTwo+TypePeak,2000,0,2000);
				TH1D *peakGreen = new TH1D(CubeIdName+faceidThree+TypePeak,CubeIdName+faceidThree+TypePeak,2000,0,2000);
				
				for (int i=0; i<n; i++)
				{
					//cout << eventID[i] << "A 3 data"<< endl;					
					amin=lowerCubeBound(a);
					amax=upperCubeBound(a);
					bmin=lowerCubeBound(b);
					bmax=upperCubeBound(b);
					cmin=lowerCubeBound(c);
					cmax=upperCubeBound(c);		
					
					for (int l=amin; l<amax+1; l++)
					{
						for (int m=bmin; m<bmax+1; m++)
						{
							for (int n=cmin; n<cmax+1; n++)
							{
								vetoA=cubeveto(energyA, l, m, n, 1, i);
								if (vetoA==false)
								goto EndOfA;
							}
						}
					}			
					EndOfA :
					//cout << eventID[i] << "B 3 data"<< endl;
					vetoB=cubeveto(energyB, a, b, c, 0.4, i);
					//cout << vetoA << "\t" << vetoB << endl;
					if (!vetoA && !vetoB)
					{
						eventnumber++;
						energyOrange -> Fill(energyB[b+i*10][a]);
						energyBlue -> Fill(energyB[b+i*10][5+c]);
						energyGreen -> Fill(energyB[5+c+i*10][a]);
						peakOrange -> Fill(peakB[b+i*10][a]);
						peakBlue -> Fill(peakB[b+i*10][5+c]);
						peakGreen -> Fill(peakB[5+c+i*10][a]);
					}	
				}

				binOrange=energyOrange->GetMaximumBin() ;
				binBlue=energyBlue->GetMaximumBin();
				binGreen=energyGreen->GetMaximumBin();

				binOrangelow=lowerBinBound(binOrange, 50);
				binOrangehigh=upperBinBound(binOrange, 50);
				binBluelow=lowerBinBound(binBlue, 50);
				binBluehigh=upperBinBound (binBlue, 50);
				binGreenlow=lowerBinBound (binGreen, 50);
				binGreenhigh=upperBinBound (binGreen, 50);

				TF1 *orangefit = new TF1("orangefit","gaus",binOrangelow*60,binOrangehigh*60);
				TF1 *bluefit = new TF1("bluefit","gaus",binBluelow*60,binBluehigh*60);
				TF1 *greenfit = new TF1("greenfit","gaus",binGreenlow*60,binGreenhigh*60);
				
				energyOrange->Fit("orangefit","R && Q");
				energyBlue->Fit("bluefit","R && Q");
				energyGreen->Fit("greenfit","R && Q");
							
				cout << "cube" << a<<b<<c<<" are analyzed after this cycle." << " " << eventnumber << " event are found."<<endl;
				cout << "Orange=" << orangefit->GetParameter(1) << endl;
				cout << "Blue=" << bluefit->GetParameter(1) << endl;
				cout << "Green=" << greenfit->GetParameter(1) << endl;

				//foutdata << "CubeID, Orange, Blue, Green, EventCount, \n";
				foutdata << cubeID << "," << orangefit->GetParameter(1) << "," << bluefit->GetParameter(1) << "," <<  greenfit->GetParameter(1) << "," << eventnumber << "," << "\n";
				
				eventnumber=0;				
				energyOrange->Write();
				energyBlue->Write();
				energyGreen->Write();
				peakOrange->Write();
				peakBlue->Write();
				peakGreen->Write();
				delete energyOrange;
				delete energyBlue;
				delete energyGreen;
				delete peakOrange;
				delete peakBlue;
				delete peakGreen;
				delete orangefit;
				delete bluefit;
				delete greenfit;				
			}
		}
	}
	//cout << "total "<< eventcount << "\t events being found" << endl;
	//foutdata.clear();
	//foutdata.close();
	//foutevent.clear();
	//foutevent.close();
	foutdata.clear();
	foutdata.close();
	f->Write();
	f->Close();
    return 0;
}


#endif
