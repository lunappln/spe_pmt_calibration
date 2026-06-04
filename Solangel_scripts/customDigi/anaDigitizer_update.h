#include "TF1.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include <fstream>
#include "TCanvas.h"

#include <iostream>
#include <cstdio>


// *************************************************
//        Variables and structure declarations
// *************************************************
Double_t ampRes = 0.244140625;// 1/4096;// amplitude resolution = 0.2 mV per bit (LSB) ;
// for sipm we run at 1GHz -> t=1/1e9 = 1ns
Double_t timeRes = 1.;// 1/5e9 samples/second -> 0.2 nsecons per bit ;


// histograms

TH2 * h2SignalShift = new TH2D("h2SignalShift", "Signal Average shifted; time (ns);Amplitude (mV);",1024*3,-timeRes*1024,timeRes*1023*2 ,4096*0.24414,-4095*0.24414,4096*0.24414);
TH2 * h2Signal = new TH2D("h2Signal", "Signal Average; time (ns);Amplitude (mV);", 1024,0,timeRes*1023 ,4096*2,-4096*0.24414,4095*0.24414);

TH1 * h1BaseLineAll = new TH1D("h1BaseLineAll", "Base Line from all events;Base line level (mV);Entries",1024,0,1024);
TH1 * h1RMSAll = new TH1D("h1RMSAll", "RMS from all events;RMS (mV);Entries",1000,0,100);
TH1 * h1charge = new TH1D("h1charge", ";Charge (fC);Entries",2000,0,1000);
TH1 * h1ampSig = new TH1D("h1ampSig", "Signal;Amplitude (mV);Entries",8000,-1000,1000);
TH1 * h1ampTr1 = new TH1D("h1ampTr1", "Trigger;Amplitude (mV);Entries",2000,-1000,1000);
TH1 * h1ThrsCrosses = new TH1D("h1ThrsCrosses", ";Number of crosses;Entries",100,0,100);
TH1 * h1Leading10 = new TH1D("h1Leading10", ";Leading time index;Entries",1024,0,1024);
TH1 * h1timeWidth = new TH1D("h1timeWidth", "; time width (ns);Entries",1024,0,1024*timeRes);

TH1 * h1thrsTime = new TH1D("h1thrsTime", "Threshold Time;Time (ns);Entries",2000,-1000,1000);

struct h1Type {
    Double_t mean;
    Double_t rms;
};

// Signal time properties
struct TimeProps {
    Int_t lead10;
    Int_t lead90;
    Int_t trail10;
    Int_t trail90;
};


// Structure to extract the information of each event
struct grType {
  UInt_t n;
  Double_t * y;// = gr->GetY();
  Double_t * x;// = gr->GetY();
  UInt_t locmax;// = TMath::LocMax(n,y);
  UInt_t locmin;// = TMath::LocMax(n,y);
  Double_t vmax;// = y[locmax];
  Double_t vmin;// = y[locmax];
  Double_t tmin;// = x[locmax];
  Double_t tmax;// = x[locmax];
  Double_t tshift;// = y[locmax];
  UInt_t loctshift;
  Bool_t UnderThr = kFALSE;
  Double_t zcross; // Zero crossing
};


// *************************************************
//        Function prototype
// *************************************************
std::ifstream::pos_type filesize(const char* filename);
h1Type GetBaseLine(float * tempBuffer);
grType GetGrProp(TGraph * gr, Double_t threshold, Int_t PulseSign , Double_t WStart, Double_t WEnd );


// *************************************************
//        Functions
// *************************************************
std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}


h1Type GetBaseLine(float * tempBuffer, const UInt_t baseLineSample=100){
  TH1D * h1BaseLine = new TH1D("h1BaseLine", "Base Line;Base line level (mV);Entries",1024,0,1024);
  for (UInt_t i = 0; i < baseLineSample; i++){ //Calculate base line
      h1BaseLine -> Fill(tempBuffer[i]*ampRes);
  }

  TGraph * gr = new TGraph();
  for (UInt_t i = 0; i < baseLineSample; i++){
    gr ->SetPoint(i, i, tempBuffer[i]*ampRes);
  }

  h1Type h1Result;
  h1Result.mean = h1BaseLine -> GetMean();
  h1Result.rms = h1BaseLine -> GetRMS();
  delete h1BaseLine;
  return h1Result;
}

grType GetGrProp(TGraph * gr, Double_t threshold = 0, Int_t PulseSign = -1, Double_t WStart =0, Double_t WEnd = 204.8){ // Graph properties
  grType graph;

  // Assign vertical (x) and horizontal (y) of the graph to pointers X and Y
  graph.n = gr ->GetN();
  graph.y = gr ->GetY();
  graph.x = gr ->GetX();

  // Locate array points for the maximum and minimum values
  graph.locmax = TMath::LocMax(graph.n,graph.y);
  graph.locmin = TMath::LocMin(graph.n,graph.y);

  // Get value of the minimum and maximum voltages and time(corresponding to the max and min voltagen)
  graph.vmax = graph.y[graph.locmax];
  graph.vmin = graph.y[graph.locmin];
  graph.tmin = graph.x[graph.locmin];
  graph.tmax = graph.x[graph.locmax];

  // Finding the location of a point over specific threshold
  for(UInt_t i = 0; i < graph.n; i++){
    if( (graph.x[i] > WStart) && (graph.x[i] < WEnd) ){
      // The fit is automatized based on a treshold and a range to analize the signal
      // Fit a pol1 funtion to get the zero crossing of a positive signals
      if(graph.y[i] >= threshold && PulseSign == 1){
        // Fit a pol1 function to the signals
        TF1 *f = new TF1("f", "[0]+[1]*x", graph.x[i], graph.x[i]+3);
        gr->Fit("f","RQ0");

        graph.zcross = -f->GetParameter(0)/f->GetParameter(1);
        // cout<< "Trigger   zero crossing time = "<< graph.zcross << "\n";
        break;
      }
      //if(graph.y[i] <= threshold && PulseSign == -1){
      if(PulseSign == -1){
        // Fit a pol1 function to the signals
        // TF1 *f = new TF1("f", "[0]+[1]*x", graph.x[i], graph.x[i]+3);
        // gr->Fit("f","R");
        //cout<< "Tmin "<< graph.tmin<< endl;
        TF1 *f = new TF1("f", "[0]+[1]*x", graph.tmin-8,graph.tmin);
        gr->Fit("f","RQ0");
        graph.zcross = -f->GetParameter(0)/f->GetParameter(1);
        // cout<< "Trigger   zero crossing time = "<< graph.zcross << "\n";
        break;
      }

      if(graph.y[i] < threshold && PulseSign == -1){
        graph.tshift = graph.x[i];
        graph.loctshift = i;
        //printf("time %f, %f, %f\r", graph.x[i], graph.y[i], PulseSign * threshold);
        break;
      }

      if(graph.y[i] > threshold && PulseSign == 1){
        graph.tshift = graph.x[i];
        graph.loctshift = i;
        //printf("time %f, %f, %f \r", graph.x[i], graph.y[i], PulseSign * threshold);
        break;
      }
    }
    if(graph.x[i] > WStart && graph.y[i] < (PulseSign * threshold)) graph.UnderThr=kTRUE;

    // TODO handle signals below the threshold
  }
  return graph;
}


double GetCharge(TGraph * gr, double  WStart=0 ,   double  WEnd = 204.8, int sign=1)
{
  // Integrate in the interval
  Double_t charge=0;
  Double_t amplitude=0;
  for(int i = 0; i < gr ->GetN(); i++){
    amplitude = sign*gr->GetPointY(i);
    if( (gr->GetPointX(i)>WStart) && (gr->GetPointX(i)<WEnd) ){
      charge += amplitude*(timeRes/50);
    }
  }
  //cout<<"Charge "<< charge<< endl;
  return charge;
}

double GetVmax(TGraph * gr, double  WStart=0 ,   double  WEnd = 204.8, int sign=1)
{
  double vMax=0;
  return vMax;
}

int GetThrsTime(TGraph * gr, Double_t threshold, int position=1){
	Double_t value;
	Int_t thrsTime=-1000;

	//if(gr->GetPointY(0)>threshold) position=-1;

	for(int i = 1; i < gr ->GetN(); i++){
		value = gr->GetPointY(i);
		if (value>threshold && position>0){
			thrsTime=i;
			break;
		}
		else if (value<threshold && position<0){
			thrsTime=i;
			break;
		}
	}
	// cout<< "Time threshold " << thrsTime << ""<< endl;
	return thrsTime;
}

int ThrsCross(TGraph * gr, Double_t threshold){
	Double_t value;
	int countNeg=0, countPos=0, above=0, minTime=10, position=1;

	if(gr->GetPointY(0)>threshold) position=-1;

	for(int i = 1; i < gr ->GetN(); i++){
		value = gr->GetPointY(i);
		if (value>threshold){
      if(above<0) above=0;
			++above;
      if(above==minTime) ++countPos;
		}
		else if (value<threshold){
			if(above>0) above=0;
      --above;
      if(above==-minTime) ++countNeg;
		}
	}

	if(position>0){
    //cout << "Threshold:" << countPos << endl;
    return countPos;
  }
    //cout << "Threshold:" << countPos << endl;
  return countNeg;
}


TimeProps GetTimes(TGraph * gr, Double_t threshold, int timeIndex){
	Double_t value;
	int position=1, maxVal=0;
  TimeProps outResults;

	if(gr->GetPointY(0)>threshold) position=-1;
  for(int i = 1; i < gr ->GetN(); i++){
		value = gr->GetPointY(i);
    value = (value-threshold)*position;
    if(value>maxVal) maxVal=value;
  }

  int i = timeIndex, first=-1, second=-1, third=-1, fourth=-1;
	for(; i < gr ->GetN(); ++i){
		value = gr->GetPointY(i);
    value = (value-threshold)*position;
		if(value<0) break;
    if(value>0.1*maxVal && first==-1) first=i;
    if(value>0.9*maxVal && second==-1) second=i;
	}
  for(--i; i > timeIndex; --i){
		value = gr->GetPointY(i);
    value = (value-threshold)*position;
		if(value<0) break;
    if(value>0.1*maxVal && fourth==-1) fourth=i;
    if(value>0.9*maxVal && third==-1) third=i;
	}
  outResults.lead10 = first;
  outResults.lead90 = second;
  outResults.trail10 = fourth;
  outResults.trail90 = third;

  cout<< "\n\r Results " <<  first << "\t" << second << "\t" << third << "\t" << fourth << endl;
	return outResults;
}
