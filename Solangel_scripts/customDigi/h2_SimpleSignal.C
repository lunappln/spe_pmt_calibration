// The program extract the signal form binary file and average them
#include <fstream>
#include <string>
#include <iostream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"
#include "TH2.h"

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

void h2_SimpleSignal(Int_t nEvents = -1, const Char_t * file = "/Users/solangel/cernbox/muonID/beamTest/data/6GeVHadrons11062023_35/wave_26.dat"){
  const UInt_t eventSz=1024;
  float  buffer[eventSz];
  FILE *ptr;

	float tres=1.; //sampling rate = 1GS/s

  // Get the size of the file
  // TODO set option for full or par of the file to be analized
  Int_t dataPointSz=sizeof buffer[0];
  cout <<  "Size of data point: " << sizeof buffer[0] << endl;
  int fSize = filesize(file);
  printf("File size %i Bytes\tEvent size (%i)\tTotal Events(%i)\n", fSize, eventSz*dataPointSz, fSize/eventSz);

  ptr = fopen(file,"rb");  // r for read, b for binary
  //cout <<  "Size of file: " << sizeof ptr << endl;

  Double_t dt[eventSz];
  Double_t amplitude;
  Double_t ampRes = 0.24414;// 1/4098;// amplitude resolution = 0.2 mV per bit (LSB) ;

  // TH2 * h2Signal = new TH2D("h2Signal", "Signal Average; time (ns);Amplitude (mV);",512,0,250, 2048,-500,100);
  TH2 * h2Signal = new TH2D("h2Signal", "Signal Average; time (ns);Amplitude (mV);", 1024,0,tres*1023 ,4096*2,-4096*0.24414,4095*0.24414);

  for (Int_t i = 0; i < eventSz; i++){
    dt[i]=i*tres;
  }

  Int_t nTotalEvents = fSize/(eventSz*dataPointSz);
   if (nEvents<0) nEvents = nTotalEvents;

  for(Int_t event=0; event < nEvents; event++){
    fread(buffer,sizeof(buffer),1,ptr); // Put one signal (1024*UInt32) in the buffer
    for (Int_t i = 0; i < eventSz; i++){
      h2Signal -> Fill(dt[i],buffer[i]*ampRes - 165);
    }
    printf("Event number %i\r",nTotalEvents);
  }
  cout << "\n\n";
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPalette(kSolar);

  TCanvas * c2 = new TCanvas("c2","c2 title",1400,800);
  h2Signal -> Draw("COLZ");
  c2 -> SetLogz();
}
