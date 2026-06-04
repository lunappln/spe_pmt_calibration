#include "anaDigitizer_update.h"
// Program to calculate amplitude and charge from the

// TODO
// * Move the edge FIT of the trigger signal to adjust properly the time alignment

void anaDigitizer_charge_v1(Int_t prcEvents = -1){

  const UInt_t eventSz=1024; // Size of each event

  TimeProps signalTimes; // Time properties structure

  //****************************************************************
  // Signal to analize inicialization
  //****************************************************************
  float  buffer[eventSz]; //Buffer to store the event
  FILE *ptr; //Pointer to the file to process
  const Char_t * file = "/Users/solangelmac/cernbox/Students/Timea/sipm/data/20230424_S1336-30xx_test0/wave_6.dat"; // File to analize

  // Parameters to analize the signal
  Double_t Wstart =200; // Start of the time expected of the signal. Begining of the time window to analize the signal.
  Double_t Wend = 600; // End of the time expected of the signal. End of the time window to analize the signal.
  Double_t thrs = 40; // Leading edge threshold
  Int_t sign = 1 ;// Sign fo the signal to analize

  // Number of entries for the analized file
  Int_t dataPointSz = sizeof buffer[0];
  cout <<  "Size of data point: " << sizeof buffer[0] <<" Bytes" << endl;
  Int_t fSize = filesize(file);
  Int_t nEvents = fSize/(eventSz*dataPointSz);
  //nEvents = 50;
  printf("File size %i Bytes\nEvent size (%i)\nTotal Events(%i): %s\r", fSize, eventSz*dataPointSz, nEvents, file);
  ptr = fopen(file,"rb");  // r for read, b for binary

  if(prcEvents>0)  nEvents = prcEvents;
  else nEvents = fSize/(eventSz*dataPointSz);


  //****************************************************************
  //****************************************************************

  //Fill in nseconds the time buffer for all the signals
  Double_t dt[eventSz];
  for (UInt_t i = 0; i < eventSz; i++) {dt[i]=i*timeRes;}
  Double_t amplitude;

  float RawAmp;

  // Plots for the signals
  TGraph* grAv  = nullptr;// For the analized signal
  TGraph* grRaw = new TGraph(eventSz);

  // create Root tree file
   gROOT->cd(); //make sure that the Tree is memory resident
   TTree *TResults = new TTree("T","test circular buffers");
   Double_t charge;
   TResults->Branch("charge",&charge,"charge/D");

  // Filling histograms and plots per each signal
  for(UInt_t event=0; event < nEvents; event++){

    cout << "Event number: " << 100*event/(nEvents+1) << "\% " << event+1 <<"/"<<nEvents << "\r";

    // Plots for the signals
    grAv = new TGraph(eventSz);
    grAv -> GetYaxis() -> SetRange(-1000,1000);
    
    // Put one signal (1024*UInt32) in the buffer
    fread(buffer,sizeof(buffer),1,ptr);


    // Producing the graphs

    // Graph from the signal input
    h1Type h1Pars = GetBaseLine(buffer,100); // Get baseline level

    for (UInt_t i = 0; i < eventSz; i++){
      RawAmp = buffer[i]; // Amplitude in ADC values
      grRaw ->SetPoint(i,dt[i],RawAmp); // Fill graph with the ADC amplitude values
      grAv ->SetPoint(i,dt[i],RawAmp*ampRes);  // Fill graph amplitude in mV values
      // All signals overlaped
      h2Signal -> Fill(dt[i],RawAmp*ampRes);
    }

    // Get signal properties
    grType sgProp = GetGrProp(grAv, thrs, sign, Wstart, Wend);
    if(sign<0) h1ampSig->Fill(sign*sgProp.vmin);
    else h1ampSig->Fill(sign*sgProp.vmax);

    // Fill histogram with the charge (no baseline subtraction)
    charge = GetCharge(grAv,Wstart, Wend, sign);
    h1charge->Fill(charge);
   
    // Baseline quality control histograms
    h1BaseLineAll->Fill(h1Pars.mean);
    h1RMSAll->Fill(h1Pars.rms);


    //keep the last graph for debuging
    if(event < nEvents-1){
      delete grAv; // Delete current graph to avoid conflict in the next iteration
      //delete grAvTrg0;
      }
    else{
      //Filling the amplitude raw values to debug and quality control
      for (UInt_t i = 0; i < eventSz; i++){
        RawAmp = buffer[i];
        grRaw ->SetPoint(i,dt[i],RawAmp);
        grAv ->SetPoint(i,dt[i],RawAmp*ampRes);
        //h2Signal -> Fill(dt[i],amplitude);
      }
    }// Keep last hsto for debuging END

     TResults->Fill(); // Store the data to the root tree file

  } // Fill histograms and plots

  TResults -> Print();
  // save Tree file
  // NOTE: to open in the classic root browser -> "root --web=off"
  TFile fout("output.root","recreate");
  fout.cd();
  TResults->Write();
  fout.Close();

  TCanvas * c1 = new TCanvas("c1","c1 title",1600,800);
  c1->Divide(3,1);
  c1->cd(1);
  grAv ->SetLineColor(kBlack);
  grAv ->SetLineWidth(1);
  grAv ->Draw("ALP*");

  //***************************************************************************************************
  //***************************************************************************************************
  // CODE FOR DEBUGGING PRUPOSSES
    h1Type h1Pars = GetBaseLine(buffer,1000);
    for (UInt_t i = 0; i < eventSz; i++){
      amplitude = buffer[i]*ampRes - h1Pars.mean;
      // Fill histo and plot with time shifted
      grAv ->SetPoint(i,dt[i],amplitude);
      //h2Signal -> Fill(dt[i],amplitude);
    }//
  grType sgProp = GetGrProp(grAv, thrs, sign, Wstart, Wend);  // Graph from the signal input
  //***************************************************************************************************
  //***************************************************************************************************



  TCanvas * cRaw = new TCanvas("cRaw","cRaw title",1400,800);
  grRaw ->SetLineColor(kBlack);
  grRaw ->SetLineWidth(1);
  grRaw ->Draw("ALP*");

  TCanvas * c2 = new TCanvas("c2","c2 title",1400,800);
  h2Signal -> Draw("COLZ");

  TCanvas * c3 = new TCanvas("c3","Base line for all events",1400,800);
  c3 -> Divide (2,1);
  c3 -> cd(1) -> SetLogy();
  h1BaseLineAll -> Draw();
  c3 -> cd(2) -> SetLogy();
  h1RMSAll -> Draw();

  TCanvas * c4 = new TCanvas("c4","charge",800,800);
  c4-> SetLogy();
  h1charge->Draw();

  TCanvas * c5 = new TCanvas("c5","Amplitude",800,800);
  h1ampSig->Draw();



}
