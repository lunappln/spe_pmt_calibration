#include "anaDigitizer_update.h"
// Program to calculate amplitude and charge from the

// TODO
// * Move the edge FIT of the trigger signal to adjust properly the time alignment

void anaDigitizer_singlePhoton_v1(Int_t prcEvents = -1){

  const UInt_t eventSz=1024; // Size of each event

  TimeProps signalTimes;

  //****************************************************************
  // Signal to analize inicialization
  //****************************************************************
  float  buffer[eventSz]; //Buffer to store the event
  FILE *ptr; //Pointer to the file to process
  const Char_t * file = "/Users/solangelmac/cernbox/Students/Vojta/Qdots/run_20230214/wave_1.dat"; // File to analize

  // Parameters to analize the signal
  Double_t Wstart =60; // Start of the time expected of the signal. Begining of the time window to analize the signal.
  Double_t Wend = 85; // End of the time expected of the signal. End of the time window to analize the signal.
  Double_t thrs = -15; // Leading edge threshold
  Int_t sign = -1 ;// Sign fo the signal to analize

  // Number of entries for the analized file
  Int_t dataPointSz = sizeof buffer[0];
  cout <<  "Size of data point: " << sizeof buffer[0] <<" Bytes" << endl;
  Int_t fSize = filesize(file);
  Int_t nEvents = fSize/(eventSz*dataPointSz);
  //nEvents = 50;
  printf("File size %i Bytes\nEvent size (%i)\nTotal Events(%i): %s\r", fSize, eventSz*dataPointSz, nEvents, file);
  ptr = fopen(file,"rb");  // r for read, b for binary


  //****************************************************************
  // TRIGGER 0 signal inicialization
  //****************************************************************
  float  bufferTrg0[eventSz]; //Buffer to store the event
  FILE *ptrTrg0; //Pointer to the file to process
  const Char_t * fileTrg0 = "/Users/solangelmac/cernbox/Students/Vojta/Qdots/run_20230214/TR_0_0.dat"; // Trigger: time reference signal
  // Parameters to analize the signal
  Double_t WstartTrg0 = 85; // Start of the time expected of the signal. Begining of the time window to analize the signal.
  Double_t WendTrg0 = 105; // End of the time expected of the signal. End of the time window to analize the signal.
  Double_t thrsTrg0 = -65; // Leading edge threshold
  Int_t signTrg0 = -1 ;// Sign fo the signal to analize

  // Number of entries for the time reference file
  dataPointSz = sizeof bufferTrg0[0];
  cout <<  "Size of data point (Trg0): " << sizeof bufferTrg0[0] <<" Bytes" << endl;
  Int_t fSizeTrg0 = filesize(fileTrg0);
  Int_t nEventsTrg0 = fSizeTrg0/(eventSz*dataPointSz);
  nEvents = nEventsTrg0;
  if(prcEvents>0)  nEvents = prcEvents;
  else nEvents = nEventsTrg0;

  printf("File size (Trg0) %i Bytes\nEvent size (%i)\nTotal Events(%i): %s\n\r", fSizeTrg0, eventSz*dataPointSz, nEventsTrg0, fileTrg0);
  ptrTrg0 = fopen(fileTrg0,"rb");  // r for read, b for binary


  //****************************************************************
  //****************************************************************

  //Fill in nseconds the time buffer for all the signals
  Double_t dt[eventSz];
  for (UInt_t i = 0; i < eventSz; i++) {dt[i]=i*timeRes;}
  Double_t amplitude;
  Double_t amplitudeTrigger;

  float RawAmp;

  // Plots for the signals
  TGraph* grAv  = nullptr;// For the analized signal
  TGraph* grRaw = new TGraph(eventSz);
  TGraph* grAvTrg0 = nullptr;// For the trigger 0 signal
  TGraph* grRawTrg0 = new TGraph(eventSz);

  // create tree file
   gROOT->cd(); //make sure that the Tree is memory resident
   TTree *TResults = new TTree("T","test circular buffers");
   Double_t charge;
   TResults->Branch("charge",&charge,"charge/D");

  // Filling histograms and plots
  for(UInt_t event=0; event < nEvents; event++){

    cout << "Event number: " << 100*event/(nEvents+1) << "\% " << event+1 <<"/"<<nEvents << "\r";

    // Plots for the signals
    grAv = new TGraph(eventSz);
    grAv -> GetYaxis() -> SetRange(-1000,1000);
    grAvTrg0 = new TGraph(eventSz);
    grAvTrg0 -> GetYaxis() -> SetRange(-1000,1000);

    // Put one signal (1024*UInt32) in the buffer
    fread(buffer,sizeof(buffer),1,ptr);
    fread(bufferTrg0,sizeof(bufferTrg0),1,ptrTrg0);

    // Producing the graphs
    // The TRIGGER 0 first to calculate the time shift from the other signal
    h1Type h1ParsTrg0 = GetBaseLine(bufferTrg0,100);
    for (UInt_t i = 0; i < eventSz; i++){
      grAvTrg0 ->SetPoint(i,dt[i],bufferTrg0[i]*ampRes-h1ParsTrg0.mean);
    }
    grType sgPropTrg0 = GetGrProp(grAvTrg0, thrsTrg0, signTrg0, WstartTrg0, WendTrg0);


    // Graph from the signal input
    h1Type h1Pars = GetBaseLine(buffer,300);
    //double ldnTRGTime = GetThrsTime(grAvTrg0,-100,-1); // leading edge time

    for (UInt_t i = 0; i < eventSz; i++){
      amplitude = buffer[i]*ampRes - h1Pars.mean;
      amplitudeTrigger = bufferTrg0[i]*ampRes-h1ParsTrg0.mean;
      
      // Fill histo and plot with time shifted
      grAv ->SetPoint(i,dt[i]-sgPropTrg0.zcross,amplitude);
      h2SignalShift -> Fill(dt[i]-sgPropTrg0.zcross,amplitude);
      h2TriggerShift -> Fill(dt[i]-sgPropTrg0.zcross,amplitudeTrigger);
     // h2TriggerShift -> Fill(dt[i]-ldnTRGTime,amplitudeTrigger);
      
      // Not shifted time histogram
      h2Signal -> Fill(dt[i],amplitude);
      h2Trigger -> Fill(dt[i],amplitudeTrigger);
    }//

    // Get signal properties
    grType sgProp = GetGrProp(grAv, thrs, sign, Wstart, Wend);
    if(sign<0) h1ampSig->Fill(sign*sgProp.vmin);
    else h1ampSig->Fill(sign*sgProp.vmax);
    charge = GetCharge(grAv,Wstart, Wend, sign);
    h1charge->Fill(charge);
    TResults->Fill();
    // Quality control histograms
    h1BaseLineAll->Fill(h1Pars.mean);
    h1RMSAll->Fill(h1Pars.rms);

    // Timeas' functions
    h1ThrsCrosses -> Fill(ThrsCross(grAv, thrs));
    signalTimes = GetTimes(grAv, thrs, GetThrsTime(grAv, thrs, sign));
    h1Leading10 -> Fill(signalTimes.lead10);
    h1timeWidth -> Fill( (signalTimes.trail10-signalTimes.lead10)*timeRes);

    h1thrsTime->Fill(GetThrsTime(grAv, thrs, sign));

    //keep the last graph for debuging
    if(event < nEvents-1){
      delete grAv;
      delete grAvTrg0;
      }
    else{
      //Filling the amplitude raw values to debug and quality control
      for (UInt_t i = 0; i < eventSz; i++){
        RawAmp = buffer[i];//&0x0000fff;
        //printf("point %d, buffer %d\n", i,buffer[i]);
        grRaw ->SetPoint(i,dt[i],RawAmp);
        grAv ->SetPoint(i,dt[i],RawAmp*ampRes- h1Pars.mean);
        grRawTrg0 ->SetPoint(i,dt[i],bufferTrg0[i]);
        grAvTrg0 ->SetPoint(i,dt[i],bufferTrg0[i]*ampRes-h1ParsTrg0.mean);
      }
    }// Keep last hsto for debuging END

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
  grType sgPropTrg0 = GetGrProp(grAvTrg0, thrsTrg0, signTrg0, WstartTrg0, WendTrg0);
  //grType sgPropTrg1 = GetGrProp(grAvTrg1, thrsTrg1, signTrg1, WstartTrg1, WendTrg1);
    h1Type h1Pars = GetBaseLine(buffer,100);
    for (UInt_t i = 0; i < eventSz; i++){
      amplitude = buffer[i]*ampRes - h1Pars.mean;
      // Fill histo and plot with time shifted
      grAv ->SetPoint(i,dt[i]-sgPropTrg0.zcross,amplitude);
      h2SignalShift -> Fill(dt[i]-sgPropTrg0.zcross,amplitude);
      // Not shifted time histogram
      h2Signal -> Fill(dt[i],amplitude);
    }//
  grType sgProp = GetGrProp(grAv, thrs, sign, Wstart, Wend);  // Graph from the signal input
  //***************************************************************************************************
  //***************************************************************************************************

  c1->cd(2);
  grAvTrg0 ->SetLineColor(kBlue);
  grAvTrg0 ->SetLineWidth(1);
  grAvTrg0 -> SetMarkerColor(kBlue);
  grAvTrg0 ->Draw("ALP*");
  c1->cd(3);

  TCanvas * cRaw = new TCanvas("cRaw","cRaw title",1400,800);
  grRaw ->SetLineColor(kBlack);
  grRaw ->SetLineWidth(1);
  grRaw ->Draw("ALP*");

  TCanvas * c2 = new TCanvas("c2","c2 title",1400,800);
  c2 -> Divide(2,1);
  c2 -> cd(1) -> SetLogz();
  h2Signal -> Draw("COLZ");
  c2 -> cd(2) ->SetLogz();
  h2SignalShift -> Draw("COLZ");

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
  //c5 -> Divide (2,1);
  //c5 -> cd(1) -> SetLogy();
  //c5 -> cd(2) -> SetLogy();
  h1ampSig->Draw();

  TCanvas * c6 = new TCanvas("c6","thrs Time",800,800);
  c6-> SetLogy();
  h1thrsTime->Draw();
  cout<<endl;

  TCanvas * c7 = new TCanvas("c7","thrs Time cross",800,800);
  h1ThrsCrosses->Draw();

  TCanvas * c8 = new TCanvas("c8","Leading edge at 10\%",800,800);
  h1Leading10->Draw();

  TCanvas * c9 = new TCanvas("c9","Time width",800,800);
  h1timeWidth->Draw();


  TCanvas * c10 = new TCanvas("c10","trigger shift",1400,800);
  c10->Divide(2,1);
  c10 -> cd(1) ->SetLogz();
  h2TriggerShift -> Draw("COLZ");
  c10 -> cd(2) ->SetLogz();
  h2Trigger-> Draw("COLZ");

}
