#include "anaDigitizer_20230611.h"
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TTree.h>
#include <algorithm>

// Program to calculate amplitude and charge from the waveform data

// TODO
// * Move the edge FIT of the trigger signal to adjust properly the time alignment

void anaDigitizer_charge_v2(Int_t prcEvents = -1) {
    const UInt_t eventSz = 1024; // Size of each event
    //const Double_t ampRes = 0.1; // Amplitude resolution in mV
    //const Double_t timeRes = 0.1; // Time resolution in ns

    TimeProps signalTimes; // Time properties structure

    //****************************************************************
    // Signal to analyze initialization
    //****************************************************************
    float buffer[eventSz]; // Buffer to store the event
    FILE *ptr; // Pointer to the file to process
    const Char_t *file = "/Users/solangel/cernbox/muonID/beamTest/data/6GeVHadrons11062023_03/wave_0.dat"; // File to analyze
    const Char_t *OutFileName = "/Users/solangel/cernbox/muonID/beamTest/data/6GeVHadrons11062023_03/pedestal.root"; // Name of the output root file

    // Parameters to analyze the signal
    Double_t Wstart = 10; // Start of the time expected of the signal. Beginning of the time window to analyze the signal.
    Double_t Wend = 210; // End of the time expected of the signal. End of the time window to analyze the signal.
    Double_t thrs = 40; // Leading edge threshold
    Int_t sign = 1; // Sign of the signal to analyze

    // Number of entries for the analyzed file
    Int_t dataPointSz = sizeof(buffer[0]);
    cout << "Size of data point: " << sizeof(buffer[0]) << " Bytes" << endl;
    Int_t fSize = filesize(file);
    Int_t nEvents= fSize/(eventSz*dataPointSz);
    //nEvents = 50;
    printf("File size %i Bytes\nEvent size (%i)\nTotal Events(%i): %s\r", fSize, eventSz*dataPointSz, nEvents, file);
    ptr = fopen(file,"rb");  // r for read, b for binary
    
    if (prcEvents > 0)
        nEvents = prcEvents;
    else
        nEvents = fSize / (eventSz * dataPointSz);

    //****************************************************************
    //****************************************************************

    // Fill in nseconds the time buffer for all the signals
    Double_t dt[eventSz];
    for (UInt_t i = 0; i < eventSz; i++) {
        dt[i] = i * timeRes;
    }
    Double_t amplitude;

    float RawAmp;

    // Plots for the signals
    TGraph *grAv = nullptr; // For the analyzed signal
    TGraph *grRaw = new TGraph(eventSz);

    // Create Root tree file
    gROOT->cd(); // Make sure that the Tree is memory resident
    TTree *TResults = new TTree("T", "test circular buffers");
    Double_t charge;
    Double_t MaxAmp;
    TResults->Branch("charge", &charge, "charge/D");
    TResults->Branch("MaxAmp", &MaxAmp, "MaxAmp/D");

    // Filling histograms and plots per each signal
    for (UInt_t event = 0; event < nEvents; event++) {

       cout << "Event number: " << 100 * event / (nEvents + 1) << "\% " << event + 1 << "/" << nEvents << "\r";

        // Plots for the signals
        grAv = new TGraph(eventSz);
        grAv->GetYaxis()->SetRange(-1000, 1000);

        // Put one signal (1024*UInt32) in the buffer
        fread(buffer, sizeof(buffer), 1, ptr);

        // Producing the graphs

        // Graph from the signal input
        h1Type h1Pars = GetBaseLine(buffer, 100); // Get baseline level

        for (UInt_t i = 0; i < eventSz; i++) {
            RawAmp = buffer[i]; // Amplitude in ADC values
            grRaw->SetPoint(i, dt[i], RawAmp); // Fill graph with the ADC amplitude values
            grAv->SetPoint(i, dt[i], RawAmp * ampRes); // Fill graph amplitude in mV values
            // All signals overlapped
            h2Signal->Fill(dt[i], amplitude);
        }

        // Get signal properties
        grType sgProp = GetGrProp(grAv, thrs, sign, Wstart, Wend);
        if (sign < 0)
            h1ampSig->Fill(sign * sgProp.vmin);
        else
            h1ampSig->Fill(sign * sgProp.vmax);

        // Fill histogram with the charge (no baseline subtraction)
        charge = GetCharge(grAv, Wstart, Wend, sign);
        h1charge->Fill(charge);


        // Baseline quality control histograms
        h1BaseLineAll->Fill(h1Pars.mean);
        h1RMSAll->Fill(h1Pars.rms);

        // Get max amplitude
        MaxAmp = GetVmax(grAv, Wstart, Wend);
        h1MaxAmp->Fill(MaxAmp );

        // Keep the last graph for debugging
        if (event < nEvents - 1) {
            delete grAv; // Delete current graph to avoid conflict in the next iteration
        } else {
            // Filling the amplitude raw values to debug and quality control
            for (UInt_t i = 0; i < eventSz; i++) {
                RawAmp = buffer[i];
                grRaw->SetPoint(i, dt[i], RawAmp);
                grAv->SetPoint(i, dt[i], RawAmp * ampRes);
                h2Signal->Fill(dt[i], amplitude);
            }
        }

        TResults->Fill(); // Store the data to the root tree file

    } // Fill histograms and plots

    TResults->Print();
    // Save Tree file
    // NOTE: to open in the classic root browser -> "root --web=off"
    TFile fout(OutFileName, "recreate");
    fout.cd();
    TResults->Write();
    fout.Close();


    //***************************************************************************************************
    //***************************************************************************************************
    // CODE FOR DEBUGGING PURPOSES
    h1Type h1Pars = GetBaseLine(buffer, 1000);
    for (UInt_t i = 0; i < eventSz; i++) {
        amplitude = buffer[i] * ampRes - h1Pars.mean;
        // Fill histogram and plot with time shifted
        grAv->SetPoint(i, dt[i], amplitude);
        //h2Signal -> Fill(dt[i],amplitude);
    }
    grType sgProp = GetGrProp(grAv, thrs, sign, Wstart, Wend); // Graph from the signal input
   
    //***************************************************************************************************
    //***************************************************************************************************
    const Int_t n=5;
    TSpectrum* spectrumA = new TSpectrum(n); // Maximum number of peaks to find
    TSpectrum* spectrumC = new TSpectrum(n); // Maximum number of peaks to find
    Int_t nPeaksA = spectrumA->Search(h1ampSig, 2, "", 40); // Search for peaks
    Int_t nPeaksC = spectrumC->Search(h1charge, 2, "", 40); // Search for peaks

    const Double_t* xPosC = spectrumC->GetPositionX(); // X positions of peaks
    const Double_t* xPosA = spectrumA->GetPositionX(); // X positions of peaks

    for(int i=0; i<n; ++i) cout << xPosA[i] <<" "<< xPosC[i]<< endl;

    Int_t ind[n];
    TMath::Sort(n,xPosC,ind,false);
    for (int i = 0; i < n; ++i)     std::cout << xPosC[i] <<" " << "INDEX   "<< ind[i] << std::endl;
    cout << "..........." << endl;
    for (int i = 0; i < n; ++i)     std::cout << xPosC[ind[i]] <<" " << "INDEX   "<< i << std::endl;

    // Fit the first two peaks with Gaussian functions

    TCanvas * c4 = new TCanvas("c4","charge",800,800);
    c4-> SetLogy();
    h1charge->Draw();


    TF1 *g1    = new TF1("g1","gaus", xPosC[0] - 3, xPosC[0] + 3);
    g1->SetLineColor(kBlue);
    h1charge->Fit(g1, "R0", "SAME");
    
    TF1 *g2    = new TF1("g2","gaus", xPosC[1] - 2, xPosC[1] + 2);
    g2->SetLineColor(kRed);
    h1charge->Fit(g2, "R0", "SAME");


    TF1 *g3    = new TF1("g3","gaus", xPosC[2] - 3, xPosC[2] + 3);
    g3->SetLineColor(kMagenta);
    h1charge->Fit(g3, "R0", "SAME");

    
    g1->SetRange( g1->GetParameter(1)-50, g1->GetParameter(1)+50);
    g2->SetRange( g2->GetParameter(1)-50, g2->GetParameter(1)+50);
    g3->SetRange( g3->GetParameter(1)-50, g3->GetParameter(1)+50);

    g1->Draw("SAME");
    g2->Draw("SAME");
    g3->Draw("SAME");

    double bin_centerx = h1charge->GetXaxis()->GetBinCenter(1);
    //double bin_centerx = h1charge->GetXaxis()->GetBinCenter(h1charge->GetXaxis()->FindFixBin(0));
    cout << "bin 0 position " << bin_centerx << endl;
    
    TCanvas * cMaxAmp = new TCanvas("cMaxAmp","Max Amplitude",800,800);
    cMaxAmp-> SetLogy();
    h1MaxAmp ->Draw();
    

} // anaDigitizer_charge_v0
