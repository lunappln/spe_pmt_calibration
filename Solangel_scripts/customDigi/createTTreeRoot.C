void createTTreeRoot(){
      // create tree file
   gROOT->cd(); //make sure that the Tree is memory resident
   TTree *TResults = new TTree("T","test circular buffers");
   Double_t charge;
   Double_t energy;
   TResults->Branch("charge",&charge,"charge/D");
   TResults->Branch("energy",&energy,"energy/D");


    for (int i = 0; i < 10; i++)
    {
        charge = i*0.5;
        energy = i*1.5;
        TResults->Fill();
    }
    
    TResults -> Print();
    // save Tree file
    // NOTE: to open in the classic root browser -> "root --web=off"
    TFile fout("output.root","recreate");
    fout.cd();
    TResults->Write();
    fout.Close();
}