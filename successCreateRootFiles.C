///original author: Rene Brun

#include "TH2.h"

 void successCreateRootFiles( ){
   auto inputFile = new TFile("c14_endofLUX_09262017.root");
   auto c14DataTree = (TTree*)inputFile->Get("filter");

   int binMax_s1area = 350, Max_s1area = 350, binMax_tlx = 10000, Max_tlx =1; 

   TH2D *S1_prompt_fraction_tlx_VS_S1_area = new TH2D("S1_prompt_fraction_tlx_VS_S1_area","S1_prompt_fraction_tlx VS S1_area", binMax_s1area, 0, Max_s1area, binMax_tlx, 0, Max_tlx);

   auto c1 = new TCanvas("c1","test",10,10,600,500);

   c14DataTree->Draw("S1_prompt_fraction_tlx:S1_area >> S1_prompt_fraction_tlx_VS_S1_area","S1_area<350 && S1_prompt_fraction_tlx <=1");

   S1_prompt_fraction_tlx_VS_S1_area = (TH2D*)gDirectory->Get("S1_prompt_fraction_tlx_VS_S1_area"); 

   ///TFile f("histos_S1_prompt_fraction_tlx_VS_S1_area_Cut_S1_area_LT300phd.root","RECREATE");
   ///TFile f("histos_S1_prompt_fraction_tlx_VS_S1_area_Cut_S1_area_LT300phd_50bin.root","RECREATE");
   TFile f("histos_S1_prompt_fraction_tlx_VS_S1_area_Cut_S1_area_LT350phd_TLX10kbin.root","RECREATE");
   ///TFile f("histos_S1_prompt_fraction_tlx_VS_S1_area_Cut_S1_area_LT300phd_100bin.root","RECREATE");
   S1_prompt_fraction_tlx_VS_S1_area->Write(); 
   /*
   Int_t binNX = 1, binNY = 1, binTmp = 0;
   Double_t binContentX = 0,binContentX95P = 0, binContentY = 0;
  /// for(binNX = 1; binNX < binMax_s1area+1; binNX++){
       ///for(binNY = 1; binNY < binMax_s1area+1 ; binNY++){
       binTmp = S1_prompt_fraction_tlx_VS_S1_area->GetBin{binNX, binNY};
       binContentX = S1_prompt_fraction_tlx_VS_S1_area->GetBinContent{binTmp};
       ///}
       binContentX95P = binContentX * 0.95;
       cout << endl << "binContentX = " << binContentX  << "binContentX95P = " << binContentX95P << endl;
   ///}
*/
}

