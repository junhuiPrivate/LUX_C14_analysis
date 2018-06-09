///original author: Rene Brun
#include "TH2.h"

   Double_t fitFunctionPol3(Double_t *x, Double_t *par) {
      return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] ; 
   }

   Double_t fitFunctionPol5(Double_t *x, Double_t *par) {
      return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0]; 
   }

   Double_t fitFunction(Double_t *x, Double_t *par) {
      return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0]
             +  par[5]*x[0]*x[0]*x[0]*x[0]*x[0] +  par[6]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
   }

   Double_t fitFunctionEvan(Double_t *x, Double_t *par) {
      return (par[0] + par[1]*x[0]) * (1 - TMath::Exp(-x[0]/par[2])) + par[3];
   }

void analysisC14Data(){
   
   TFile *file = new TFile("histos_S1_prompt_fraction_tlx_VS_S1_area_Cut_S1_area_LT350phd_TLX10kbin.root"); 
   TH2D *S1_prompt_fraction_tlx_VS_S1_area;
   file->GetObject("S1_prompt_fraction_tlx_VS_S1_area;1", S1_prompt_fraction_tlx_VS_S1_area); 
   ///The ";1" should not be ignored.
   
   int binMax_s1area = 350, Max_s1area = binMax_s1area, binMax_tlx = 10000, Max_tlx =1; 
   Int_t binNX = 0, binNY = 0, binTmp = 0, binNY_5Percent = 0;
   Int_t S1phdUpCut = 350;
   Double_t binContentX = 0,binContentX95P = 0, binContentY = 0,binContentX_5Percent = 0, binContentXTmp = 0, binContentXAcum = 0;
   const Int_t n1 = binMax_s1area ;///301
   Double_t x1[n1], y1[n1], accept[n1];
   
   for(binNX = 1; binNX < binMax_s1area+1; binNX++){
       binContentX = 0;
       for(binNY = 1; binNY < binMax_tlx+1 ; binNY++){
       binContentX += S1_prompt_fraction_tlx_VS_S1_area->GetBinContent(binNX, binNY);
       }
       binContentX95P = binContentX * 0.995;
       binContentX_5Percent = binContentX * 0.005;

       binContentXAcum = 0;
       for(binNY = 1; binNY < binMax_tlx+1 ; binNY++){
       binContentXAcum += S1_prompt_fraction_tlx_VS_S1_area->GetBinContent(binNX, binNY);
          if (binContentXAcum >= binContentX_5Percent){
             binTmp = binNY;
             break ;
          }
       }
       x1[binNX] = binNX;
       y1[binNX] = binTmp* double(Max_tlx) / double(binMax_tlx);
       ///cout << x1[binNX] << "   " << y1[binNX] << endl;
       ///cout << "x1[binNX] = " << x1[binNX] << ", y1[binNX] = " << y1[binNX] << " ~~~ ... EndLoopX ..."<< endl<< endl;
       ///cout << "binTmp = " << binTmp << ", binContentX_5Percent = " << binContentX_5Percent << " ~~~ ... EndLoopX ..."<< endl<< endl;
   }
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
   TGraph *gr1 = new TGraph(S1phdUpCut,x1,y1);
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(7);
   gr1->SetLineColor(kGreen);
   gr1->SetLineStyle(2);
   gr1->SetLineWidth(4);

   TAxis* Xaxis = gr1->GetXaxis();
   TAxis* Yaxis = gr1->GetYaxis();

   c1->Divide(1,2);
   c1->cd(1);
   gr1->Draw("AP");

   gStyle->SetOptFit(1);

   Double_t par[4], evaluedValues = 0., acceptance = 0.;
   TF1 *fitFcnEvan = new TF1("fitFcnEvan",fitFunctionEvan,1,S1phdUpCut,4);
   fitFcnEvan->SetParameters(0.507021473,1.46005822e-04,8.93312493,0.268578463);
   gr1->Fit("fitFcnEvan","R");
   fitFcnEvan->GetParameters(&par[0]);

   TF1 *fitFcnPol5H3 = new TF1("fitFcnPol5H3",fitFunctionPol5,1,S1phdUpCut-300,5);
   fitFcnPol5H3->SetParameters(3.77155251e-01,3.88637498e-02,-1.63849839e-03,3.19166994e-05,-2.31131979e-07);
   fitFcnPol5H3->SetFillColor(19);
   fitFcnPol5H3->SetFillStyle(0);
   fitFcnPol5H3->SetLineColor(7);
   fitFcnPol5H3->Draw("Same");


   ///The following loop is to get the acceptance 

   for(binNX = 1; binNX < S1phdUpCut; binNX++){
       binContentX = 0;
       for(binNY = 1; binNY < binMax_tlx+1 ; binNY++){
       binContentX += S1_prompt_fraction_tlx_VS_S1_area->GetBinContent(binNX, binNY);
       }
       binContentX95P = binContentX * 0.995;
       binContentX_5Percent = binContentX * 0.005;

       binContentXAcum = 0;
       for(binNY = 1; binNY < binMax_tlx+1 ; binNY++){
       binContentXAcum += S1_prompt_fraction_tlx_VS_S1_area->GetBinContent(binNX, binNY);
          if (binContentXAcum >= binContentX_5Percent){
             binTmp = binNY;
             break ;
          }
       }

       Double_t binTmpNormalized = binTmp* double(Max_tlx) / double(binMax_tlx);
       ///evaluedValues = fitFcn->Eval(binNX);
       evaluedValues = fitFcnEvan->Eval(binNX);

       if(evaluedValues >= binTmpNormalized) acceptance = 1.;
       else acceptance = evaluedValues / binTmpNormalized  ;
       ///cout << "bin # = " << binNX <<  ", acc = " << acceptance  << endl;
       x1[binNX] = binNX;
       accept[binNX] = acceptance;
   }



   gPad->Modified();
   gPad->Update(); 

   TGraph *gr2 = new TGraph(S1phdUpCut,x1,accept);
   gr2->SetMarkerColor(kRed);
   gr2->SetMarkerStyle(3);
   gr2->SetLineColor(kGreen);
   gr2->SetLineStyle(2);
   gr2->SetLineWidth(4);

   c1->cd(2);
   gr2->Draw("AP");

   TPaveStats *Graph_st = (TPaveStats*)(gr1->GetListOfFunctions()->FindObject("stats"));
   Graph_st->SetOptFit(1);
   gPad->Modified();
   gPad->Update(); // make sure it's (re)drawn
   Graph_st->SetX1NDC(0.2);
   Graph_st->SetX2NDC(0.6);
   Graph_st->SetY1NDC(0.2);
   Graph_st->SetY2NDC(0.6);

   c1->Update();
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();





   ///Find the bin on Y axis which has accumulated bin content over 5% of binContentX. 
   ///auto c1 = new TCanvas("c1","test",10,10,600,500);
   ///S1_prompt_fraction_tlx_VS_S1_area->Draw();

   /// The following are reading histograms from a ROOT file    
/*   
   auto inputFile = new TFile("c14_endofLUX_09262017.root");
   auto c14DataTree = (TTree*)inputFile->Get("filter");

   int binMax_s1area = 300, Max_s1area = 300, binMax_tlx =50, Max_tlx =1; 

   TH2F *S1_prompt_fraction_tlx_VS_S1_area = new TH2F("S1_prompt_fraction_tlx_VS_S1_area","S1_prompt_fraction_tlx VS S1_area", binMax_s1area, 0, Max_s1area, binMax_tlx, 0, Max_tlx);

   auto c1 = new TCanvas("c1","test",10,10,600,500);

   c14DataTree->Draw("S1_prompt_fraction_tlx:S1_area >> S1_prompt_fraction_tlx_VS_S1_area","S1_area<300 && S1_prompt_fraction_tlx <=1");

   S1_prompt_fraction_tlx_VS_S1_area = (TH2F*)gDirectory->Get("S1_prompt_fraction_tlx_VS_S1_area"); 

   TFile f("Histos_S1_prompt_fraction_tlx_VS_S1_area_Cut_S1_area_LT300phd.root","RECREATE");
   S1_prompt_fraction_tlx_VS_S1_area->Write(); 
*/
 }

