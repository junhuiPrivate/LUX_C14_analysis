///original author: Rene Brun
#include "TH2.h"

   Double_t fitFunctionPol2(Double_t *x, Double_t *par) {
      return par[0] + par[1]*x[0] ;
   }

   Double_t fitFunctionPol3(Double_t *x, Double_t *par) {
      return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] ; 
   }

   Double_t fitFunctionPol5(Double_t *x, Double_t *par) {
      return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0] 
             + par[5]*x[0]*x[0]*x[0]*x[0]*x[0] ;
   }

   Double_t fitFunction(Double_t *x, Double_t *par) {
      return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0]
             +  par[5]*x[0]*x[0]*x[0]*x[0]*x[0] +  par[6]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
   }

   Double_t fitFunctionEvan(Double_t *x, Double_t *par) {
      return (par[0] + par[1]*x[0]) * (1 - TMath::Exp(-x[0]/par[2])) + par[3];
   }

void analyze_S1_max_peak_area_VS_S1_area(){
   
   ///TFile *file = new TFile("histos_S1_prompt_fraction_tlx_VS_S1_area_Cut_S1_area_LT350phd_TLX10kbin.root"); 
   TFile *file = new TFile("histos_S1_max_peak_area_VS_S1_area_Cut_S1_area_LT350phd_S1MaxPeakArea10Kbin.root"); 
   TH2D *S1_max_peak_area_VS_S1_area;
   file->GetObject("S1_max_peak_area_VS_S1_area;1", S1_max_peak_area_VS_S1_area); 
   ///The ";1" should not be ignored.
   
   int binMax_s1area = 350, Max_s1area = binMax_s1area, binMax_tlx = 10000, Max_tlx =100; 
   Int_t binNX = 0, binNY = 0, binTmp = 0, binNY_5Percent = 0;
   Int_t S1phdUpCut = 350;
   Double_t binContentX = 0,binContentX99P = 0, binContentY = 0,binContentX_5Percent = 0, binContentXTmp = 0, binContentXAcum = 0;
   const Int_t n1 = binMax_s1area ;
   Double_t x1[n1], y1[n1], accept[n1];
   
   for(binNX = 1; binNX < binMax_s1area+1; binNX++){
       binContentX = 0;
       for(binNY = 1; binNY < binMax_tlx+1 ; binNY++){
       binContentX += S1_max_peak_area_VS_S1_area->GetBinContent(binNX, binNY);
       }
       binContentX99P = binContentX * 0.99;
       ///binContentX_5Percent = binContentX * 0.005;

       binContentXAcum = 0;
       for(binNY = 1; binNY < binMax_tlx+1 ; binNY++){
       binContentXAcum += S1_max_peak_area_VS_S1_area->GetBinContent(binNX, binNY);
          if (binContentXAcum >= binContentX99P){
             binTmp = binNY;
             break ;
          }
       }
       x1[binNX] = binNX;
       y1[binNX] = binTmp* double(Max_tlx) / double(binMax_tlx);
       ///if(binNX > 35 && binNX < 45) cout << "binNX = " << binNX << ", x1[binNX] =  " <<  x1[binNX] << ", y1[binNX] = " << y1[binNX] << endl;
       ///cout << "x1[binNX] = " << x1[binNX] << ", y1[binNX] = " << y1[binNX] << " ~~~ ... EndLoopX ..."<< endl<< endl;
       ///cout << "binTmp = " << binTmp << ", binContentX_5Percent = " << binContentX_5Percent << " ~~~ ... EndLoopX ..."<< endl<< endl;
   }
   TCanvas *c1 = new TCanvas("c1","Graph Example",200,10,700,500);
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

   ///Double_t par[4], evaluedValues = 0., acceptance = 0.;
   Double_t par[2], evaluedValues = 0., acceptance = 0.;
   TF1 *fitFcnPol2 = new TF1("fitFcnPol2",fitFunctionPol2,1,S1phdUpCut,2);
   gr1->Fit("fitFcnPol2","R");
   fitFcnPol2->GetParameters(&par[0]);

   TF1 *fitFcnPol2H3 = new TF1("fitFcnPol2H3",fitFunctionPol2,1,S1phdUpCut-300,2);
   fitFcnPol2H3->SetParameters(2.94551845,0.10207864);
   fitFcnPol2H3->SetFillColor(19);
   fitFcnPol2H3->SetFillStyle(0);
   fitFcnPol2H3->SetLineColor(7);
   fitFcnPol2H3->Draw("Same");

   ///The following loop is to get the acceptance 

   for(binNX = 1; binNX < S1phdUpCut; binNX++){
       binContentX = 0;
       for(binNY = 1; binNY < binMax_tlx+1 ; binNY++){
       binContentX += S1_max_peak_area_VS_S1_area->GetBinContent(binNX, binNY);
       }
       binContentX99P = binContentX * 0.99;
       ///binContentX_5Percent = binContentX * 0.005;

       binContentXAcum = 0;
       for(binNY = 1; binNY < binMax_tlx+1 ; binNY++){
       binContentXAcum += S1_max_peak_area_VS_S1_area->GetBinContent(binNX, binNY);
          if (binContentXAcum >= binContentX99P){
             binTmp = binNY;
             break ;
          }
       }

       Double_t binTmpNormalized = binTmp* double(Max_tlx) / double(binMax_tlx);
       ///evaluedValues = fitFcn->Eval(binNX);
       evaluedValues = fitFcnPol2->Eval(binNX);

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

 }

