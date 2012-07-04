#include <cstdlib>
#include <cmath>
#include <iostream>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>
#include "MultiCanvas.h"
//#include "plot.C"
using namespace std;

const double pi=acos(-1.);
const double pi2=2*pi -1;

//! pt binning
double ptbins[] ={80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400};
const int bins  = sizeof(ptbins)/sizeof(Double_t) - 1;
const int nbins = bins;

//! actual nj  are 7;
//const char *calgo[knj] = {"icPu5","ak2PF","ak3PF","ak4PF","akPu2PF","akPu3PF","akPu4PF"}; 
const int knj=3;
const char *calgo[knj]= {"akPu2PF","akPu3PF","akPu4PF"}; 

const int ncen=7;
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","pp"};

const char *fopt="MEWLRQ0+"; int iFit=0; 
const int knpx=2000;
float fitmin=0.00;
float fitmax=5.00;

double xmin=80.;
double xmax=400.;
int maxEntry=10;

int GetPtBin(float /*pt*/);
void MakeHist(TH1 */*hist*/,int /*istat*/);
void MakeHistRMS(TH1 */*hRMS*/,float /*max*/,float /*min*/);
void MakeHistMean(TH1 */*Mean*/,float /*max*/,float /*min*/);

void FillMeanSigma(int /*ip*/,TH1 */*h1*/,TH1 */*ArM*/,TH1 */*RMS*/,TH1 */*Mean*/,TH1 */*Sigma*/);

int SmearResponse(int rfit=0)
{
  const char *reta = "eta2.0";
  bool iSave=true;

  float ketacut=1.6;
  if(strcmp(reta,"eta2.0")==0)ketacut=2;
  bool iSigma=false;


  if(rfit==0){fitmin=0.01;fitmax=2.75;}   
  else if(rfit==1){fitmin=0.00;fitmax=1.25;}
  else if(rfit==2){fitmin=0.00;fitmax=1.50;}
  else if(rfit==3){fitmin=0.00;fitmax=1.75;}
  else if(rfit==4){fitmin=0.00;fitmax=2.50;}
  else if(rfit==5){fitmin=0.00;fitmax=3.50;} 
  else if(rfit==6){fitmin=0.50;fitmax=1.25;}
  else if(rfit==7){fitmin=0.50;fitmax=1.50;}
  else if(rfit==8){fitmin=0.50;fitmax=1.75;}
  else if(rfit==9){fitmin=0.50;fitmax=2.50;}
  else if(rfit==10){fitmin=0.50;fitmax=3.50;}
  else if(rfit==11){fitmin=0.75;fitmax=1.25;}
  else if(rfit==12){fitmin=0.75;fitmax=1.50;}
  else if(rfit==13){fitmin=0.75;fitmax=1.75;}
  else if(rfit==14){fitmin=0.75;fitmax=2.50;}
  else if(rfit==15){fitmin=0.75;fitmax=3.50;}
  else if(rfit==16){fitmin=1.00;fitmax=1.25;}
  else if(rfit==17){fitmin=1.00;fitmax=1.50;}
  else if(rfit==18){fitmin=1.00;fitmax=1.75;}
  else if(rfit==19){fitmin=1.00;fitmax=2.50;}
  else if(rfit==20){fitmin=1.00;fitmax=3.50;}

  //! Input files 
  TFile *fin_pbpb = new TFile("input/Convolution/pbpb/Response_newHiForest_DJ_merged_pbpb_2012.root","r");
  //! 0   : pp w/o smearing
  //! 1-6 : pp to match pbpb smearing
  TFile *fin_pp     = new TFile("input/Convolution/ppsmeared/Response_newHiForest_DJ_merged_pp_2012.root","r");
  
  cout<<"\t"<<endl
      <<"rfit : "<<rfit<<"\t fitmin : "<<fitmin<<"\t fitmax : "<<fitmax<<endl
      <<"Input file name  pbpb : "<<fin_pbpb->GetName()<<endl
      <<"Input file name  pp   : "<<fin_pp->GetName()<<endl
      <<"# of pt bins : "<<bins<<endl
      <<"\t"<<endl;

  //return 0;

  //const char *mcent[ncen] = {"010"  ,"1020"  ,"2030"  ,"3050"  ,"5070"  ,"7090"  ,"pp"};
  const int   icol [ncen] = {     2,       2,       2,       2,      2,       2,    1};

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetErrorX(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetPadBorderSize(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetPalette(1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetOptFit(1);
  gStyle->SetPadBorderMode(0);


  //! PbPb Resposnse
  TH2F *hratiocorrrefpt[knj][ncen];

  TH1F *hratiocorrrefpt1D[knj][ncen][nbins];  
  TH1F *hMean[knj][ncen], *hSigma[knj][ncen], *hRMS[knj][ncen], *hArM[knj][ncen];

  //! pp Response
  TH2F *hratiocorrrefpt_pp[knj][ncen];
  TH1F *hratiocorrrefpt1D_pp[knj][ncen][nbins];  
  TH1F *hMean_pp [knj][ncen], *hArM_pp[knj][ncen], *hSigma_pp [knj][ncen], *hRMS_pp[knj][ncen];

  //! Ratio of smeared pp and pbpb
  TH1F *hRatioS[knj][ncen],  *hRatioUS[knj][ncen];
  TH1F *hRatioM[knj][ncen], *hRatioUM[knj][ncen];
  for(int nj=0;nj<knj;nj++){
    cout<<"nj : "<<nj<<Form("\t %s",calgo[nj])<<endl;

    //! Heavy-ion starts here and also read the smeared pp from heavy-ion 
    //! resolution matching
    //! Centrality dependence

    for(int icen=0;icen<ncen;icen++){

      //! pp /////////////////////////////
      hMean_pp [nj][icen] = new TH1F(Form("hMean_pp%d_%d",nj,icen),Form("<reco p_{T}/gen p_{T}> %s pp %d",calgo[nj],icen),nbins,ptbins);
      MakeHist(hMean_pp[nj][icen],0);
      hArM_pp [nj][icen] = new TH1F(Form("hArM_pp%d_%d",nj,icen),Form("<reco p_{T}/gen p_{T}> %s pp %d",calgo[nj],icen),nbins,ptbins);
      MakeHist(hArM_pp[nj][icen],0);
      hSigma_pp [nj][icen] = new TH1F(Form("hSigma_pp%d_%d",nj,icen),Form("#sigma(reco p_{T}/gen p_{T}) %s pp %d",calgo[nj],icen),nbins,ptbins);
      MakeHist(hSigma_pp[nj][icen],0);    
      hRMS_pp [nj][icen] = new TH1F(Form("hRMS_pp%d_%d",nj,icen),Form("RMS(reco p_{T}/gen p_{T}) %s pp %d",calgo[nj],icen),nbins,ptbins);
      MakeHist(hRMS_pp[nj][icen],0);        

      hratiocorrrefpt_pp[nj][icen]  = (TH2F*)fin_pp->Get(Form("hratiocorrrefpt%d_%d",nj+1,icen));
      hratiocorrrefpt_pp[nj][icen]->SetName(Form("hratiocorrrefpt_pp%d_%d",nj,icen));
      
      for(int ip=0;ip<hratiocorrrefpt_pp[nj][icen]->GetNbinsX();ip++){
	hratiocorrrefpt1D_pp[nj][icen][ip]  = (TH1F*)hratiocorrrefpt_pp [nj][icen]->ProjectionY(Form("hratiocorrrefpt1D_pp%d_%d_%d",nj,icen,ip),ip+1,ip+1);
	MakeHist(hratiocorrrefpt1D_pp[nj][icen][ip],1);
	FillMeanSigma(ip,hratiocorrrefpt1D_pp[nj][icen][ip],hArM_pp[nj][icen],hRMS_pp[nj][icen],hMean_pp[nj][icen],hSigma_pp[nj][icen]);
      }//! ip bin


      ///////////////////////////////////
      if(icen>ncen-2)continue;

      hMean [nj][icen] = new TH1F(Form("hMean%d_%d",nj,icen),Form("<reco p_{T}/gen p_{T}> %s pbpb %d",calgo[nj],icen),nbins,ptbins);
      MakeHist(hMean[nj][icen],0);
      hArM [nj][icen] = new TH1F(Form("hArM%d_%d",nj,icen),Form("<reco p_{T}/gen p_{T}> %s pbpb %d",calgo[nj],icen),nbins,ptbins);
      MakeHist(hArM[nj][icen],0);
      hSigma [nj][icen] = new TH1F(Form("hSigma%d_%d",nj,icen),Form("#sigma(reco p_{T}/gen p_{T}) %s pbpb %d",calgo[nj],icen),nbins,ptbins);
      MakeHist(hSigma[nj][icen],0);    
      hRMS [nj][icen] = new TH1F(Form("hRMS%d_%d",nj,icen),Form("RMS(reco p_{T}/gen p_{T}) %s pbpb %d",calgo[nj],icen),nbins,ptbins);
      MakeHist(hRMS[nj][icen],0);        
      
      cout<<"icen : "<<icen <<"\t ccent "<<ccent[icen]<<endl;
      hratiocorrrefpt[nj][icen] = (TH2F*)fin_pbpb->Get(Form("hratiocorrrefpt%d_%d",nj+4,icen));
      hratiocorrrefpt[nj][icen]->SetName(Form("hratiocorrrefpt_pbpb%d_%d",nj,icen));
      
      for(int ip=0;ip<hratiocorrrefpt[nj][icen]->GetNbinsX();ip++){
	hratiocorrrefpt1D[nj][icen][ip]  = (TH1F*)hratiocorrrefpt [nj][icen]->ProjectionY(Form("hratiocorrrefpt1D%d_%d_%d",nj,icen,ip),ip+1,ip+1);
	MakeHist(hratiocorrrefpt1D[nj][icen][ip],1);
	FillMeanSigma(ip,hratiocorrrefpt1D[nj][icen][ip],hArM[nj][icen],hRMS[nj][icen],hMean[nj][icen],hSigma[nj][icen]);
      }

    }//! icen loop ends
  }//! nj loop ends
  //return 0;

  //! Get the ratio
  for(int nj=0;nj<knj;nj++){
    for(int icen=0;icen<ncen-1;icen++){
      if(iSigma){

 	hRatioS[nj][icen] = (TH1F*)hSigma[nj][icen]->Clone(Form("hRatioS%d_%d",nj,icen));
	hRatioS[nj][icen]->SetName(Form("hRatioS%d_%d",nj,icen));
	hRatioS[nj][icen]->Divide(hSigma_pp[nj][icen]);

 	hRatioUS[nj][icen] = (TH1F*)hSigma[nj][icen]->Clone(Form("hRatioUS%d_%d",nj,icen));
	hRatioUS[nj][icen]->SetName(Form("hRatioUS%d_%d",nj,icen));
	hRatioUS[nj][icen]->Divide(hSigma_pp[nj][ncen-1]);

	hRatioM[nj][icen] = (TH1F*)hMean[nj][icen]->Clone(Form("hRatioM%d_%d",nj,icen));
	hRatioM[nj][icen]->SetName(Form("hRatioM%d_%d",nj,icen));
	hRatioM[nj][icen]->Divide(hMean_pp[nj][icen]);

	hRatioUM[nj][icen] = (TH1F*)hMean[nj][icen]->Clone(Form("hRatioUM%d_%d",nj,icen));
	hRatioUM[nj][icen]->SetName(Form("hRatioUM%d_%d",nj,icen));
	hRatioUM[nj][icen]->Divide(hMean_pp[nj][ncen-1]);

	
      }else{

	hRatioS[nj][icen] = (TH1F*)hRMS[nj][icen]->Clone(Form("hRatioS%d_%d",nj,icen));
	hRatioS[nj][icen]->SetName(Form("hRatioS%d_%d",nj,icen));
	hRatioS[nj][icen]->Divide(hRMS_pp[nj][icen]);

	hRatioUS[nj][icen] = (TH1F*)hRMS[nj][icen]->Clone(Form("hRatioUS%d_%d",nj,icen));
	hRatioUS[nj][icen]->SetName(Form("hRatioUS%d_%d",nj,icen));
	hRatioUS[nj][icen]->Divide(hRMS_pp[nj][ncen-1]);

	hRatioM[nj][icen] = (TH1F*)hArM[nj][icen]->Clone(Form("hRatioM%d_%d",nj,icen));
	hRatioM[nj][icen]->SetName(Form("hRatioM%d_%d",nj,icen));
	hRatioM[nj][icen]->Divide(hArM_pp[nj][icen]);

	hRatioUM[nj][icen] = (TH1F*)hArM[nj][icen]->Clone(Form("hRatioUM%d_%d",nj,icen));
	hRatioUM[nj][icen]->SetName(Form("hRatioUM%d_%d",nj,icen));
	hRatioUM[nj][icen]->Divide(hArM_pp[nj][ncen-1]);

      }
    }
  }

  /*
  for(int nj=0;nj<knj;nj++){
    for(int icen=0;icen<ncen-1;icen++){
      for(int ix=1;ix<=hRatioS[nj][icen]->GetNbinsX();ix++){
	hRatioS[nj][icen]->SetBinError(ix,0);
      }
    }
  }
  */
  /*
  TFile *fout = new TFile("Ratio_Smearedppandpbpb.root","RECREATE");
  fout->cd();
  for(int nj=0;nj<knj;nj++){
    for(int icen=0;icen<ncen-1;icen++){
      hRatioS[nj][icen]->SetName(Form("hRatioS_%s_%s",calgo[nj],ccent[icen]));
      hRatioS[nj][icen]->SetTitle(Form("Ratio Reco jet p_{T} PbPb / pp smeared MC %s %s",calgo[nj],ccent[icen]));
      hRatioS[nj][icen]->Write();
      hRatioUS[nj][icen]->SetName(Form("hRatioUS_%s_%s",calgo[nj],ccent[icen]));
      hRatioUS[nj][icen]->SetTitle(Form("Ratio Reco jet p_{T} pp smeared / pp MC %s %s",calgo[nj],ccent[icen]));
      hRatioUS[nj][icen]->Write();
    }  
  }
  fout->Close();
  return 0;
  */

  TText *text=0;
  int ipad=0;

 ipad=0;
 TCanvas *c99[ncen-1];
 for(int ic=0;ic<ncen-1;ic++){
   c99[ic] = new TCanvas(Form("c99_%d",ic),Form("%s %s Fitting plots",calgo[1],ccent[ic]),98,46,1215,930);
   c99[ic]->Divide(4,4,0,0);
   ipad=0;
   for(int ip=0;ip<nbins;ip++){      
     if( hratiocorrrefpt1D[1][ic][ip]->Integral()==0)continue;
     c99[ic]->cd(++ipad);
     if(ipad%4==0)gPad->SetRightMargin(0.02);
     gPad->SetBottomMargin(0.15);
     //if(ipad>20)gPad->SetBottomMargin(0.15);
     gPad->SetLogy();
     
     hratiocorrrefpt1D[1][ic][ip]->SetMaximum(14.634);
     hratiocorrrefpt1D[1][ic][ip]->SetMinimum(0.001);
     hratiocorrrefpt1D[1][ic][ip]->SetTitle(0);
     hratiocorrrefpt1D[1][ic][ip]->SetStats(0);
     hratiocorrrefpt1D[1][ic][ip]->GetXaxis()->SetTitleFont(42);
     hratiocorrrefpt1D[1][ic][ip]->GetXaxis()->SetLabelFont(42);
     hratiocorrrefpt1D[1][ic][ip]->GetXaxis()->SetLabelSize(0.08);
     hratiocorrrefpt1D[1][ic][ip]->GetXaxis()->SetTitleSize(0.07);
     
     hratiocorrrefpt1D[1][ic][ip]->GetYaxis()->SetTitle("");
     hratiocorrrefpt1D[1][ic][ip]->GetYaxis()->SetTitleFont(42);
     hratiocorrrefpt1D[1][ic][ip]->GetYaxis()->SetLabelFont(42);
     hratiocorrrefpt1D[1][ic][ip]->GetYaxis()->SetLabelSize(0.08);
     
     hratiocorrrefpt1D[1][ic][ip]->Draw("p");  
     
     hratiocorrrefpt1D_pp[1][ic][ip]->SetMarkerStyle(29);
     hratiocorrrefpt1D_pp[1][ic][ip]->SetMarkerColor(2);
     hratiocorrrefpt1D_pp[1][ic][ip]->SetLineColor(2);
     hratiocorrrefpt1D_pp[1][ic][ip]->Draw("psame");     
     
     hratiocorrrefpt1D_pp[1][ncen-1][ip]->SetMarkerColor(1);
     hratiocorrrefpt1D_pp[1][ncen-1][ip]->SetLineColor(1);
     hratiocorrrefpt1D_pp[1][ncen-1][ip]->Draw("histsame");     
     
     
     c99[ic]->Update();
     
     /*
       TPaveStats *ps = (TPaveStats*)  hratiocorrrefpt1D[1][ic][ip]->GetListOfFunctions()->FindObject("stats");
       ps->SetX1NDC(0.44);
       ps->SetY1NDC(0.25);       
       ps->SetX2NDC(0.95);
       ps->SetY2NDC(0.79);
       ps->SetTextFont(42);
       ps->Draw();
     */ 
     TPaveText *pt   = new TPaveText(0.4784302,0.8368631,0.8946255,0.9521693,"brNDC");
     pt->SetBorderSize(0);
     pt->SetFillColor(10);
     pt->SetTextFont(42);
     TText *text2 = pt->AddText(Form("%0.0f < p_{T} (GeV/c) < %0.0f",ptbins[ip], ptbins[ip+1]));
     text2->SetTextSize(0.07);
     pt->Draw();
     
     
     //! PbPb
     TPaveText *ptm1   = new TPaveText(0.45,0.65,0.95,0.85,"brNDC");
     ptm1->SetBorderSize(0);
     ptm1->SetFillColor(10);
     ptm1->SetTextFont(42);
     TText *tc1 = ptm1->AddText(Form("%s",ccent[ic]));
     tc1->SetTextSize(0.06);
     TText *tm1 = ptm1->AddText(Form("Mean = %0.4f #pm %0.04f",hratiocorrrefpt1D[1][ic][ip]->GetMean(),hratiocorrrefpt1D[1][ic][ip]->GetMeanError()));
     tm1->SetTextSize(0.06);
     TText *ts1 = ptm1->AddText(Form("RMS = %0.4f #pm %0.04f",hratiocorrrefpt1D[1][ic][ip]->GetRMS(),hratiocorrrefpt1D[1][ic][ip]->GetRMSError()));
     ts1->SetTextSize(0.06);
     ptm1->Draw();
     
     //! pp smeared
     TPaveText *ptm2   = new TPaveText(0.45,0.40,0.95,0.65,"brNDC");
     ptm2->SetBorderSize(0);
     ptm2->SetFillColor(10);
     ptm2->SetTextFont(42);
     TText *tc2 = ptm2->AddText("pp smeared");
     tc2->SetTextSize(0.06);
     tc2->SetTextColor(2);
     TText *tm2 = ptm2->AddText(Form("Mean = %0.4f #pm %0.04f",hratiocorrrefpt1D_pp[1][ic][ip]->GetMean(),hratiocorrrefpt1D_pp[1][ic][ip]->GetMeanError()));
     tm2->SetTextSize(0.06);
     tm2->SetTextColor(2);
     TText *ts2 = ptm2->AddText(Form("RMS = %0.4f #pm %0.04f",hratiocorrrefpt1D_pp[1][ic][ip]->GetRMS(),hratiocorrrefpt1D_pp[1][ic][ip]->GetRMSError()));
     ts2->SetTextSize(0.06);
     ts2->SetTextColor(2);
     ptm2->Draw();
     
     //! pp
     TPaveText *ptm3   = new TPaveText(0.45,0.15,0.95,0.40,"brNDC");
     ptm3->SetBorderSize(0);
     ptm3->SetFillColor(10);
     ptm3->SetTextFont(42);
     TText *tc3 = ptm3->AddText("pp");
     tc3->SetTextSize(0.06);
     tc3->SetTextColor(2);
     TText *tm3 = ptm3->AddText(Form("Mean = %0.4f #pm %0.04f",hratiocorrrefpt1D_pp[1][ncen-1][ip]->GetMean(),hratiocorrrefpt1D_pp[1][ncen-1][ip]->GetMeanError()));
     tm3->SetTextSize(0.06);
     tm3->SetTextColor(2);
     TText *ts3 = ptm3->AddText(Form("RMS = %0.4f #pm %0.04f",hratiocorrrefpt1D_pp[1][ncen-1][ip]->GetRMS(),hratiocorrrefpt1D_pp[1][ncen-1][ip]->GetRMSError()));
     ts3->SetTextSize(0.06);
     ts3->SetTextColor(2);
     ptm3->Draw();
     
     if(ipad==1){
       TPaveText *pt1 = new TPaveText(0.12,0.82,0.40,0.96,"brNDC");
       pt1->SetBorderSize(0);
       pt1->SetFillColor(10);
       pt1->SetTextFont(42);
       TText *text1 = pt1->AddText(calgo[1]);
       text1->SetTextSize(0.1);	
       pt1->Draw();
     }
   }
   c99[ic]->Update();
   if(iSave){
     c99[ic]->SaveAs(Form("AN/Smearing/Smearing_Distributions_%s_%s.png",calgo[1],ccent[ic]));
     c99[ic]->SaveAs(Form("AN/Smearing/Smearing_Distributions_%s_%s.pdf",calgo[1],ccent[ic]));
     c99[ic]->SaveAs(Form("AN/Smearing/Smearing_Distributions_%s_%s.C",calgo[1],ccent[ic]));
     c99[ic]->SaveAs(Form("AN/Smearing/Smearing_Distributions_%s_%s.eps",calgo[1],ccent[ic]));
   }
 }

  TPaveText *pt[ncen];
  ipad=0;
  TCanvas *c3 = new TCanvas("c3","sigma",15,131,1885,546);
  //TCanvas *c3 = new TCanvas("c3","sigma",316,23,1349,323);
  makeMultiPanelCanvas(c3,ncen-1,2,0.0,0.0,0.22,0.22,0.02);

  hSigma_pp[1][ncen-1]->SetLineColor(1);
  hRMS_pp  [1][ncen-1]->SetLineColor(1);

  //! here starts the centrality  
  for(int icen=ncen-2;icen>=0;icen--){
    c3->cd(++ipad); 
    //gPad->SetLogx();

    if(iSigma){
      MakeHistRMS(hSigma[1][icen],0.23,0.001);

      hSigma[1][icen]->Draw("p");

      hSigma_pp[1][icen]->SetMarkerStyle(30);
      hSigma_pp[1][icen]->SetMarkerColor(icol[icen]);
      hSigma_pp[1][icen]->SetLineColor(icol[icen]);
      hSigma_pp[1][icen]->SetMarkerSize(1.3);
      hSigma_pp[1][icen]->SetLineStyle(1);
      hSigma_pp[1][icen]->Draw("psame");    

      hSigma_pp[1][ncen-1]->SetLineColor(1);
      hSigma_pp[1][ncen-1]->SetLineStyle(1);
      hSigma_pp[1][ncen-1]->Draw("histsame");    

    }else{

      MakeHistRMS(hRMS[1][icen],0.23,0.001);
      hRMS[1][icen]->SetMarkerStyle(20);
      hRMS[1][icen]->SetMarkerColor(1);
      hRMS[1][icen]->SetLineColor(1);
      hRMS[1][icen]->SetMarkerSize(1.0);
      hRMS[1][icen]->Draw("p");

      hRMS_pp[1][icen]->SetMarkerStyle(30);
      hRMS_pp[1][icen]->SetMarkerColor(icol[icen]);
      hRMS_pp[1][icen]->SetLineColor(icol[icen]);
      hRMS_pp[1][icen]->SetMarkerSize(1.3);
      hRMS_pp[1][icen]->SetLineStyle(1);
      hRMS_pp[1][icen]->Draw("psame");    

      hRMS_pp[1][ncen-1]->SetLineColor(1);
      hRMS_pp[1][ncen-1]->SetLineStyle(1);
      hRMS_pp[1][ncen-1]->Draw("histsame");    
    }

    pt[icen]   = new TPaveText(0.269672,0.04089943,0.4756715,0.1469053,"brNDC");
    pt[icen]->SetBorderSize(0);
    pt[icen]->SetFillColor(10);
    pt[icen]->SetTextFont(42);
    text = pt[icen]->AddText(Form("%s",ccent[icen]));
    text->SetTextSize(0.09);
    pt[icen]->Draw();
    
    if(ipad==1){
      TPaveText *p1 = new TPaveText(0.45,0.65,0.75,0.85,"brNDC");
      p1->SetBorderSize(0);
      p1->SetFillColor(10);
      p1->SetTextFont(42);
      TText *text1 = p1->AddText("CMS Preliminary");
      TText *text2 = p1->AddText("Anti-k_{T}, PF, R =0.3");
      text1->SetTextSize(0.08);
      text2->SetTextSize(0.08);
      p1->Draw();
    }
    else if(ipad==2){
      TPaveText *p2 = new TPaveText(0.45,0.65,0.75,0.85,"brNDC");
      p2->SetBorderSize(0);
      p2->SetFillColor(10);
      p2->SetTextFont(42);
      TText *text3 = p2->AddText("PYTHIA+HYDJET 1.8");
      TText *text4 = p2->AddText(Form("|#eta_{jet}|<%0.0f",ketacut));
      text3->SetTextSize(0.08);
      text4->SetTextSize(0.08);
      p2->Draw();

    }else if(ipad==3){
      TLegend *l5 = new TLegend(0.23,0.60,0.65,0.96,NULL,"BRNDC");
      l5->SetHeader("");
      l5->SetBorderSize(0);
      l5->SetTextFont(42);
      l5->SetTextSize(0.07);
      l5->SetLineColor(1);
      l5->SetLineStyle(1);
      l5->SetLineWidth(1);
      l5->SetFillColor(10);
      l5->SetFillStyle(1001);
      l5->SetHeader("");			  
      
      if(iSigma){
	l5->AddEntry(hSigma[1][0] ,"PbPb","p"); 
	l5->AddEntry(hSigma_pp[1][0] ,"pp (smeared)","p"); 
	l5->AddEntry(hSigma_pp[1][ncen-1] ,"pp","l"); 
      }else{
	l5->AddEntry(hRMS[1][0] ,"PbPb","p"); 
	l5->AddEntry(hRMS_pp[1][0] ,"pp (smeared)","p"); 
	l5->AddEntry(hRMS_pp[1][ncen-1] ,"pp","l"); 
      }
      l5->Draw();
    }

    c3->cd(ipad+(ncen-1));
    //gPad->SetLogx();
    TLine *l2 = new TLine(xmin,1,xmax+250,1);
    l2->SetLineWidth(1);
    l2->SetLineStyle(2);

    
    hRatioS[1][icen]->SetTitle("");
    hRatioS[1][icen]->GetXaxis()->SetRangeUser(xmin,xmax);
    hRatioS[1][icen]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    hRatioS[1][icen]->GetXaxis()->CenterTitle(true);
    hRatioS[1][icen]->GetXaxis()->SetMoreLogLabels();
    hRatioS[1][icen]->GetXaxis()->SetNoExponent();
    hRatioS[1][icen]->GetXaxis()->SetNdivisions(507);
    hRatioS[1][icen]->GetXaxis()->SetLabelFont(42);
    hRatioS[1][icen]->GetXaxis()->SetLabelOffset(0.01);
    hRatioS[1][icen]->GetXaxis()->SetLabelSize(0.07);
    hRatioS[1][icen]->GetXaxis()->SetTitleSize(0.07);
    hRatioS[1][icen]->GetXaxis()->SetTitleOffset(1.18);
    hRatioS[1][icen]->GetXaxis()->SetTitleFont(42);
    hRatioS[1][icen]->GetYaxis()->SetTitle("Ratio");
    hRatioS[1][icen]->GetYaxis()->CenterTitle(true);
    hRatioS[1][icen]->GetYaxis()->SetNdivisions(507);
    hRatioS[1][icen]->GetYaxis()->SetLabelFont(42);
    hRatioS[1][icen]->GetYaxis()->SetLabelOffset(0.01);
    hRatioS[1][icen]->GetYaxis()->SetLabelSize(0.07);
    hRatioS[1][icen]->GetYaxis()->SetTitleSize(0.07);
    hRatioS[1][icen]->GetYaxis()->SetTitleOffset(1.44);
    hRatioS[1][icen]->GetYaxis()->SetDecimals(true);

    //hRatioS[1][icen]->SetMaximum(1.1893);
    //hRatioS[1][icen]->SetMinimum(0.8993);

    hRatioS[1][icen]->SetMaximum(1.3893);
    hRatioS[1][icen]->SetMinimum(0.5993);


    //hRatioS[1][icen]->SetMaximum(1.0593);
    //hRatioS[1][icen]->SetMinimum(0.943);
    hRatioS[1][icen]->SetMarkerStyle(30);
    hRatioS[1][icen]->SetMarkerSize(1.0);
    hRatioS[1][icen]->SetMarkerColor(2);
    hRatioS[1][icen]->SetLineColor(2);
    hRatioS[1][icen]->Draw("p");

    hRatioUS[1][icen]->SetMarkerStyle(24);
    hRatioUS[1][icen]->SetMarkerSize(0.9);
    hRatioUS[1][icen]->Draw("psame");

    if(ipad==1){
      TLegend *l6 = new TLegend(0.2508922,0.7662123,0.851988,0.9672997,NULL,"BRNDC");
      l6->SetHeader("");
      l6->SetBorderSize(0);
      l6->SetTextFont(42);
      l6->SetTextSize(0.07);
      l6->SetLineColor(1);
      l6->SetLineStyle(1);
      l6->SetLineWidth(1);
      l6->SetFillColor(10);
      l6->SetFillStyle(1001);
      l6->SetHeader("");			  
      l6->AddEntry(hRatioUS[1][icen],"PbPb / pp","p"); 
      l6->AddEntry(hRatioS[1][icen] ,"PbPb / pp (smeared)","p"); 
      l6->Draw();
    }
    l2->Draw();
  }
  if(iSave){
    c3->SaveAs(Form("AN/Smearing/Smearing_Resolutions_%s.png",calgo[1]));
    c3->SaveAs(Form("AN/Smearing/Smearing_Resolutions_%s.pdf",calgo[1]));
    c3->SaveAs(Form("AN/Smearing/Smearing_Resolutions_%s.C",calgo[1]));
    c3->SaveAs(Form("AN/Smearing/Smearing_Resolutions_%s.eps",calgo[1]));
  }


  //! Mean
  ipad=0;
  TCanvas *c4 = new TCanvas("c4","Mean",15,131,1885,546);
  makeMultiPanelCanvas(c4,ncen-1,2,0.0,0.0,0.22,0.22,0.02);

  hMean_pp[1][ncen-1]->SetLineColor(1);
  hArM_pp  [1][ncen-1]->SetLineColor(1);

  //! here starts the centrality  
  for(int icen=ncen-2;icen>=0;icen--){
    c4->cd(++ipad); 
    //gPad->SetLogx();

    if(iSigma){
      MakeHistMean(hMean[1][icen],1.046,0.934);
      
      hMean[1][icen]->Draw("p");

      hMean_pp[1][icen]->SetMarkerStyle(30);
      hMean_pp[1][icen]->SetMarkerColor(icol[icen]);
      hMean_pp[1][icen]->SetLineColor(icol[icen]);
      hMean_pp[1][icen]->SetMarkerSize(1.3);
      hMean_pp[1][icen]->SetLineStyle(1);
      hMean_pp[1][icen]->Draw("psame");    

      hMean_pp[1][ncen-1]->SetLineColor(1);
      hMean_pp[1][ncen-1]->SetLineStyle(1);
      hMean_pp[1][ncen-1]->Draw("histsame");    
      
    }else{

      MakeHistMean(hArM[1][icen],1.046,0.934);
      hArM[1][icen]->SetMarkerStyle(20);
      hArM[1][icen]->SetMarkerColor(1);
      hArM[1][icen]->SetLineColor(1);
      hArM[1][icen]->SetMarkerSize(1.0);
      hArM[1][icen]->Draw("p");

      hArM_pp[1][icen]->SetMarkerStyle(30);
      hArM_pp[1][icen]->SetMarkerColor(icol[icen]);
      hArM_pp[1][icen]->SetLineColor(icol[icen]);
      hArM_pp[1][icen]->SetMarkerSize(1.3);
      hArM_pp[1][icen]->SetLineStyle(1);
      hArM_pp[1][icen]->Draw("psame");    

      hArM_pp[1][ncen-1]->SetLineColor(1);
      hArM_pp[1][ncen-1]->SetLineStyle(1);
      hArM_pp[1][ncen-1]->Draw("histsame");    
    }

    pt[icen]   = new TPaveText(0.269672,0.04089943,0.4756715,0.1469053,"brNDC");
    pt[icen]->SetBorderSize(0);
    pt[icen]->SetFillColor(10);
    pt[icen]->SetTextFont(42);
    text = pt[icen]->AddText(Form("%s",ccent[icen]));
    text->SetTextSize(0.09);
    pt[icen]->Draw();
    
    if(ipad==1){
      TPaveText *p1 = new TPaveText(0.45,0.65,0.75,0.85,"brNDC");
      p1->SetBorderSize(0);
      p1->SetFillColor(10);
      p1->SetTextFont(42);
      TText *text1 = p1->AddText("CMS Preliminary");
      TText *text2 = p1->AddText("Anti-k_{T}, PF, R =0.3");
      text1->SetTextSize(0.08);
      text2->SetTextSize(0.08);
      p1->Draw();
    }
    else if(ipad==2){
      TPaveText *p2 = new TPaveText(0.45,0.65,0.75,0.85,"brNDC");
      p2->SetBorderSize(0);
      p2->SetFillColor(10);
      p2->SetTextFont(42);
      TText *text3 = p2->AddText("PYTHIA+HYDJET 1.8");
      TText *text4 = p2->AddText(Form("|#eta_{jet}|<%0.0f",ketacut));
      text3->SetTextSize(0.08);
      text4->SetTextSize(0.08);
      p2->Draw();

    }else if(ipad==3){
      TLegend *l5 = new TLegend(0.2275405,0.6732997,0.6497149,0.960796,NULL,"BRNDC");
      l5->SetHeader("");
      l5->SetBorderSize(0);
      l5->SetTextFont(42);
      l5->SetTextSize(0.07);
      l5->SetLineColor(1);
      l5->SetLineStyle(1);
      l5->SetLineWidth(1);
      l5->SetFillColor(10);
      l5->SetFillStyle(1001);
      l5->SetHeader("");			  
      
      if(iSigma){
	l5->AddEntry(hMean[1][0] ,"PbPb","p"); 
	l5->AddEntry(hMean_pp[1][0] ,"pp (smeared)","p"); 
	l5->AddEntry(hMean_pp[1][ncen-1] ,"pp","l"); 
      }else{
	l5->AddEntry(hArM[1][0] ,"PbPb","p"); 
	l5->AddEntry(hArM_pp[1][0] ,"pp (smeared)","p"); 
	l5->AddEntry(hArM_pp[1][ncen-1] ,"pp","l"); 
      }
      l5->Draw();
    }

    c4->cd(ipad+(ncen-1));
    //gPad->SetLogx();
    TLine *l2 = new TLine(xmin,1,xmax+250,1);
    l2->SetLineWidth(1);
    l2->SetLineStyle(2);

    
    hRatioM[1][icen]->SetTitle("");
    hRatioM[1][icen]->GetXaxis()->SetRangeUser(xmin,xmax);
    hRatioM[1][icen]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    hRatioM[1][icen]->GetXaxis()->CenterTitle(true);
    hRatioM[1][icen]->GetXaxis()->SetMoreLogLabels();
    hRatioM[1][icen]->GetXaxis()->SetNoExponent();
    hRatioM[1][icen]->GetXaxis()->SetNdivisions(507);
    hRatioM[1][icen]->GetXaxis()->SetLabelFont(42);
    hRatioM[1][icen]->GetXaxis()->SetLabelOffset(0.01);
    hRatioM[1][icen]->GetXaxis()->SetLabelSize(0.07);
    hRatioM[1][icen]->GetXaxis()->SetTitleSize(0.07);
    hRatioM[1][icen]->GetXaxis()->SetTitleOffset(1.18);
    hRatioM[1][icen]->GetXaxis()->SetTitleFont(42);
    hRatioM[1][icen]->GetYaxis()->SetTitle("Ratio");
    hRatioM[1][icen]->GetYaxis()->CenterTitle(true);
    hRatioM[1][icen]->GetYaxis()->SetNdivisions(507);
    hRatioM[1][icen]->GetYaxis()->SetLabelFont(42);
    hRatioM[1][icen]->GetYaxis()->SetLabelOffset(0.01);
    hRatioM[1][icen]->GetYaxis()->SetLabelSize(0.07);
    hRatioM[1][icen]->GetYaxis()->SetTitleSize(0.07);
    hRatioM[1][icen]->GetYaxis()->SetTitleOffset(1.44);
    hRatioM[1][icen]->GetYaxis()->SetDecimals(true);

    //hRatioM[1][icen]->SetMaximum(1.1893);
    //hRatioM[1][icen]->SetMinimum(0.8993);

    hRatioM[1][icen]->SetMaximum(1.0593);
    hRatioM[1][icen]->SetMinimum(0.955);


    //hRatioM[1][icen]->SetMaximum(1.0593);
    //hRatioM[1][icen]->SetMinimum(0.943);
    hRatioM[1][icen]->SetMarkerStyle(30);
    hRatioM[1][icen]->SetMarkerSize(1.0);
    hRatioM[1][icen]->SetMarkerColor(2);
    hRatioM[1][icen]->SetLineColor(2);
    hRatioM[1][icen]->Draw("p");

    hRatioUM[1][icen]->SetMarkerStyle(24);
    hRatioUM[1][icen]->SetMarkerSize(0.9);
    hRatioUM[1][icen]->Draw("psame");

    if(ipad==1){
      TLegend *l6 = new TLegend(0.2665729,0.7350091,0.6847265,0.9603656,NULL,"BRNDC");
      l6->SetHeader("");
      l6->SetBorderSize(0);
      l6->SetTextFont(42);
      l6->SetTextSize(0.07);
      l6->SetLineColor(1);
      l6->SetLineStyle(1);
      l6->SetLineWidth(1);
      l6->SetFillColor(10);
      l6->SetFillStyle(1001);
      l6->SetHeader("");			  
      l6->AddEntry(hRatioUM[1][icen],"PbPb / pp","p"); 
      l6->AddEntry(hRatioM[1][icen] ,"PbPb / pp (smeared)","p"); 
      l6->Draw();
    }
    l2->Draw();
  }
  if(iSave){
    c4->SaveAs(Form("AN/Smearing/Smearing_Means_%s.png",calgo[1]));
    c4->SaveAs(Form("AN/Smearing/Smearing_Means_%s.pdf",calgo[1]));
    c4->SaveAs(Form("AN/Smearing/Smearing_Means_%s.C",calgo[1]));
    c4->SaveAs(Form("AN/Smearing/Smearing_Means_%s.eps",calgo[1]));
  }







  return 0;
}
void MakeHist(TH1 *histo,int istat)
{
  histo->SetStats(istat);
  histo->SetMarkerStyle(24);
  histo->SetMarkerColor(1);
  histo->SetLineColor(1);
  histo->SetLineStyle(1);
  histo->GetXaxis()->SetTitle("p_{T}^{GenJet} (GeV/c)");
  histo->GetXaxis()->CenterTitle(true);
  histo->GetYaxis()->SetTitle("<p_{T}^{RecoJet}/p_{T}^{GenJet}>");
  histo->GetYaxis()->CenterTitle(true);
}
void FillMeanSigma(int ip,TH1 *h1F,TH1 *hArM,TH1 *hRMS,TH1 *hMean,TH1 *hSigma)
{

  TF1 *f1 = new TF1("f1","(([0]/(2*pi*[1]*[1]))*exp(-1.*((x-[2])*(x-[2])/(2*[1]*[1]))))",0,5);
  f1->SetParameters(1,0.1,1);
  f1->SetParNames("A","#sigma","mean");
  f1->SetParLimits(0,0,20);     //! A
  f1->SetParLimits(1,0.,10.0);  //! Sigma
  f1->SetParLimits(2,fitmin,fitmax); //! Mean
  f1->SetLineWidth(1);
  f1->SetNpx(knpx);
  
  float mm=0,ss=0,p0=0;  
  if(h1F->GetEntries()<maxEntry){
    h1F->Scale(0.);
    hArM  ->SetBinContent(ip+1,-9999);
    hArM  ->SetBinError  (ip+1,0);
    hRMS  ->SetBinContent(ip+1,-9999);
    hRMS  ->SetBinError  (ip+1,0);
  }
  if(h1F->Integral()>0){
    h1F->Scale(1./h1F->Integral());
    if(iFit==0){
      h1F->Fit("gaus",fopt,"",fitmin,fitmax);
      TF1* f2 = (TF1*)h1F->GetFunction("gaus");
      f2->SetLineWidth(1);
      f2->SetLineStyle(2);
      f2->SetNpx(knpx);
      hMean ->SetBinContent(ip+1,f2->GetParameter(1));
      hSigma->SetBinContent(ip+1,f2->GetParameter(2));
      
      if(strcmp(fopt,"MLLRQ+")==0){
	hMean ->SetBinError  (ip+1,h1F->GetMeanError());
	hSigma->SetBinError  (ip+1,h1F->GetRMSError());
      }else{
	hMean ->SetBinError(ip+1,f2->GetParError(1));
	hSigma->SetBinError(ip+1,f2->GetParError(2));
      }
    }else{
      mm = h1F->GetMean();
      ss = h1F->GetRMS();
      p0 = h1F->GetMaximum();
      f1->SetParameters(p0,ss,mm);
      f1->SetParLimits(0,0,2*p0);
      f1->SetParLimits(1,0,2*ss);
      f1->SetParLimits(2,fitmin,fitmax);
      //f1->SetParLimits(2,mm-2.5*ss,mm+2.5*ss);
      
      h1F->Fit("f1",fopt,"",fitmin,fitmax);
      hMean ->SetBinContent(ip+1,f1->GetParameter(2));
      hSigma->SetBinContent(ip+1,f1->GetParameter(1));
      
      //hMean ->SetBinError  (ip+1,f1->GetParError(2));
      //hSigma->SetBinError  (ip+1,f1->GetParError(1));
      
      if(strcmp(fopt,"MLLRQ+")==0){
	hMean ->SetBinError  (ip+1,h1F->GetMeanError());
	hSigma->SetBinError  (ip+1,h1F->GetRMSError());
      }else{
	hMean ->SetBinError  (ip+1,f1->GetParError(2));
	hSigma->SetBinError  (ip+1,f1->GetParError(1));
      }
    }
    hArM  ->SetBinContent(ip+1,h1F->GetMean());
    hArM  ->SetBinError  (ip+1,h1F->GetMeanError());
    hRMS  ->SetBinContent(ip+1,h1F->GetRMS());
    hRMS  ->SetBinError  (ip+1,h1F->GetRMSError());
  }
}
void MakeHistRMS(TH1 *h1,float ymax,float ymin)
{

  h1->SetTitle("");
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelOffset(0.01);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitle("#sigma (RecoJet p_{T} / GenJet p_{T})");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelOffset(0.01);
  h1->GetYaxis()->SetLabelSize(0.09);
  h1->GetYaxis()->SetTitleSize(0.09);
  h1->GetYaxis()->SetTitleOffset(1.12);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetDecimals(true);

}
void MakeHistMean(TH1 *h1,float ymax,float ymin)
{
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->SetTitle("");
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetLabelOffset(0.005);
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(1.50);
  h1->GetYaxis()->SetLabelSize(0.07);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelFont(42);
}
int GetPtBin(float pt)
{
  for(int ix=0;ix<bins;ix++){
    if(pt>=ptbins[ix] && pt<ptbins[ix+1]){
      return ix;
    }
  }
  return -1;
}
