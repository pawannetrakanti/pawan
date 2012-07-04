#include <TROOT.h>
#include <TMath.h>
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

using namespace std;
double xmin=100.;
double xmax=300.;

void GetPileup(TH2 */*h1*/,TH1 */*hmean*/,TH1 */*hrms*/);
int PileupAn()
{
  float ketacut=2.0;
  bool iSave=true;

  const int ncen=6;
  const char *ksp  [ncen] = {"pbpb","pbpb" ,"pbpb"  ,"pbpb"  ,"pbpb"  ,"pbpb"  };
  const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%"};
  Color_t dcol[ncen] = {kRed,kBlue,kMagenta,kCyan,kGreen+1,kBlue+4};
  int     dsty[ncen] = {24,25,26,28,30,32};
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
  //gStyle->SetOptFit(1);
  gStyle->SetPadBorderMode(0);
  //gStyle->SetOptStat("neMR");
  gStyle->SetOptStat("0");

  //! data pt binning     
  double ptbins[] ={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
  const int bins = sizeof(ptbins)/sizeof(double) - 1;

  TFile *fin[2];
  fin[0]  = new TFile("/home/pawan/CMS/Analysis/Raa/Output/HLT_HIJet80_v1_HiForest2_pbpb_pt80_2012.root","r"); //! data
  fin[1]  = new TFile("input/Convolution/pbpb/Response_newHiForest_DJ_merged_pbpb_2012.root","r"); //! MC

  const int knj=3;
  const char *calgo[knj]= {
    "akPu2PF","akPu3PF","akPu4PF"
  };             

  const char *meta[2] = {"|#eta|<1.3","1.3<|#eta|<2.0"};
  //! pu
  TH2F *hjetptpu_data[knj][ncen];
  TH1F *hpmean_data[knj][ncen];
  TH1F *hprms_data [knj][ncen];

  TH2F *hjetptpu_etab_data[knj][ncen][2];
  TH1F *hpmean_etab_data[knj][ncen][2];
  TH1F *hprms_etab_data [knj][ncen][2];

  TH2F *hjetptpu_mc[knj][ncen];
  TH1F *hpmean_mc[knj][ncen];
  TH1F *hprms_mc [knj][ncen];

  TH2F *hjetptpu_etab_mc[knj][ncen][2];
  TH1F *hpmean_etab_mc[knj][ncen][2];
  TH1F *hprms_etab_mc [knj][ncen][2];

  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){
      
      //! PbPb
      //! Pile up subtraction histograms
      hjetptpu_data  [nj][ic] = (TH2F*)fin[0]->Get(Form("hjetptpu%d_%d",nj,ic));
      hjetptpu_data  [nj][ic]->SetName(Form("hjetptpu_%s_%s_%d",ksp[ic],calgo[nj],ic));
      
      //! pileup mean from MC
      hpmean_data[nj][ic] = new TH1F(Form("hpmean_data_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hpmea_data_n%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hpmean_data[nj][ic]->Sumw2();
      hprms_data [nj][ic] = new TH1F(Form("hprms_data_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hprms_data_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hprms_data [nj][ic]->Sumw2();

      hjetptpu_mc  [nj][ic] = (TH2F*)fin[1]->Get(Form("hjetptpu%d_%d",nj+4,ic));
      hjetptpu_mc  [nj][ic]->SetName(Form("hjetptpu_%s_%s_%d",ksp[ic],calgo[nj],ic));

      //! pileup mean from MC
      hpmean_mc[nj][ic] = new TH1F(Form("hpmean_mc_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hpmea_mc_n%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hpmean_mc[nj][ic]->Sumw2();
      hprms_mc[nj][ic] = new TH1F(Form("hprms_mc_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hprms_mc_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hprms_mc[nj][ic]->Sumw2();

      GetPileup(hjetptpu_data[nj][ic],hpmean_data[nj][ic],hprms_data [nj][ic]);
      GetPileup(hjetptpu_mc[nj][ic],hpmean_mc[nj][ic],hprms_mc[nj][ic]);

      //! eta dependence of pile up subtraction
      for(int ie=0;ie<2;ie++){
	hjetptpu_etab_data  [nj][ic][ie] = (TH2F*)fin[0]->Get(Form("hjetptpu_etab%d_%d_%d",nj,ic,ie));
	hjetptpu_etab_data  [nj][ic][ie]->SetName(Form("hjetptpu_etab_data_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie));

	//! pileup mean 
	hpmean_etab_data[nj][ic][ie] = new TH1F(Form("hpmean_etab_data_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie),Form("hpmean_etab_data_%s_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic],meta[ie]),bins,ptbins);
	hpmean_etab_data[nj][ic][ie]->Sumw2();
	hprms_etab_data[nj][ic][ie] = new TH1F(Form("hprms_etab_data_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie),Form("hprms_etab_data_%s_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic],meta[ie]),bins,ptbins);
	hprms_etab_data[nj][ic][ie]->Sumw2();

	GetPileup(hjetptpu_etab_data  [nj][ic][ie],hpmean_etab_data[nj][ic][ie],hprms_etab_data[nj][ic][ie]);
	
	hjetptpu_etab_mc  [nj][ic][ie] = (TH2F*)fin[1]->Get(Form("hjetptpu_etab%d_%d_%d",nj,ic,ie));
	hjetptpu_etab_mc  [nj][ic][ie]->SetName(Form("hjetptpu_etab_mc_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie));

	//! pileup mean 
	hpmean_etab_mc[nj][ic][ie] = new TH1F(Form("hpmean_etab_mc_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie),Form("hpmean_etab_mc_%s_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic],meta[ie]),bins,ptbins);
	hpmean_etab_mc[nj][ic][ie]->Sumw2();
	hprms_etab_mc[nj][ic][ie] = new TH1F(Form("hprms_etab_mc_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie),Form("hprms_etab_mc_%s_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic],meta[ie]),bins,ptbins);
	hprms_etab_mc[nj][ic][ie]->Sumw2();

	GetPileup(hjetptpu_etab_mc[nj][ic][ie],hpmean_etab_mc[nj][ic][ie],hprms_etab_mc[nj][ic][ie]);
      }
    }
  }

  int ipad=0;
  int imin=1;
  int imax=bins;

  TCanvas *c80 = new TCanvas("c80","<pu> pT akpu3pf",211,393,1271,549);
  c80->Divide(2,1);
  c80->cd(++ipad);
  gPad->SetTopMargin(0.03);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.035);
  gPad->SetBottomMargin(0.15);

  TLegend *l1 = new TLegend(0.1978591,0.6282855,0.4609478,0.922069,"BRNDC");
  l1->SetHeader("PbPb");
  l1->SetFillColor(10);
  l1->SetLineColor(10);
  l1->SetTextFont(42);
  l1->SetTextSize(0.04);

  for(int ic=0;ic<ncen;ic++){

    hpmean_data[1][ic]->SetTitle("");
    hpmean_data[1][ic]->SetMaximum(54);
    hpmean_data[1][ic]->SetMinimum(-0.8);
    hpmean_data[1][ic]->GetXaxis()->SetRange(imin,imax);
    hpmean_data[1][ic]->GetXaxis()->SetNdivisions(509);
    hpmean_data[1][ic]->GetXaxis()->SetMoreLogLabels();
    hpmean_data[1][ic]->GetXaxis()->SetNoExponent();
    hpmean_data[1][ic]->GetXaxis()->SetLabelFont(42);      
    hpmean_data[1][ic]->GetXaxis()->SetLabelSize(0.05);
    hpmean_data[1][ic]->GetXaxis()->SetTitleFont(42);      
    hpmean_data[1][ic]->GetXaxis()->SetTitleSize(0.06);
    hpmean_data[1][ic]->GetXaxis()->SetTitleOffset(1.00);
    hpmean_data[1][ic]->GetXaxis()->SetLabelOffset(0.003);
    hpmean_data[1][ic]->GetXaxis()->SetTitle("p^{Jet}_{T} (GeV/c)");
    hpmean_data[1][ic]->GetXaxis()->CenterTitle(true);
    hpmean_data[1][ic]->GetYaxis()->SetTitle("<Background p_{T}> (GeV/c)");
    hpmean_data[1][ic]->GetYaxis()->CenterTitle(true);
    hpmean_data[1][ic]->GetYaxis()->SetTitleSize(0.06);
    hpmean_data[1][ic]->GetYaxis()->SetTitleFont(42);      
    hpmean_data[1][ic]->GetYaxis()->SetTitleOffset(1.09);
    hpmean_data[1][ic]->GetYaxis()->SetLabelSize(0.05);
    hpmean_data[1][ic]->GetYaxis()->SetLabelFont(42);      
    hpmean_data[1][ic]->GetYaxis()->SetLabelOffset(0.006);    
    
    hpmean_data[1][ic]->SetLineColor(dcol[ic]);
    hpmean_data[1][ic]->SetMarkerColor(dcol[ic]);
    hpmean_data[1][ic]->SetMarkerStyle(dsty[ic]);
    hpmean_data[1][ic]->SetMarkerSize(1.2);
    
    if(ic==0)hpmean_data[1][ic]->Draw("p");
    else hpmean_data[1][ic]->Draw("psame");
    hpmean_mc[1][ic]->SetLineColor(dcol[ic]);
    hpmean_mc[1][ic]->Draw("ehistsame");
    l1->AddEntry(hpmean_data[1][ic],Form("%s ",ccent[ic]),"p");
  }
  l1->AddEntry(hpmean_mc[1][0],"PYTHIA (Z2) + HYDJET 1.8","l");
  l1->Draw();

  gPad->Update();
  TPaveText *pt0 = new TPaveText(0.56,0.75,0.92,0.86,"brNDC");
  pt0->SetBorderSize(0);
  pt0->SetFillColor(10);
  pt0->SetTextFont(42);
  TText *text0 = pt0->AddText("Anti-k_{T} (R = 0.3)");
  text0->SetTextSize(0.045);
  TText *text1 = pt0->AddText("Particle Flow Jets");
  text1->SetTextSize(0.045);
  pt0->Draw();

  TPaveText *pt1   = new TPaveText(0.6303111,0.6811665,0.8194061,0.749716,"brNDC");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(10);
  pt1->SetTextFont(42);
  TText *text2 = pt1->AddText(Form("|#eta_{jet}| < %0.0f",ketacut));  text2->SetTextSize(0.05);
  pt1->Draw();
  
  TPaveText *pt2   = new TPaveText(0.21,0.86,0.34,0.94,"brNDC");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(10);
  pt2->SetTextFont(42);
  //TText *text3 = pt2->AddText("PbPb");  text3->SetTextSize(0.04);
  //pt2->Draw();

  TPaveText *pt3   = new TPaveText(0.559606,0.8711465,0.6895061,0.9416546,"brNDC");
  pt3->SetBorderSize(0);
  pt3->SetFillColor(10);
  pt3->SetTextFont(42);
  TText *text4 = pt3->AddText("CMS");  text4->SetTextSize(0.045);
  pt3->Draw();

  c80->cd(++ipad);
  gPad->SetTopMargin(0.03);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.035);
  gPad->SetBottomMargin(0.15);

  TLegend *l2 = new TLegend(0.2005208,0.6137865,0.4631076,0.9074292,"BRNDC");
  l2->SetHeader("PbPb");
  l2->SetFillColor(10);
  l2->SetLineColor(10);
  l2->SetTextFont(42);
  l2->SetTextSize(0.04);
  
  for(int ic=0;ic<ncen;ic++){  
    hprms_data[1][ic]->SetTitle("");
    hprms_data[1][ic]->SetMaximum(18);
    hprms_data[1][ic]->SetMinimum(-0.2);
    hprms_data[1][ic]->GetXaxis()->SetRange(imin,imax);
    hprms_data[1][ic]->GetXaxis()->SetNdivisions(509);
    hprms_data[1][ic]->GetXaxis()->SetMoreLogLabels();
    hprms_data[1][ic]->GetXaxis()->SetNoExponent();
    hprms_data[1][ic]->GetXaxis()->SetLabelFont(42);      
    hprms_data[1][ic]->GetXaxis()->SetLabelSize(0.05);
    hprms_data[1][ic]->GetXaxis()->SetTitleFont(42);      
    hprms_data[1][ic]->GetXaxis()->SetTitleSize(0.06);
    hprms_data[1][ic]->GetXaxis()->SetTitleOffset(1.00);
    hprms_data[1][ic]->GetXaxis()->SetLabelOffset(0.003);
    hprms_data[1][ic]->GetXaxis()->SetTitle("p^{Jet}_{T} (GeV/c)");
    hprms_data[1][ic]->GetXaxis()->CenterTitle(true);
    hprms_data[1][ic]->GetYaxis()->SetTitle("#sigma( Background p_{T} ) (GeV/c)");
    hprms_data[1][ic]->GetYaxis()->CenterTitle(true);
    hprms_data[1][ic]->GetYaxis()->SetNdivisions(509);
    hprms_data[1][ic]->GetYaxis()->SetTitleSize(0.06);
    hprms_data[1][ic]->GetYaxis()->SetTitleFont(42);      
    hprms_data[1][ic]->GetYaxis()->SetTitleOffset(1.09);
    hprms_data[1][ic]->GetYaxis()->SetLabelSize(0.05);
    hprms_data[1][ic]->GetYaxis()->SetLabelFont(42);      
    hprms_data[1][ic]->GetYaxis()->SetLabelOffset(0.006);    
    
    hprms_data[1][ic]->SetLineColor(dcol[ic]);
    hprms_data[1][ic]->SetMarkerColor(dcol[ic]);
    hprms_data[1][ic]->SetMarkerStyle(dsty[ic]);
    hprms_data[1][ic]->SetMarkerSize(1.2);
    
    if(ic==0)hprms_data[1][ic]->Draw("p");
    else hprms_data[1][ic]->Draw("psame");
    hprms_mc[1][ic]->SetLineColor(dcol[ic]);
    hprms_mc[1][ic]->Draw("ehistsame");
    l2->AddEntry(hprms_data[1][ic],Form("%s ",ccent[ic]),"p");
  }
  l2->AddEntry(hprms_mc[1][0],"PYTHIA (Z2) + HYDJET 1.8","l");
  l2->Draw();

  gPad->Update();

  TPaveText *pt00 = new TPaveText(0.56,0.75,0.92,0.86,"brNDC");
  pt00->SetBorderSize(0);
  pt00->SetFillColor(10);
  pt00->SetTextFont(42);
  TText *text00 = pt00->AddText("Anti-k_{T} (R = 0.3)");
  text00->SetTextSize(0.045);
  TText *text11 = pt00->AddText("Particle Flow Jets");
  text11->SetTextSize(0.045);
  pt00->Draw();

  TPaveText *pt11   = new TPaveText(0.6303111,0.6811665,0.8194061,0.749716,"brNDC");
  pt11->SetBorderSize(0);
  pt11->SetFillColor(10);
  pt11->SetTextFont(42);
  TText *text12 = pt11->AddText(Form("|#eta_{jet}| < %0.0f",ketacut));  text12->SetTextSize(0.05);
  pt11->Draw();
  
  TPaveText *pt22   = new TPaveText(0.21,0.86,0.34,0.94,"brNDC");
  pt22->SetBorderSize(0);
  pt22->SetFillColor(10);
  pt22->SetTextFont(42);
  //TText *text33 = pt22->AddText("PbPb");  text33->SetTextSize(0.055);
  //pt22->Draw();

  TPaveText *pt33   = new TPaveText(0.559606,0.8711465,0.6895061,0.9416546,"brNDC");
  pt33->SetBorderSize(0);
  pt33->SetFillColor(10);
  pt33->SetTextFont(42);
  TText *text44 = pt33->AddText("CMS");  text44->SetTextSize(0.045);
  pt33->Draw();


  if(iSave){
    c80->SaveAs("AN/Background_pT_datamc.C");
    c80->SaveAs("AN/Background_pT_datamc.eps");
    c80->SaveAs("AN/Background_pT_datamc.pdf");
    c80->SaveAs("AN/Background_pT_datamc.gif");
  }


  return 0;
}
void GetPileup(TH2 *h2,TH1 *hmean, TH1*hrms)
{
  for(int ix=1;ix<=h2->GetNbinsX();ix++){
    double binwx = h2->GetXaxis()->GetBinWidth(ix);
    for(int iy=1;iy<=h2->GetNbinsY();iy++){
      double binwy = h2->GetYaxis()->GetBinWidth(iy);
      double binc  = h2->GetBinContent(ix,iy)/binwx/binwy;
      double bine  = h2->GetBinError(ix,iy)/binwx/binwy;
      h2->SetBinContent(ix,iy,binc);
      h2->SetBinError(ix,iy,bine);
    }//! iy

    //! pile up pt distribution
    TH1F *hpu = (TH1F*)h2->ProjectionY("hpu",ix,ix);
    hpu->SetName("hpu");
    hpu->Scale(binwx);
    
    hmean->SetBinContent(ix,hpu->GetMean());
    hmean->SetBinError(ix,hpu->GetMeanError());
      
    hrms->SetBinContent(ix,hpu->GetRMS());
    hrms->SetBinError(ix,hpu->GetRMSError());
    delete hpu;
  }//! ix
}
