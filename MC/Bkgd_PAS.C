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

using namespace std;
double xmin=100.;
double xmax=300.;

void drawText(const char *text, float xp, float yp, int size);
void putCMSPrel(double x, double y, double size);
void drawText(const char *text, float xp, float yp, int size);
void drawText2(const char *text, float xp, float yp, int size);
//TH1F *functionHist(TF1 */*f*/, TH1F* /*h*/,char */*fHistname*/);
void makeMultiPanelCanvasWithGap(TCanvas*& canv,
				 const Int_t columns,
				 const Int_t rows,
				 const Float_t leftOffset,
				 const Float_t bottomOffset,
				 const Float_t leftMargin,
				 const Float_t bottomMargin,
				 const Float_t edge, const Float_t asyoffset); 

void GetPileup(TH2 */*h1*/,TH1 */*hmean*/,TH1 */*hrms*/,int /*iFit*/);
int Bkgd_PAS()
{
  int iFit=0;
  float ketacut=2.0;
  bool iSave=false;

  const int ncen=6;
  const char *ksp  [ncen] = {"pbpb","pbpb" ,"pbpb"  ,"pbpb"  ,"pbpb"  ,"pbpb"  };
  const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%"};
  const char *mcent[ncen] = {"0-5","5-10","10-30","30-50","50-70","70-90"};

  Color_t dcol[ncen] = {kRed,kBlue,kBlack,kAzure+7,kGreen+3,kViolet-5};
  int     dsty[ncen] = {20,24,21,25,33,27};

  //Color_t dcol[ncen] = {kBlack,kBlack,kBlack,kBlack,kBlack,kBlack};
  //Color_t dcol[ncen] = {kRed,kBlue,kRed+1,kBlue+1,kRed+1,kBlue+1};
  //int     dsty[ncen] = {24,20,25,21,27,33};

  //Color_t dcol[ncen] = {kRed+4,kRed+2,kRed,kRed-3,kRed-6,kRed-9};
  //Color_t dcol[ncen] = {kRed-4,kRed+2,kRed-3,kRed+3,kRed-6,kRed+4};
  //int     dsty[ncen] = {20,24,21,25,33,27};

  const float npart [ncen]={381.29,329.41 ,224.28  ,108.12  ,42.04    ,11.43};
  const float enpart[ncen]={2.56  ,2.99   ,5.33    ,5.89    ,4.36     ,1.61 };
  const float rnpart[ncen]={19.2  ,22.5   ,45.9    ,27.1    ,14.4     ,5.73 };

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
  const int nbins = bins;
  cout<<"# of pt bins : "<<nbins<<endl;

  TFile *fin[2];
  //fin[0]  = new TFile("/home/pawan/CMS/Analysis/Raa/Output/OrgUsed/HLT_HIJet80_v1_HiForest2_pbpb_pt80_2012.root","r"); //! data
  fin[0]  = new TFile("/home/pawan/CMS/Analysis/Raa/Output/HLT_HIJet80_v1_HiForest2_pbpb_pt80_2012.root","r"); //! data
  fin[1]  = new TFile("../input/pbpb/PFjets/genmatched/Response_newHiForest_DJ_PFjets_merged_pbpb_2012.root","r");//! MC
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


  //! Background for jets
  TH2F *hjetptbkgd_data[knj][ncen];
  TH1F *hbmean_data[knj][ncen];
  TH1F *hbrms_data [knj][ncen];

  TH2F *hjetptbkgd_mc[knj][ncen];
  TH1F *hbmean_mc[knj][ncen];
  TH1F *hbrms_mc [knj][ncen];

  TH1F *hjet1D_data[knj][ncen][nbins];
  TH1F *hjet1D_mc[knj][ncen][nbins];

  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){
      
      //! PbPb
      //! Pile up subtraction histograms
      hjetptpu_data  [nj][ic] = (TH2F*)fin[0]->Get(Form("hjetptpu%d_%d",nj,ic));
      hjetptpu_data  [nj][ic]->SetName(Form("hjetptpu_%s_%s_%d",ksp[ic],calgo[nj],ic));
      
      //! pileup mean from MC
      hpmean_data[nj][ic] = new TH1F(Form("hpmean_data_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hpmean_data_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hpmean_data[nj][ic]->Sumw2();
      hprms_data [nj][ic] = new TH1F(Form("hprms_data_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hprms_data_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hprms_data [nj][ic]->Sumw2();

      hjetptpu_mc  [nj][ic] = (TH2F*)fin[1]->Get(Form("hjetptpu_genm%d_%d",nj+4,ic));
      hjetptpu_mc  [nj][ic]->SetName(Form("hjetptpu_%s_%s_%d",ksp[ic],calgo[nj],ic));

      //! pileup mean from MC
      hpmean_mc[nj][ic] = new TH1F(Form("hpmean_mc_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hpmean_mc_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hpmean_mc[nj][ic]->Sumw2();
      hprms_mc[nj][ic] = new TH1F(Form("hprms_mc_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hprms_mc_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hprms_mc[nj][ic]->Sumw2();

      GetPileup(hjetptpu_data[nj][ic],hpmean_data[nj][ic],hprms_data [nj][ic],0);
      GetPileup(hjetptpu_mc[nj][ic],hpmean_mc[nj][ic],hprms_mc[nj][ic],0);

      //! Jet background
      //! PbPb data
      hjetptbkgd_data  [nj][ic] = (TH2F*)fin[0]->Get(Form("hjetptbkgd%d_%d",nj,ic));
      hjetptbkgd_data  [nj][ic]->SetName(Form("hjetptbkgd_%s_%s_%d",ksp[ic],calgo[nj],ic));

      //! pileup mean from data
      hbmean_data[nj][ic] = new TH1F(Form("hbmean_data_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hbmean_data_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hbmean_data[nj][ic]->Sumw2();
      hbrms_data [nj][ic] = new TH1F(Form("hbrms_data_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hbrms_data_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hbrms_data [nj][ic]->Sumw2();


      //! PbPb MC
      hjetptbkgd_mc  [nj][ic] = (TH2F*)fin[1]->Get(Form("hjetptbkgd_genm%d_%d",nj+4,ic));
      hjetptbkgd_mc  [nj][ic]->SetName(Form("hjetptbkgd_%s_%s_%d",ksp[ic],calgo[nj],ic));

      //! pileup mean from data
      hbmean_mc[nj][ic] = new TH1F(Form("hbmean_mc_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hbmean_mc_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hbmean_mc[nj][ic]->Sumw2();
      hbrms_mc [nj][ic] = new TH1F(Form("hbrms_mc_%s_%s_%d",ksp[ic],calgo[nj],ic),Form("hbrms_mc_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic]),bins,ptbins);
      hbrms_mc [nj][ic]->Sumw2();

      GetPileup(hjetptbkgd_data[nj][ic],hbmean_data[nj][ic],hbrms_data[nj][ic],iFit);
      GetPileup(hjetptbkgd_mc  [nj][ic],hbmean_mc  [nj][ic],hbrms_mc  [nj][ic],iFit);

      for(int ix=1;ix<=hjetptbkgd_data[nj][ic]->GetNbinsX();ix++){
	hjet1D_data[nj][ic][ix-1] = (TH1F*)hjetptbkgd_data[nj][ic]->ProjectionY(Form("hjet1D_data_%d_%d_%d",nj,ic,ix-1),ix,ix);
	hjet1D_data[nj][ic][ix-1]->Scale(1./hjet1D_data[nj][ic][ix-1]->Integral());
	hjet1D_data[nj][ic][ix-1]->SetMarkerColor(1);
	hjet1D_data[nj][ic][ix-1]->SetLineColor(1);
	hjet1D_data[nj][ic][ix-1]->SetMarkerStyle(20);

	hjet1D_mc[nj][ic][ix-1]   = (TH1F*)hjetptbkgd_mc[nj][ic]  ->ProjectionY(Form("hjet1D_mc_%d_%d_%d",nj,ic,ix-1),ix,ix);
	hjet1D_mc[nj][ic][ix-1]->Scale(1./hjet1D_mc[nj][ic][ix-1]->Integral());
	hjet1D_mc[nj][ic][ix-1]->SetLineColor(1);
	hjet1D_mc[nj][ic][ix-1]->SetLineStyle(2);

	//if(nj==1)cout<<"\t ic : "<<ic<<"\t pT : "<< hjetptbkgd_mc[nj][ic]->GetBinCenter(ix)<<"\t RMS : "<<hjet1D_mc[nj][ic][ix-1]->GetRMS()<<"\t RMS Error : "<<hjet1D_mc[nj][ic][ix-1]->GetRMSError()<<endl;
      }

      //! eta dependence of pile up subtraction
      for(int ie=0;ie<2;ie++){
	hjetptpu_etab_data  [nj][ic][ie] = (TH2F*)fin[0]->Get(Form("hjetptpu_etab%d_%d_%d",nj,ic,ie));
	hjetptpu_etab_data  [nj][ic][ie]->SetName(Form("hjetptpu_etab_data_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie));

	//! pileup mean 
	hpmean_etab_data[nj][ic][ie] = new TH1F(Form("hpmean_etab_data_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie),Form("hpmean_etab_data_%s_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic],meta[ie]),bins,ptbins);
	hpmean_etab_data[nj][ic][ie]->Sumw2();
	hprms_etab_data[nj][ic][ie] = new TH1F(Form("hprms_etab_data_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie),Form("hprms_etab_data_%s_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic],meta[ie]),bins,ptbins);
	hprms_etab_data[nj][ic][ie]->Sumw2();

	GetPileup(hjetptpu_etab_data  [nj][ic][ie],hpmean_etab_data[nj][ic][ie],hprms_etab_data[nj][ic][ie],0);
	
	hjetptpu_etab_mc  [nj][ic][ie] = (TH2F*)fin[1]->Get(Form("hjetptpu_etab_genm%d_%d_%d",nj+4,ic,ie));
	hjetptpu_etab_mc  [nj][ic][ie]->SetName(Form("hjetptpu_etab_mc_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie));

	//! pileup mean 
	hpmean_etab_mc[nj][ic][ie] = new TH1F(Form("hpmean_etab_mc_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie),Form("hpmean_etab_mc_%s_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic],meta[ie]),bins,ptbins);
	hpmean_etab_mc[nj][ic][ie]->Sumw2();
	hprms_etab_mc[nj][ic][ie] = new TH1F(Form("hprms_etab_mc_%s_%s_%s_%d",ksp[ic],calgo[nj],ccent[ic],ie),Form("hprms_etab_mc_%s_%s_%s_%s",ksp[ic],calgo[nj],ccent[ic],meta[ie]),bins,ptbins);
	hprms_etab_mc[nj][ic][ie]->Sumw2();

	GetPileup(hjetptpu_etab_mc[nj][ic][ie],hpmean_etab_mc[nj][ic][ie],hprms_etab_mc[nj][ic][ie],0);
      }
    }
  }

  TH1F *hBkgd_mc[nbins], *hBkgd_data[nbins];
  //! as a function of NPart
  TGraphErrors *grBkgdNP_mc[nbins], *grBkgdNP_data[nbins];

  for(int ix=0;ix<nbins;ix++){
    hBkgd_mc[ix] = new TH1F(Form("hBkgd_mc_%d",ix),"",ncen,-0.5,ncen-0.5);
    hBkgd_mc[ix]->SetMarkerColor(1);
    hBkgd_mc[ix]->SetLineColor(1);
    hBkgd_mc[ix]->SetMarkerStyle(24);
    hBkgd_mc[ix]->SetMarkerSize(1.2);


    hBkgd_data[ix] = new TH1F(Form("hBkgd_data_%d",ix),"",ncen,-0.5,ncen-0.5);
    hBkgd_data[ix]->SetMarkerColor(1);
    hBkgd_data[ix]->SetLineColor(1);
    hBkgd_data[ix]->SetMarkerStyle(20);
    hBkgd_data[ix]->SetMarkerSize(1.2);
    //hBkgd_data[ix]->GetXaxis()->LabelsOption("d");
    for(int ic=0;ic<ncen;ic++){
      hBkgd_mc  [ix]->GetXaxis()->SetBinLabel(ic+1,ccent[(ncen-ic)-1]);
      hBkgd_data[ix]->GetXaxis()->SetBinLabel(ic+1,ccent[(ncen-ic)-1]);
    }


  }
  //hBkgd[0]->Draw();
  //return 0;

  for(int ic=0;ic<ncen;ic++){
    for(int ix=1;ix<=hbrms_mc[1][ic]->GetNbinsX();ix++){
      hBkgd_mc[ix-1]->SetBinContent((ncen-ic),hbrms_mc[1][ic]->GetBinContent(ix));
      hBkgd_mc[ix-1]->SetBinError  ((ncen-ic),hbrms_mc[1][ic]->GetBinError(ix));
      hBkgd_data[ix-1]->SetBinContent((ncen-ic),hbrms_data[1][ic]->GetBinContent(ix));
      hBkgd_data[ix-1]->SetBinError  ((ncen-ic),hbrms_data[1][ic]->GetBinError(ix));
    }
  }


  for(int ic=ncen-1;ic>=0;ic--){
    cout<<ccent[ic]<<endl;
    for(int ix=1;ix<=hbrms_data[1][ic]->GetNbinsX();ix++){
      if(ic==1 || ic==2)cout<<"pT : "<<hbrms_data[1][ic]->GetBinLowEdge(ix)<<"\t RMS : "<<hbrms_data[1][ic]->GetBinContent(ix)<<"\t RMS  ERROR : "<<hbrms_data[1][ic]->GetBinError(ix)<<"\t Ind RMS : "<<hjet1D_data[1][ic][ix-1]->GetRMS()<<"\t Ind RMS Errors : "<<hjet1D_data[1][ic][ix-1]->GetRMSError()<<endl;
      //if(ic==1 || ic==2)cout<<"pT : "<<hbrms_mc[1][ic]->GetBinLowEdge(ix)<<"\t RMS : "<<hbrms_mc[1][ic]->GetBinContent(ix)<<"\t RMS  ERROR : "<<hbrms_mc[1][ic]->GetBinError(ix)<<"\t Ind RMS : "<<hjet1D_mc[1][ic][ix-1]->GetRMS()<<"\t Ind RMS Errors : "<<hjet1D_mc[1][ic][ix-1]->GetRMSError()<<endl;
    }
    cout<<endl;
  }

  for(int ix=0;ix<nbins;ix++){
    grBkgdNP_mc[ix]   = new TGraphErrors(hBkgd_mc[ix]);
    grBkgdNP_data[ix] = new TGraphErrors(hBkgd_data[ix]);
    
    for(int ic=0;ic<ncen;ic++){
      grBkgdNP_mc[ix]->SetPoint     (ic,npart [ic],hBkgd_mc[ix]->GetBinContent(ncen-ic));
      grBkgdNP_mc[ix]->SetPointError(ic,enpart[ic],hBkgd_mc[ix]->GetBinError(ncen-ic));

      grBkgdNP_data[ix]->SetPoint(ic,npart[ic],hBkgd_data[ix]->GetBinContent(ncen-ic));
      grBkgdNP_data[ix]->SetPointError(ic,enpart[ic],hBkgd_data[ix]->GetBinError(ncen-ic));
    }
  }
  

  const char *cpad[] = {"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)"};

  int ipad=0;
  int imin=1;
  int imax=bins;

  TCanvas *c84  = new TCanvas("c84","Background jet pt RMS vs Centrality",1100,770);
  makeMultiPanelCanvasWithGap(c84,3,2,0.01,0.01,0.16,0.2,0.04,0.04); 
  ipad=0;

  //TLegend *l09 = new TLegend(0.1802593,0.645195,0.5831155,0.8404753,"BRNDC");  //! centrality
  TLegend *l09 = new TLegend(0.4730181,0.2326311,0.8758743,0.4279113,NULL,"brNDC"); //! NPart
  l09->SetTextSize(0.05);
  l09->SetHeader("");
  l09->SetFillColor(10);
  l09->SetLineColor(10);
  l09->SetTextFont(42);
  l09->AddEntry(hBkgd_data[0],"PbPb Data","p");  
  l09->AddEntry(hBkgd_mc  [0],"PYTHIA+HYDJET","p");  

  for(int ix=0;ix<6;ix++){
    c84->cd(++ipad);
    //hDum0->Draw("");
    hBkgd_data[ix]->SetMaximum(12.45);     
    hBkgd_data[ix]->SetMinimum(0);     
    //hBkgd_data[ix]->GetYaxis()->SetNdivisions(610); 
    hBkgd_data[ix]->GetYaxis()->SetTitle("RMS ( Background p_{T} ) (GeV/c)");
    hBkgd_data[ix]->GetYaxis()->SetLabelFont(43);
    hBkgd_data[ix]->GetYaxis()->SetTitleFont(43);
    hBkgd_data[ix]->GetYaxis()->SetLabelSize(20);
    hBkgd_data[ix]->GetYaxis()->SetTitleSize(22);
    hBkgd_data[ix]->GetYaxis()->SetTitleOffset(2.6);

    hBkgd_data[ix]->GetXaxis()->SetTitle("Centrality");
    hBkgd_data[ix]->GetXaxis()->SetLabelFont(43);
    hBkgd_data[ix]->GetXaxis()->SetTitleFont(43);
    hBkgd_data[ix]->GetXaxis()->SetLabelSize(17);
    hBkgd_data[ix]->GetXaxis()->SetTitleSize(22);
    hBkgd_data[ix]->GetXaxis()->SetTitleOffset(3.1);
    hBkgd_data[ix]->GetXaxis()->SetNoExponent();
    hBkgd_data[ix]->GetXaxis()->SetMoreLogLabels();
    hBkgd_data[ix]->GetXaxis()->SetTitleOffset(3.5);
    hBkgd_data[ix]->GetXaxis()->LabelsOption("d");  
    //hBkgd_data[ix]->Draw("p");
    //hBkgd_mc  [ix]->Draw("psame");
    


    grBkgdNP_data[ix]->GetYaxis()->SetTitle("RMS ( Background p_{T} ) (GeV/c)");
    grBkgdNP_data[ix]->GetYaxis()->SetLabelFont(43);
    grBkgdNP_data[ix]->GetYaxis()->SetTitleFont(43);
    grBkgdNP_data[ix]->GetYaxis()->SetLabelSize(20);
    grBkgdNP_data[ix]->GetYaxis()->SetTitleSize(22);
    grBkgdNP_data[ix]->GetYaxis()->SetTitleOffset(2.6);
    grBkgdNP_data[ix]->GetXaxis()->SetTitle("N_{part}");
    grBkgdNP_data[ix]->GetXaxis()->SetLabelFont(43);
    grBkgdNP_data[ix]->GetXaxis()->SetTitleFont(43);
    grBkgdNP_data[ix]->GetXaxis()->SetLabelSize(17);
    grBkgdNP_data[ix]->GetXaxis()->SetTitleSize(22);
    grBkgdNP_data[ix]->GetXaxis()->SetTitleOffset(3.1);
    grBkgdNP_data[ix]->GetXaxis()->SetNoExponent();
    grBkgdNP_data[ix]->GetXaxis()->SetMoreLogLabels();
    grBkgdNP_data[ix]->GetXaxis()->SetTitleOffset(2.1);
    grBkgdNP_data[ix]->SetMaximum(12.45);     
    grBkgdNP_data[ix]->SetMinimum(0);     

    grBkgdNP_data[ix]->Draw("ap");
    grBkgdNP_mc  [ix]->Draw("psame");

    if(ipad==1 || ipad==4){
    }
    else{
      grBkgdNP_data[ix]->GetYaxis()->SetLabelSize(0);
      grBkgdNP_data[ix]->GetYaxis()->SetTitleSize(0);
      grBkgdNP_data[ix]->GetYaxis()->SetTitle("");
    }    
    gPad->Update();
    drawText2(Form("%0.0f< p^{Jet}_{T} (GeV/c) < %0.0f",ptbins[ix],ptbins[ix+1]),0.23,0.86,18);  
    if(ipad==1)l09->Draw();
  }
  //return 0;

  
  //! PAS plot
  TH1D *hDum, *hDum_1;   
  float rymin = 0.00001;
  float rymax = 20.386;
  hDum = GetDummyHist(0.0,155,rymin,rymax,"< Background p_{T} > (GeV/c)","Normalized Yield (a. u.)",false);
  hDum->GetYaxis()->SetNdivisions(610); 
  hDum->GetYaxis()->SetLabelFont(43);
  hDum->GetYaxis()->SetTitleFont(43);
  hDum->GetYaxis()->SetLabelSize(20);
  hDum->GetYaxis()->SetTitleSize(22);
  hDum->GetYaxis()->SetTitleOffset(2.6);
  hDum->GetXaxis()->SetLabelFont(43);
  hDum->GetXaxis()->SetTitleFont(43);
  hDum->GetXaxis()->SetLabelSize(20);
  hDum->GetXaxis()->SetTitleSize(22);
  hDum->GetXaxis()->SetTitleOffset(3.1);
  hDum->GetXaxis()->SetNoExponent();
  hDum->GetXaxis()->SetMoreLogLabels();

  // to show in up and down
  hDum->GetXaxis()->SetTitleOffset(2.4);
  hDum_1 = (TH1D*) hDum->Clone("hDum_1");
  hDum_1->SetAxisRange(0.00001,20.386,"Y");
  ipad=0;
  TCanvas *c83  = new TCanvas("c83","Background jet pt",1100,770);
  makeMultiPanelCanvasWithGap(c83,3,2,0.01,0.01,0.16,0.2,0.04,0.04); 

  TLegend *l07 = new TLegend(0.1799019,0.6671984,0.5820431,0.9174872,"BRNDC");
  l07->SetTextSize(0.08);
  l07->SetHeader("");
  l07->SetFillColor(10);
  l07->SetLineColor(10);
  l07->SetTextFont(42);
  TH1F *hDumCl[ncen];
  for(int ic=ncen-1;ic>=0;ic--){
    c83->cd(++ipad);
    gPad->SetLogy();
    //hDum->Draw("");
    hDumCl[ic] = (TH1F*)hDum->Clone(Form("hDumCl%d",ic));

    if(ipad==1 || ipad==4){
    }else{
      hDumCl[ic]->GetYaxis()->SetLabelSize(0);
      hDumCl[ic]->GetYaxis()->SetTitleSize(0);
    }

    hDumCl[ic]->Draw("");
    drawText(ccent[ic],0.648173,0.8459761,22);

    if(ipad==1){
      drawText("CMS Preliminary",0.23,0.80,19);  
      drawText2("Anti-k_{T} Particle Flow Jets R = 0.3",0.23,0.73,19);  
      drawText2("100 < p^{Jet}_{T} (GeV/c) < 110, |#eta_{Jet}| < 2",0.23,0.66,18);  
    }

    if(ipad==3)drawText("#int L^{PbPb} dt = 129 #mub^{-1}",0.15,0.75,21);
    if(ipad==1 || ipad==4)drawText(cpad[ipad-1],0.2,0.86,24);
    else drawText(cpad[ipad-1],0.08,0.86,24);

    hjet1D_data[1][ic][0]->SetMarkerStyle(20);
    hjet1D_data[1][ic][0]->SetMarkerSize(1.2);
    hjet1D_data[1][ic][0]->SetMarkerColor(1);
    hjet1D_data[1][ic][0]->SetLineColor(1);
    hjet1D_data[1][ic][0]->Draw("psame");
    

    hjet1D_mc[1][ic][0]->SetFillStyle(3004);
    hjet1D_mc[1][ic][0]->SetFillColor(2);
    hjet1D_mc[1][ic][0]->SetLineColor(2);
    hjet1D_mc[1][ic][0]->Draw("histsame");
    gPad->Update();
  }
  c83->cd(2);
  l07->AddEntry(hjet1D_data[1][0][0],"PbPb Data","p");  
  l07->AddEntry(hjet1D_mc  [1][0][0],"PYTHIA+HYDJET","f");  
  l07->Draw();
  c83->Update();
  c83->SaveAs("Plots/BackgroundpT_DataMc_PAS.C");
  c83->SaveAs("Plots/BackgroundpT_DataMc_PAS.eps");
  c83->SaveAs("Plots/BackgroundpT_DataMc_PAS.pdf");
  c83->SaveAs("Plots/BackgroundpT_DataMc_PAS.gif");
  //return 0;


  // Jet background 
  ipad=0;
  imin=1;
  imax=bins;

  TCanvas *c81 = new TCanvas("c81","<bkgd> pT akpu3pf",211,393,1271,549);
  c81->Divide(2,1);
  c81->cd(++ipad);
  gPad->SetTopMargin(0.03);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.035);
  gPad->SetBottomMargin(0.15);
  
  TLegend *l3 = new TLegend(0.1879933,0.5499432,0.451082,0.8574366,"BRNDC");
  l3->SetHeader("PbPb data");
  l3->SetFillColor(10);
  l3->SetLineColor(10);
  l3->SetTextFont(42);
  l3->SetTextSize(0.04);

  for(int ic=0;ic<ncen;ic++){

    hbmean_data[1][ic]->SetTitle("");
    //hbmean_data[1][ic]->SetMaximum(94);

    hbmean_data[1][ic]->SetMaximum(114);
    hbmean_data[1][ic]->SetMinimum(-0.8);

    hbmean_data[1][ic]->GetXaxis()->SetRange(imin,imax);
    hbmean_data[1][ic]->GetXaxis()->SetNdivisions(509);
    hbmean_data[1][ic]->GetXaxis()->SetMoreLogLabels();
    hbmean_data[1][ic]->GetXaxis()->SetNoExponent();
    hbmean_data[1][ic]->GetXaxis()->SetLabelFont(42);      
    hbmean_data[1][ic]->GetXaxis()->SetLabelSize(0.05);
    hbmean_data[1][ic]->GetXaxis()->SetTitleFont(42);      
    hbmean_data[1][ic]->GetXaxis()->SetTitleSize(0.06);
    hbmean_data[1][ic]->GetXaxis()->SetTitleOffset(1.00);
    hbmean_data[1][ic]->GetXaxis()->SetLabelOffset(0.003);
    hbmean_data[1][ic]->GetXaxis()->SetTitle("p^{Jet}_{T} (GeV/c)");
    hbmean_data[1][ic]->GetXaxis()->CenterTitle(true);
    hbmean_data[1][ic]->GetYaxis()->SetTitle("< Background p_{T} > (GeV/c)");
    hbmean_data[1][ic]->GetYaxis()->CenterTitle(true);
    hbmean_data[1][ic]->GetYaxis()->SetTitleSize(0.06);
    hbmean_data[1][ic]->GetYaxis()->SetTitleFont(42);      
    hbmean_data[1][ic]->GetYaxis()->SetTitleOffset(1.09);
    hbmean_data[1][ic]->GetYaxis()->SetLabelSize(0.05);
    hbmean_data[1][ic]->GetYaxis()->SetLabelFont(42);      
    hbmean_data[1][ic]->GetYaxis()->SetLabelOffset(0.006);    
    
    hbmean_data[1][ic]->SetLineColor(dcol[ic]);
    hbmean_data[1][ic]->SetMarkerColor(dcol[ic]);
    hbmean_data[1][ic]->SetMarkerStyle(dsty[ic]);
    if(ic==4 || ic==5)hbmean_data[1][ic]->SetMarkerSize(2.0);
    else hbmean_data[1][ic]->SetMarkerSize(1.2);

    if(ic==0)hbmean_data[1][ic]->Draw("p");
    else hbmean_data[1][ic]->Draw("psame");
    hbmean_mc[1][ic]->SetLineColor(dcol[ic]);
    hbmean_mc[1][ic]->Draw("histsame");
    l3->AddEntry(hbmean_data[1][ic],Form("%s ",ccent[ic]),"p");
  }
  l3->AddEntry(hbmean_mc[1][0],"PYTHIA + HYDJET","l");
  l3->Draw();

  gPad->Update();
  TPaveText *pt90 = new TPaveText(0.4823237,0.7830115,0.8424263,0.8926907,"brNDC");
  pt90->SetBorderSize(0);
  pt90->SetFillColor(10);
  pt90->SetTextFont(42);
  TText *text90 = pt90->AddText("Anti-k_{T} Particle Flow Jets R = 0.3");
  text90->SetTextSize(0.045);
  //TText *text91 = pt90->AddText("Particle Flow Jets");
  //text91->SetTextSize(0.045);
  pt90->Draw();

  TPaveText *pt81   = new TPaveText(0.5382301,0.7438403,0.727325,0.8123898,"brNDC");
  pt81->SetBorderSize(0);
  pt81->SetFillColor(10);
  pt81->SetTextFont(42);
  TText *text82 = pt81->AddText(Form("|#eta_{Jet}| < %0.0f",ketacut));  text82->SetTextSize(0.05);
  pt81->Draw();
  
  TPaveText *pt92   = new TPaveText(0.21,0.86,0.34,0.94,"brNDC");
  pt92->SetBorderSize(0);
  pt92->SetFillColor(10);
  pt92->SetTextFont(42);
  //TText *text93 = pt92->AddText("PbPb");  text93->SetTextSize(0.04);
  //pt92->Draw();

  TPaveText *pt93   = new TPaveText(0.559606,0.8711465,0.6895061,0.9416546,"brNDC");
  pt93->SetBorderSize(0);
  pt93->SetFillColor(10);
  pt93->SetTextFont(62);
  TText *text94 = pt93->AddText("CMS Preliminary");  text94->SetTextSize(0.045);
  pt93->Draw();

  drawText(cpad[ipad-1],0.2,0.90,24);
  drawText("#int L^{PbPb} dt = 129 #mub^{-1}",0.56,0.68,20);
  gPad->Update();

  c81->cd(++ipad);
  gPad->SetTopMargin(0.03);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.035);
  gPad->SetBottomMargin(0.15);

  TLegend *l4 = new TLegend(0.1921041,0.5303576,0.4535484,0.851561,"BRNDC");
  l4->SetHeader("PbPb data");
  l4->SetFillColor(10);
  l4->SetLineColor(10);
  l4->SetTextFont(42);
  l4->SetTextSize(0.04);
  
  for(int ic=0;ic<ncen;ic++){  
    hbrms_data[1][ic]->SetTitle("");
    //hbrms_data[1][ic]->SetMaximum(21);
    hbrms_data[1][ic]->SetMaximum(24);
    hbrms_data[1][ic]->SetMinimum(-0.2);
    hbrms_data[1][ic]->GetXaxis()->SetRange(imin,imax);
    hbrms_data[1][ic]->GetXaxis()->SetNdivisions(509);
    hbrms_data[1][ic]->GetXaxis()->SetMoreLogLabels();
    hbrms_data[1][ic]->GetXaxis()->SetNoExponent();
    hbrms_data[1][ic]->GetXaxis()->SetLabelFont(42);      
    hbrms_data[1][ic]->GetXaxis()->SetLabelSize(0.05);
    hbrms_data[1][ic]->GetXaxis()->SetTitleFont(42);      
    hbrms_data[1][ic]->GetXaxis()->SetTitleSize(0.06);
    hbrms_data[1][ic]->GetXaxis()->SetTitleOffset(1.00);
    hbrms_data[1][ic]->GetXaxis()->SetLabelOffset(0.003);
    hbrms_data[1][ic]->GetXaxis()->SetTitle("p^{Jet}_{T} (GeV/c)");
    hbrms_data[1][ic]->GetXaxis()->CenterTitle(true);
    hbrms_data[1][ic]->GetYaxis()->SetTitle("RMS ( Background p_{T} ) (GeV/c)");
    hbrms_data[1][ic]->GetYaxis()->CenterTitle(true);
    hbrms_data[1][ic]->GetYaxis()->SetNdivisions(509);
    hbrms_data[1][ic]->GetYaxis()->SetTitleSize(0.06);
    hbrms_data[1][ic]->GetYaxis()->SetTitleFont(42);      
    hbrms_data[1][ic]->GetYaxis()->SetTitleOffset(1.09);
    hbrms_data[1][ic]->GetYaxis()->SetLabelSize(0.05);
    hbrms_data[1][ic]->GetYaxis()->SetLabelFont(42);      
    hbrms_data[1][ic]->GetYaxis()->SetLabelOffset(0.006);    
    
    hbrms_data[1][ic]->SetLineColor(dcol[ic]);
    hbrms_data[1][ic]->SetMarkerColor(dcol[ic]);
    hbrms_data[1][ic]->SetMarkerStyle(dsty[ic]);
    //hbrms_data[1][ic]->SetMarkerSize(1.2);
    if(ic==4 || ic==5)hbrms_data[1][ic]->SetMarkerSize(2.0);
    else hbrms_data[1][ic]->SetMarkerSize(1.2);
    
    if(ic==0)hbrms_data[1][ic]->Draw("p");
    else hbrms_data[1][ic]->Draw("psame");
    hbrms_mc[1][ic]->SetLineColor(dcol[ic]);
    hbrms_mc[1][ic]->Draw("histsame");
    l4->AddEntry(hbrms_data[1][ic],Form("%s ",ccent[ic]),"p");
  }
  l4->AddEntry(hbrms_mc[1][0],"PYTHIA + HYDJET","l");
  l4->Draw();
  gPad->Update();


  TPaveText *pt900 = new TPaveText(0.4823237,0.7830115,0.8424263,0.8926907,"brNDC");
  pt900->SetBorderSize(0);
  pt900->SetFillColor(10);
  pt900->SetTextFont(42);
  TText *text900 = pt900->AddText("Anti-k_{T} Particle Flow Jets R = 0.3");
  text900->SetTextSize(0.045);
  //TText *text911 = pt900->AddText("Particle Flow Jets");
  //text911->SetTextSize(0.045);
  pt900->Draw();

  TPaveText *pt911   = new TPaveText(0.5382301,0.7438403,0.727325,0.8123898,"brNDC");
  pt911->SetBorderSize(0);
  pt911->SetFillColor(10);
  pt911->SetTextFont(42);
  TText *text912 = pt911->AddText(Form("|#eta_{Jet}| < %0.0f",ketacut));  text912->SetTextSize(0.05);
  pt911->Draw();
  
  TPaveText *pt922   = new TPaveText(0.21,0.86,0.34,0.94,"brNDC");
  pt922->SetBorderSize(0);
  pt922->SetFillColor(10);
  pt922->SetTextFont(42);
  //TText *text933 = pt922->AddText("PbPb");  text933->SetTextSize(0.055);
  //pt922->Draw();

  TPaveText *pt933   = new TPaveText(0.559606,0.8711465,0.6895061,0.9416546,"brNDC");
  pt933->SetBorderSize(0);
  pt933->SetFillColor(10);
  pt933->SetTextFont(62);
  TText *text944 = pt933->AddText("CMS Preliminary");  text944->SetTextSize(0.045);
  pt933->Draw();

  drawText(cpad[ipad-1],0.18,0.90,24);
  drawText("#int L^{PbPb} dt = 129 #mub^{-1}",0.56,0.68,20);
  gStyle->SetErrorX(0);
  gPad->Update();

  c81->SaveAs("Plots/Background_pT_SumRaw_datamc.C");
  c81->SaveAs("Plots/Background_pT_SumRaw_datamc.eps");
  c81->SaveAs("Plots/Background_pT_SumRaw_datamc.pdf");
  c81->SaveAs("Plots/Background_pT_SumRaw_datamc.gif");

  //! Re arrange the plot
  return 0;



  TCanvas *c82[ncen];
  TLegend *l82 = new TLegend(0.1835252,0.7500871,0.4461073,0.8894808,"BRNDC");
  l82->SetHeader("PbPb");
  l82->SetFillColor(10);
  l82->SetLineColor(10);
  l82->SetTextFont(42);
  l82->SetTextSize(0.05);
  l82->AddEntry(hjet1D_data[1][0][0],"Data","p");
  l82->AddEntry(hjet1D_mc[1][0][0],"MC","l");

  for(int ic=0;ic<ncen;ic++){
    c82[ic] = new TCanvas(Form("c82_%d",ic),Form("%s Ind akpu3pf",ccent[ic]),158,52,1351,664);
    c82[ic]->Divide(4,2,0,0);
    ipad=0;
    cout<<ccent[ic]<<endl;
    for(int ip=0;ip<nbins-2;ip++){
      c82[ic]->cd(++ipad);
      gPad->SetLogy();

      if(ipad%4==0)gPad->SetRightMargin(0.03);
      //gPad->SetTopMargin(0.03);
      //gPad->SetLeftMargin(0.15);
      //gPad->SetRightMargin(0.035);
      //gPad->SetBottomMargin(0.15);

      hjet1D_data[1][ic][ip]->GetXaxis()->SetNdivisions(509);
      hjet1D_data[1][ic][ip]->GetXaxis()->SetMoreLogLabels();
      hjet1D_data[1][ic][ip]->GetXaxis()->SetNoExponent();
      hjet1D_data[1][ic][ip]->GetXaxis()->SetLabelFont(42);      
      hjet1D_data[1][ic][ip]->GetXaxis()->SetLabelSize(0.05);
      hjet1D_data[1][ic][ip]->GetXaxis()->SetTitleFont(42);      
      hjet1D_data[1][ic][ip]->GetXaxis()->SetTitleSize(0.05);
      hjet1D_data[1][ic][ip]->GetXaxis()->SetTitleOffset(0.93);
      hjet1D_data[1][ic][ip]->GetXaxis()->SetLabelOffset(0.003);
      hjet1D_data[1][ic][ip]->GetXaxis()->SetTitle("<Background p_{T}> (GeV/c)");
      hjet1D_data[1][ic][ip]->GetXaxis()->CenterTitle(true);
      hjet1D_data[1][ic][ip]->GetYaxis()->SetTitle("");
      hjet1D_data[1][ic][ip]->GetYaxis()->CenterTitle(true);
      hjet1D_data[1][ic][ip]->GetYaxis()->SetTitleSize(0.06);
      hjet1D_data[1][ic][ip]->GetYaxis()->SetTitleFont(42);      
      hjet1D_data[1][ic][ip]->GetYaxis()->SetTitleOffset(1.09);
      hjet1D_data[1][ic][ip]->GetYaxis()->SetLabelSize(0.05);
      hjet1D_data[1][ic][ip]->GetYaxis()->SetLabelFont(42);      
      hjet1D_data[1][ic][ip]->GetYaxis()->SetLabelOffset(0.006);    

      hjet1D_data[1][ic][ip]->GetXaxis()->SetRangeUser(0,120);
      hjet1D_data[1][ic][ip]->SetTitle("");
      hjet1D_data[1][ic][ip]->SetMaximum(20.386);
      hjet1D_data[1][ic][ip]->SetMinimum(0.00005678);

      hjet1D_data[1][ic][ip]->Fit("gaus","RQ0","",0,120);
      TF1 *fd = (TF1*)hjet1D_data[1][ic][ip]->GetFunction("gaus");
      double arm_data   = hjet1D_data[1][ic][ip]->GetMean();
      double mean_data  = fd->GetParameter(1);
      double rms_data   = hjet1D_data[1][ic][ip]->GetRMS();
      double sigma_data = fd->GetParameter(2);
      delete fd;

      hjet1D_data[1][ic][ip]->Draw("p");

      hjet1D_mc[1][ic][ip]->SetFillStyle(3004);
      hjet1D_mc[1][ic][ip]->SetFillColor(4);
      hjet1D_mc[1][ic][ip]->SetLineColor(4);

      hjet1D_mc[1][ic][ip]->Fit("gaus","RQ0","",0,120);
      TF1 *fm = (TF1*)hjet1D_mc[1][ic][ip]->GetFunction("gaus");
      double arm_mc   = hjet1D_mc[1][ic][ip]->GetMean();
      double mean_mc  = fm->GetParameter(1);
      double rms_mc   = hjet1D_mc[1][ic][ip]->GetRMS();
      double sigma_mc = fm->GetParameter(2);
      delete fm;
      hjet1D_mc[1][ic][ip]->Draw("histsame");
      
      //cout<<"\t \t Data :  "<<ptbins[ip]<<"\t ArM : "<<arm_data<<"\t Mean : "<<mean_data<<"\t  RMS : "<<rms_data<<"\t  Sigma : "<<sigma_data<<endl
      //  <<"\t \t MC   :  "<<ptbins[ip]<<"\t ArM : "<<arm_mc  <<"\t Mean : "<<mean_mc  <<"\t  RMS : "<<rms_mc  <<"\t  Sigma : "<<sigma_mc<<endl
      //  <<endl;

      if(ipad==1){
	TPaveText *pt1   = new TPaveText(0.6303111,0.6811665,0.8194061,0.749716,"brNDC");
	pt1->SetBorderSize(0);
	pt1->SetFillColor(10);
	pt1->SetTextFont(42);
	TText *text2 = pt1->AddText(Form("%s",ccent[ic]));  text2->SetTextSize(0.07);
	pt1->Draw();
	l82->Draw();

	TPaveText *pt3   = new TPaveText(0.559606,0.8711465,0.6895061,0.9416546,"brNDC");
	pt3->SetBorderSize(0);
	pt3->SetFillColor(10);
	pt3->SetTextFont(62);
	TText *text4 = pt3->AddText("CMS Preliminary");  text4->SetTextSize(0.045);
	pt3->Draw();

      }

      TPaveText *pt   = new TPaveText(0.2936403,0.8511475,0.7086895,0.9661473,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(10);
      pt->SetTextFont(42);
      TText *text = pt->AddText(Form("%0.0f < p_{T} (GeV/c) < %0.0f",ptbins[ip], ptbins[ip+1]));
      text->SetTextSize(0.07);
      pt->Draw();


    }
    if(iSave)c82[ic]->SaveAs(Form("AN/Bkgd/AvBackground_%s.pdf",mcent[ic]));
  }
  //return 0;

  ipad=0;
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
  TText *text2 = pt1->AddText(Form("|#eta_{Jet}| < %0.0f",ketacut));  text2->SetTextSize(0.05);
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
  pt3->SetTextFont(62);
  TText *text4 = pt3->AddText("CMS Preliminary");  text4->SetTextSize(0.045);
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
    hprms_data[1][ic]->GetYaxis()->SetTitle("RMS ( Background p_{T} ) (GeV/c)");
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
  TText *text12 = pt11->AddText(Form("|#eta_{Jet}| < %0.0f",ketacut));  text12->SetTextSize(0.05);
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
  pt33->SetTextFont(62);
  TText *text44 = pt33->AddText("CMS");  text44->SetTextSize(0.045);
  pt33->Draw();


  
  if(iSave){
    c80->SaveAs("AN/Bkgd/Background_pT_datamc.C");
    c80->SaveAs("AN/Bkgd/Background_pT_datamc.eps");
    c80->SaveAs("AN/Bkgd/Background_pT_datamc.pdf");
    c80->SaveAs("AN/Bkgd/Background_pT_datamc.gif");
  }



  return 0;
}
void GetPileup(TH2 *h2,TH1 *hmean, TH1*hrms,int iFit)
{
  for(int ix=1;ix<=h2->GetNbinsX();ix++){
    
    /*
    double binwx = h2->GetXaxis()->GetBinWidth(ix);
    for(int iy=1;iy<=h2->GetNbinsY();iy++){
      double binwy = h2->GetYaxis()->GetBinWidth(iy);
      double binc  = h2->GetBinContent(ix,iy);
      double bine  = h2->GetBinError(ix,iy);
      h2->SetBinContent(ix,iy,binc);
      h2->SetBinError(ix,iy,bine);
      }//! iy
    */

    //! pile up pt distribution
    TH1F *hpu = (TH1F*)h2->ProjectionY("hpu",ix,ix);
    hpu->SetName("hpu");
    //hpu->Scale(binwx);

    if(iFit){
      hpu->Fit("gaus","RQ0","",0,250);
      TF1 *fgaus = (TF1*)hpu->GetFunction("gaus");
      hmean->SetBinContent(ix,fgaus->GetParameter(1));
      hmean->SetBinError  (ix,fgaus->GetParError(1));
      hrms->SetBinContent (ix,fgaus->GetParameter(2));
      hrms->SetBinError   (ix,fgaus->GetParError(2));

    }else{
      hmean->SetBinContent(ix,hpu->GetMean());
      hmean->SetBinError  (ix,hpu->GetMeanError());
      hrms->SetBinContent(ix,hpu->GetRMS());
      hrms->SetBinError  (ix,hpu->GetRMSError());
    }

    delete hpu;

  }//! ix
}
void putCMSPrel(double x, double y, double size){
  TLatex *tex=0;
  tex = new TLatex(x,y,"CMS Preliminary");
  tex->SetTextSize(size);
  tex->SetLineWidth(2);
  tex->SetNDC();
  tex->Draw();
}
void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}

void drawText2(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
void makeMultiPanelCanvasWithGap(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge, const Float_t asyoffset) {
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];

  Float_t PadWidth =
     (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
                       (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
     (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
                         (1.0/(1.0-edge))+(Float_t)rows-2.0);

   //PadHeight = 0.5*PadHeight;

   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
     Xlow[i] = Xup[0] + (i-1)*PadWidth;
     Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
     if(i==rows-2){
       Ylow[i] = Yup[rows-1] + ct*PadHeight;
       Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
     }else{
       Ylow[i] = Yup[rows-1] + ct*PadHeight;
       Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
       //Yup[i] = 0.2*Yup[i];
     }
     ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
     for(Int_t j=0;j<rows;j++) {
       canv->cd();
       padName = Form("p_%d_%d",i,j);
       //pad[i][j] = new TPad(padName.Data(),padName.Data(),
       //Xlow[i],Ylow[j],Xup[i],Yup[j]);
       // this is hacked version to create aysmmetric pads around low 
       if(j==0){
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
       }else{
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
       }


       if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
       else pad[i][j]->SetLeftMargin(0);

       if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
       else pad[i][j]->SetRightMargin(0);

       if(j==0) pad[i][j]->SetTopMargin(edge);
       //else pad[i][j]->SetTopMargin(0.01);
       else pad[i][j]->SetTopMargin(0.02);

       //if(j==0) pad[i][j]->SetTopMargin(edge*3.5);
       //else pad[i][j]->SetTopMargin(0.0);

       if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
       else pad[i][j]->SetBottomMargin(0.15);

       pad[i][j]->Draw();
       pad[i][j]->cd();
       pad[i][j]->SetNumber(columns*j+i+1);
     }
   }
}
