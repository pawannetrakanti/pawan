
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

using namespace std;

const double pi=acos(-1.);
const double pi2=2*pi -1;

int bins =500;
float minval = 0;
float maxval = 1000;
float binw   = 2;

const int knj=4;
const char *calgo[knj]= {"akPu2PF","akPu3PF","akPu4PF","akPu5PF"}; 
const char *algn[knj] = {"Anti-k_{T}, PF, R = 0.2","Anti-k_{T}, PF, R = 0.3","Anti-k_{T}, PF, R = 0.4","Anti-k_{T}, PF, R = 0.5"};


//! Note :  for pp we are using the jet algorithm with out pu subtraction ie. ak3PF, ak4PF and ak5PF

const int ncen=7;
//!                          0      1        2        3        4        5      6
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","pp"};
const char *ksp  [ncen] = {"pbpb","pbpb" ,"pbpb"  ,"pbpb"  ,"pbpb"  ,"pbpb"  ,"pp"};

//const char *fopt="MLRQ+"; 
int iFit=0; 
const char *fopt="RQ+";
const int knpx=2000;
float fitmin=0.00;
float fitmax=5.00;

double xmin=50.;
double xmax=400.;
double xfitmin=60;
double xfitmax=360;
int maxEntry=5;

int GetPtBin(float /*pt*/);
void MakeHist(TH1 */*hist*/,int /*istat*/);
void MakeHistRMS(TH1 */*hRMS*/,float /*max*/,float /*min*/);
void MakeHistMean(TH1 */*Mean*/,float /*max*/,float /*min*/);
void CleanHist(TH1 */*h1D*/,float /*lxrange*/,float /*hxrange*/);
void FillMeanSigma(int /*ip*/,TH1 */*h1*/,TH1 */*ArM*/,TH1 */*RMS*/,TH1 */*Mean*/,TH1 */*Sigma*/);


void LoadStyle();
void drawText(const char */*text*/, float /*xp*/, float /*yp*/, int /*size*/);
void drawText2(const char */*text*/, float /*xp*/, float /*yp*/, int /*size*/);
void makeMultiPanelCanvas(TCanvas*& /*canv*/,
                          const Int_t /*columns*/,
                          const Int_t /*rows*/,
                          const Float_t /*leftOffset*/,
                          const Float_t /*bottomOffset*/,
                          const Float_t /*leftMargin*/,
                          const Float_t /*bottomMargin*/,
                          const Float_t /*edge*/,const Float_t /*asyoffset*/); 


int JERJES(const char *reta = "eta2")
{

  LoadStyle();

  int rfit=0;
  int statop=1;
  int kRebin=1;
  //const char *reta = "eta3.0";
  //bool iSave=false;

  float ketacut=2.0;
  if(strcmp(reta,"eta3.0")==0)ketacut=3;
  bool iSigma=false;

  if(rfit==0){fitmin=0.01;fitmax=2.00;}   
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

  if(kRebin){
    if(bins%kRebin!=0){
      cout<<"Cannot be divided in these bins chose another combination : "<<endl;
      return 0;
    }
    bins /= kRebin;
    binw = (maxval  - minval)/(1.0*bins);
    cout<<"kRebin : "<<kRebin<<"\t bins : "<<bins<<"\t binw  : "<<binw<<endl;
  }
 
  //! Input files 
  TFile *fin_pbpb = new TFile(Form("output/test/%s/MC_pbpb_2013_merged.root",reta),"r");
  //! 0   : pp w/o 
  TFile *fin_pp   = new TFile(Form("output/test/%s/MC_pp_2013_merged.root",reta),"r");
  
  cout<<"\t"<<endl
      <<"rfit : "<<rfit<<"\t fitmin : "<<fitmin<<"\t fitmax : "<<fitmax<<endl
      <<"Input file name  pbpb : "<<fin_pbpb->GetName()<<endl
      <<"Input file name  pp   : "<<fin_pp->GetName()<<endl
      <<"\t"<<endl;


  //! used for fit
  double ptbins[] = {31,41,61,81,101,121,141,161,181,201,221,241,261,281,291,321};
  //double ptbins[] = {50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,340,400,440,600};
  const int b1  = sizeof(ptbins)/sizeof(Double_t) - 1;
  const int nbins = b1;
  cout<<"# of pt bins : "<<nbins<<endl;
  int maxr=3;

  int ipad=0;

  //! 0-5 PbPb Resposnse and 6th is pp response
  //! Response corrected pT  (jtpt) / gen pT vs gen pT
  TH2F *hratiocorrrefpt[knj][ncen];
  TH1F *hratiocorrrefpt1D[knj][ncen][nbins];  
  TH1F *hMean[knj][ncen], *hSigma[knj][ncen], *hRMS[knj][ncen], *hArM[knj][ncen];

  //! Response raw pT / gen pT  vs gen pT
  TH2F *hratiorawrefpt[knj][ncen];
  TH1F *hratiorawrefpt1D[knj][ncen][nbins];  
  TH1F *hMean_r[knj][ncen], *hSigma_r[knj][ncen], *hRMS_r[knj][ncen], *hArM_r[knj][ncen];

  //! Ratio of Fit and geometric mean/RMS
  //TH1F *hRatio_Mean[knj][ncen], *hRatio_RMS[knj][ncen];


  for(int nj=0;nj<knj;nj++){
    //cout<<"nj : "<<nj<<Form("\t %s",calgo[nj])<<endl;
    
    for(int icen=0;icen<ncen;icen++){

      if(icen < ncen-1){
	if(nj==1)cout<<"icen : "<<ccent[icen]<<endl;

	//! PbPb Histograms
	//! corrected jet pT
	hMean [nj][icen] = new TH1F(Form("hMean%d_%d",nj,icen),Form("<reco p_{T}/gen p_{T}> %s pbpb %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hMean[nj][icen],statop);
	hArM [nj][icen] = new TH1F(Form("hArM%d_%d",nj,icen),Form("<reco p_{T}/gen p_{T}> %s pbpb %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hArM[nj][icen],statop);
	hSigma [nj][icen] = new TH1F(Form("hSigma%d_%d",nj,icen),Form("#sigma(reco p_{T}/gen p_{T}) %s pbpb %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hSigma[nj][icen],statop);    
	hRMS [nj][icen] = new TH1F(Form("hRMS%d_%d",nj,icen),Form("RMS(reco p_{T}/gen p_{T}) %s pbpb %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hRMS[nj][icen],statop);        
	
	//! Raw jet pt
	hMean_r [nj][icen] = new TH1F(Form("hMean_r%d_%d",nj,icen),Form("<raw p_{T}/gen p_{T}> %s pbpb %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hMean_r[nj][icen],statop);
	hArM_r [nj][icen] = new TH1F(Form("hArM_r%d_%d",nj,icen),Form("<raw p_{T}/gen p_{T}> %s pbpb %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hArM_r[nj][icen],statop);
	hSigma_r [nj][icen] = new TH1F(Form("hSigma_r%d_%d",nj,icen),Form("#sigma(raw p_{T}/gen p_{T}) %s pbpb %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hSigma_r[nj][icen],statop);    
	hRMS_r [nj][icen] = new TH1F(Form("hRMS_r%d_%d",nj,icen),Form("RMS(raw p_{T}/gen p_{T}) %s pbpb %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hRMS_r[nj][icen],statop);        

	//! PbPb  starts here
	hratiocorrrefpt[nj][icen] = (TH2F*)fin_pbpb->Get(Form("hrescrpt_genm%d_%d",nj+knj,icen));
	hratiocorrrefpt[nj][icen]->SetName(Form("hratiocorrrefpt_pbpb%d_%d",nj,icen));
	
	hratiorawrefpt[nj][icen] = (TH2F*)fin_pbpb->Get(Form("hresrrpt_genm%d_%d",nj+knj,icen));
	hratiorawrefpt[nj][icen]->SetName(Form("hratiorawrefpt_pbpb%d_%d",nj,icen));
	
	for(int ip=0;ip<nbins;ip++){
	  int lbin = (int)(ptbins[ip]   - minval)/binw +1;
	  int hbin = (int)(ptbins[ip+1] - minval)/binw +1;
	  hratiocorrrefpt1D[nj][icen][ip]  = (TH1F*)hratiocorrrefpt [nj][icen]->ProjectionY(Form("hratiocorrrefpt1D%d_%d_%d",nj,icen,ip),lbin,hbin,"e");
	  if(hratiocorrrefpt1D[nj][icen][ip]->GetEntries()<maxEntry)continue;
	  hratiocorrrefpt1D[nj][icen][ip]->Rebin(2);
	  FillMeanSigma(ip,hratiocorrrefpt1D[nj][icen][ip],
			hArM[nj][icen],hRMS[nj][icen],
			hMean[nj][icen],hSigma[nj][icen]);
	  
	  //! Raw / Gen
	  hratiorawrefpt1D[nj][icen][ip]  = (TH1F*)hratiorawrefpt [nj][icen]->ProjectionY(Form("hratiorawrefpt1D%d_%d_%d",nj,icen,ip),lbin,hbin,"e");
	  if(hratiorawrefpt1D[nj][icen][ip]->GetEntries()<maxEntry)continue;
	  hratiorawrefpt1D[nj][icen][ip]->Rebin(2);
	  FillMeanSigma(ip,hratiorawrefpt1D[nj][icen][ip],
			hArM_r[nj][icen],hRMS_r[nj][icen],
			hMean_r[nj][icen],hSigma_r[nj][icen]);
	}

      	//! PbPb ends
      }else{

	if(nj==1)cout<<"icen : "<<ccent[icen]<<endl;
	//! pp /////////////////////////////
	hMean [nj][icen] = new TH1F(Form("hMean%d_%d",nj,icen),Form("<reco p_{T}/gen p_{T}> %s pp %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hMean[nj][icen],statop);
	hArM [nj][icen] = new TH1F(Form("hArM%d_%d",nj,icen),Form("<reco p_{T}/gen p_{T}> %s pp %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hArM[nj][icen],statop);
	hSigma [nj][icen] = new TH1F(Form("hSigma%d_%d",nj,icen),Form("#sigma(reco p_{T}/gen p_{T}) %s pp %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hSigma[nj][icen],statop);    
	hRMS [nj][icen] = new TH1F(Form("hRMS%d_%d",nj,icen),Form("RMS(reco p_{T}/gen p_{T}) %s pp %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hRMS[nj][icen],statop);        
	
	hMean_r [nj][icen] = new TH1F(Form("hMean_r%d_%d",nj,icen),Form("<raw p_{T}/gen p_{T}> %s pp %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hMean_r[nj][icen],statop);
	hArM_r [nj][icen] = new TH1F(Form("hArM_r%d_%d",nj,icen),Form("<raw p_{T}/gen p_{T}> %s pp %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hArM_r[nj][icen],statop);
	hSigma_r [nj][icen] = new TH1F(Form("hSigma_r%d_%d",nj,icen),Form("#sigma(raw p_{T}/gen p_{T}) %s pp %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hSigma_r[nj][icen],statop);    
	hRMS_r [nj][icen] = new TH1F(Form("hRMS_r%d_%d",nj,icen),Form("RMS(raw p_{T}/gen p_{T}) %s pp %d",calgo[nj],icen),nbins,ptbins);
	MakeHist(hRMS_r[nj][icen],statop);        
      
	//! Histograms for the response calculation
	hratiocorrrefpt[nj][icen]  = (TH2F*)fin_pp->Get(Form("hrescrpt_genm%d_%d",nj,0));
	hratiocorrrefpt[nj][icen]->SetName(Form("hratiocorrrefpt_pp%d_%d",nj,icen));
	
	hratiorawrefpt[nj][icen]  = (TH2F*)fin_pp->Get(Form("hresrrpt_genm%d_%d",nj,0));
	hratiorawrefpt[nj][icen]->SetName(Form("hratiorawrefpt_pp%d_%d",nj,icen));
	
	/*
	//! Pile up subtracted algorithms
	hratiocorrrefpt[nj][icen]  = (TH2F*)fin_pp->Get(Form("hrescrpt%d_%d",nj+knj,icen));
	hratiocorrrefpt[nj][icen]->SetName(Form("hratiocorrrefpt_pp%d_%d",nj,icen));
	
	hratiorawrefpt[nj][icen]  = (TH2F*)fin_pp->Get(Form("hresrrpt%d_%d",nj+knj,icen));
	hratiorawrefpt[nj][icen]->SetName(Form("hratiorawrefpt_pp%d_%d",nj,icen));
	*/
	
	
	for(int ip=0;ip<nbins;ip++){
	  int lbin = (int)(ptbins[ip]   - minval)/binw +1;
	  int hbin = (int)(ptbins[ip+1] - minval)/binw +1;
	  hratiocorrrefpt1D[nj][icen][ip]  = (TH1F*)hratiocorrrefpt [nj][icen]->ProjectionY(Form("hratiocorrrefpt1D_pp%d_%d_%d",nj,icen,ip),lbin,hbin,"e");
	  if(hratiocorrrefpt1D[nj][icen][ip]->GetEntries()<maxEntry)continue;
	  hratiocorrrefpt1D[nj][icen][ip]->Rebin(2);
	  FillMeanSigma(ip,hratiocorrrefpt1D[nj][icen][ip],
			hArM[nj][icen],hRMS[nj][icen],
			hMean[nj][icen],hSigma[nj][icen]);
	  
	  //! Raw/Gen
	  hratiorawrefpt1D[nj][icen][ip]  = (TH1F*)hratiorawrefpt [nj][icen]->ProjectionY(Form("hratiorawrefpt1D_pp%d_%d_%d",nj,icen,ip),lbin,hbin,"e");
	  if(hratiorawrefpt1D[nj][icen][ip]->GetEntries()<maxEntry)continue;
	  hratiorawrefpt1D[nj][icen][ip]->Rebin(2);
	  FillMeanSigma(ip,hratiorawrefpt1D[nj][icen][ip],
			hArM_r[nj][icen],hRMS_r[nj][icen],
			hMean_r[nj][icen],hSigma_r[nj][icen]);
	}//! ip bin
      }//! pp ends    
    }//! icen loop ends
  }//! nj loop ends

  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){
      
      if(ic<ncen-1){
	MakeHistMean(hArM  [nj][ic],1.058,0.955);
	MakeHistRMS (hRMS  [nj][ic],0.43,0.001);
	MakeHistMean(hMean [nj][ic],1.058,0.955);
	MakeHistRMS (hSigma[nj][ic],0.43,0.001);
	
	MakeHistMean(hArM_r  [nj][ic],1.058,0.955);
	MakeHistRMS (hRMS_r  [nj][ic],0.43,0.001);
	MakeHistMean(hMean_r [nj][ic],1.058,0.955);
	MakeHistRMS (hSigma_r[nj][ic],0.43,0.001);
      }else{
	MakeHistMean(hArM  [nj][ic],1.058,0.955);
	MakeHistRMS (hRMS  [nj][ic],0.43,0.001);
	MakeHistMean(hMean [nj][ic],1.058,0.955);
	MakeHistRMS (hSigma[nj][ic],0.43,0.001);
	
	MakeHistMean(hArM_r  [nj][ic],1.058,0.955);
	MakeHistRMS (hRMS_r  [nj][ic],0.43,0.001);
	MakeHistMean(hMean_r [nj][ic],1.058,0.955);
	MakeHistRMS (hSigma_r[nj][ic],0.43,0.001);
      }
    }      
  }


  TCanvas *c3[knj];
  TLegend *l3[knj];
  int maxc=ncen-1;
  maxr=2;
  for(int nj=0;nj<knj;nj++){
    //if(nj!=1)continue;
    c3[nj] = new TCanvas(Form("c3_%d",nj),Form("%d JES JER",nj),97,118,1739,479);
    makeMultiPanelCanvas(c3[nj],maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
    //std::cout<<std::endl;
    ipad=0;

    l3[nj] = new TLegend(0.2453033,0.6247655,0.4838159,0.9399625,NULL,"BRNDC");
    l3[nj]->SetHeader("");
    l3[nj]->SetBorderSize(0);
    l3[nj]->SetTextFont(42);
    l3[nj]->SetTextSize(0.09);
    l3[nj]->SetLineColor(1);
    l3[nj]->SetLineStyle(1);
    l3[nj]->SetLineWidth(1);
    l3[nj]->SetFillColor(10);
    l3[nj]->SetFillStyle(1001);
    l3[nj]->SetHeader("");			  

    for(int ic=ncen-2;ic>=0;ic--){

      hRMS[nj][ic]->SetMarkerStyle(24);
      hRMS[nj][ic]->SetMarkerColor(1);
      hRMS[nj][ic]->SetLineColor(1);
      hRMS[nj][ic]->SetMarkerSize(1.3);

      hRMS[nj][ncen-1]->SetMarkerStyle(20);
      hRMS[nj][ncen-1]->SetMarkerColor(1);
      hRMS[nj][ncen-1]->SetLineColor(1);
      hRMS[nj][ncen-1]->SetMarkerSize(1.1);

      hArM[nj][ic]->SetMarkerStyle(24);
      hArM[nj][ic]->SetMarkerColor(1);
      hArM[nj][ic]->SetLineColor(1);
      hArM[nj][ic]->SetMarkerSize(1.3);

      hArM[nj][ncen-1]->SetMarkerStyle(20);
      hArM[nj][ncen-1]->SetMarkerColor(1);
      hArM[nj][ncen-1]->SetLineColor(1);
      hArM[nj][ncen-1]->SetMarkerSize(1.1);

      hArM_r[nj][ic]->SetMarkerStyle(25);
      hArM_r[nj][ic]->SetMarkerColor(1);
      hArM_r[nj][ic]->SetLineColor(1);
      hArM_r[nj][ic]->SetMarkerSize(1.3);

      hArM_r[nj][ncen-1]->SetMarkerStyle(21);
      hArM_r[nj][ncen-1]->SetMarkerColor(1);
      hArM_r[nj][ncen-1]->SetLineColor(1);
      hArM_r[nj][ncen-1]->SetMarkerSize(1.1);
      
      c3[nj]->cd(++ipad);
      //gPad->SetLogx();
      hRMS[nj][ic]->Draw("p");
      hRMS[nj][ncen-1]->Draw("psame");

      if(ipad==6){
	l3[nj]->AddEntry(hRMS[nj][ncen-1],"pp","p"); 
	//l3[nj]->AddEntry(hRMS[nj][ncen-1],"pp w pu subtraction","p"); 
	l3[nj]->AddEntry(hRMS[nj][0],"PbPb","p"); 
	l3[nj]->Draw();
      }
      if(ipad==1){
	TPaveText *pt3 = new TPaveText(0.4193406,0.7098186,0.7083456,0.8899312,"brNDC");
	pt3->SetBorderSize(0);
	pt3->SetFillColor(10);
	pt3->SetTextFont(42);
	TText *text3 = pt3->AddText("CMS Preliminary");
	text3->SetTextSize(0.09);
	TText *text4 = pt3->AddText("PYTHIA Z2 + HYDJET 1.8");
	text4->SetTextSize(0.09);
	pt3->Draw();
      }
      
      TPaveText *pt = new TPaveText(0.2663379,0.02939336,0.4986753,0.1544715,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(10);
      pt->SetTextFont(42);
      TText *text = pt->AddText(Form("%s",ccent[ic]));
      text->SetTextSize(0.09);
      pt->Draw();

      if(ipad==2){
	TPaveText *pt1 = new TPaveText(0.1346399,0.6797999,0.9011061,0.8749218,"brNDC");
	pt1->SetBorderSize(0);
	pt1->SetFillColor(10);
	pt1->SetTextFont(42);
	//TText *text1 = pt1->AddText(Form("%s",calgo[nj]));
	TText *text1 = pt1->AddText(Form("%s",algn[nj]));
	text1->SetTextSize(0.09);
	TText *text2 = pt1->AddText(Form("|#eta| < %0.0f",ketacut));
	text2->SetTextSize(0.09);
	pt1->Draw();
      }


      c3[nj]->cd(ipad+(ncen-1));

      //cout<<"pad : "<<(ipad+(ncen-1))<<"\t "<<calgo[nj]<<"\t centrality : "<<ccent[ic]<<"\t name : "<<hArM[nj][ic]->GetName()<<endl;
      //gPad->SetLogx();

      hArM[nj][ic]->Draw("p");
      hArM_r[nj][ic]->Draw("psame");      
      hArM[nj][ncen-1]->Draw("psame");
      hArM_r[nj][ncen-1]->Draw("psame");


      if(ipad==6){
	TLatex *tex1 = new TLatex(175,0.68,"Raw");
	tex1->SetTextFont(42);
	tex1->SetTextSize(0.08);
	tex1->SetLineWidth(2);
	tex1->Draw();
	TLatex *tex2 = new TLatex(174,0.63,"Corrected");
	tex2->SetTextFont(42);
	tex2->SetTextSize(0.08);
	tex2->SetLineWidth(2);
	tex2->Draw();
	TMarker *marker1 = new TMarker(162,0.65,25);
	marker1->SetMarkerStyle(25);
	marker1->Draw();
	TMarker *marker2 = new TMarker(139,0.65,21);
	marker2->SetMarkerStyle(21);
	marker2->Draw();
	TMarker *marker3 = new TMarker(162,0.70,24);
	marker3->SetMarkerStyle(24);
	marker3->Draw();
	TMarker *marker4 = new TMarker(139,0.70,20);
	marker4->SetMarkerStyle(20);
	marker4->Draw();
      }
      TLine *line = new TLine(xmin,1,xmax,1);
      line->SetLineWidth(1);
      line->SetLineStyle(2);
      line->Draw();
    }//! icen
    /*
    if(iSave){
      c3[nj]->SaveAs(Form("AN/JERJES/JERJES_%s.png",calgo[nj]));
      c3[nj]->SaveAs(Form("AN/JERJES/JERJES_%s.pdf",calgo[nj]));
      c3[nj]->SaveAs(Form("AN/JERJES/JERJES_%s.C",calgo[nj]));
      c3[nj]->SaveAs(Form("AN/JERJES/JERJES_%s.eps",calgo[nj]));
    }
    */
  }//! algo
  return 0;

  /*
  ipad=0;
  TCanvas *c99[knj][ncen];
  maxr=3;
  for(int nj=1;nj<2;nj++){
    for(int ic=0;ic<ncen;ic++){
      c99[nj][ic] = new TCanvas(Form("c99_%d_%d",nj,ic),Form("%s Fitting plots %s",calgo[nj],ccent[ic]),100,102,1399,942);
      c99[nj][ic]->Divide(6,maxr,0,0);
      ipad=0;
      
      for(int ip=0;ip<nbins;ip++){      
	c99[nj][ic]->cd(++ipad);
	if(ipad%6==0)gPad->SetRightMargin(0.02);
	gPad->SetBottomMargin(0.15);
	gPad->SetLogy();
	

	if(ic<ncen-1){
	  hratiocorrrefpt1D[nj][ic][ip]->SetMaximum(25.634);
	  hratiocorrrefpt1D[nj][ic][ip]->SetMinimum(1e-09);
	  hratiocorrrefpt1D[nj][ic][ip]->SetTitle(0);
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitleFont(42);
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetLabelFont(42);
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetLabelSize(0.08);
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitleSize(0.07);
	
	  hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetTitle("");
	  hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetTitleFont(42);
	  hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetLabelFont(42);
	  hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetLabelSize(0.08);

	  hratiocorrrefpt1D[nj][ic][ip]->SetMarkerStyle(24);
	  hratiocorrrefpt1D[nj][ic][ip]->SetMarkerColor(1);
	  hratiocorrrefpt1D[nj][ic][ip]->SetLineColor(1);
	  hratiocorrrefpt1D[nj][ic][ip]->SetMarkerSize(0.9);
	  hratiocorrrefpt1D[nj][ic][ip]->Draw("p");  

	  hratiorawrefpt1D[nj][ic][ip]->SetMarkerStyle(24);
	  hratiorawrefpt1D[nj][ic][ip]->SetMarkerColor(4);
	  hratiorawrefpt1D[nj][ic][ip]->SetLineColor(4);
	  hratiorawrefpt1D[nj][ic][ip]->SetMarkerSize(0.9);
	  //hratiorawrefpt1D[nj][ic][ip]->Draw("psame");  
	
	  c99[nj][ic]->Update();
	  TPaveStats *ps = (TPaveStats*)  hratiocorrrefpt1D[nj][ic][ip]->GetListOfFunctions()->FindObject("stats");
	  ps->SetX1NDC(0.50);
	  ps->SetY1NDC(0.41);       
	  ps->SetX2NDC(0.95);
	  ps->SetY2NDC(0.79);
	  ps->SetTextFont(42);
	  ps->Draw();

	}else{
	  
	  hratiocorrrefpt1D[nj][ic][ip]->SetMaximum(25.634);
	  hratiocorrrefpt1D[nj][ic][ip]->SetMinimum(1e-09);
	  hratiocorrrefpt1D[nj][ic][ip]->SetTitle(0);
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitleFont(42);
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetLabelFont(42);
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetLabelSize(0.08);
	  hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitleSize(0.07);
	
	  hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetTitle("");
	  hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetTitleFont(42);
	  hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetLabelFont(42);
	  hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetLabelSize(0.08);

	  hratiocorrrefpt1D[nj][ic][ip]->SetMarkerStyle(24);
	  hratiocorrrefpt1D[nj][ic][ip]->SetMarkerColor(1);
	  hratiocorrrefpt1D[nj][ic][ip]->SetLineColor(1);
	  hratiocorrrefpt1D[nj][ic][ip]->SetMarkerSize(0.9);
	  hratiocorrrefpt1D[nj][ic][ip]->Draw("p");  
	
	  c99[nj][ic]->Update();
	  TPaveStats *ps = (TPaveStats*)  hratiocorrrefpt1D[nj][ic][ip]->GetListOfFunctions()->FindObject("stats");
	  ps->SetX1NDC(0.50);
	  ps->SetY1NDC(0.41);       
	  ps->SetX2NDC(0.95);
	  ps->SetY2NDC(0.79);
	  ps->SetTextFont(42);
	  ps->Draw();
	}	

	
	TPaveText *pt   = new TPaveText(0.4524683,0.8914759,0.7023389,0.9597512,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(10);
	pt->SetTextFont(42);
	TText *text = pt->AddText(Form("%0.0f < p_{T} (GeV/c) < %0.0f",ptbins[ip], ptbins[ip+1]));
	text->SetTextSize(0.07);
	pt->Draw();
	
	if(ipad==1){
	  TPaveText *pt1 = new TPaveText(0.6044166,0.2194909,0.8644171,0.3668644,"brNDC");
	  pt1->SetBorderSize(0);
	  pt1->SetFillColor(10);
	  pt1->SetTextFont(42);
	  TText *text1 = 0;
	  if(ic==ncen-1)text1 = pt1->AddText("ak3PF");
	  else text1 = pt1->AddText(calgo[nj]);
	  text1->SetTextSize(0.09);	
	  pt1->Draw();
	  
	  TPaveText *pt2 = new TPaveText(0.60104,0.8160025,0.8475339,0.9142515,"brNDC");
	  pt2->SetBorderSize(0);
	  pt2->SetFillColor(10);
	  pt->SetTextFont(42);
	  TText *text2 = pt2->AddText(Form("%s",ccent[ic]));
	  text2->SetTextSize(0.08);
	  pt2->Draw();
	}
      }
    }
  }
  return 0;
  */  

  /*
  //! Difference between Gaussian fits mean/sigma and ArM/RMS
  cout<<"Sigma and RMS ....." <<endl;
  ipad=0;
  ialg=1;
  TCanvas *c7 = new TCanvas("c7","Sigma and RMSFit resolution",15,131,1885,546);
  makeMultiPanelCanvas(c7,ncen,2,0.0,0.0,0.22,0.22,0.02,0);
  //c7->SetGridx();
  TLegend *l7 = new TLegend(0.2453033,0.6247655,0.4838159,0.9399625,NULL,"BRNDC");
  l7->SetHeader("");
  l7->SetBorderSize(0);
  l7->SetTextFont(42);
  l7->SetTextSize(0.09);
  l7->SetLineColor(1);
  l7->SetLineStyle(1);
  l7->SetLineWidth(1);
  l7->SetFillColor(10);
  l7->SetFillStyle(1001);
  l7->SetHeader("");			  
  l7->AddEntry(hRMS[ialg][0],"RMS","p");
  l7->AddEntry(hSigma[ialg][0],"Gaussian Sigma","p");

  for(int icen=ncen-1;icen>=0;icen--){
    c7->cd(++ipad); 
    //gPad->SetGridx();
    //gPad->SetLogx();

    if(ipad==1){//! pp

      hSigma[ialg][icen]->SetStats(0);
      hSigma[ialg][icen]->SetMarkerStyle(20);
      hSigma[ialg][icen]->SetMarkerColor(1);
      hSigma[ialg][icen]->SetLineColor(1);
      hSigma[ialg][icen]->SetMarkerSize(1.0);

      hRMS[ialg][icen]->SetStats(0);      
      hRMS[ialg][icen]->SetMarkerStyle(24);
      hRMS[ialg][icen]->SetMarkerColor(1);
      hRMS[ialg][icen]->SetLineColor(1);
      hRMS[ialg][icen]->SetMarkerSize(1.2);
      
      hRMS[ialg][icen]->Draw("p");
      hSigma[ialg][icen]->Draw("psame");
      hRatio_RMS[ialg][icen] = (TH1F*)hRMS[ialg][icen]->Clone(Form("hRatio_RMS%d_%d",1,icen));
      hRatio_RMS[ialg][icen]->Divide(hSigma[ialg][icen]);

    }else{//! pbpb

      hSigma[ialg][icen]->SetStats(0);      
      hSigma[ialg][icen]->SetMarkerStyle(20);
      hSigma[ialg][icen]->SetMarkerColor(1);
      hSigma[ialg][icen]->SetLineColor(1);
      hSigma[ialg][icen]->SetMarkerSize(1.0);

      hRMS[ialg][icen]->SetStats(0);            
      hRMS[ialg][icen]->SetMarkerStyle(24);
      hRMS[ialg][icen]->SetMarkerColor(1);
      hRMS[ialg][icen]->SetLineColor(1);
      hRMS[ialg][icen]->SetMarkerSize(1.2);
      
      hRMS[ialg][icen]->Draw("p");
      hSigma[ialg][icen]->Draw("psame");
      


      hRatio_RMS[ialg][icen] = (TH1F*)hRMS[ialg][icen]->Clone(Form("hRatio_RMS%d_%d",1,icen));
      hRatio_RMS[ialg][icen]->Divide(hSigma[ialg][icen]);
    }
  
    gPad->Update();
    
    TPaveText *pt1 = new TPaveText(0.2846445,0.06345905,0.4451356,0.2028512,"brNDC");
    pt1->SetBorderSize(0);
    pt1->SetFillColor(10);
    pt1->SetTextFont(42);
    TText *text1 = pt1->AddText(Form("%s",ccent[icen]));
    text1->SetTextSize(0.08);
    pt1->Draw();

    if(ipad==1)l7->Draw();
    
    gPad->Update();

    c7->cd(ipad+ncen);
    //gPad->SetLogx();
    //gPad->SetGridx();
    
    TLine *l2 = new TLine(xmin,1,xmax+250,1);
    l2->SetLineWidth(1);
    l2->SetLineStyle(2);

    hRatio_RMS[ialg][icen]->SetStats(0);
    hRatio_RMS[ialg][icen]->SetMaximum(1.546);
    hRatio_RMS[ialg][icen]->SetMinimum(0.846);
    hRatio_RMS[ialg][icen]->GetYaxis()->SetTitle("RMS/Sigma");
    hRatio_RMS[ialg][icen]->Draw("p");
    l2->Draw();
  }

  // c7->SaveAs("AN/JER/Fits_RMSSigma_JER.png");
  // c7->SaveAs("AN/JER/Fits_RMSSigma_JER.C");
  // c7->SaveAs("AN/JER/Fits_RMSSigma_JER.eps");
  // c7->SaveAs("AN/JER/Fits_RMsSigma_JER.pdf");


  cout<<endl;
  cout<<"ArM and Mean ....." <<endl;
  //! Mean
  ipad=0;
  TCanvas *c8 = new TCanvas("c8","ArM and Mean",15,131,1885,546);
  makeMultiPanelCanvas(c8,ncen,2,0.0,0.0,0.22,0.22,0.02,0);

  TLegend *l8 = new TLegend(0.2453033,0.6247655,0.4838159,0.9399625,NULL,"BRNDC");
  l8->SetHeader("");
  l8->SetBorderSize(0);
  l8->SetTextFont(42);
  l8->SetTextSize(0.09);
  l8->SetLineColor(1);
  l8->SetLineStyle(1);
  l8->SetLineWidth(1);
  l8->SetFillColor(10);
  l8->SetFillStyle(1001);
  l8->SetHeader("");			  
  l8->AddEntry(hArM[ialg][0],"Ar Mean","p");
  l8->AddEntry(hMean[ialg][0],"Gaussian Mean","p");

  for(int icen=ncen-1;icen>=0;icen--){
    
    c8->cd(++ipad); 

    if(ipad==1){//! pp
      hMean[ialg][icen]->SetStats(0);
      hMean[ialg][icen]->SetMarkerStyle(20);
      hMean[ialg][icen]->SetMarkerColor(1);
      hMean[ialg][icen]->SetLineColor(1);
      hMean[ialg][icen]->SetMarkerSize(1.0);
      
      hArM[ialg][icen]->SetStats(0);
      hArM[ialg][icen]->SetMarkerStyle(24);
      hArM[ialg][icen]->SetMarkerColor(1);
      hArM[ialg][icen]->SetLineColor(1);
      hArM[ialg][icen]->SetMarkerSize(1.2);

      hArM[ialg][icen]->Draw("p");
      hMean[ialg][icen]->Draw("psame");
      
      hRatio_Mean[ialg][icen] = (TH1F*)hArM[ialg][icen]->Clone(Form("hRatio_ArM%d_%d",ialg,icen));
      hRatio_Mean[ialg][icen]->Divide(hMean[ialg][icen]);
      
    }else{

      hMean[ialg][icen]->SetStats(0);
      hMean[ialg][icen]->SetMarkerStyle(20);
      hMean[ialg][icen]->SetMarkerColor(1);
      hMean[ialg][icen]->SetLineColor(1);
      hMean[ialg][icen]->SetMarkerSize(1.0);

      hArM[ialg][icen]->SetStats(0);
      hArM[ialg][icen]->SetMarkerStyle(24);
      hArM[ialg][icen]->SetMarkerColor(1);
      hArM[ialg][icen]->SetLineColor(1);
      hArM[ialg][icen]->SetMarkerSize(1.2);

      hArM[ialg][icen]->Draw("p");      
      hMean[ialg][icen]->Draw("psame");

      hRatio_Mean[ialg][icen] = (TH1F*)hArM[ialg][icen]->Clone(Form("hRatio_Mean%d_%d",ialg,icen));
      hRatio_Mean[ialg][icen]->Divide(hMean[ialg][icen]);
    }      

    TPaveText *pt1 = new TPaveText(0.2846445,0.06345905,0.4451356,0.2028512,"brNDC");
    pt1->SetBorderSize(0);
    pt1->SetFillColor(10);
    pt1->SetTextFont(42);
    TText *text1 = pt1->AddText(Form("%s",ccent[icen]));
    text1->SetTextSize(0.08);
    pt1->Draw();

    if(ipad==1)l8->Draw();

    c8->cd(ipad+ncen);
    TLine *l2 = new TLine(xmin,1,xmax+250,1);
    l2->SetLineWidth(1);
    l2->SetLineStyle(2);

    hRatio_Mean[ialg][icen]->SetStats(0);
    hRatio_Mean[ialg][icen]->SetMaximum(1.046);
    hRatio_Mean[ialg][icen]->SetMinimum(0.946);
    hRatio_Mean[ialg][icen]->GetYaxis()->SetTitle("ArM/Mean");
    hRatio_Mean[ialg][icen]->Draw("p");
    l2->Draw();
  }
  // c8->SaveAs("AN/JER/Fits_RMSSigma_JES.png");
  // c8->SaveAs("AN/JER/Fits_RMSSigma_JES.C");
  // c8->SaveAs("AN/JER/Fits_RMSSigma_JES.eps");
  // c8->SaveAs("AN/JER/Fits_RMsSigma_JES.pdf");
  */
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
      //h1F->Fit("gaus",fopt,"",fitmin,fitmax);
      h1F->Fit("gaus","MLLRQ0+","",fitmin,fitmax);
      TF1* f2 = (TF1*)h1F->GetFunction("gaus");
      f2->SetLineWidth(1);
      f2->SetLineStyle(2);
      f2->SetNpx(knpx);
      hMean ->SetBinContent(ip+1,f2->GetParameter(1));
      hSigma->SetBinContent(ip+1,f2->GetParameter(2));
      
      if(strcmp(fopt,"MLRQ+")==0){
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
      
      //h1F->Fit("f1",fopt,"",fitmin,fitmax);
      h1F->Fit("f1","MLLRQ0+","",fitmin,fitmax);
      hMean ->SetBinContent(ip+1,f1->GetParameter(2));
      hSigma->SetBinContent(ip+1,f1->GetParameter(1));
      
      if(strcmp(fopt,"MLRQ+")==0){
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
  int ibin=-1;
  ibin = (int)(pt + minval)/binw;
  return ibin;
}
void CleanHist(TH1 *h1F,float lxval, float hxval)
{
  for(int ix=1;ix<=h1F->GetNbinsX();ix++){
    double val = h1F->GetBinCenter(ix);
    if(val<lxval || val>hxval){
      h1F->SetBinContent(ix,0);
      h1F->SetBinError(ix,0);
    }
  }
}
void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge,const Float_t asyoffset) {
  if (canv==0) {
    cout<<"makeMultiPanelCanvas :  Got null canvas."<<endl;
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
  
  Xlow[0] = leftOffset - asyoffset;
  Xup[0]  = leftOffset + PadWidth/(1.0-leftMargin);
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
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  
  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
			   Xlow[i],Ylow[j],Xup[i],Yup[j]);
      
      // this is hacked version to create aysmmetric pads around low 
      if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
      else pad[i][j]->SetLeftMargin(0);
      
      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);
      
      if(j==0) pad[i][j]->SetTopMargin(edge);
      else pad[i][j]->SetTopMargin(0);
      
      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else pad[i][j]->SetBottomMargin(0);
      
      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
    }
  }
}
void LoadStyle()
{
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
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
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
