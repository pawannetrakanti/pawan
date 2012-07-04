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
#include <TCanvas.h>
#include <TMarker.h>
#include <TString.h>
#include <TStopwatch.h>
#include "MultiCanvas.h"

using namespace std;


const double pi=acos(-1.);
const double pi2=2*pi -1;

//! pt binning
//  double ptbins[] ={10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,200.,220.,240.,260.,280.,300.,325.,350.,375.,400.,450.,500.,825.};
/*
double ptbins[] ={30.,40.,50.,60.,70.,80.,90.,
		  100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,
		  200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,
		  300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,
		  400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,
		  500.};
*/

//double ptbins[] ={80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 220, 240, 260, 300, 340, 400};
double ptbins[] ={80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400};
const int bins  = sizeof(ptbins)/sizeof(Double_t) - 1;
const int nbins = bins;

//! data pt binning     
double ptbins_data[] ={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
const int dbins = sizeof(ptbins_data)/sizeof(double) - 1;



const int knj=7;
const char *calgo[knj]= {
  "icPu5",
  "ak2PF","ak3PF","ak4PF",
  "akPu2PF","akPu3PF","akPu4PF"
};             
const int maxe=5;
const char *ceta[maxe] = {"0<|#eta|<0.4",
			  "0.4<|#eta|<0.8",
			  "0.8<|#eta|<1.2",
			  "1.2<|#eta|<1.6",
			  "1.6<|#eta|<2.0"
};
const char *meta[2] = {"|#eta|<1.3","1.3<|#eta|<2.0"};
const int ncen=7;
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","pp"};

const int maxEntry=50;
float kgenptcutL=80;
float kgenptcutH=399;

const char *fopt="MEWLRQ+"; int iFit=0; 
const int knpx=2000;
float fitmin=0.00;
float fitmax=5.00;



int GetPtBin(float /*pt*/);
void MakeHist(TH1 */*hist*/,int /*istat*/);
void FillMeanSigma(int /*ip*/,TH1 */*h1*/,TH1 */*ArM*/,TH1 */*RMS*/,TH1 */*Mean*/,TH1 */*Sigma*/);

int PlotResponse(int cent=0,const char *ksp="pbpb",const char *algo="pu",int rfit=0)
{
  std::cout<<"# of pt bins : " <<nbins<<std::endl;
  std::cout<<std::endl;

  
  int iBinMin = -1;
  iBinMin = GetPtBin(kgenptcutL);
  int iBinMax = -1;
  iBinMax = GetPtBin(kgenptcutH);

  if(iBinMin<0 || iBinMax<0){
    cout<<"Something wrong ... please check the kgenpt cut : "<<endl;
    exit(0);
  }
  

  //int ksmooth=0;
  const char *reta = "eta2.0";
  float ketacut=1.6;
  if(strcmp(reta,"eta2.0")==0)ketacut=2;
  bool iSave=true;
  bool iSigma=false;
  const char *cdR="dR03";
  //const char *cdphi="2pi3";

  int minnj=0;
  int maxnj=7; 
  //double lowy =0.945;
  //double highy=1.045;

  //! Sigma
  double lowy1 =0.0;
  double highy1=0.845;

  //! Mean
  double lowy =0.045;
  double highy=2.245;

  double min=30.;
  double max=250.;
  //double max=300.;


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
  TFile *fin=0;
  //if(strcmp(ksp,"pp")==0)fin   = new TFile("input/NewCent/coarseptbins/withvtxcut/Response_newHiForest_DJ_merged_pp_2012.root","r");  
  //else fin   = new TFile("input/NewCent/coarseptbins/withvtxcut/Response_newHiForest_DJ_merged_pbpb_2012.root","r");

  if(strcmp(ksp,"pp")==0)fin   = new TFile("input/Convolution/pp/Response_newHiForest_DJ_merged_pp_2012.root","r");  
  else fin   = new TFile("input/Convolution/pbpb/Response_newHiForest_DJ_merged_pbpb_2012.root","r");

  std::cout<<"\t"<<std::endl
	   <<"rfit : "<<rfit<<"\t fitmin : "<<fitmin<<"\t fitmax : "<<fitmax<<std::endl
	   <<"Input file name   : "<<fin->GetName()<<std::endl
	   <<"\t"<<std::endl;

  //return 0;


  /*
  const int knj=25;
  const char *calgo[knj]= {
    "icPu5Calo", 
    "ak1PF","ak2PF","ak3PF","ak4PF","ak5PF","ak6PF",
    "ak1Calo","ak2Calo","ak3Calo","ak4Calo","ak5Calo","ak6Calo",
    "akPu1PF","akPu2PF","akPu3PF","akPu4PF","akPu5PF","akPu6PF",
    "akPu1Calo","akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo","akPu6Calo"
  };                                     
  */

  //for(int nj=0;nj<knj;nj++)std::cout<<"nj : " <<nj<<Form("\t %s",calgo[nj])<<"\t type : "<<mtyp[nj]<<std::endl;



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
  gStyle->SetOptStat("neMR");

  //! pp Response
  TH2F *hratiocorrrefpt[knj];
  TH2F *hratiorawrefpt[knj];
  TH1F *hratiocorrrefpt1D[knj][nbins];  
  TH1F *hratiorawrefpt1D [knj][nbins];  
  TH1F *hMean [knj], *hArM[knj], *hSigma [knj], *hRMS[knj];
  TH1F *hMean_r[knj], *hSigma_r[knj];
  TH1F *hArM_r[knj] , *hRMS_r[knj];

  TH2F *hratiocorrrefpt_lead[knj],*hratiocorrrefpt_slead[knj],*hratiocorrrefpt_remain[knj],*hratiocorrrefpt_genm[knj];
  TH1F *hratiocorrrefpt1D_lead[knj][nbins],*hratiocorrrefpt1D_slead[knj][nbins],*hratiocorrrefpt1D_remain[knj][nbins],*hratiocorrrefpt1D_genm[knj][nbins];      
  TH1F *hMean_lead [knj], *hArM_lead[knj], *hSigma_lead [knj], *hRMS_lead[knj];
  TH1F *hMean_slead [knj], *hArM_slead[knj], *hSigma_slead [knj], *hRMS_slead[knj];
  TH1F *hMean_remain [knj], *hArM_remain[knj], *hSigma_remain [knj], *hRMS_remain[knj];
  TH1F *hMean_genm [knj], *hArM_genm[knj], *hSigma_genm [knj], *hRMS_genm[knj];


  //! Integrate response
  TH1F *hratiocorrrefpt1D_int[knj];


  //! eta dependent
  TH2F *hpteta[knj][maxe];
  TH1F *hetad[knj][maxe][nbins];
  TH1F *hMepp  [knj][maxe], *hSepp  [knj][maxe];
  TH1F *hArMepp[knj][maxe], *hRMSepp[knj][maxe];
  
  //! in two eta regions |eta|<1.3 and 1.3<|eta|<2.0  
  TH2F *hratiocorrrefpt_eta[knj][2];
  TH1F *hratiocorrrefpt1D_eta[knj][2][nbins];  
  TH1F *hMean_eta [knj][2], *hSigma_eta [knj][2];
  TH1F *hArM_eta [knj][2], *hRMS_eta [knj][2];



  //! Gen Reco plots
  TH2F *hgenjrecoj[knj];

  //! pu
  TH2F *hjetptpu[knj];
  TH1F *hpmean[knj];
  TH1F *hprms [knj];

  TH2F *hjetptpu_etab[knj][2];
  TH1F *hpmean_etab[knj][2];
  TH1F *hprms_etab [knj][2];


  for(int nj=0;nj<knj;nj++){
    //std::cout<<"nj : "<<nj<<Form("\t %s",calgo[nj])<<std::endl;

    //int ksty=25; //! raw 
    //int msty=24; //! corrected

    //int cent=0;
    
    //! pp /////////////////////////////
    hMean [nj] = new TH1F(Form("hMean_%s%d_%d",ksp,nj,cent),Form("<reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hMean[nj],0);
    hArM [nj] = new TH1F(Form("hArM_%s%d_%d",ksp,nj,cent),Form("<reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hArM[nj],0);
    hSigma [nj] = new TH1F(Form("hSigma_%s%d_%d",ksp,nj,cent),Form("#sigma(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hSigma[nj],0);    
    hRMS [nj] = new TH1F(Form("hRMS_%s%d_%d",ksp,nj,cent),Form("RMS(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hRMS[nj],0);        
    
    hMean_lead [nj] = new TH1F(Form("hMean_lead_%s%d_%d",ksp,nj,cent),Form("Leading <reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hMean_lead[nj],0);
    hArM_lead [nj] = new TH1F(Form("hArM_lead_%s%d_%d",ksp,nj,cent),Form("Leading <reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hArM_lead[nj],0);
    hSigma_lead [nj] = new TH1F(Form("hSigma_lead_%s%d_%d",ksp,nj,cent),Form("Leading #sigma(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hSigma_lead[nj],0);    
    hRMS_lead [nj] = new TH1F(Form("hRMS_lead_%s%d_%d",ksp,nj,cent),Form("Leading RMS(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hRMS_lead[nj],0);        

    hMean_slead [nj] = new TH1F(Form("hMean_slead_%s%d_%d",ksp,nj,cent),Form("Sub-Leading <reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hMean_slead[nj],0);
    hArM_slead [nj] = new TH1F(Form("hArM_slead_%s%d_%d",ksp,nj,cent),Form("Sub-Leading <reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hArM_slead[nj],0);
    hSigma_slead [nj] = new TH1F(Form("hSigma_slead_%s%d_%d",ksp,nj,cent),Form("Sub-Leading #sigma(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hSigma_slead[nj],0);    
    hRMS_slead [nj] = new TH1F(Form("hRMS_slead_%s%d_%d",ksp,nj,cent),Form("Sub-Leading RMS(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hRMS_slead[nj],0);        

    hMean_remain [nj] = new TH1F(Form("hMean_remain_%s%d_%d",ksp,nj,cent),Form("Remaining <reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hMean_remain[nj],0);
    hArM_remain [nj] = new TH1F(Form("hArM_remain_%s%d_%d",ksp,nj,cent),Form("Remaining <reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hArM_remain[nj],0);
    hSigma_remain [nj] = new TH1F(Form("hSigma_remain_%s%d_%d",ksp,nj,cent),Form("Remaining #sigma(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hSigma_remain[nj],0);    
    hRMS_remain [nj] = new TH1F(Form("hRMS_remain_%s%d_%d",ksp,nj,cent),Form("Remaining RMS(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hRMS_remain[nj],0);        

    hMean_genm [nj] = new TH1F(Form("hMean_genm_%s%d_%d",ksp,nj,cent),Form("Gen Matched <reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hMean_genm[nj],0);
    hArM_genm [nj] = new TH1F(Form("hArM_genm_%s%d_%d",ksp,nj,cent),Form("Gen Matched <reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hArM_genm[nj],0);
    hSigma_genm [nj] = new TH1F(Form("hSigma_genm_%s%d_%d",ksp,nj,cent),Form("Gen Matched #sigma(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hSigma_genm[nj],0);    
    hRMS_genm [nj] = new TH1F(Form("hRMS_genm_%s%d_%d",ksp,nj,cent),Form("Gen Matched RMS(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hRMS_genm[nj],0);        


    hMean_r [nj] = new TH1F(Form("hMean_r%s%d_%d",ksp,nj,cent),Form("<raw p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hMean_r[nj],0);        
    hSigma_r [nj] = new TH1F(Form("hSigma_r%s%d_%d",ksp,nj,cent),Form("#sigma(raw p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),bins,ptbins);
    MakeHist(hSigma_r[nj],0);        
    hArM_r [nj] = new TH1F(Form("hArM_r%s%d_%d",ksp,nj,cent),Form("<raw p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),nbins,ptbins);
    MakeHist(hArM_r[nj],0);        
    hRMS_r [nj] = new TH1F(Form("hRMS_r%s%d_%d",ksp,nj,cent),Form("#sigma(raw p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),bins,ptbins);
    MakeHist(hRMS_r[nj],0);            
    
    
    //! Gen Reco scatter plots
    hgenjrecoj[nj] = (TH2F*)fin->Get(Form("hgenjrecoj%d_%d",nj,cent));
    hgenjrecoj[nj]->SetName(Form("hgenjrecoj_%s%d_%d",ksp,nj,cent));

    //! Ratio plots
    hratiocorrrefpt[nj]  = (TH2F*)fin->Get(Form("hratiocorrrefpt%d_%d",nj,cent));
    hratiocorrrefpt[nj]->SetName(Form("hratiocorrrefpt_%s%d_%d",ksp,nj,cent));
    
    hratiorawrefpt [nj]  = (TH2F*)fin->Get(Form("hratiorawrefpt%d_%d",nj,cent));
    hratiorawrefpt [nj]->SetName(Form("hratiorawrefpt_%s%d_%d",ksp,nj,cent));
    
    //! Leading jet 
    hratiocorrrefpt_lead[nj]  = (TH2F*)fin->Get(Form("hratiocorrrefpt_lead%d_%d",nj,cent));
    hratiocorrrefpt_lead[nj]->SetName(Form("hratiocorrrefpt_lead_%s%d_%d",ksp,nj,cent));

    //! Leading jet 
    hratiocorrrefpt_slead[nj]  = (TH2F*)fin->Get(Form("hratiocorrrefpt_slead%d_%d",nj,cent));
    hratiocorrrefpt_slead[nj]->SetName(Form("hratiocorrrefpt_slead_%s%d_%d",ksp,nj,cent));

    //! Remaing jet 
    hratiocorrrefpt_remain[nj]  = (TH2F*)fin->Get(Form("hratiocorrrefpt_remain%d_%d",nj,cent));
    hratiocorrrefpt_remain[nj]->SetName(Form("hratiocorrrefpt_remain_%s%d_%d",ksp,nj,cent));

    //! genmatch jet 
    hratiocorrrefpt_genm[nj]  = (TH2F*)fin->Get(Form("hratiocorrrefpt_genm%d_%d",nj,cent));
    hratiocorrrefpt_genm[nj]->SetName(Form("hratiocorrrefpt_genm_%s%d_%d",ksp,nj,cent));


    //! Integrated from 80 to 400 GeV
    hratiocorrrefpt1D_int[nj]  = (TH1F*)hratiocorrrefpt [nj]->ProjectionY(Form("hratiocorrrefpt1D_int%d_%d",nj,cent),iBinMin+1,iBinMax+1);
    MakeHist(hratiocorrrefpt1D_int[nj],1);

    for(int ip=0;ip<hratiocorrrefpt[nj]->GetNbinsX();ip++){
      hratiocorrrefpt1D[nj][ip]  = (TH1F*)hratiocorrrefpt [nj]->ProjectionY(Form("hratiocorrrefpt1D_%s%d_%d_%d",ksp,nj,cent,ip),ip+1,ip+1);
      MakeHist(hratiocorrrefpt1D[nj][ip],1);
      FillMeanSigma(ip,hratiocorrrefpt1D[nj][ip],hArM[nj],hRMS[nj],hMean[nj],hSigma[nj]);

      //! raw
      hratiorawrefpt1D[nj][ip]  = (TH1F*)hratiorawrefpt [nj]->ProjectionY(Form("hratiorawrefpt1D_%s%d_%d_%d",ksp,nj,cent,ip),ip+1,ip+1);
      MakeHist(hratiorawrefpt1D[nj][ip],1);
      FillMeanSigma(ip,hratiocorrrefpt1D[nj][ip],hArM_r[nj],hRMS_r[nj],hMean_r[nj],hSigma_r[nj]);

      hratiocorrrefpt1D_lead[nj][ip]  = (TH1F*)hratiocorrrefpt_lead [nj]->ProjectionY(Form("hratiocorrrefpt1D_lead_%s%d_%d_%d",ksp,nj,cent,ip),ip+1,ip+1);
      MakeHist(hratiocorrrefpt1D_lead[nj][ip],1);
      FillMeanSigma(ip,hratiocorrrefpt1D_lead[nj][ip],hArM_lead[nj],hRMS_lead[nj],hMean_lead[nj],hSigma_lead[nj]);

      hratiocorrrefpt1D_slead[nj][ip]  = (TH1F*)hratiocorrrefpt_slead [nj]->ProjectionY(Form("hratiocorrrefpt1D_slead_%s%d_%d_%d",ksp,nj,cent,ip),ip+1,ip+1);
      MakeHist(hratiocorrrefpt1D_slead[nj][ip],1);
      FillMeanSigma(ip,hratiocorrrefpt1D_slead[nj][ip],hArM_slead[nj],hRMS_slead[nj],hMean_slead[nj],hSigma_slead[nj]);

      hratiocorrrefpt1D_remain[nj][ip]  = (TH1F*)hratiocorrrefpt_remain [nj]->ProjectionY(Form("hratiocorrrefpt1D_remain_%s%d_%d_%d",ksp,nj,cent,ip),ip+1,ip+1);
      MakeHist(hratiocorrrefpt1D_remain[nj][ip],1);
      FillMeanSigma(ip,hratiocorrrefpt1D_remain[nj][ip],hArM_remain[nj],hRMS_remain[nj],hMean_remain[nj],hSigma_remain[nj]);

      hratiocorrrefpt1D_genm[nj][ip]  = (TH1F*)hratiocorrrefpt_genm [nj]->ProjectionY(Form("hratiocorrrefpt1D_genm_%s%d_%d_%d",ksp,nj,cent,ip),ip+1,ip+1);
      MakeHist(hratiocorrrefpt1D_genm[nj][ip],1);
      FillMeanSigma(ip,hratiocorrrefpt1D_genm[nj][ip],hArM_genm[nj],hRMS_genm[nj],hMean_genm[nj],hSigma_genm[nj]);
    }//! ip bin
      
    //! eta dependence in finer bins pp
    for(int im=0;im<maxe;im++){
      hpteta[nj][im] =(TH2F*)fin->Get(Form("hpteta%d_%d_%d",nj,cent,im));
      hpteta[nj][im]->SetName(Form("hpteta_%s%d_%d_%d",ksp,nj,cent,im));

      hMepp [nj][im] = new TH1F(Form("hMe%s%d_%d_%d",ksp,nj,cent,im),Form("<reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),bins,ptbins);
      MakeHist(hMepp[nj][im],0);            	
      hSepp [nj][im] = new TH1F(Form("hSe%s%d_%d_%d",ksp,nj,cent,im),Form("#sigma(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),bins,ptbins);
      MakeHist(hSepp[nj][im],0); 
      hArMepp [nj][im] = new TH1F(Form("hArMe%s%d_%d_%d",ksp,nj,cent,im),Form("<reco p_{T}/gen p_{T}> %s %s %d",ksp,calgo[nj],cent),bins,ptbins);
      MakeHist(hArMepp[nj][im],0);    
      hRMSepp [nj][im] = new TH1F(Form("hRMSe%s%d_%d_%d",ksp,nj,cent,im),Form("#sigma(reco p_{T}/gen p_{T}) %s %s %d",ksp,calgo[nj],cent),bins,ptbins);
      MakeHist(hRMSepp[nj][im],0);
	
      for(int il=0;il<hpteta[nj][im]->GetNbinsX();il++){
	hetad [nj][im][il] = (TH1F*)hpteta[nj][im]->ProjectionY(Form("hetad_%s%d_%d_%d_%d",ksp,nj,cent,im,il),il+1,il+1);
	hetad [nj][im][il]->SetName(Form("hetad_%s%d_%d_%d_%d",ksp,nj,cent,im,il));
	MakeHist(hetad [nj][im][il],1);
	FillMeanSigma(il,hetad [nj][im][il],hArMepp[nj][im],hRMSepp[nj][im],hMepp[nj][im],hSepp[nj][im]);
      }
    }//! finer eta bins

    //! Only to check the JEC for |eta|<1.3 and 1.3<|eta|<2.0 regions
    for(int ie=0;ie<2;ie++){
      hratiocorrrefpt_eta[nj][ie] = (TH2F*)fin->Get(Form("hratiocorrrefpt_eta%d_%d_%d",nj,cent,ie));
      hratiocorrrefpt_eta[nj][ie]->SetName(Form("hratiocorrrefpt_%s_eta%d_%d_%d",ksp,nj,cent,ie));

      hMean_eta [nj][ie] = new TH1F(Form("hMean_eta_%s%d_%d_%d",ksp,nj,cent,ie),Form("<reco p_{T}/gen p_{T}> %s %s %s",calgo[nj],ksp,ccent[cent]),bins,ptbins);
      MakeHist(hMean_eta [nj][ie],0);

      hSigma_eta [nj][ie] = new TH1F(Form("hSigma_eta_%s%d_%d_%d",ksp,nj,cent,ie),Form("#sigma(reco p_{T}/gen p_{T}) %s %s %s",calgo[nj],ksp,ccent[cent]),bins,ptbins);
      MakeHist(hSigma_eta [nj][ie],0);

      hArM_eta [nj][ie] = new TH1F(Form("hArM_eta_%s%d_%d_%d",ksp,nj,cent,ie),Form("<reco p_{T}/gen p_{T}> %s %s %s",calgo[nj],ksp,ccent[cent]),bins,ptbins);
      MakeHist(hArM_eta [nj][ie],0);	

      hRMS_eta [nj][ie] = new TH1F(Form("hRMS_eta_%s%d_%d_%d",ksp,nj,cent,ie),Form("#sigma(reco p_{T}/gen p_{T}) %s %s %s",calgo[nj],ksp,ccent[cent]),bins,ptbins);
      MakeHist(hRMS_eta [nj][ie],0);
      
      for(int ip=0;ip<hratiocorrrefpt_eta[nj][ie]->GetNbinsX();ip++){
	hratiocorrrefpt1D_eta[nj][ie][ip]  = (TH1F*)hratiocorrrefpt_eta [nj][ie]->ProjectionY(Form("hratiocorrrefpt1D_eta_%s%d_%d_%d_%d",ksp,nj,cent,ie,ip),ip+1,ip+1);
	MakeHist(hratiocorrrefpt1D_eta[nj][ie][ip],1);
	FillMeanSigma(ip,hratiocorrrefpt1D_eta[nj][ie][ip],hArM_eta[nj][ie],hRMS_eta[nj][ie],hMean_eta [nj][ie],hSigma_eta[nj][ie]);	
      }//! ip loop
    }//! ie |eta|<1.3 and 1.3<|eta|<2.0
    
    //! Pile up subtraction histograms
    hjetptpu  [nj] = (TH2F*)fin->Get(Form("hjetptpu%d_%d",nj,cent));
    hjetptpu  [nj]->SetName(Form("hjetptpu_%s_%s_%s",ksp,calgo[nj],ccent[cent]));

    //! pileup mean from MC
    hpmean[nj] = new TH1F(Form("hpmean%s_%s_%s",ksp,calgo[nj],ccent[cent]),Form("hpmean%s_%s_%s",ksp,calgo[nj],ccent[cent]),dbins,ptbins_data);
    hpmean[nj]->Sumw2();
    hprms[nj] = new TH1F(Form("hprms%s_%s_%s",ksp,calgo[nj],ccent[cent]),Form("hprms%s_%s_%s",ksp,calgo[nj],ccent[cent]),dbins,ptbins_data);
    hprms[nj]->Sumw2();

    for(int ix=1;ix<=hjetptpu[nj]->GetNbinsX();ix++){
      double binwx = hjetptpu[nj]->GetXaxis()->GetBinWidth(ix);

      for(int iy=1;iy<=hjetptpu[nj]->GetNbinsY();iy++){
	double binwy = hjetptpu[nj]->GetYaxis()->GetBinWidth(iy);
	double binc  = hjetptpu[nj]->GetBinContent(ix,iy)/binwx/binwy;
	double bine  = hjetptpu[nj]->GetBinError(ix,iy)/binwx/binwy;
	hjetptpu[nj]->SetBinContent(ix,iy,binc);
	hjetptpu[nj]->SetBinError(ix,iy,bine);
      }//! iy

      //! pile up pt distribution
      TH1F *htpileup = (TH1F*)hjetptpu[nj]->ProjectionY("htpileup",ix,ix);
      htpileup->SetName("htpileup");
      htpileup->Scale(binwx);
      
      hpmean[nj]->SetBinContent(ix,htpileup->GetMean());
      hpmean[nj]->SetBinError(ix,htpileup->GetMeanError());
      
      hprms[nj]->SetBinContent(ix,htpileup->GetRMS());
      hprms[nj]->SetBinError(ix,htpileup->GetRMSError());
      delete htpileup;
    }//! ix
    
    //! eta dependence of pile up subtraction
    for(int ie=0;ie<2;ie++){
      hjetptpu_etab  [nj][ie] = (TH2F*)fin->Get(Form("hjetptpu_etab%d_%d_%d",nj,cent,ie));
      hjetptpu_etab  [nj][ie]->SetName(Form("hjetptpu_etab_%s_%s_%s_%d",ksp,calgo[nj],ccent[cent],ie));

      //! pileup mean from MC
      hpmean_etab[nj][ie] = new TH1F(Form("hpmean_etab%s_%s_%s_%d",ksp,calgo[nj],ccent[cent],ie),Form("hpmean_etab_%s_%s_%s_%d",ksp,calgo[nj],ccent[cent],ie),dbins,ptbins_data);
      hpmean_etab[nj][ie]->Sumw2();
      hprms_etab[nj][ie] = new TH1F(Form("hprms_etab%s_%s_%s_%d",ksp,calgo[nj],ccent[cent],ie),Form("hprms_etab_%s_%s_%s_%d",ksp,calgo[nj],ccent[cent],ie),dbins,ptbins_data);
      hprms_etab[nj][ie]->Sumw2();

      for(int ix=1;ix<=hjetptpu_etab[nj][ie]->GetNbinsX();ix++){
	double binwx = hjetptpu_etab[nj][ie]->GetXaxis()->GetBinWidth(ix);
	
	for(int iy=1;iy<=hjetptpu_etab[nj][ie]->GetNbinsY();iy++){
	  double binwy = hjetptpu_etab[nj][ie]->GetYaxis()->GetBinWidth(iy);
	  double binc  = hjetptpu_etab[nj][ie]->GetBinContent(ix,iy)/binwx/binwy;
	  double bine  = hjetptpu_etab[nj][ie]->GetBinError(ix,iy)/binwx/binwy;
	  hjetptpu_etab[nj][ie]->SetBinContent(ix,iy,binc);
	  hjetptpu_etab[nj][ie]->SetBinError(ix,iy,bine);
	}//! iy

      	//! pile up pt distribution
	TH1F *htpileup = (TH1F*)hjetptpu[nj]->ProjectionY("htpileup",ix,ix);
	htpileup->SetName("htpileup");
	htpileup->Scale(binwx);
	
	hpmean_etab[nj][ie]->SetBinContent(ix,htpileup->GetMean());
	hpmean_etab[nj][ie]->SetBinError(ix,htpileup->GetMeanError());
	
	hprms_etab[nj][ie]->SetBinContent(ix,htpileup->GetRMS());
	hprms_etab[nj][ie]->SetBinError(ix,htpileup->GetRMSError());
	delete htpileup;
      }
    }//! ie eta loop for pu
  }//! nj loop ends


  //! smooth the histograms 
  /*
  for(int nj=0;nj<knj;nj++){
    hSigma[nj]->Smooth(ksmooth);
    hRMS  [nj]->Smooth(ksmooth);
    hMean [nj]->Smooth(ksmooth);      
    hArM  [nj]->Smooth(ksmooth);      
    for(int im=0;im<maxe;im++){
      hMepp  [nj][im]->Smooth(ksmooth);
      hSepp  [nj][im]->Smooth(ksmooth);
      hRMSepp[nj][im]->Smooth(ksmooth);
      hArMepp[nj][im]->Smooth(ksmooth);
    }
    for(int ie=0;ie<2;ie++){
      hMean_eta  [nj][ie]->Smooth(ksmooth);
      hSigma_eta [nj][ie]->Smooth(ksmooth);
      hRMS_eta   [nj][ie]->Smooth(ksmooth);
      hArM_eta   [nj][ie]->Smooth(ksmooth);
    }
  }
  */

  int ipad=0;
  int sty[6] = {24,25,26,28,30,32};
  int col[6] = { 2, 4, 6, 8,48,44};
  int ic=0;
  int maxc=3;
  int maxr=2;

  //hratiocorrrefpt1D[4][15]->Draw("p");
  //return 0;

  ipad=0;
  TCanvas *c99[knj];
  for(int nj=5;nj<6;nj++){
    c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%s Fitting plots",calgo[nj]),98,46,1215,930);
    c99[nj]->Divide(4,4,0,0);
    ipad=0;
    for(int ip=0;ip<nbins;ip++){      
      //if( hratiocorrrefpt1D[nj][ip]->Integral()==0)continue;
      c99[nj]->cd(++ipad);
      if(ipad%4==0)gPad->SetRightMargin(0.02);
      gPad->SetBottomMargin(0.15);
      //if(ipad>20)gPad->SetBottomMargin(0.15);

      hratiocorrrefpt1D[nj][ip]->SetMaximum(0.634);
      hratiocorrrefpt1D[nj][ip]->SetMinimum(-0.008);
      hratiocorrrefpt1D[nj][ip]->SetTitle(0);
      hratiocorrrefpt1D[nj][ip]->GetXaxis()->SetTitleFont(42);
      hratiocorrrefpt1D[nj][ip]->GetXaxis()->SetLabelFont(42);
      hratiocorrrefpt1D[nj][ip]->GetXaxis()->SetLabelSize(0.08);
      hratiocorrrefpt1D[nj][ip]->GetXaxis()->SetTitleSize(0.07);

      hratiocorrrefpt1D[nj][ip]->GetYaxis()->SetTitle("");
      hratiocorrrefpt1D[nj][ip]->GetYaxis()->SetTitleFont(42);
      hratiocorrrefpt1D[nj][ip]->GetYaxis()->SetLabelFont(42);
      hratiocorrrefpt1D[nj][ip]->GetYaxis()->SetLabelSize(0.08);
      hratiocorrrefpt1D[nj][ip]->Draw("p");  

      c99[nj]->Update();

      TPaveStats *ps = (TPaveStats*)  hratiocorrrefpt1D[nj][ip]->GetListOfFunctions()->FindObject("stats");
      ps->SetX1NDC(0.44);
      ps->SetY1NDC(0.25);       
      ps->SetX2NDC(0.95);
      ps->SetY2NDC(0.79);
      ps->SetTextFont(42);
      ps->Draw();

      
      TPaveText *pt   = new TPaveText(0.4784302,0.8368631,0.8946255,0.9521693,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(10);
      pt->SetTextFont(42);
      TText *text = pt->AddText(Form("%0.0f < p_{T} (GeV/c) < %0.0f",ptbins[ip], ptbins[ip+1]));
      text->SetTextSize(0.07);
      pt->Draw();

      if(ipad==1){
	TPaveText *pt1 = new TPaveText(0.12,0.65,0.42,0.96,"brNDC");
	pt1->SetBorderSize(0);
	pt1->SetFillColor(10);
	pt1->SetTextFont(42);
	TText *text1 = pt1->AddText(calgo[nj]);
	text1->SetTextSize(0.14);	
	pt1->Draw();
      }
    }
    c99[nj]->Update();
    //if(iSave)c99[nj]->SaveAs(Form("plots/Fitting_%s_%s.png",cdR,calgo[nj]));
  }
  //return 0;
  

  TFile *fout = new TFile(Form("rootfiles/JERJES/Histo_%s_%s.root",ksp,ccent[cent]),"RECREATE");
  fout->cd();
  for(int nj=0;nj<knj;nj++){


    hSigma[nj]->SetName(Form("hSigma_%s_%s_%d",ksp,calgo[nj],cent));
    hSigma[nj]->Write();
    hRMS[nj]->SetName(Form("hRMS_%s_%s_%d",ksp,calgo[nj],cent));
    hRMS  [nj]->Write();
    hMean[nj]->SetName(Form("hMean_%s_%s_%d",ksp,calgo[nj],cent));
    hMean [nj]->Write();
    hArM[nj]->SetName(Form("hArM_%s_%s_%d",ksp,calgo[nj],cent));
    hArM  [nj]->Write();

    hSigma_lead[nj]->SetName(Form("hSigma_lead_%s_%s_%d",ksp,calgo[nj],cent));
    hSigma_lead[nj]->Write();
    hRMS_lead[nj]->SetName(Form("hRMS_lead_%s_%s_%d",ksp,calgo[nj],cent));
    hRMS_lead  [nj]->Write();
    hMean_lead[nj]->SetName(Form("hMean_lead_%s_%s_%d",ksp,calgo[nj],cent));
    hMean_lead [nj]->Write();
    hArM_lead[nj]->SetName(Form("hArM_lead_%s_%s_%d",ksp,calgo[nj],cent));
    hArM_lead  [nj]->Write();

    hSigma_slead[nj]->SetName(Form("hSigma_slead_%s_%s_%d",ksp,calgo[nj],cent));
    hSigma_slead[nj]->Write();
    hRMS_slead[nj]->SetName(Form("hRMS_slead_%s_%s_%d",ksp,calgo[nj],cent));
    hRMS_slead  [nj]->Write();
    hMean_slead[nj]->SetName(Form("hMean_slead_%s_%s_%d",ksp,calgo[nj],cent));
    hMean_slead [nj]->Write();
    hArM_slead[nj]->SetName(Form("hArM_slead_%s_%s_%d",ksp,calgo[nj],cent));
    hArM_slead  [nj]->Write();

    hSigma_remain[nj]->SetName(Form("hSigma_remain_%s_%s_%d",ksp,calgo[nj],cent));
    hSigma_remain[nj]->Write();
    hRMS_remain[nj]->SetName(Form("hRMS_remain_%s_%s_%d",ksp,calgo[nj],cent));
    hRMS_remain  [nj]->Write();
    hMean_remain[nj]->SetName(Form("hMean_remain_%s_%s_%d",ksp,calgo[nj],cent));
    hMean_remain [nj]->Write();
    hArM_remain[nj]->SetName(Form("hArM_remain_%s_%s_%d",ksp,calgo[nj],cent));
    hArM_remain  [nj]->Write();

    hSigma_genm[nj]->SetName(Form("hSigma_genm_%s_%s_%d",ksp,calgo[nj],cent));
    hSigma_genm[nj]->Write();
    hRMS_genm[nj]->SetName(Form("hRMS_genm_%s_%s_%d",ksp,calgo[nj],cent));
    hRMS_genm  [nj]->Write();
    hMean_genm[nj]->SetName(Form("hMean_genm_%s_%s_%d",ksp,calgo[nj],cent));
    hMean_genm [nj]->Write();
    hArM_genm[nj]->SetName(Form("hArM_genm_%s_%s_%d",ksp,calgo[nj],cent));
    hArM_genm  [nj]->Write();


    hArM_r[nj]->SetName(Form("hArM_r_%s_%s_%d",ksp,calgo[nj],cent));
    hArM_r  [nj]->Write();
    hRMS_r[nj]->SetName(Form("hRMS_r_%s_%s_%d",ksp,calgo[nj],cent));
    hRMS_r  [nj]->Write();
    hMean_r[nj]->SetName(Form("hMean_r_%s_%s_%d",ksp,calgo[nj],cent));
    hMean_r [nj]->Write();
    hSigma_r[nj]->SetName(Form("hSigma_r_%s_%s_%d",ksp,calgo[nj],cent));
    hSigma_r[nj]->Write();


//
//    for(int im=0;im<maxe;im++){
//    hMepp  [nj][im]->Write();
//    hSepp  [nj][im]->Write();
//    hRMSepp[nj][im]->Write();
//    hArMepp[nj][im]->Write();
//    }
//
    for(int ie=0;ie<2;ie++){
      hMean_eta [nj][ie]->SetName(Form("hMean_eta_%s_%s_%d_%d",ksp,calgo[nj],cent,ie));      
      hMean_eta [nj][ie]->Write();
      hSigma_eta[nj][ie]->SetName(Form("hSigma_eta_%s_%s_%d_%d",ksp,calgo[nj],cent,ie));      
      hSigma_eta[nj][ie]->Write();
      hRMS_eta  [nj][ie]->SetName(Form("hRMS_eta_%s_%s_%d_%d",ksp,calgo[nj],cent,ie));      
      hRMS_eta  [nj][ie]->Write();
      hArM_eta  [nj][ie]->SetName(Form("hArM_eta_%s_%s_%d_%d",ksp,calgo[nj],cent,ie));      
      hArM_eta  [nj][ie]->Write();
    }

    hpmean[nj]->SetName(Form("hpmean_%s_%s_%d",ksp,calgo[nj],cent));
    hpmean[nj]->Write();
    hprms[nj]->SetName(Form("hprms_%s_%s_%d",ksp,calgo[nj],cent));
    hprms [nj]->Write();
    for(int ie=0;ie<2;ie++){
      hpmean_etab[nj][ie]->SetName(Form("hpmean_eta_%s_%s_%d_%d",ksp,calgo[nj],cent,ie));      
      hpmean_etab[nj][ie]->Write();
      hprms_etab [nj][ie]->SetName(Form("hprms_eta_%s_%s_%d_%d",ksp,calgo[nj],cent,ie));      
      hprms_etab [nj][ie]->Write();
    }

    hratiocorrrefpt[nj]->SetName(Form("hratiocorrrefpt_%s_%s_%d",ksp,calgo[nj],cent));   
    hratiocorrrefpt[nj]->Write();
    hratiocorrrefpt1D_int[nj]->SetName(Form("hratiocorrrefpt1D_int_%s_%s_%d",ksp,calgo[nj],cent));   
    hratiocorrrefpt1D_int[nj]->Write();
    hgenjrecoj[nj]->SetName(Form("hgenjrecoj_%s_%s_%d",ksp,calgo[nj],cent));   
    hgenjrecoj[nj]->Write();
  }
  fout->Close();
  return 0;


  
  TLine *line = new TLine(min,1.0,800,1.0);
  line->SetLineWidth(1);
  line->SetLineStyle(2);

  TCanvas *c3 = new TCanvas("c3","pp sigma and mean",92,95,1143,689);
  makeMultiPanelCanvas(c3,maxc,maxr,0.0,0.0,0.22,0.22,0.02);
  std::cout<<std::endl;
  ipad=0;  
  TLegend *l3[maxc];
  int ij=-1;
  ic=0;
  for(int nj=1;nj<7;nj++){
    
    if((nj-1)%3==0){
      c3->cd(++ipad);
      //gPad->SetLogx();
      
      ic=0;

      //l3[++ij] = new TLegend(0.72,0.44,0.97,0.94,NULL,"BRNDC");
      l3[++ij] = new TLegend(0.6452845,0.3598836,0.8942526,0.7199491,NULL,"BRNDC");
      l3[ij]->SetHeader("");
      l3[ij]->SetBorderSize(0);
      l3[ij]->SetTextFont(42);
      l3[ij]->SetTextSize(0.05);
      l3[ij]->SetLineColor(1);
      l3[ij]->SetLineStyle(1);
      l3[ij]->SetLineWidth(1);
      l3[ij]->SetFillColor(10);
      l3[ij]->SetFillStyle(1001);
      l3[ij]->SetHeader("");			  

      if(iSigma){
	hSigma  [nj]->SetMarkerStyle(sty[ic]);
	hSigma  [nj]->SetMarkerColor(col[ic]);
	hSigma  [nj]->SetLineColor(col[ic]);
	hSigma  [nj]->SetMarkerSize(1.3);
	hSigma  [nj]->Draw("p");
	l3[ij]->AddEntry(hSigma[nj] ,calgo[nj],"p"); 

	//std::cout<<"nj : "<<nj<<"\t ipad : "<< ipad<<"\t im : "<<ij<<std::endl;
      }else{
	hRMS  [nj]->SetMaximum(highy1);
	hRMS  [nj]->SetMinimum(lowy1);	
	hRMS  [nj]->SetMarkerStyle(sty[ic]);
	hRMS  [nj]->SetMarkerColor(col[ic]);
	hRMS  [nj]->SetLineColor(col[ic]);
	hRMS  [nj]->SetMarkerSize(1.3);
	hRMS  [nj]->Draw("p");
	l3[ij]->AddEntry(hRMS[nj] ,calgo[nj],"p"); 
      }
      l3[ij]->Draw();

      c3->cd(ipad+maxc);
      //gPad->SetLogx();
      if(iSigma){
	hMean  [nj]->SetMarkerStyle(sty[ic]);
	hMean  [nj]->SetMarkerColor(col[ic]);
	hMean  [nj]->SetLineColor(col[ic]);
	hMean  [nj]->SetMarkerSize(1.3);
	hMean  [nj]->Draw("p");
      }else{
	hArM  [nj]->SetMaximum(highy);
	hArM  [nj]->SetMinimum(lowy);
	hArM  [nj]->SetMarkerStyle(sty[ic]);
	hArM  [nj]->SetMarkerColor(col[ic]);
	hArM  [nj]->SetLineColor(col[ic]);
	hArM  [nj]->SetMarkerSize(1.3);
	hArM  [nj]->Draw("p");
      }
      line->Draw();
      ic++;

    }else{

      c3->cd(ipad);
      //gPad->SetLogx();
      if(iSigma){
	hSigma  [nj]->SetMarkerStyle(sty[ic]);
	hSigma  [nj]->SetMarkerColor(col[ic]);
	hSigma  [nj]->SetLineColor(col[ic]);
	hSigma  [nj]->SetMarkerSize(1.3);
	hSigma  [nj]->Draw("psame");
	l3[ij]->AddEntry(hSigma[nj] ,calgo[nj],"p"); 
      }else{
	hRMS[nj]->SetMaximum(highy1);
	hRMS[nj]->SetMinimum(lowy1);
	hRMS  [nj]->SetMarkerStyle(sty[ic]);
	hRMS  [nj]->SetMarkerColor(col[ic]);
	hRMS  [nj]->SetLineColor(col[ic]);
	hRMS  [nj]->SetMarkerSize(1.3);
	hRMS  [nj]->Draw("psame");
	l3[ij]->AddEntry(hRMS[nj] ,calgo[nj],"p"); 
      }
      //std::cout<<"\t \t nj : "<<nj<<"\t ipad : "<< ipad<<std::endl;

      c3->cd(ipad+maxc);
      //gPad->SetLogx();
      if(iSigma){

	hMean[nj]->SetMarkerStyle(sty[ic]);
	hMean[nj]->SetMarkerColor(col[ic]);
	hMean[nj]->SetLineColor(col[ic]);
	hMean[nj]->SetMarkerSize(1.3);
	hMean[nj]->Draw("psame");

      }else{
	hArM[nj]->SetMaximum(highy);
	hArM[nj]->SetMinimum(lowy);
	hArM[nj]->SetMarkerStyle(sty[ic]);
	hArM[nj]->SetMarkerColor(col[ic]);
	hArM[nj]->SetLineColor(col[ic]);
	hArM[nj]->SetMarkerSize(1.3);
	hArM[nj]->Draw("psame");
      }
      line->Draw();
      ic++;
    }
  }
  c3->cd(1);
  TPaveText *pt3 = new TPaveText(0.3935117,0.7644594,0.68073,0.9419683,"brNDC");
  pt3->SetBorderSize(0);
  pt3->SetFillColor(10);
  pt3->SetTextFont(42);
  TText *text3 = pt3->AddText("CMS Simulation");
  text3->SetTextSize(0.08);
  TText *text4 = 0;
  if(strcmp(ksp,"pp")==0)text4 = pt3->AddText("PYTHIA Z2");
  else text4 = pt3->AddText("PYTHIA Z2+HYDJET1.8");
  text4->SetTextSize(0.08);
  pt3->Draw();

  c3->cd(2);
  TPaveText *pt4 = new TPaveText(0.3935117,0.7644594,0.68073,0.9419683,"brNDC");
  pt4->SetBorderSize(0);
  pt4->SetFillColor(10);
  pt4->SetTextFont(42);
  TText *text5 = 0;
  if(strcmp(ksp,"pp")==0)text5 = pt4->AddText(Form("%s, |#eta|<%0.0f",ccent[cent],ketacut));
  else text5 = pt4->AddText(Form("PbPb %s, |#eta|<%0.0f",ccent[cent],ketacut));
  text5->SetTextSize(0.08);
  pt4->Draw();

  //return 0;

  c3->cd(maxc);
  //gPad->SetLogx();

  hSigma[0]->SetMaximum(highy1);
  hSigma[0]->SetMinimum(lowy1);
  hSigma[0]->SetTitle("");
  hSigma[0]->GetXaxis()->SetRangeUser(min,max);
  hSigma[0]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  hSigma[0]->GetXaxis()->CenterTitle(true);
  hSigma[0]->GetXaxis()->SetTitleFont(42);
  hSigma[0]->GetXaxis()->SetLabelFont(42);
  hSigma[0]->GetXaxis()->SetMoreLogLabels();
  hSigma[0]->GetXaxis()->SetNoExponent();
  hSigma[0]->GetXaxis()->SetTitleSize(0.07);
  hSigma[0]->GetXaxis()->SetTitleOffset(1.18);
  hSigma[0]->GetXaxis()->SetLabelSize(0.07);
  hSigma[0]->GetXaxis()->SetNdivisions(507);
  hSigma[0]->GetYaxis()->SetTitle("#sigma (RecoJet p_{T} / GenJet p_{T})");    
  hSigma[0]->GetYaxis()->SetTitleSize(0.07);
  hSigma[0]->GetYaxis()->SetTitleOffset(1.26);
  hSigma[0]->GetYaxis()->SetLabelSize(0.08);
  hSigma[0]->GetYaxis()->SetNdivisions(507);
  hSigma[0]->GetYaxis()->SetTitleFont(42);
  hSigma[0]->GetYaxis()->SetLabelFont(42);
  hSigma[0]->SetMarkerStyle(20);
  hSigma[0]->SetMarkerColor(1);
  hSigma[0]->SetLineColor(1);
  hSigma[0]->SetMarkerSize(1.2);

  hRMS[0]->SetMaximum(highy1);
  hRMS[0]->SetMinimum(lowy1);
  hRMS[0]->SetTitle("");
  hRMS[0]->GetXaxis()->SetRangeUser(min,max);
  hRMS[0]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  hRMS[0]->GetXaxis()->CenterTitle(true);
  hRMS[0]->GetXaxis()->SetTitleFont(42);
  hRMS[0]->GetXaxis()->SetLabelFont(42);
  hRMS[0]->GetXaxis()->SetMoreLogLabels();
  hRMS[0]->GetXaxis()->SetNoExponent();
  hRMS[0]->GetXaxis()->SetTitleSize(0.07);
  hRMS[0]->GetXaxis()->SetTitleOffset(1.18);
  hRMS[0]->GetXaxis()->SetLabelSize(0.07);
  hRMS[0]->GetXaxis()->SetNdivisions(507);
  hRMS[0]->GetYaxis()->SetTitle("#sigma (RecoJet p_{T} / GenJet p_{T})");    
  hRMS[0]->GetYaxis()->SetTitleSize(0.07);
  hRMS[0]->GetYaxis()->SetTitleOffset(1.26);
  hRMS[0]->GetYaxis()->SetLabelSize(0.08);
  hRMS[0]->GetYaxis()->SetNdivisions(507);
  hRMS[0]->GetYaxis()->SetTitleFont(42);
  hRMS[0]->GetYaxis()->SetLabelFont(42);

  hRMS[0]->SetMarkerStyle(20);
  hRMS[0]->SetMarkerColor(1);
  hRMS[0]->SetLineColor(1);
  hRMS[0]->SetMarkerSize(1.2);

  l3[++ij] = new TLegend(0.60,0.45,0.90,0.95,NULL,"BRNDC");
  l3[ij]->SetHeader("");
  l3[ij]->SetBorderSize(0);
  l3[ij]->SetTextFont(42);
  l3[ij]->SetTextSize(0.05);
  l3[ij]->SetLineColor(1);
  l3[ij]->SetLineStyle(1);
  l3[ij]->SetLineWidth(1);
  l3[ij]->SetFillColor(10);
  l3[ij]->SetFillStyle(1001);
  l3[ij]->SetHeader("");			  

  if(iSigma){
    hSigma [0]->Draw("p");
    l3[ij]->AddEntry(hSigma[0] ,calgo[0],"p"); 
  }else{
    hRMS [0]->Draw("p");
    l3[ij]->AddEntry(hRMS[0] ,calgo[0],"p"); 
  }
  l3[ij]->Draw();


  c3->cd(2*maxc);
  //gPad->SetLogx();
  hArM[0]->SetMaximum(highy);
  hArM[0]->SetMinimum(lowy);
  hArM[0]->SetTitle("");
  hArM[0]->GetXaxis()->SetRangeUser(min,max);
  hArM[0]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  hArM[0]->GetXaxis()->CenterTitle(true);
  hArM[0]->GetXaxis()->SetMoreLogLabels();
  hArM[0]->GetXaxis()->SetNoExponent();
  hArM[0]->GetXaxis()->SetTitleFont(42);
  hArM[0]->GetXaxis()->SetLabelFont(42);
  hArM[0]->GetXaxis()->SetTitleSize(0.06);
  hArM[0]->GetXaxis()->SetTitleOffset(1.15);
  hArM[0]->GetXaxis()->SetLabelSize(0.06);
  hArM[0]->GetXaxis()->SetLabelOffset(0.005);
  hArM[0]->GetXaxis()->SetNdivisions(507);
  hArM[0]->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
  hArM[0]->GetYaxis()->SetTitleSize(0.06);
  hArM[0]->GetYaxis()->SetTitleOffset(1.50);
  hArM[0]->GetYaxis()->SetLabelSize(0.06);
  hArM[0]->GetYaxis()->SetNdivisions(507);
  hArM[0]->GetYaxis()->SetDecimals(true);
  hArM[0]->GetYaxis()->SetTitleFont(42);
  hArM[0]->GetYaxis()->SetLabelFont(42);
  hArM [0]->SetMarkerStyle(20);
  hArM [0]->SetMarkerColor(4);
  hArM [0]->SetLineColor(4);
  hArM [0]->SetMarkerSize(1.2);

  hMean[0]->SetMaximum(highy);
  hMean[0]->SetMinimum(lowy);
  hMean[0]->SetTitle("");
  hMean[0]->GetXaxis()->SetRangeUser(min,max);
  hMean[0]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  hMean[0]->GetXaxis()->CenterTitle(true);
  hMean[0]->GetXaxis()->SetMoreLogLabels();
  hMean[0]->GetXaxis()->SetNoExponent();
  hMean[0]->GetXaxis()->SetTitleFont(42);
  hMean[0]->GetXaxis()->SetLabelFont(42);
  hMean[0]->GetXaxis()->SetTitleSize(0.06);
  hMean[0]->GetXaxis()->SetTitleOffset(1.15);
  hMean[0]->GetXaxis()->SetLabelSize(0.06);
  hMean[0]->GetXaxis()->SetLabelOffset(0.005);
  hMean[0]->GetXaxis()->SetNdivisions(507);
  hMean[0]->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
  hMean[0]->GetYaxis()->SetTitleSize(0.06);
  hMean[0]->GetYaxis()->SetTitleOffset(1.50);
  hMean[0]->GetYaxis()->SetLabelSize(0.06);
  hMean[0]->GetYaxis()->SetNdivisions(507);
  hMean[0]->GetYaxis()->SetDecimals(true);
  hMean[0]->GetYaxis()->SetTitleFont(42);
  hMean[0]->GetYaxis()->SetLabelFont(42);
  hMean[0]->SetMarkerStyle(20);
  hMean[0]->SetMarkerColor(4);
  hMean[0]->SetLineColor(4);
  hMean[0]->SetMarkerSize(1.2);

  if(iSigma)hMean[0]->Draw("p");  
  else hArM[0]->Draw("p");
  line->Draw();

  if(iSave)c3->SaveAs(Form("plots/JEC_%s_%s_%s_2012Forest.png",algo,ccent[cent],cdR));

  fout->Close();
  return 0;




  //TLine *line = new TLine(min,1.0,max+350,1.0);

  for(int nj=0;nj<knj;nj++){
    hprms[nj]->SetStats(0);
    hprms[nj]->SetTitle("");

    //! 0-5%
    //hprms[nj]->SetMaximum(28.65);
    //hprms[nj]->SetMinimum(0.01);

    //! pp
    //hprms[nj]->SetMaximum(6.65);
    hprms[nj]->SetMaximum(3.65);
    hprms[nj]->SetMinimum(0);

    hprms[nj]->GetXaxis()->SetRangeUser(min,max);
    hprms[nj]->GetXaxis()->SetTitle("Reco Jet p_{T} (GeV/c)");
    hprms[nj]->GetXaxis()->CenterTitle(true);
    hprms[nj]->GetXaxis()->SetTitleFont(42);
    hprms[nj]->GetXaxis()->SetLabelFont(42);
    hprms[nj]->GetXaxis()->SetMoreLogLabels();
    hprms[nj]->GetXaxis()->SetNoExponent();
    hprms[nj]->GetXaxis()->SetTitleSize(0.08);
    hprms[nj]->GetXaxis()->SetTitleOffset(1.15);
    hprms[nj]->GetXaxis()->SetLabelSize(0.08);
    hprms[nj]->GetXaxis()->SetNdivisions(507);
    hprms[nj]->GetYaxis()->SetTitle("#sigma (pu)");    
    hprms[nj]->GetYaxis()->SetTitleSize(0.08);
    hprms[nj]->GetYaxis()->SetTitleOffset(1.26);
    hprms[nj]->GetYaxis()->SetLabelSize(0.08);
    hprms[nj]->GetYaxis()->SetNdivisions(507);
    hprms[nj]->GetYaxis()->SetTitleFont(42);
    hprms[nj]->GetYaxis()->SetLabelFont(42);

    //hpmean[nj]->SetMaximum(100.05);
    //hpmean[nj]->SetMinimum(0);

    //! pp
    //hpmean[nj]->SetMaximum(11.05);
    hpmean[nj]->SetMaximum(4.05);
    hpmean[nj]->SetMinimum(0);

    hpmean[nj]->SetStats(0);
    hpmean[nj]->SetTitle("");
    hpmean[nj]->GetXaxis()->SetRangeUser(min,max);
    hpmean[nj]->GetXaxis()->SetTitle("Reco Jet p_{T} (GeV/c)");
    hpmean[nj]->GetXaxis()->CenterTitle(true);
    hpmean[nj]->GetXaxis()->SetMoreLogLabels();
    hpmean[nj]->GetXaxis()->SetNoExponent();
    hpmean[nj]->GetXaxis()->SetTitleFont(42);
    hpmean[nj]->GetXaxis()->SetLabelFont(42);
    hpmean[nj]->GetXaxis()->SetTitleSize(0.08);
    hpmean[nj]->GetXaxis()->SetTitleOffset(1.15);
    hpmean[nj]->GetXaxis()->SetLabelSize(0.08);
    hpmean[nj]->GetXaxis()->SetLabelOffset(0.005);
    hpmean[nj]->GetXaxis()->SetNdivisions(507);
    hpmean[nj]->GetYaxis()->SetTitle("<pu>");
    hpmean[nj]->GetYaxis()->SetTitleSize(0.08);
    hpmean[nj]->GetYaxis()->SetTitleOffset(1.50);
    hpmean[nj]->GetYaxis()->SetLabelSize(0.08);
    hpmean[nj]->GetYaxis()->SetNdivisions(507);
    hpmean[nj]->GetYaxis()->SetDecimals(true);
    hpmean[nj]->GetYaxis()->SetTitleFont(42);
    hpmean[nj]->GetYaxis()->SetLabelFont(42);
  }


  TCanvas *c10 = new TCanvas("c10","pu ",92,95,1143,689);
  makeMultiPanelCanvas(c10,maxc,maxr,0.0,0.0,0.22,0.22,0.02);
  std::cout<<std::endl;
  ipad=0;  
  TLegend *l10[maxc];
  ij=-1;
  ic=0;
  for(int nj=minnj;nj<maxnj;nj++){
    
    if((nj-1)%6==0){
      c10->cd(++ipad);
      gPad->SetLogx();
      
      ic=0;

      l10[++ij] = new TLegend(0.2555564,0.6330544,0.5043282,0.9916947,NULL,"BRNDC");
      l10[ij]->SetHeader("");
      l10[ij]->SetBorderSize(0);
      l10[ij]->SetTextFont(42);
      l10[ij]->SetTextSize(0.05);
      l10[ij]->SetLineColor(1);
      l10[ij]->SetLineStyle(1);
      l10[ij]->SetLineWidth(1);
      l10[ij]->SetFillColor(10);
      l10[ij]->SetFillStyle(1001);
      l10[ij]->SetHeader("");			  


      hpmean  [nj]->SetMarkerStyle(sty[ic]);
      hpmean  [nj]->SetMarkerColor(col[ic]);
      hpmean  [nj]->SetLineColor(col[ic]);
      hpmean  [nj]->SetMarkerSize(1.3);
      hpmean  [nj]->Draw("p");
      l10[ij]->AddEntry(hpmean[nj] ,calgo[nj],"p"); 
      
      //std::cout<<"nj : "<<nj<<"\t ipad : "<< ipad<<"\t im : "<<ij<<std::endl;


      c10->cd(ipad+maxc);
      gPad->SetLogx();

      hprms  [nj]->SetMarkerStyle(sty[ic]);
      hprms  [nj]->SetMarkerColor(col[ic]);
      hprms  [nj]->SetLineColor(col[ic]);
      hprms  [nj]->SetMarkerSize(1.3);
      hprms  [nj]->Draw("p");

      l10[ij]->Draw();
      ic++;

    }else{

      c10->cd(ipad);
      gPad->SetLogx();
      hpmean[nj]->SetMarkerStyle(sty[ic]);
      hpmean[nj]->SetMarkerColor(col[ic]);
      hpmean[nj]->SetLineColor(col[ic]);
      hpmean[nj]->SetMarkerSize(1.3);
      hpmean[nj]->Draw("psame");
      //line->Draw();
      l10[ij]->AddEntry(hpmean[nj] ,calgo[nj],"p"); 

      //std::cout<<"\t \t nj : "<<nj<<"\t ipad : "<< ipad<<std::endl;

      c10->cd(ipad+maxc);
      gPad->SetLogx();
      hprms  [nj]->SetMarkerStyle(sty[ic]);
      hprms  [nj]->SetMarkerColor(col[ic]);
      hprms  [nj]->SetLineColor(col[ic]);
      hprms  [nj]->SetMarkerSize(1.3);
      hprms  [nj]->Draw("psame");
      ic++;
    }
  }
  c10->cd(1);
  TPaveText *pt10 = new TPaveText(0.3935117,0.7644594,0.68073,0.9419683,"brNDC");
  pt10->SetBorderSize(0);
  pt10->SetFillColor(10);
  pt10->SetTextFont(42);
  TText *text10 = pt10->AddText("CMS Simulation");
  text10->SetTextSize(0.08);
  TText *text11 = 0;
  if(strcmp(ksp,"pp")==0)text11 = pt10->AddText("PYTHIA Z2");
  else text11 = pt10->AddText("PYTHIA Z2+HYDJET1.8");
  text11->SetTextSize(0.08);
  pt10->Draw();

  c10->cd(2);
  TPaveText *pt12 = new TPaveText(0.3935117,0.7644594,0.68073,0.9419683,"brNDC");
  pt12->SetBorderSize(0);
  pt12->SetFillColor(10);
  pt12->SetTextFont(42);
  TText *text12 = 0;
  if(strcmp(ksp,"pp")==0)text12 = pt12->AddText(Form("%s, |#eta|<%0.0f",ccent[cent],ketacut));
  else text12 = pt12->AddText(Form("PbPb %s, |#eta|<%0.0f",ccent[cent],ketacut));
  text12->SetTextSize(0.08);
  pt12->Draw();

  //return 0;

  c10->cd(maxc);
  gPad->SetLogx();

  hpmean[0]->SetMarkerStyle(20);
  hpmean[0]->SetMarkerColor(4);
  hpmean[0]->SetLineColor(4);
  hpmean[0]->SetMarkerSize(1.2);
  hpmean[0]->Draw("p");  

  l10[++ij] = new TLegend(0.2555564,0.6330544,0.5043282,0.9916947,NULL,"BRNDC");
  l10[ij]->SetHeader("");
  l10[ij]->SetBorderSize(0);
  l10[ij]->SetTextFont(42);
  l10[ij]->SetTextSize(0.05);
  l10[ij]->SetLineColor(1);
  l10[ij]->SetLineStyle(1);
  l10[ij]->SetLineWidth(1);
  l10[ij]->SetFillColor(10);
  l10[ij]->SetFillStyle(1001);
  l10[ij]->SetHeader("");			  
  l10[ij]->AddEntry(hpmean[0] ,calgo[0],"p"); 



  c10->cd(2*maxc);
  gPad->SetLogx();

  hprms[0]->SetMarkerStyle(20);
  hprms[0]->SetMarkerColor(1);
  hprms[0]->SetLineColor(1);
  hprms[0]->SetMarkerSize(1.2);
  hprms[0]->Draw("p");
  l10[ij]->Draw();

  //line->Draw();

  //return 0;



  for(int nj=1;nj<knj;nj++){
    hRMS[nj]->SetTitle("");
    hRMS[nj]->GetXaxis()->SetRangeUser(min,max);
    hRMS[nj]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    hRMS[nj]->GetXaxis()->CenterTitle(true);
    hRMS[nj]->GetXaxis()->SetTitleFont(42);
    hRMS[nj]->GetXaxis()->SetLabelFont(42);
    hRMS[nj]->GetXaxis()->SetMoreLogLabels();
    hRMS[nj]->GetXaxis()->SetNoExponent();
    hRMS[nj]->GetXaxis()->SetTitleSize(0.06);
    hRMS[nj]->GetXaxis()->SetTitleOffset(1.15);
    hRMS[nj]->GetXaxis()->SetLabelSize(0.06);
    hRMS[nj]->GetXaxis()->SetNdivisions(507);
    hRMS[nj]->GetYaxis()->SetTitle("#sigma (RecoJet p_{T} / GenJet p_{T})");    
    hRMS[nj]->GetYaxis()->SetTitleSize(0.07);
    hRMS[nj]->GetYaxis()->SetTitleOffset(1.26);
    hRMS[nj]->GetYaxis()->SetLabelSize(0.08);
    hRMS[nj]->GetYaxis()->SetNdivisions(507);
    hRMS[nj]->GetYaxis()->SetTitleFont(42);
    hRMS[nj]->GetYaxis()->SetLabelFont(42);
    hRMS[nj]->SetMaximum(highy1);
    hRMS[nj]->SetMinimum(lowy1);


    hSigma[nj]->SetTitle("");
    hSigma[nj]->GetXaxis()->SetRangeUser(min,max);
    hSigma[nj]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    hSigma[nj]->GetXaxis()->CenterTitle(true);
    hSigma[nj]->GetXaxis()->SetTitleFont(42);
    hSigma[nj]->GetXaxis()->SetLabelFont(42);
    hSigma[nj]->GetXaxis()->SetMoreLogLabels();
    hSigma[nj]->GetXaxis()->SetNoExponent();
    hSigma[nj]->GetXaxis()->SetTitleSize(0.06);
    hSigma[nj]->GetXaxis()->SetTitleOffset(1.15);
    hSigma[nj]->GetXaxis()->SetLabelSize(0.06);
    hSigma[nj]->GetXaxis()->SetNdivisions(507);
    hSigma[nj]->GetYaxis()->SetTitle("#sigma (RecoJet p_{T} / GenJet p_{T})");    
    hSigma[nj]->GetYaxis()->SetTitleSize(0.07);
    hSigma[nj]->GetYaxis()->SetTitleOffset(1.26);
    hSigma[nj]->GetYaxis()->SetLabelSize(0.08);
    hSigma[nj]->GetYaxis()->SetNdivisions(507);
    hSigma[nj]->GetYaxis()->SetTitleFont(42);
    hSigma[nj]->GetYaxis()->SetLabelFont(42);
    hSigma[nj]->SetMaximum(highy1);
    hSigma[nj]->SetMinimum(lowy1);

    //hArM[nj]->SetMaximum(1.045);
    //hArM[nj]->SetMinimum(0.945);
    hArM[nj]->SetMaximum(highy);
    hArM[nj]->SetMinimum(lowy);
    hArM[nj]->SetTitle("");
    hArM[nj]->GetXaxis()->SetRangeUser(min,max);
    hArM[nj]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    hArM[nj]->GetXaxis()->CenterTitle(true);
    hArM[nj]->GetXaxis()->SetMoreLogLabels();
    hArM[nj]->GetXaxis()->SetNoExponent();
    hArM[nj]->GetXaxis()->SetTitleFont(42);
    hArM[nj]->GetXaxis()->SetLabelFont(42);
    hArM[nj]->GetXaxis()->SetTitleSize(0.06);
    hArM[nj]->GetXaxis()->SetTitleOffset(1.15);
    hArM[nj]->GetXaxis()->SetLabelSize(0.06);
    hArM[nj]->GetXaxis()->SetLabelOffset(0.005);
    hArM[nj]->GetXaxis()->SetNdivisions(507);
    hArM[nj]->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
    hArM[nj]->GetYaxis()->SetTitleSize(0.06);
    hArM[nj]->GetYaxis()->SetTitleOffset(1.50);
    hArM[nj]->GetYaxis()->SetLabelSize(0.06);
    hArM[nj]->GetYaxis()->SetNdivisions(507);
    hArM[nj]->GetYaxis()->SetDecimals(true);
    hArM[nj]->GetYaxis()->SetTitleFont(42);
    hArM[nj]->GetYaxis()->SetLabelFont(42);
    

    //hMean[nj]->SetMaximum(1.045);
    //hMean[nj]->SetMinimum(0.945);
    hMean[nj]->SetMaximum(highy);
    hMean[nj]->SetMinimum(lowy);
    hMean[nj]->SetTitle("");
    hMean[nj]->GetXaxis()->SetRangeUser(min,max);
    hMean[nj]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    hMean[nj]->GetXaxis()->CenterTitle(true);
    hMean[nj]->GetXaxis()->SetMoreLogLabels();
    hMean[nj]->GetXaxis()->SetNoExponent();
    hMean[nj]->GetXaxis()->SetTitleFont(42);
    hMean[nj]->GetXaxis()->SetLabelFont(42);
    hMean[nj]->GetXaxis()->SetTitleSize(0.06);
    hMean[nj]->GetXaxis()->SetTitleOffset(1.15);
    hMean[nj]->GetXaxis()->SetLabelSize(0.06);
    hMean[nj]->GetXaxis()->SetLabelOffset(0.005);
    hMean[nj]->GetXaxis()->SetNdivisions(507);
    hMean[nj]->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
    hMean[nj]->GetYaxis()->SetTitleSize(0.06);
    hMean[nj]->GetYaxis()->SetTitleOffset(1.50);
    hMean[nj]->GetYaxis()->SetLabelSize(0.06);
    hMean[nj]->GetYaxis()->SetNdivisions(507);
    hMean[nj]->GetYaxis()->SetDecimals(true);
    hMean[nj]->GetYaxis()->SetTitleFont(42);
    hMean[nj]->GetYaxis()->SetLabelFont(42);
  }  



  /*
  const int maxa=knj/6;
  const char *cann[maxa] = {"akPF","akCalo","akPuPF","akPuCalo"};
  //float minS[maxa] = {0};
  //float maxS[maxa] = {1.6};

  std::cout<<std::endl;
  std::cout<<"maxa : "<<maxa<<std::endl;

  TCanvas *c5[maxa], *c6[maxa], *c7[maxa];
  ic=-1;
  for(int nj=1;nj<knj;nj++){
    if((nj-1)%6==0){
      ++ic;

      //c5[ic] = new TCanvas(Form("c5_%d",ic),Form("%d eta mean",ic),152,29,1691,1020);
      //makeMultiPanelCanvas(c5[ic],6,maxe,0.0,0.0,0.22,0.22,0.02);

      //c6[ic] = new TCanvas(Form("c6_%d",ic),Form("%d eta sigma",ic),152,29,1691,1020);
      //makeMultiPanelCanvas(c6[ic],6,maxe,0.0,0.0,0.22,0.22,0.02);
      cout<<"ic : "<<ic<<" nj : "<<nj<<endl;

      c7[ic] = new TCanvas(Form("c7_%d",ic),Form("%s Sigma-mean eta ",cann[ic]),157,52,1664,521);
      makeMultiPanelCanvas(c7[ic],6,2,0.0,0.0,0.22,0.22,0.02);
    }	
  }
  //return 0;


  for(int nj=0;nj<knj;nj++){
    for(int ie=0;ie<2;ie++){
	hArM_eta[nj][ie]->SetTitle("");
	hArM_eta[nj][ie]->SetMinimum(lowy);
	hArM_eta[nj][ie]->SetMaximum(highy);
	hArM_eta[nj][ie]->SetMarkerColor(col[ie]);
	hArM_eta[nj][ie]->SetLineColor(col[ie]);
	hArM_eta[nj][ie]->SetMarkerStyle(sty[ie]);
	hArM_eta[nj][ie]->SetMarkerSize(1.1);
	hArM_eta[nj][ie]->GetXaxis()->SetRangeUser(min,max);
	hArM_eta[nj][ie]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	hArM_eta[nj][ie]->GetXaxis()->CenterTitle(true);
	hArM_eta[nj][ie]->GetXaxis()->SetMoreLogLabels();
	hArM_eta[nj][ie]->GetXaxis()->SetNoExponent();
	hArM_eta[nj][ie]->GetXaxis()->SetTitleFont(42);
	hArM_eta[nj][ie]->GetXaxis()->SetLabelFont(42);
	hArM_eta[nj][ie]->GetXaxis()->SetTitleSize(0.08);
	hArM_eta[nj][ie]->GetXaxis()->SetTitleOffset(1.15);
	hArM_eta[nj][ie]->GetXaxis()->SetLabelSize(0.08);
	hArM_eta[nj][ie]->GetXaxis()->SetLabelOffset(0.005);
	hArM_eta[nj][ie]->GetXaxis()->SetNdivisions(507);
	hArM_eta[nj][ie]->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
	hArM_eta[nj][ie]->GetYaxis()->SetTitleSize(0.08);
	hArM_eta[nj][ie]->GetYaxis()->SetTitleOffset(1.25);
	hArM_eta[nj][ie]->GetYaxis()->SetLabelSize(0.08);
	hArM_eta[nj][ie]->GetYaxis()->SetNdivisions(507);
	hArM_eta[nj][ie]->GetYaxis()->SetDecimals(true);
	hArM_eta[nj][ie]->GetYaxis()->SetTitleFont(42);
	hArM_eta[nj][ie]->GetYaxis()->SetLabelFont(42);
	
	hMean_eta[nj][ie]->SetTitle("");	
	hMean_eta[nj][ie]->SetMinimum(lowy);
	hMean_eta[nj][ie]->SetMaximum(highy);
	hMean_eta[nj][ie]->SetMarkerColor(col[ie]);
	hMean_eta[nj][ie]->SetLineColor(col[ie]);
	hMean_eta[nj][ie]->SetMarkerStyle(sty[ie]);
	hMean_eta[nj][ie]->SetMarkerSize(1.1);
	hMean_eta[nj][ie]->GetXaxis()->SetRangeUser(min,max);
	hMean_eta[nj][ie]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	hMean_eta[nj][ie]->GetXaxis()->CenterTitle(true);
	hMean_eta[nj][ie]->GetXaxis()->SetMoreLogLabels();
	hMean_eta[nj][ie]->GetXaxis()->SetNoExponent();
	hMean_eta[nj][ie]->GetXaxis()->SetTitleFont(42);
	hMean_eta[nj][ie]->GetXaxis()->SetLabelFont(42);
	hMean_eta[nj][ie]->GetXaxis()->SetTitleSize(0.08);
	hMean_eta[nj][ie]->GetXaxis()->SetTitleOffset(1.15);
	hMean_eta[nj][ie]->GetXaxis()->SetLabelSize(0.08);
	hMean_eta[nj][ie]->GetXaxis()->SetLabelOffset(0.005);
	hMean_eta[nj][ie]->GetXaxis()->SetNdivisions(507);
	hMean_eta[nj][ie]->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
	hMean_eta[nj][ie]->GetYaxis()->SetTitleSize(0.08);
	hMean_eta[nj][ie]->GetYaxis()->SetTitleOffset(1.25);
	hMean_eta[nj][ie]->GetYaxis()->SetLabelSize(0.08);
	hMean_eta[nj][ie]->GetYaxis()->SetNdivisions(507);
	hMean_eta[nj][ie]->GetYaxis()->SetDecimals(true);
	hMean_eta[nj][ie]->GetYaxis()->SetTitleFont(42);
	hMean_eta[nj][ie]->GetYaxis()->SetLabelFont(42);

	hRMS_eta[nj][ie]->SetTitle("");	
	hRMS_eta[nj][ie]->SetMinimum(lowy1);
	hRMS_eta[nj][ie]->SetMaximum(highy1);
	hRMS_eta[nj][ie]->SetMarkerColor(col[ie]);
	hRMS_eta[nj][ie]->SetLineColor(col[ie]);
	hRMS_eta[nj][ie]->SetMarkerStyle(sty[ie]);
	hRMS_eta[nj][ie]->SetMarkerSize(1.1);
	hRMS_eta[nj][ie]->GetXaxis()->SetRangeUser(min,max);
	hRMS_eta[nj][ie]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	hRMS_eta[nj][ie]->GetXaxis()->CenterTitle(true);
	hRMS_eta[nj][ie]->GetXaxis()->SetMoreLogLabels();
	hRMS_eta[nj][ie]->GetXaxis()->SetNoExponent();
	hRMS_eta[nj][ie]->GetXaxis()->SetTitleFont(42);
	hRMS_eta[nj][ie]->GetXaxis()->SetLabelFont(42);
	hRMS_eta[nj][ie]->GetXaxis()->SetTitleSize(0.08);
	hRMS_eta[nj][ie]->GetXaxis()->SetTitleOffset(1.15);
	hRMS_eta[nj][ie]->GetXaxis()->SetLabelSize(0.08);
	hRMS_eta[nj][ie]->GetXaxis()->SetLabelOffset(0.005);
	hRMS_eta[nj][ie]->GetXaxis()->SetNdivisions(507);
	hRMS_eta[nj][ie]->GetYaxis()->SetTitle("#sigma(<RecoJet p_{T} / GenJet p_{T}>)");
	hRMS_eta[nj][ie]->GetYaxis()->SetTitleSize(0.08);
	hRMS_eta[nj][ie]->GetYaxis()->SetTitleOffset(1.25);
	hRMS_eta[nj][ie]->GetYaxis()->SetLabelSize(0.08);
	hRMS_eta[nj][ie]->GetYaxis()->SetNdivisions(507);
	hRMS_eta[nj][ie]->GetYaxis()->SetDecimals(true);
	hRMS_eta[nj][ie]->GetYaxis()->SetTitleFont(42);
	hRMS_eta[nj][ie]->GetYaxis()->SetLabelFont(42);
	
	hSigma_eta[nj][ie]->SetTitle("");
	hSigma_eta[nj][ie]->SetMinimum(lowy1);
	hSigma_eta[nj][ie]->SetMaximum(highy1);
	hSigma_eta[nj][ie]->SetMarkerColor(col[ie]);
	hSigma_eta[nj][ie]->SetLineColor(col[ie]);
	hSigma_eta[nj][ie]->SetMarkerStyle(sty[ie]);
	hSigma_eta[nj][ie]->SetMarkerSize(1.1);
	hSigma_eta[nj][ie]->GetXaxis()->SetRangeUser(min,max);
	hSigma_eta[nj][ie]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	hSigma_eta[nj][ie]->GetXaxis()->CenterTitle(true);
	hSigma_eta[nj][ie]->GetXaxis()->SetMoreLogLabels();
	hSigma_eta[nj][ie]->GetXaxis()->SetNoExponent();
	hSigma_eta[nj][ie]->GetXaxis()->SetTitleFont(42);
	hSigma_eta[nj][ie]->GetXaxis()->SetLabelFont(42);
	hSigma_eta[nj][ie]->GetXaxis()->SetTitleSize(0.08);
	hSigma_eta[nj][ie]->GetXaxis()->SetTitleOffset(1.15);
	hSigma_eta[nj][ie]->GetXaxis()->SetLabelSize(0.08);
	hSigma_eta[nj][ie]->GetXaxis()->SetLabelOffset(0.005);
	hSigma_eta[nj][ie]->GetXaxis()->SetNdivisions(507);
	hSigma_eta[nj][ie]->GetYaxis()->SetTitle("#sigma(<RecoJet p_{T} / GenJet p_{T}>)");
	hSigma_eta[nj][ie]->GetYaxis()->SetTitleSize(0.08);
	hSigma_eta[nj][ie]->GetYaxis()->SetTitleOffset(1.25);
	hSigma_eta[nj][ie]->GetYaxis()->SetLabelSize(0.08);
	hSigma_eta[nj][ie]->GetYaxis()->SetNdivisions(507);
	hSigma_eta[nj][ie]->GetYaxis()->SetDecimals(true);
	hSigma_eta[nj][ie]->GetYaxis()->SetTitleFont(42);
	hSigma_eta[nj][ie]->GetYaxis()->SetLabelFont(42);
    }
  }

  cout<<endl;
  cout<<"eta dependence for Barel and Forward regions only"<<endl;
  cout<<endl;
  TLegend *l7 = new TLegend(0.6189319,0.5743486,0.918033,0.7803089,NULL,"BRNDC");
  l7->SetHeader("");
  l7->SetBorderSize(0);
  l7->SetTextFont(42);
  l7->SetTextSize(0.07);
  l7->SetLineColor(1);
  l7->SetLineStyle(1);
  l7->SetLineWidth(1);
  l7->SetFillColor(10);
  l7->SetFillStyle(1001);
  l7->SetHeader("");			  


  ipad=0;
  for(int i=0;i<maxa/2;i++){
    std::cout<<"i : "<<i<<std::endl;

    ipad=0;
    for(int nj=1;nj<7;nj++){
      
      if(i==0 && nj==1){
	l7->AddEntry(hSigma_eta[i*6+nj][0] ,meta[0],"p"); 
	l7->AddEntry(hSigma_eta[i*6+nj][1] ,meta[1],"p"); 
      }

      c7[i]->cd(++ipad);
      std::cout<<"\t  \t ***** c7 "<<i<<" nj : "<<nj<<"\t"<<calgo[i*6+nj]<<"\t pad Sigma : "<<(i*6+nj)<<std::endl;
      gPad->SetLogx();

      if(iSigma){
	hSigma_eta[i*6+nj][0]->Draw("p");
	hSigma_eta[i*6+nj][1]->Draw("psame");
      }
      else {
	hRMS_eta[i*6+nj][0]->Draw("p");
	hRMS_eta[i*6+nj][1]->Draw("psame");
      }

      if(ipad==1){
	l7->Draw();
	TPaveText *pt6 = new TPaveText(0.3464834,0.7436938,0.4057113,0.9221927,"brNDC");
	pt6->SetBorderSize(0);
	pt6->SetFillColor(10);
	pt6->SetTextFont(42);
	TText *text6 = 0;
	if(strcmp(ksp,"pp")==0)text6 = pt6->AddText(ccent[cent]);
	else text6 = pt6->AddText(Form("PbPb %s",ccent[cent]));
	text6->SetTextSize(0.1);
	pt6->Draw();
      }
      
      TPaveText *pt71 = new TPaveText(0.6323381,0.7904674,0.9324801,0.9423785,"brNDC");
      pt71->SetBorderSize(0);
      pt71->SetFillColor(10);
      pt71->SetTextFont(42);
      TText *text71 = pt71->AddText(calgo[i*6+nj]);
      text71->SetTextSize(0.08);	
      //TText *text72 = pt71->AddText(meta[0]);
      //text72->SetTextSize(0.08);	
      //TText *text73 = pt71->AddText(meta[1]);
      //text73->SetTextSize(0.08);	
      pt71->Draw();
      

      c7[i]->cd(ipad+6);
      gPad->SetLogx();
      std::cout<<"\t  \t ***** c7 "<<i<<" nj : "<<nj<<"\t"<<calgo[i*6+nj]<<"\t pad Mean : "<<(i*6+nj)<<std::endl;
      if(iSigma){
	hMean_eta[i*6+nj][0]->Draw("p");
	hMean_eta[i*6+nj][1]->Draw("psame");
      }
      else {
	hArM_eta[i*6+nj][0]->Draw("p");
	hArM_eta[i*6+nj][1]->Draw("psame");
      }
      line->Draw();

      //! Pu ones
      c7[i+2]->cd(ipad);
      gPad->SetLogx();
      std::cout<<"\t  \t ***** c7 "<<(i+2)<<" nj : "<<nj<<"\t"<<calgo[(i+2)*6+nj]<<"\t pad Sigma : "<<(i*6+nj)<<std::endl;
      if(iSigma){
	hSigma_eta[(i+2)*6+nj][0]->Draw("p");
	hSigma_eta[(i+2)*6+nj][1]->Draw("psame");
      }	else {
	hRMS_eta[(i+2)*6+nj][0]->Draw("p");
	hRMS_eta[(i+2)*6+nj][1]->Draw("psame");
      }

      if(ipad==1){ 
	l7->Draw();
	TPaveText *pt6 = new TPaveText(0.3464834,0.7436938,0.4057113,0.9221927,"brNDC");
	pt6->SetBorderSize(0);
	pt6->SetFillColor(10);
	pt6->SetTextFont(42);
	TText *text6 = 0;
	if(strcmp(ksp,"pp")==0)text6 = pt6->AddText(ccent[cent]);
	else text6 = pt6->AddText(Form("PbPb %s",ccent[cent]));
	text6->SetTextSize(0.08);
	pt6->Draw();
      }

      TPaveText *pt72 = new TPaveText(0.6323381,0.7904674,0.9324801,0.9423785,"brNDC");
      pt72->SetBorderSize(0);
      pt72->SetFillColor(10);
      pt72->SetTextFont(42);
      TText *text73 = pt72->AddText(calgo[(i+2)*6+nj]);
      text73->SetTextSize(0.08);	
      //TText *text74 = pt72->AddText(meta[ie]);
      //text74->SetTextSize(0.08);	
      pt72->Draw();
      
      
      c7[i+2]->cd(ipad+6);
      gPad->SetLogx();	
      std::cout<<"\t  \t ***** c7 "<<(i+2)<<" nj : "<<nj<<"\t"<<calgo[(i+2)*6+nj]<<"\t pad Mean : "<<(i*6+nj)<<std::endl;
      if(iSigma){
	hMean_eta[(i+2)*6+nj][0]->Draw("p");
	hMean_eta[(i+2)*6+nj][1]->Draw("psame");
      }
      else {
	hArM_eta[(i+2)*6+nj][0]->Draw("p");
	hArM_eta[(i+2)*6+nj][1]->Draw("psame");
      }
      line->Draw();

    }//! nj
    std::cout<<std::endl;
  }//! i canvas loop

  if(iSave){
    if(strcmp(algo,"pu")==0){
      c7[2]->SaveAs(Form("plots/EtaDependent_JERJES_akPuPF_%s.png",ccent[cent]));
      c7[3]->SaveAs(Form("plots/EtaDependent_JERJES_akPuCalo_%s.png",ccent[cent]));
    }
    else{
      c7[0]->SaveAs(Form("plots/EtaDependent_JERJES_akPF_%s.png",ccent[cent]));
      c7[1]->SaveAs(Form("plots/EtaDependent_JERJES_akCalo_%s.png",ccent[cent]));
    } 
  }
  return 0;


  ipad=0;
  for(int i=0;i<maxa;i++){
    std::cout<<"i : "<<i<<std::endl;
    if(i==1)return 0;
    for(int ie=0;ie<maxe;ie++){
      std::cout<<"\t  ie : "<<ceta[ie]<<std::endl;
      for(int nj=1;nj<7;nj++){

	c5[i]->cd(ie*6+nj);
	std::cout<<"\t  \t nj : "<<calgo[i*6+nj]<<"\t"<<ie*6+nj<<std::endl;
	gPad->SetLogx();

	hArMepp[i*6+nj][ie]->SetTitle("");
	hArMepp[i*6+nj][ie]->SetMinimum(0.745);
	hArMepp[i*6+nj][ie]->SetMaximum(1.145);

	hArMepp[i*6+nj][ie]->SetMarkerColor(1);
	hArMepp[i*6+nj][ie]->SetLineColor(1);
	hArMepp[i*6+nj][ie]->SetMarkerStyle(20);
	hArMepp[i*6+nj][ie]->SetMarkerSize(1.1);

	hArMepp[i*6+nj][ie]->GetXaxis()->SetRangeUser(min,max);
	hArMepp[i*6+nj][ie]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	hArMepp[i*6+nj][ie]->GetXaxis()->CenterTitle(true);
	hArMepp[i*6+nj][ie]->GetXaxis()->SetMoreLogLabels();
	hArMepp[i*6+nj][ie]->GetXaxis()->SetNoExponent();
	hArMepp[i*6+nj][ie]->GetXaxis()->SetTitleFont(42);
	hArMepp[i*6+nj][ie]->GetXaxis()->SetLabelFont(42);
	hArMepp[i*6+nj][ie]->GetXaxis()->SetTitleSize(0.08);
	hArMepp[i*6+nj][ie]->GetXaxis()->SetTitleOffset(1.15);
	hArMepp[i*6+nj][ie]->GetXaxis()->SetLabelSize(0.08);
	hArMepp[i*6+nj][ie]->GetXaxis()->SetLabelOffset(0.005);
	hArMepp[i*6+nj][ie]->GetXaxis()->SetNdivisions(507);
	hArMepp[i*6+nj][ie]->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
	hArMepp[i*6+nj][ie]->GetYaxis()->SetTitleSize(0.08);
	hArMepp[i*6+nj][ie]->GetYaxis()->SetTitleOffset(1.25);
	hArMepp[i*6+nj][ie]->GetYaxis()->SetLabelSize(0.08);
	hArMepp[i*6+nj][ie]->GetYaxis()->SetNdivisions(507);
	hArMepp[i*6+nj][ie]->GetYaxis()->SetDecimals(true);
	hArMepp[i*6+nj][ie]->GetYaxis()->SetTitleFont(42);
	hArMepp[i*6+nj][ie]->GetYaxis()->SetLabelFont(42);


	hMepp[i*6+nj][ie]->SetMinimum(0.745);
	hMepp[i*6+nj][ie]->SetMaximum(1.145);
	hMepp[i*6+nj][ie]->SetMarkerColor(1);
	hMepp[i*6+nj][ie]->SetLineColor(1);
	hMepp[i*6+nj][ie]->SetMarkerStyle(20);
	hMepp[i*6+nj][ie]->SetMarkerSize(1.1);
	hMepp[i*6+nj][ie]->GetXaxis()->SetRangeUser(min,max);
	hMepp[i*6+nj][ie]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	hMepp[i*6+nj][ie]->GetXaxis()->CenterTitle(true);
	hMepp[i*6+nj][ie]->GetXaxis()->SetMoreLogLabels();
	hMepp[i*6+nj][ie]->GetXaxis()->SetNoExponent();
	hMepp[i*6+nj][ie]->GetXaxis()->SetTitleFont(42);
	hMepp[i*6+nj][ie]->GetXaxis()->SetLabelFont(42);
	hMepp[i*6+nj][ie]->GetXaxis()->SetTitleSize(0.08);
	hMepp[i*6+nj][ie]->GetXaxis()->SetTitleOffset(1.15);
	hMepp[i*6+nj][ie]->GetXaxis()->SetLabelSize(0.08);
	hMepp[i*6+nj][ie]->GetXaxis()->SetLabelOffset(0.005);
	hMepp[i*6+nj][ie]->GetXaxis()->SetNdivisions(507);
	hMepp[i*6+nj][ie]->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
	hMepp[i*6+nj][ie]->GetYaxis()->SetTitleSize(0.08);
	hMepp[i*6+nj][ie]->GetYaxis()->SetTitleOffset(1.25);
	hMepp[i*6+nj][ie]->GetYaxis()->SetLabelSize(0.08);
	hMepp[i*6+nj][ie]->GetYaxis()->SetNdivisions(507);
	hMepp[i*6+nj][ie]->GetYaxis()->SetDecimals(true);
	hMepp[i*6+nj][ie]->GetYaxis()->SetTitleFont(42);
	hMepp[i*6+nj][ie]->GetYaxis()->SetLabelFont(42);
	hMepp[i*6+nj][ie]->SetTitle("");

	if(iSigma)hMepp[i*6+nj][ie]->Draw();
	else hArMepp[i*6+nj][ie]->Draw();

	TPaveText *pt1 = new TPaveText(0.6323381,0.7904674,0.9324801,0.9423785,"brNDC");
	pt1->SetBorderSize(0);
	pt1->SetFillColor(10);
	pt1->SetTextFont(42);
	TText *text1 = pt1->AddText(calgo[i*6+nj]);
	text1->SetTextSize(0.08);	
	TText *text2 = pt1->AddText(ceta[ie]);
	text2->SetTextSize(0.08);	
	pt1->Draw();
	
	if(ie*6+nj>=25){
	  text1->SetTextSize(0.07);
	  text2->SetTextSize(0.07);
	  hArMepp[i*6+nj][ie]->GetYaxis()->SetTitleSize(0.065);
	  hArMepp[i*6+nj][ie]->GetYaxis()->SetLabelSize(0.065);
	  hArMepp[i*6+nj][ie]->GetYaxis()->SetTitleOffset(1.45);

	  hMepp[i*6+nj][ie]->GetYaxis()->SetTitleSize(0.065);
	  hMepp[i*6+nj][ie]->GetYaxis()->SetLabelSize(0.065);
	  hMepp[i*6+nj][ie]->GetYaxis()->SetTitleOffset(1.45);
	}
	line->Draw();

	c6[i]->cd(ie*6+nj);
	gPad->SetLogx();
	hRMSepp[i*6+nj][ie]->SetMinimum(0.005);
	hRMSepp[i*6+nj][ie]->SetMaximum(0.845);
	
	hRMSepp[i*6+nj][ie]->SetMarkerColor(1);
	hRMSepp[i*6+nj][ie]->SetLineColor(1);
	hRMSepp[i*6+nj][ie]->SetMarkerStyle(20);
	hRMSepp[i*6+nj][ie]->SetMarkerSize(1.1);
	
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetRangeUser(min,max);
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	hRMSepp[i*6+nj][ie]->GetXaxis()->CenterTitle(true);
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetMoreLogLabels();
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetNoExponent();
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetTitleFont(42);
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetLabelFont(42);
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetTitleSize(0.08);
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetTitleOffset(1.15);
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetLabelSize(0.08);
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetLabelOffset(0.005);
	hRMSepp[i*6+nj][ie]->GetXaxis()->SetNdivisions(507);
	hRMSepp[i*6+nj][ie]->GetYaxis()->SetTitle("#sigma(<RecoJet p_{T} / GenJet p_{T}>)");
	hRMSepp[i*6+nj][ie]->GetYaxis()->SetTitleSize(0.08);
	hRMSepp[i*6+nj][ie]->GetYaxis()->SetTitleOffset(1.25);
	hRMSepp[i*6+nj][ie]->GetYaxis()->SetLabelSize(0.08);
	hRMSepp[i*6+nj][ie]->GetYaxis()->SetNdivisions(507);
	hRMSepp[i*6+nj][ie]->GetYaxis()->SetDecimals(true);
	hRMSepp[i*6+nj][ie]->GetYaxis()->SetTitleFont(42);
	hRMSepp[i*6+nj][ie]->GetYaxis()->SetLabelFont(42);

	hRMSepp[i*6+nj][ie]->SetTitle("");


	//
	hSepp[i*6+nj][ie]->SetMinimum(0.005);
	hSepp[i*6+nj][ie]->SetMaximum(0.845);
	hSepp[i*6+nj][ie]->SetMarkerColor(1);
	hSepp[i*6+nj][ie]->SetLineColor(1);
	hSepp[i*6+nj][ie]->SetMarkerStyle(20);
	hSepp[i*6+nj][ie]->SetMarkerSize(1.1);
	hSepp[i*6+nj][ie]->GetXaxis()->SetRangeUser(min,max);
	hSepp[i*6+nj][ie]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	hSepp[i*6+nj][ie]->GetXaxis()->CenterTitle(true);
	hSepp[i*6+nj][ie]->GetXaxis()->SetMoreLogLabels();
	hSepp[i*6+nj][ie]->GetXaxis()->SetNoExponent();
	hSepp[i*6+nj][ie]->GetXaxis()->SetTitleFont(42);
	hSepp[i*6+nj][ie]->GetXaxis()->SetLabelFont(42);
	hSepp[i*6+nj][ie]->GetXaxis()->SetTitleSize(0.08);
	hSepp[i*6+nj][ie]->GetXaxis()->SetTitleOffset(1.15);
	hSepp[i*6+nj][ie]->GetXaxis()->SetLabelSize(0.08);
	hSepp[i*6+nj][ie]->GetXaxis()->SetLabelOffset(0.005);
	hSepp[i*6+nj][ie]->GetXaxis()->SetNdivisions(507);
	hSepp[i*6+nj][ie]->GetYaxis()->SetTitle("#sigma(<RecoJet p_{T} / GenJet p_{T}>)");
	hSepp[i*6+nj][ie]->GetYaxis()->SetTitleSize(0.08);
	hSepp[i*6+nj][ie]->GetYaxis()->SetTitleOffset(1.25);
	hSepp[i*6+nj][ie]->GetYaxis()->SetLabelSize(0.08);
	hSepp[i*6+nj][ie]->GetYaxis()->SetNdivisions(507);
	hSepp[i*6+nj][ie]->GetYaxis()->SetDecimals(true);
	hSepp[i*6+nj][ie]->GetYaxis()->SetTitleFont(42);
	hSepp[i*6+nj][ie]->GetYaxis()->SetLabelFont(42);
	hSepp[i*6+nj][ie]->SetTitle("");

	if(iSigma)hSepp[i*6+nj][ie]->Draw();
	else hRMSepp[i*6+nj][ie]->Draw();

	TPaveText *pt2 = new TPaveText(0.6323381,0.7904674,0.9324801,0.9423785,"brNDC");
	pt2->SetBorderSize(0);
	pt2->SetFillColor(10);
	pt2->SetTextFont(42);
	TText *text6 = pt2->AddText(calgo[i*6+nj]);
	text6->SetTextSize(0.08);	
	TText *text7 = pt2->AddText(ceta[ie]);
	text7->SetTextSize(0.08);	
	pt2->Draw();
	
	if(ie*6+nj>=25){
	  text6->SetTextSize(0.065);
	  text7->SetTextSize(0.065);
	  hRMSepp[i*6+nj][ie]->GetYaxis()->SetTitleSize(0.065);
	  hRMSepp[i*6+nj][ie]->GetYaxis()->SetLabelSize(0.065);
	  hRMSepp[i*6+nj][ie]->GetYaxis()->SetTitleOffset(1.45);

	  hSepp[i*6+nj][ie]->GetYaxis()->SetTitleSize(0.08);
	  hSepp[i*6+nj][ie]->GetYaxis()->SetLabelSize(0.08);
	  hSepp[i*6+nj][ie]->GetYaxis()->SetTitleOffset(1.45);
	}
      }//! nj
    }//! ie
  }//! i canvas loop
  //return 0;
  
  if(iSave){
    for(int i=0;i<maxa;i++){
      c5[i]->SaveAs(Form("plots/etaDependence_JES_%s_%s_%d.png",algo,cdR,i));
      c6[i]->SaveAs(Form("plots/etaDependence_JER_%s_%s_%d.png",algo,cdR,i));
    }  
  }


  ipad=0;
  TCanvas *c1 = new TCanvas("c","Gen:Reco jet pT",136,175,1752,577);
  c1->Divide(6,2,0,0.1);
  for(int nj=1;nj<knj;nj++){

    c1->cd(++ipad);

    if((nj-1)%6==0)gPad->SetLeftMargin(0.15);
    if(ipad%6==0)gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.15);

    hgenjrecoj[nj]->SetStats(0);
    //hgenjrecoj[nj]->GetXaxis()->SetRangeUser(0.,214);
    //hgenjrecoj[nj]->GetYaxis()->SetRangeUser(0.,214);

    hgenjrecoj[nj]->GetXaxis()->SetRangeUser(0.,128);
    hgenjrecoj[nj]->GetYaxis()->SetRangeUser(0.,128);

    hgenjrecoj[nj]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    hgenjrecoj[nj]->GetXaxis()->CenterTitle(true);
    hgenjrecoj[nj]->GetXaxis()->SetTitleFont(42);
    hgenjrecoj[nj]->GetXaxis()->SetLabelFont(42);
    hgenjrecoj[nj]->GetXaxis()->SetMoreLogLabels();
    hgenjrecoj[nj]->GetXaxis()->SetNoExponent();
    hgenjrecoj[nj]->GetXaxis()->SetTitleSize(0.07);
    hgenjrecoj[nj]->GetXaxis()->SetTitleOffset(1.0);
    hgenjrecoj[nj]->GetXaxis()->SetLabelSize(0.07);
    hgenjrecoj[nj]->GetXaxis()->SetNdivisions(507);
    hgenjrecoj[nj]->GetYaxis()->SetTitle("RecoJet p_{T} (GeV/c)");    
    hgenjrecoj[nj]->GetYaxis()->SetTitleSize(0.07);
    hgenjrecoj[nj]->GetYaxis()->SetTitleOffset(1.05);
    hgenjrecoj[nj]->GetYaxis()->SetLabelSize(0.07);
    hgenjrecoj[nj]->GetYaxis()->SetNdivisions(507);
    hgenjrecoj[nj]->GetYaxis()->SetTitleFont(42);
    hgenjrecoj[nj]->GetYaxis()->SetLabelFont(42);
    hgenjrecoj[nj]->SetTitle("");
    
    hgenjrecoj[nj]->Draw("col");
   
    TPaveText *pt   = new TPaveText(0.15,0.88,0.45,0.98,"brNDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(10);
    pt->SetTextFont(42);
    TText *text = pt->AddText(calgo[nj]);
    text->SetTextSize(0.09);
    pt->Draw();
  }
  if(iSave)c1->SaveAs(Form("plots/GenReco_%s_JEC.png",cdR));
  */
  /*
  if(iSave){
    for(int i=0;i<maxa;i++){
      c7[i]->SaveAs(Form("plots/etaDependence_TwoBins_JES_%s_%s_%d.png",algo,cdR,i));
      c8[i]->SaveAs(Form("plots/etaDependence_TwoBins_JER_%s_%d.png",algo,cdR,i));
    }  
  }
  */

  //if(iSave)c10->SaveAs(Form("plots/JEC_%s_%s_%s_2012Forest.png",algo,ccent[cent],cdR));
  return 0;
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
