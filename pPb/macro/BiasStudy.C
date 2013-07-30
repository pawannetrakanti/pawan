// From the prompt reco
#include "../HiForest/V2/hiForest.h"
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

const double pi=acos(-1.);
//const double pi2=2*pi -1;
const float ketacut=3.;
const double kptrecocut=30.;
const double kalpha=0.8;

const int knj = 4;
const char *calgo[knj] = {"ak3PF"    ,"akPu3PF", "ak3Calo", "akPu3Calo"};
//const char *calgo[knj] = {"ak3Calo"    ,"akPu3Calo"};
const double kDelR[knj]  = {0.3,0.3,0.3,0.3};


void LoadLib();
void ShutoffBranches(HiForest */*hi*/,const char */*ksp*/);
void FindLeadSubLeadJets(Jets */*mJets*/, int */*ljet*/);
double delphi(double /*phi1*/, double /*phi2*/);


int GetNtracksBin(int /*hiNtracks*/);
int GetNPixBin (int /*hiNPix*/);
int GetNPixelTracksBin(int /*hiNPixelTracks*/);
int GetHFBin   (double /*hiHF*/);
int GetAnaBin  (double /*hiHFplusEta4*/);
int GetNewAnaBin(double /*hiHFparam*/);
int GetHFPE4Bin(double /*hiHFplusEta4*/);
int GetHFME4Bin(double /*hiHFplusEta4*/);
int GetHFPBin  (double /*hiHFplus*/);
int GetHFMBin  (double /*hiHFminus*/);
int GetZDCBin  (double /*hiZDCplus*/);
int GetZDCPBin (double /*hiZDCplus*/);
int GetZDCMBin (double /*hiZDCminus*/);


const int nmult =66; //! # of bins for different centralities in groups of 5
const int nbins =5;  //! # of bins in centrality
const int ngrp  =14; //! # of groups for centrality
const int lbins =6;  //! lead-jet pT bins
const char *ctrig[lbins] = {"p_{T}^{Jet}>80 GeV/c",
                            "p_{T}^{Jet}>100 GeV/c",
                            "p_{T}^{Jet}>120 GeV/c",
                            "p_{T}^{Jet}>150 GeV/c",
                            "p_{T}^{Jet}>200 GeV/c",
                            "p_{T}^{Jet}>300 GeV/c",
};

const char *cassoc[lbins] = {"p_{T}^{Jet}>30 GeV/c",
			     "p_{T}^{Jet}>40 GeV/c",
			     "p_{T}^{Jet}>50 GeV/c",
			     "p_{T}^{Jet}>60 GeV/c",
			     "p_{T}^{Jet}>70 GeV/c",
			     "p_{T}^{Jet}>80 GeV/c",
};
const char *clead[lbins] = {"80GeV",
                            "100GeV",
                            "120GeV",
                            "150GeV",
                            "200GeV",
                            "300GeV"
};
const char *cslead[lbins] = {"30GeV",
                             "40GeV",
                             "50GeV",
                             "60GeV",
                             "70GeV",
                             "80GeV"
};


int bpixhit  [nbins] = {661 , 527 , 434 , 357 , 0}; //! "PixelHits"
int bpixtr   [nbins] = {59  , 45  , 36  , 29  , 0}; //! "PixelTracks"
int btrks    [nbins] = {82  , 64  , 52  , 42  , 0}; //! "Tracks"

double bhf   [nbins] = {72.2, 57.3, 46.8, 38.4, 0}; //! "HFtowers" (hiHF)
double bhfp  [nbins] = {46.6, 34.9, 27.1, 21.1, 0}; //! "HFtowersPlus" (hiHFplus)
double bhfm  [nbins] = {28.7, 23.3, 19.4, 16.2, 0}; //! "HFtowersMinus" (hiHFminus)

double bana      [nbins] = {30.0, 20.0, 15.0, 10.0, 0}; //! "HFtowersPlusTrunc" (hiHFplusEta4) used in highpT analysis 
double bana_new  [nbins] = {40.0, 30.0, 25.0, 20.0, 0}; //! "HFSum" (hiHFSum = 1.17 *hiHFplus + 6) used in highpT analysis new parametrization

double bhfpe4[nbins] = {19.3, 14.4, 11.1, 8.56, 0}; //! "HFtowersPlusTrunc" (hiHFplusEta4)
double bhfme4[nbins] = {11.0, 8.87, 7.36, 6.12, 0}; //! "HFtowersMinusTrunc" (hiHFminusEta4)

double bzdc  [nbins] = {9085.01, 8230.39, 7514.89, 6808.38, 0}; //! "ZDC"
double bzdcp [nbins] = {84.7248, 37.2032, 21.323, 11.1515 , 0}; //! "ZDCplus"
double bzdcm [nbins] = {9052.98, 8197.01, 7478.34, 6764.27, 0}; //! "ZDCminus"                                                                                                                                      

/*                                                                                                                                                                                                                  
double bhf   [nbins] = {72.2, 57.3, 46.8, 38.4, 30.9, 0}; //! "HFtowers" (hiHF)                                                                                                                                     
double bhfp  [nbins] = {46.6, 34.9, 27.1, 21.1, 16.1, 0}; //! "HFtowersPlus" (hiHFplus)                                                                                                                             
double bhfm  [nbins] = {28.7, 23.3, 19.4, 16.2, 13.3, 0}; //! "HFtowersMinus" (hiHFminus)                                                                                                                           
double bhfpe4[nbins] = {19.3, 14.4, 11.1, 8.56, 6.49, 0}; //! "HFtowersPlusTrunc" (hiHFplusEta4)                                                                                                                    
double bhfme4[nbins] = {11.0, 8.87, 7.36, 6.12, 5.02, 0}; //! "HFtowersMinusTrunc" (hiHFminusEta4)                                                                                                                  
int bpixhit  [nbins] = {661 , 527 , 434 , 357 , 290 , 0}; //! "PixelHits"                                                                                                                                           
int bpixtr   [nbins] = {59  , 45  , 36  , 29  , 22  , 0}; //! "PixelTracks"                                                                                                                                         
int btrks    [nbins] = {82  , 64  , 52  , 42  , 33  , 0}; //! "Tracks"                                                                                                                                              
//"ZDC" {9085.01, 8230.39, 7514.89, 6808.38, 6031.04, 0}                                                                                                                                                            
//"ZDCplus" {84.7248, 37.2032, 21.323, 11.1515, 3.25617, 0}                                                                                                                                                         
//"ZDCminus" {9052.98, 8197.01, 7478.34, 6764.27, 5973.12, 0}                                                                                                                                                       
*/


TStopwatch timer;
int BiasStudy(int leadb=2,int sleadb=0,const char *ksp="ppb")
{

  timer.Start();

  LoadLib();

  TString inname="";
  if(strcmp(ksp,"ppb")==0 || strcmp(ksp,"pbp")==0 ){
    //! Prompt reco for Jet80
    //!  cmsPfn /store/caf/user/velicanu/PA2013_merged_HiForest/pPb_hiForest2_1_15_test.root //! to get the path from CAF
    
    inname="root://eoscms//eos/cms/store/caf/user/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";
    
    //! Minbias
    //inname="root://eoscms//eos/cms/store/group/phys_heavyions/kjung/MinBiasUPCForest_v71/hiForest_MinBiasUPC_v71_Merged1.root";
    //inname="root://eoscms//eos/cms/store/group/phys_heavyions/kjung/MinBiasUPCForest_v71/hiForest_MinBiasUPC_v71_Merged2.root";
  }else  if(strcmp(ksp,"pbpb")==0){
    inname="root://eoscms//eos/cms/store/caf/user/yjlee/PbPHiForest2_PbPbPAHighPtJet80_cent50-100_pprereco.root";
  }else if(strcmp(ksp,"pp")==0){
    //inname="root://eoscms//eos/cms/store/caf/user/yjlee/pp/hiForest2_pp_ppreco_415_90percent.root";
    inname="root://eoscms//eos/cms/store/caf/user/yjlee/pp/PP2013_HiForest_PromptReco_HighPtJet_forestv77.root";
  }



  //! Load Lib
  //gSystem->Load("../HiForest/V2/hiForest_h.so");


  //! Define the input file and HiForest
  //! CMSSW_5_3_3
  HiForest *c = 0;
  if(strcmp(ksp,"ppb")==0 || strcmp(ksp,"pbp")==0)c = new HiForest(inname,Form("Forest%s",ksp),cPPb);
  else if(strcmp(ksp,"pp")==0)c = new HiForest(inname,Form("Forest%s",ksp),cPP);
  else c = new HiForest(inname,Form("Forest%s",ksp),cPbPb);

  cout<<"Loaded the hiforest tree : "<<c->GetName()<<endl;
  ShutoffBranches(c,ksp);



  //! Output file
  //! HIHighPt
  TFile *fout = new TFile(Form("output/promptreco/BiasStudy/%s/BiasStudy_%s_Lead%s_SubLead%s_v77.root",ksp,ksp,clead[leadb],cslead[sleadb]),"RECREATE");  

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Running for %s ",ksp)<<std::endl;
  std::cout<<Form("pT  cut for %0.3f ",kptrecocut)<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"My hiForest Tree : " <<c->GetName()<<"\t Entries "<<c->GetEntries()<<std::endl;
  std::cout<<"Output file  "<<fout->GetName()<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;

  //! shut off jet trees
  c->hasIcPu5JetTree=0;
  c->hasAk1CaloJetTree=0;
  c->hasAk2CaloJetTree=0;
  c->hasAk4CaloJetTree=0;
  //c->hasAk3CaloJetTree=0;
  c->hasAk5CaloJetTree=0;
  c->hasAk6CaloJetTree=0;

  c->hasAkPu1CaloJetTree=0;
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  //c->hasAkPu3CaloJetTree=0;
  c->hasAkPu5CaloJetTree=0;
  c->hasAkPu6CaloJetTree=0;

  c->hasAk1PFJetTree=0;
  c->hasAk2PFJetTree=0;
  //c->hasAk3PFJetTree=0;
  c->hasAk4PFJetTree=0;
  c->hasAk5PFJetTree=0;
  c->hasAk6PFJetTree=0;

  c->hasAkPu1PFJetTree=0;
  c->hasAkPu2PFJetTree=0;
  //c->hasAkPu3PFJetTree=0;
  c->hasAkPu4PFJetTree=0;
  c->hasAkPu5PFJetTree=0;
  c->hasAkPu6PFJetTree=0;

  c->hasTrackTree=0;

  //! For jets
  Jets *mJets=0;

  //! Define histograms here
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  TProfile::SetDefaultSumw2();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1D *hEvt = new TH1D("hEvt","# of events ",4,0,4);
  TH1D *hVz  = new TH1D("hVz","vertex z distribution",80,-20,20);

  TH2F *hHFplusEta4NTracks = new TH2F("hHFplusEta4NTracks","HFplusEta4 vs NTracks",800,0,200,1000,0,1000);
  TH2F *hHFplusEta4Npix    = new TH2F("hHFplusEta4Npix","HFplusEta4 vs Npix",800,0,200,5000,0,5000);
  TH2F *hHFplusEta4NpixelTracks    = new TH2F("hHFplusEta4NpixelTracks","HFplusEta4 vs Npix",800,0,200,1000,0,1000);

  TH2F *hHFplusEta4HF = new TH2F("hHFplusEta4HF","HFplusEta4 vs HF",800,0,200,1600,0,400);
  TH2F *hHFplusEta4HFminusEta4 = new TH2F("hHFplusEta4HFminusEta4","HFplusEta4 vs HFminusEta4",800,0,200,800,0,200);
  TH2F *hHFplusEta4HFminus = new TH2F("hHFplusEta4HFminus","HFplusEta4 vs HFminus",800,0,200,800,0,200);
  TH2F *hHFplusEta4HFplus = new TH2F("hHFplusEta4HFplus","HFplusEta4 vs HFplus",800,0,200,800,0,200);

  TH2F *hHFplusEta4ZDC = new TH2F("hHFplusEta4ZDC","HFplusEta4 vs ZDC",800,0,200,5000,0,45000);
  TH2F *hHFplusEta4ZDCplus = new TH2F("hHFplusEta4ZDCplus","HFplusEta4 vs ZDCplus",800,0,200,5000,0,45000);
  TH2F *hHFplusEta4ZDCminus = new TH2F("hHFplusEta4ZDCminus","HFplusEta4 vs ZDCminus",800,0,200,5000,0,45000);

  //! HFSum                                                                                                                                                                                                         
  TH2F *hHFSumNTracks = new TH2F("hHFSumNTracks","HFplusEta4 vs NTracks",800,0,200,1000,0,1000);
  TH2F *hHFSumNpix    = new TH2F("hHFSumNpix","HFplusEta4 vs Npix",800,0,200,5000,0,5000);
  TH2F *hHFSumNpixelTracks    = new TH2F("hHFSumNpixelTracks","HFplusEta4 vs Npix",800,0,200,1000,0,1000);

  TH2F *hHFSumHF = new TH2F("hHFSumHF","HFSum vs HF",800,0,200,1600,0,400);
  TH2F *hHFSumHFplusEta4 = new TH2F("hHFSumHFplusEta4","HFSum vs HFplusEta4",800,0,200,800,0,200);
  TH2F *hHFSumHFminusEta4 = new TH2F("hHFSumHFminusEta4","HFSum vs HFminusEta4",800,0,200,800,0,200);
  TH2F *hHFSumHFminus = new TH2F("hHFSumHFminus","HFSum vs HFminus",800,0,200,800,0,200);
  TH2F *hHFSumHFplus = new TH2F("hHFSumHFplus","HFSum vs HFplus",800,0,200,800,0,200);

  TH2F *hHFSumZDC = new TH2F("hHFSumZDC","HFSum vs ZDC",800,0,200,5000,0,45000);
  TH2F *hHFSumZDCplus = new TH2F("hHFSumZDCplus","HFSum vs ZDCplus",800,0,200,5000,0,45000);
  TH2F *hHFSumZDCminus = new TH2F("hHFSumZDCminus","HFSum vs ZDCminus",800,0,200,5000,0,45000);
  
  TH3D *hHFSum = new TH3D("hHFSum","Lead pT bin, hf bin and HFsum",nbins,-0.5,nbins-0.5,30,80,680,200,0,120);
  TH1D *hHFSum_LeadJet[lbins];
  for(int i=0;i<lbins;i++){
    hHFSum_LeadJet[i] = new TH1D(Form("hHFSum_LeadJet%d",i),Form("HFSum %s",ctrig[i]),800,0,200);
  }

  //! Reco jet                                                                                                                                                                                                      
  TH3D *hjetpt [knj][3];
  TH3D *hjeteta[knj][3];
  TH3D *hjetphi[knj][3];

  TH3D *hdjeta [knj][2];
  TH3D *hpt2pt1[knj][2];
  TH3D *hdeta  [knj][2];
  TH3D *hdphi  [knj][2];

  TH3D *hypt2pt1   [knj][2][lbins*lbins];
  TH3D *hdedp      [knj][2][lbins*lbins];
  TH3D *hdijept2pt1[knj][2][lbins*lbins];

  TH3D *hJetEnergyBalance[knj][nbins+1];
  TH3D *hJetEnergyRes[knj][nbins+1];

  TH2D *hetaprobe[knj];
  TH2D *hetabarel[knj];
  TH3D *hjes[knj];
  TH3D *hptbptp[knj];


  for(int nj=0;nj<knj;nj++){
    hetaprobe[nj] = new TH2D(Form("hetaprobe%d",nj),Form("hetaprobe%d",nj),24,-1.*ketacut,1.*ketacut,40,80,480);
    hetabarel[nj] = new TH2D(Form("hetabarel%d",nj),Form("hetabarel%d",nj),24,-1.*ketacut,1.*ketacut,40,80,480);
    hjes     [nj] = new TH3D(Form("hjes%d",nj),Form("hjes%d",nj),24,-1.*ketacut,1.*ketacut,40,80,480,50,0,5);
    hptbptp  [nj] = new TH3D(Form("hptbptp%d",nj),Form("hptbptp%d",nj),24,-1.*ketacut,1.*ketacut,40,80,480,40,80,480); 

    for(int ic=0;ic<nbins+1;ic++){
      hJetEnergyRes    [nj][ic] = new TH3D(Form("hJetEnergyRes%d_%d",nj,ic),Form("hJetEnergyRes%d_%d",nj,ic),4,-0.5,3.5,25,50,550,50,-1.00,1.00);
      hJetEnergyBalance[nj][ic] = new TH3D(Form("hJetEnergyBalance%d_%d",nj,ic),Form("hJetEnergyBalance%d_%d",nj,ic),12,-1.*ketacut,ketacut,25,50,550,50,-1.00,1.00);
    }
    for(int j=0;j<3;j++){

      //std::cout<<Form("Histograms for algo %s",calgo[nj])<<std::endl;
      hjetpt [nj][j] = new TH3D(Form("hjetpt%d_%d" ,nj,j) ,Form("Jet pt sepctra %s %d",calgo[nj],j),nmult,-0.5,nmult-0.5,lbins*lbins,-0.5,lbins*lbins-0.5,500,0,1000);
      hjeteta[nj][j] = new TH3D(Form("hjeteta%d_%d",nj,j),Form("Jet eta %s %d",calgo[nj],j),nmult,-0.5,nmult-0.5,lbins*lbins,-0.5,lbins*lbins-0.5,12,-1.*ketacut,ketacut);
      hjetphi[nj][j] = new TH3D(Form("hjetphi%d_%d",nj,j),Form("Jet phi %s %d",calgo[nj],j),nmult,-0.5,nmult-0.5,lbins*lbins,-0.5,lbins*lbins-0.5,30,-pi,pi);

      if(j<2){
	for(int il=0;il<lbins;il++){
	  for(int sil=0;sil<lbins;sil++){
	    int ibin  = il*lbins + sil;
	    hypt2pt1   [nj][j][ibin] = new TH3D(Form("hypt2pt1%d_%d_%d",nj,j,ibin),Form("y^{*} pT2/pT1 %s %d %s %s",calgo[nj],j,ctrig[il],cassoc[sil]),nmult,-0.5,nmult-0.5,24,-2.*ketacut,2*ketacut,10,0,1);
	    hdijept2pt1[nj][j][ibin] = new TH3D(Form("hdijept2pt1%d_%d_%d",nj,j,ibin),Form("(eta1+eta2)/2 vs pT2/pT1 %s %d %s %s",calgo[nj],j,ctrig[il],cassoc[sil]),nmult,-0.5,nmult-0.5,12,-1.*ketacut,ketacut,10,0,1);
	    hdedp      [nj][j][ibin] = new TH3D(Form("hdedp%d_%d_%d",nj,j,ibin),Form("deta-dphi %s %d %s %s",calgo[nj],j,ctrig[il],cassoc[sil]),nmult,-0.5,nmult-0.5,24,-2.*ketacut,2*ketacut,30,0,pi);
	  }
	}
	hpt2pt1 [nj][j] = new TH3D(Form("hpt2pt1%d_%d",nj,j),Form("Jet pT2 / pT1 lead pT bins %s %d ",calgo[nj],j),nmult,-0.5,nmult-0.5,lbins*lbins,-0.5,lbins*lbins-0.5,10,0,1);
        hdphi   [nj][j] = new TH3D(Form("hdphi%d_%d"  ,nj,j),Form("Jet dphi lead pT bins %s %d",calgo[nj],j),nmult,-0.5,nmult-0.5,lbins*lbins,-0.5,lbins*lbins-0.5,30,0,pi);
        hdjeta  [nj][j] = new TH3D(Form("hdjeta%d_%d" ,nj,j),Form("di-jet lead pT bins eta %s %d",calgo[nj],j),nmult,-0.5,nmult-0.5,lbins*lbins,-0.5,lbins*lbins-0.5,12,-1.*ketacut,ketacut);
        hdeta   [nj][j] = new TH3D(Form("hdeta%d_%d"  ,nj,j),Form("Jet deta lead pT bins %s %d",calgo[nj],j),nmult,-0.5,nmult-0.5,lbins*lbins,-0.5,lbins*lbins-0.5,24,-2.*ketacut,2*ketacut);

      }
      //std::cout<<Form("Initialized the histograms %s",calgo[nj]) <<std::endl;
    }
  }
  std::cout<<"\t"<<std::endl;
  /////////////////////////////////////////////////////////////////////////////////////////


  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;
  hEvt->Fill(2,nentries);
  double wxs=1;

  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
    //for (Long64_t ievt=0; ievt<50000;ievt++) {//! event loop
    //! load the hiForest event
    c->GetEntry(ievt);

    
    //! only selecting pPb runs
    if(c->evt.run>211256  && strcmp(ksp,"ppb")==0)continue;
    if(c->evt.run<=211256 && strcmp(ksp,"pbp")==0)continue;

    //! events with Single vertex
    bool evSel = false;

    if(strcmp(ksp,"ppb")==0 || strcmp(ksp,"pbp")==0){
      evSel        = c->skim.pHBHENoiseFilter  && c->skim.pPAcollisionEventSelectionPA && c->skim.pVertexFilterCutGplus 
	&& (c->hlt.HLT_PAJet80_NoJetID_v1 || c->hlt.HLT_PAJet100_NoJetID_v1)
	/*&& c->hlt.HLT_PAJet100_NoJetID_v1*/
	/*&& c->hlt.HLT_PAZeroBiasPixel_SingleTrack_v1*/
	;
    }else if(strcmp(ksp,"pbpb")==0){
      evSel = c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_HIJet80_v1;
    }else if(strcmp(ksp,"pp")==0){
      //evSel = c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_Jet60_v1;
      evSel   = c->skim.pHBHENoiseFilter && (c->hlt.HLT_PAJet80_NoJetID_v1 || c->hlt.HLT_PAJet100_NoJetID_v1) && c->skim.phfPosFilter1 && c->skim.phfNegFilter1;
    }
    //cout<<"Getting the trigger :" <<endl;

    //! select good events
    bool goodEvent = false;
    goodEvent = evSel && fabs(c->evt.vz)<15.; 
    if(!goodEvent)continue;

    //! Event variables
    double vz            = c->evt.vz;    
    int RunNo            = c->evt.run;

    int hiNtracks        = c->evt.hiNtracks;
    int hiNpix           = c->evt.hiNpix;
    int hiNpixelTracks   = c->evt.hiNpixelTracks;

    double hiHF           = c->evt.hiHF;
    double hiHFplusEta4   = c->evt.hiHFplusEta4;
    double hiHFminusEta4  = c->evt.hiHFminusEta4;
    double hiHFplus       = c->evt.hiHFplus;
    double hiHFminus      = c->evt.hiHFminus;

    //! Parametrization of hfPlus and hfsum
    //double hiHFSum = 1.17 * hiHFplus + 6; //! Doga's parametrization
    double hiHFSum   = hiHFplusEta4 + hiHFminusEta4; //! same as hiHF but wih different centrality selections

    double hiZDC           = c->evt.hiZDC;
    double hiZDCplus       = c->evt.hiZDCplus;
    double hiZDCminus      = c->evt.hiZDCminus;

    //! bins
    int bin[ngrp];
    for(int i=0; i<ngrp; i++)bin[i] = -1;
    int msel=-1;
    bin[++msel]  = GetNtracksBin(hiNtracks);  //! 0
    bin[++msel]  = GetNPixBin(hiNpix);        //! 1 
    bin[++msel]  = GetNPixelTracksBin(hiNpixelTracks); //! 2
    bin[++msel]  = GetHFBin(hiHF);                     //! 3
    bin[++msel]  = GetAnaBin(hiHFplusEta4);            //! 4
    bin[++msel]  = GetNewAnaBin(hiHFSum);              //! 5
    bin[++msel]  = GetHFPE4Bin(hiHFplusEta4);          //! 6
    bin[++msel]  = GetHFME4Bin(hiHFminusEta4);         //! 7
    bin[++msel]  = GetHFPBin(hiHFplus);                //! 8
    bin[++msel]  = GetHFMBin(hiHFminus);               //! 9

    //! ZDC bad runs
    bool badRun = RunNo==210658 || RunNo==210676 || RunNo==210635 || RunNo==210638;
    if(!badRun){
      bin[++msel]  = GetZDCBin(hiZDC);  //! 10
      bin[++msel]  = GetZDCPBin(hiZDCplus);//! 11
      bin[++msel]  = GetZDCMBin(hiZDCminus);//! 12
    }else{
      bin[++msel]  = -1;
      bin[++msel]  = -1;
      bin[++msel]  = -1;
    }
    if(hiHFSum<100)bin[++msel]  = 0; //! 13  inclusive (0-100%) 

    if(ievt%10000==0)std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<RunNo<<"\t HF : "<<hiHF<<std::endl;

    for(int nj=0;nj<knj;nj++){
      
      if(nj==0)mJets = &(c->ak3PF);
      else if(nj==1)mJets = &(c->akPu3PF);
      else if(nj==2)mJets = &(c->ak3Calo);
      else if(nj==3)mJets = &(c->akPu3Calo);

      int *ljet = new int[3];
      FindLeadSubLeadJets(mJets,ljet);

      //! Jet energy resolution and scale
      //if(ljet[0]>=0 && ljet[1]>=0 && ljet[2]<0 && mJets->jtpt[ljet[0]]>80. && mJets->jtpt[ljet[1]]>80. && (mJets->trackMax[ljet[0]]>4. || mJets->trackMax[ljet[1]]>4.)){
      if(ljet[0]>=0 && ljet[1]>=0 && mJets->jtpt[ljet[0]]>80. && mJets->jtpt[ljet[1]]>80. && (mJets->trackMax[ljet[0]]>4. || mJets->trackMax[ljet[1]]>4.)){

	double ptave    = (mJets->jtpt[ljet[0]] + mJets->jtpt[ljet[1]])/2.;
	double dphi  = delphi(mJets->jtphi [ljet[1]],mJets->jtphi [ljet[0]]);

	int mstat=1;
	if(ljet[2]>=0){
	  if((mJets->jtpt[ljet[2]]/ptave)  > kalpha)mstat=0;
	}

	if(dphi >= 2.7 && mstat){

          double rn1 = gRandom->Rndm();
          double rn2 = gRandom->Rndm();

	  if((fabs(mJets->jteta[ljet[0]])<1.3 || fabs(mJets->jteta[ljet[1]])<1.3)){ //! balance
	    int probe=-1, barel=-1;
	    if(fabs(mJets->jteta[ljet[0]])<1.3 && fabs(mJets->jteta[ljet[1]])>1.3){barel=0;probe=1;}
            else if(fabs(mJets->jteta[ljet[1]])<1.3 && fabs(mJets->jteta[ljet[0]])>1.3){barel=1;probe=0;}
            else if(fabs(mJets->jteta[ljet[0]])<1.3 && fabs(mJets->jteta[ljet[1]])<1.3){
              if(rn1 > rn2){
                barel=0;probe=1;
              }else{
                barel=1;probe=0;
              }
            }
	    if(probe>=0 && barel>=0){
	      double ptprobe  = mJets->jtpt[ljet[probe]];
              double ptbarel  = mJets->jtpt[ljet[barel]];
	      double etaprobe = mJets->jteta[ljet[probe]];
	      double etabarel = mJets->jteta[ljet[barel]];
	      double B = (mJets->jtpt[ljet[probe]] - mJets->jtpt[ljet[barel]])/ptave;
	      
	      if(bin[5]>=0 && bin[5]<nbins){
		hetaprobe[nj]->Fill(etaprobe,ptprobe);
                hetabarel[nj]->Fill(etabarel,ptbarel);
                hjes[nj]->Fill(etaprobe,ptave,(ptprobe/ptbarel));
		hptbptp[nj]->Fill(etaprobe,ptprobe,ptbarel);
		hJetEnergyBalance[nj][bin[5]]->Fill(etaprobe,ptave,B);
		hJetEnergyBalance[nj][nbins]->Fill(etaprobe,ptave,B);	  
	      }
	    }
	  }//! balance 
	  
	  double A=-9999;
	  if(rn1 > rn2){
            A = (mJets->jtpt[ljet[0]] - mJets->jtpt[ljet[1]])/(mJets->jtpt[ljet[0]] + mJets->jtpt[ljet[1]]);
          }else{
            A = (mJets->jtpt[ljet[1]] - mJets->jtpt[ljet[0]])/(mJets->jtpt[ljet[1]] + mJets->jtpt[ljet[0]]);
          }

	  int etab=-1;
	  if((fabs(mJets->jteta[ljet[0]])<1.3 && fabs(mJets->jteta[ljet[1]])<1.3))etab=0;
	  if((fabs(mJets->jteta[ljet[0]])<1.3 && fabs(mJets->jteta[ljet[1]])>1.3))etab=1;
	  if((fabs(mJets->jteta[ljet[0]])>1.3 && fabs(mJets->jteta[ljet[1]])<1.3))etab=2;
	  if((fabs(mJets->jteta[ljet[0]])>1.3 && fabs(mJets->jteta[ljet[1]])>1.3))etab=3;

	  if(etab>=0 && fabs(A)<0.4){
	    hJetEnergyRes[nj][bin[5]]->Fill(etab,ptave,A);
	    hJetEnergyRes[nj][nbins]->Fill(etab,ptave,A);
	  }
	}
      }
      delete [] ljet;
      continue;

      int lead_jet [lbins]={0}; //! 6 bins         
      int slead_jet[lbins]={0}; //! 6 bins         
      //! {>80, >100, > 120, > 150, > 200, > 300}
                                   
      double ljetpt  = -999;
      double sljetpt = -999;
      if(ljet[0]<0)continue;
      else if(ljet[0] >=0 && ljet[1] >=0){

        ljetpt  = mJets->jtpt[ljet[0]];
	sljetpt = mJets->jtpt[ljet[1]];

        if(ljetpt > 80 )lead_jet[0]=1;
        if(ljetpt > 100)lead_jet[1]=1;
        if(ljetpt > 120)lead_jet[2]=1;
        if(ljetpt > 150)lead_jet[3]=1;
        if(ljetpt > 200)lead_jet[4]=1;
        if(ljetpt > 300)lead_jet[5]=1;

        if(sljetpt > 30)slead_jet[0]=1;
        if(sljetpt > 40)slead_jet[1]=1;
        if(sljetpt > 50)slead_jet[2]=1;
        if(sljetpt > 60)slead_jet[3]=1;
        if(sljetpt > 70)slead_jet[4]=1;
        if(sljetpt > 80)slead_jet[5]=1;

	if(nj==0){
	  hHFSum->Fill(bin[5],ljetpt,hiHFSum);
	  for(int il=0;il<lbins;il++){
	    if(lead_jet[il]){
	      hHFSum_LeadJet[il]->Fill(hiHFSum);
	    }
	  }
	}

	double dphi  = delphi(mJets->jtphi [ljet[1]],mJets->jtphi [ljet[0]]);
	
	if(lead_jet[leadb]){
	  if(slead_jet[sleadb]){
	    if(mJets->trackMax[ljet[0]]>4. || mJets->trackMax[ljet[1]]>4.){
	      
	      int jetbin = leadb*lbins + sleadb; 
	      
	      for(int i=0;i<ngrp;i++){
		
		if(bin[i]<0 || bin[i]>=nbins)continue;
		int histb =  bin[i] + i*nbins;
		
		hdphi  [nj][0]->Fill(histb,jetbin,dphi,wxs);
		hdedp  [nj][0][jetbin]->Fill(histb,(mJets->jteta[ljet[1]] - mJets->jteta[ljet[0]]),dphi,wxs);
		
		if ( dphi >= 2*pi/3. ){
		  hpt2pt1[nj][0]->Fill(histb,jetbin,(mJets->jtpt [ljet[1]] / mJets->jtpt [ljet[0]]),wxs);
		  hdeta  [nj][0]->Fill(histb,jetbin,(mJets->jteta[ljet[1]] - mJets->jteta[ljet[0]]),wxs);
		  hdjeta [nj][0]->Fill(histb,jetbin,(mJets->jteta[ljet[1]] + mJets->jteta[ljet[0]])/2.,wxs);
		  hypt2pt1   [nj][0][jetbin]->Fill(histb,(mJets->jteta[ljet[1]] - mJets->jteta[ljet[0]])   ,(mJets->jtpt[ljet[1]] / mJets->jtpt[ljet[0]]),wxs);
		  hdijept2pt1[nj][0][jetbin]->Fill(histb,(mJets->jteta[ljet[1]] + mJets->jteta[ljet[0]])/2.,(mJets->jtpt[ljet[1]] / mJets->jtpt[ljet[0]]),wxs);
		  
		  for(int j=0;j<2;j++){
		    hjetpt [nj][j]->Fill(histb,jetbin,mJets->jtpt [ljet[j]],wxs);
		    hjetphi[nj][j]->Fill(histb,jetbin,mJets->jtphi[ljet[j]],wxs);
		    hjeteta[nj][j]->Fill(histb,jetbin,mJets->jteta[ljet[j]],wxs);
		  }
		}
	      }
	    }//! trackMax cut
	  }//! sub-leading jet bin

	  //! Inclusive away-side jet
	  for(int ij=0; ij<mJets->nref; ij++){
	    if(ljet[0] == ij)continue;
	    //! jet selction cut 
	    
	    bool seljet = false;
	    if(strcmp(ksp,"ppb")==0 || strcmp(ksp,"pbp")==0){
	      seljet = mJets->jtpt[ij]>kptrecocut && mJets->rawpt[ij]>15. && fabs(mJets->jteta[ij])<ketacut && (mJets->trackMax[ij] > 4. || mJets->trackMax[ljet[0]] > 4.);
	    }else {
	      seljet = mJets->jtpt[ij]>kptrecocut && mJets->rawpt[ij]>15. && fabs(mJets->jteta[ij])<ketacut && (mJets->trackMax[ij] > 4. || mJets->trackMax[ljet[0]] > 4.);
	      //seljet = mJets->jtpt[ij]>kptrecocut && mJets->rawpt[ij]>15. && fabs(mJets->jteta[ij])<ketacut && (mJets->trackMax[ij]/mJets->jtpt[ij]) > 0.01;
	    }
	    if(!seljet)continue;
	    
	    double jetpt  = mJets->jtpt[ij];
	    double jetphi = mJets->jtphi[ij];
	    double jeteta = mJets->jteta[ij];
	    
	    int sljetbin[lbins] = {-1,-1,-1,-1,-1,-1};
	    if(jetpt > 30)sljetbin[0]=1;
	    if(jetpt > 40)sljetbin[1]=1;
	    if(jetpt > 50)sljetbin[2]=1;
	    if(jetpt > 60)sljetbin[3]=1;
	    if(jetpt > 70)sljetbin[4]=1;
	    if(jetpt > 80)sljetbin[5]=1;
	    
	    double dphi1    = delphi(jetphi,mJets->jtphi[ljet[0]]);
	    double deta1    = jeteta - mJets->jteta[ljet[0]];
	    
	    if(sljetbin[sleadb]){
	      int incljetb = leadb*lbins + sleadb;
	      
	      for(int i=0;i<ngrp;i++){
		if(bin[i]<0 || bin[i]>=nbins)continue;
		int histb =  bin[i] + i*nbins;
		
		hdphi  [nj][1]->Fill(histb,incljetb,dphi1,wxs);
		hdedp  [nj][1][incljetb]->Fill(histb,deta1,dphi1);
		
		if ( dphi1 >= 2*pi/3. ){
		  hpt2pt1[nj][1] ->Fill(histb,incljetb,(jetpt/mJets->jtpt[ljet[0]]),wxs);
		  hdeta  [nj][1] ->Fill(histb,incljetb,(jeteta - mJets->jteta[ljet[0]]),wxs);
		  hdjeta [nj][1] ->Fill(histb,incljetb,(jeteta + mJets->jteta[ljet[0]])/2.,wxs);
		  hypt2pt1   [nj][1][incljetb]->Fill(histb,(jeteta - mJets->jteta[ljet[0]])   ,(jetpt / mJets->jtpt[ljet[0]]),wxs);
		  hdijept2pt1[nj][1][incljetb]->Fill(histb,(jeteta + mJets->jteta[ljet[0]])/2.,(jetpt / mJets->jtpt[ljet[0]]),wxs);
		  
		  hjetpt [nj][2]->Fill(histb,incljetb,jetpt,wxs);
		  hjetphi[nj][2]->Fill(histb,incljetb,jetphi,wxs);
		  hjeteta[nj][2]->Fill(histb,incljetb,jeteta,wxs);
		}
	      }//! ngrp loop
	    }//! sljetbin
	  }//! nref loop
	}//! lead-jet bin
      }//! lead jet condition              
      //delete [] ljet;
    }//! nj loop
    
    //! Global distribution
    hVz->Fill(vz);
    hHFplusEta4NTracks->Fill(hiHFplusEta4,hiNtracks);
    hHFplusEta4Npix   ->Fill(hiHFplusEta4,hiNpix);
    hHFplusEta4NpixelTracks->Fill(hiHFplusEta4,hiNpixelTracks);
    hHFplusEta4HF          ->Fill(hiHFplusEta4,hiHF);
    hHFplusEta4HFminusEta4 ->Fill(hiHFplusEta4,hiHFminusEta4);
    hHFplusEta4HFminus->Fill(hiHFplusEta4,hiHFminus);
    hHFplusEta4HFplus->Fill(hiHFplusEta4,hiHFplus);
    hHFplusEta4ZDC->Fill(hiHFplusEta4,hiZDC);
    hHFplusEta4ZDCplus->Fill(hiHFplusEta4,hiZDCplus);
    hHFplusEta4ZDCminus->Fill(hiHFplusEta4,hiZDCminus);

    hHFSumNTracks     ->Fill(hiHFSum,hiNtracks,wxs);
    hHFSumNpix        ->Fill(hiHFSum,hiNpix,wxs);
    hHFSumNpixelTracks->Fill(hiHFSum,hiNpixelTracks,wxs);
    hHFSumHF          ->Fill(hiHFSum,hiHF,wxs);
    hHFSumHFplusEta4  ->Fill(hiHFSum,hiHFplusEta4,wxs);
    hHFSumHFminusEta4 ->Fill(hiHFSum,hiHFminusEta4,wxs);
    hHFSumHFminus     ->Fill(hiHFSum,hiHFminus,wxs);
    hHFSumHFplus      ->Fill(hiHFSum,hiHFplus,wxs);
    hHFSumZDC         ->Fill(hiHFSum,hiZDC,wxs);
    hHFSumZDCplus     ->Fill(hiHFSum,hiZDCplus,wxs);
    hHFSumZDCminus    ->Fill(hiHFSum,hiZDCminus,wxs);

  }//! event loop ends


  //! Write to output file
  fout->cd();
  fout->Write();
  fout->Close();

  //! Check
  timer.Stop();
  double rtime  = timer.RealTime();
  double ctime  = timer.CpuTime();

  std::cout<<"\t"<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;

  return 0;
}
double delphi(double phi1, double phi2)
{
  double dphi = phi2 - phi1;
  dphi = fabs(atan2(sin(dphi),cos(dphi)));
  return dphi;
}
int GetNtracksBin(int ntrks)
{
  int ibin=-1;
  if(ntrks>=btrks[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(ntrks<btrks[i] && ntrks>=btrks[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetNPixBin (int npix)
{
  int ibin=-1;
  if(npix>=bpixhit[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(npix<bpixhit[i] && npix>=bpixhit[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetNPixelTracksBin(int npixtrk)
{
  int ibin=-1;
  if(npixtrk>=bpixtr[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(npixtrk<bpixtr[i] && npixtrk>=bpixtr[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetHFBin   (double hf)
{
  int ibin=-1;
  if(hf>=bhf[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(hf<bhf[i] && hf>=bhf[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetAnaBin(double hf)
{
  int ibin=-1;
  if(hf>=bana[0] && hf<70){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(hf<bana[i] && hf>=bana[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetNewAnaBin(double hf)
{
  int ibin=-1;
  if(hf>=bana_new[0] && hf<100){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(hf<bana_new[i] && hf>=bana_new[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}

int GetHFPE4Bin(double hf)
{
  int ibin=-1;
  if(hf>=bhfpe4[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(hf<bhfpe4[i] && hf>=bhfpe4[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetHFME4Bin(double hf)
{
  int ibin=-1;
  if(hf>=bhfme4[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(hf<bhfme4[i] && hf>=bhfme4[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetHFPBin(double hf)
{
  int ibin=-1;
  if(hf>=bhfp[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(hf<bhfp[i] && hf>=bhfp[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetHFMBin(double hf)
{
  int ibin=-1;
  if(hf>=bhfm[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(hf<bhfm[i] && hf>=bhfm[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetZDCBin(double zdc)
{
  int ibin=-1;
  if(zdc>=bzdc[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(zdc<bzdc[i] && zdc>=bzdc[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetZDCPBin(double zdc)
{
  int ibin=-1;
  if(zdc>=bzdcp[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(zdc<bzdcp[i] && zdc>=bzdcp[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}
int GetZDCMBin(double zdc)
{
  int ibin=-1;
  if(zdc>=bzdcm[0]){
    ibin=0;
    return ibin;
  }
  for(int i=0;i<nbins-1;i++){
    if(zdc<bzdcm[i] && zdc>=bzdcm[i+1]){
      ibin=i+1;
      break;
    }
  }
  return ibin;
}


void FindLeadSubLeadJets(Jets *jetc, int *ljet)
{
  ljet[0]=-1; ljet[1]=-2; ljet[2]=-3;
  
  double tempt=-9;
  //! Get the leading jet
  for(int ij=0; ij<jetc->nref; ij++){
    //if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<80 || jetc->rawpt[ij]<15 || jetc->trackMax[ij]<4.)continue;
    if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<50 || jetc->rawpt[ij]<15)continue;
    double jetpt = jetc->jtpt[ij];
    if(jetpt > tempt){
      tempt = jetpt;
      ljet[0] = ij;
    }
  }

  if(ljet[0]>=0){
    // Subleading
    tempt=-9;
    for(int ij=0; ij<jetc->nref; ij++){
      if(ij==ljet[0])continue;
      //if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<30 || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15 || jetc->trackMax[ij]<4.)continue;
      if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<30 || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15)continue;
      double jetpt = jetc->jtpt[ij];
      //double dphi  = delphi(mJets->jtphi[ljet[1]],mJets->jtphi[ljet[0]]);
      //double dphi  = jetc->jtphi[ij] - jetc->jtphi[ljet[0]];
      //if (dphi > pi ) dphi = dphi - 2 * pi;
      //if (dphi < -pi) dphi = dphi + 2 * pi;
      //if(dphi < 2*pi/3.)continue;
      if (jetpt > tempt){
	tempt = jetpt;
	ljet[1] = ij;
      }
    }

    if(ljet[1]>=0){
      tempt=-9;
      for(int ij=0; ij<jetc->nref; ij++){
	if(ij==ljet[0] || ij==ljet[1])continue;
        if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<30 || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15)continue;
        float jetpt = jetc->jtpt[ij];
        if (jetpt > tempt){
          tempt = jetpt;
          ljet[2] = ij;
        }
      }
    }
  }
}
void LoadLib()
{
  gSystem->Load("../HiForest/V2/hiForest_h.so");
}
void ShutoffBranches(HiForest *hi,const char *ksp)
{

  //! added by pawan
  //! Select only the branches you want to use for the analysis
  //! This increases the speed for running

  //! For Hlt
  hi->hltTree->SetBranchStatus("*",0,0);
  if(strcmp(ksp,"ppb")==0 || strcmp(ksp,"pbp")==0){//! pPb
    hi->hltTree->SetBranchStatus("HLT_PAJet40_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet60_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet80_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet100_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAZeroBiasPixel_SingleTrack_v1",1,0);
  }else if(strcmp(ksp,"pbpb")==0){  //! PbPb
    hi->hltTree->SetBranchStatus("HLT_HIJet80_v1",1,0);
  }
  else if(strcmp(ksp,"pp")==0){  //! pp
    //hi->hltTree->SetBranchStatus("HLT_Jet60_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet80_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet100_NoJetID_v1",1,0);
  }

  //! for Skim Tree
  hi->skimTree->SetBranchStatus("*",0,0);
  hi->skimTree->SetBranchStatus("pHBHENoiseFilter",1,0);
  if(strcmp(ksp,"ppb")==0 || strcmp(ksp,"pbp")==0){
    hi->skimTree->SetBranchStatus("pVertexFilterCutGplus",1,0);
    hi->skimTree->SetBranchStatus("pPAcollisionEventSelectionPA",1,0);
  }else if(strcmp(ksp,"pbpb")==0){
    hi->skimTree->SetBranchStatus("pcollisionEventSelection",1,0);
  }
  else if(strcmp(ksp,"pp")==0){
    hi->skimTree->SetBranchStatus("pcollisionEventSelection",1,0);
    hi->skimTree->SetBranchStatus("phfPosFilter1",1,0);
    hi->skimTree->SetBranchStatus("phfNegFilter1",1,0);
  }




  //! Evt tree
  hi->evtTree->SetBranchStatus("*",0,0);
  hi->evtTree->SetBranchStatus("run",1,0);
  hi->evtTree->SetBranchStatus("vz",1,0);
  hi->evtTree->SetBranchStatus("hiHF",1,0);//! HF
  hi->evtTree->SetBranchStatus("hiHFplusEta4",1,0);
  hi->evtTree->SetBranchStatus("hiHFminusEta4",1,0);
  hi->evtTree->SetBranchStatus("hiHFplus",1,0);
  hi->evtTree->SetBranchStatus("hiHFminus",1,0);
  hi->evtTree->SetBranchStatus("hiNtracks",1,0);
  hi->evtTree->SetBranchStatus("hiNpix",1,0);
  hi->evtTree->SetBranchStatus("hiNpixelTracks",1,0);
  hi->evtTree->SetBranchStatus("hiZDC",1,0);
  hi->evtTree->SetBranchStatus("hiZDCplus",1,0);
  hi->evtTree->SetBranchStatus("hiZDCminus",1,0);


  /*
  //! Track tree
  hi->trackTree->SetBranchStatus("*",0,0);
  hi->trackTree->SetBranchStatus("nTrk",1,0);
  hi->trackTree->SetBranchStatus("trkPt",1,0);
  hi->trackTree->SetBranchStatus("trkEta",1,0);
  hi->trackTree->SetBranchStatus("trkPhi",1,0);
  hi->trackTree->SetBranchStatus("highPurity",1,0);
  hi->trackTree->SetBranchStatus("trkDz1",1,0);
  hi->trackTree->SetBranchStatus("trkDzError1",1,0);
  hi->trackTree->SetBranchStatus("trkDxy1",1,0);
  hi->trackTree->SetBranchStatus("trkDxyError1",1,0);
  */


  hi->ak3PFJetTree->SetBranchStatus("*",0,0);
  hi->ak3PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak3PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak3PFJetTree->SetBranchStatus("trackMax",1,0);

  hi->akPu3PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu3PFJetTree->SetBranchStatus("nref",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("trackMax",1,0);

  hi->ak3CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak3CaloJetTree->SetBranchStatus("nref",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("trackMax",1,0);


  hi->akPu3CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu3CaloJetTree->SetBranchStatus("nref",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("trackMax",1,0);



  /*
  //! PF jet tree
  hi->ak2PFJetTree->SetBranchStatus("*",0,0);
  hi->ak2PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak2PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak2PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak2PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak2PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak2PFJetTree->SetBranchStatus("trackMax",1,0);
  */


  /*
  hi->ak4PFJetTree->SetBranchStatus("*",0,0);
  hi->ak4PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak4PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak4PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak4PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak4PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak4PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  /*
  hi->ak5PFJetTree->SetBranchStatus("*",0,0);
  hi->ak5PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak5PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak5PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak5PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak5PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak5PFJetTree->SetBranchStatus("trackMax",1,0);
  */

  //
  /*
  hi->akPu2PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu2PFJetTree->SetBranchStatus("nref",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("trackMax",1,0);
  */


  /*
  hi->akPu4PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu4PFJetTree->SetBranchStatus("nref",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  /*
  hi->akPu5PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu5PFJetTree->SetBranchStatus("nref",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  //


  //! Calo jet trees
  /*
  hi->ak2CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak2CaloJetTree->SetBranchStatus("nref",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("trackMax",1,0);


  hi->ak4CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak4CaloJetTree->SetBranchStatus("nref",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("trackMax",1,0);

  hi->ak5CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak5CaloJetTree->SetBranchStatus("nref",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("trackMax",1,0);
  */

  ////
  /*
  hi->akPu2CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu2CaloJetTree->SetBranchStatus("nref",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("trackMax",1,0);


  hi->akPu4CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu4CaloJetTree->SetBranchStatus("nref",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("trackMax",1,0);

  hi->akPu5CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu5CaloJetTree->SetBranchStatus("nref",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("trackMax",1,0);
  */
  ///
  std::cout<<"Loaded all tree variables "<<std::endl;

}
