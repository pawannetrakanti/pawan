#include "/afs/cern.ch/user/p/pawan/scratch0/CMSSW_6_2_0/src/work/pPb/HiForest/V3/hiForest.h"
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


#ifdef _MAKECINT_
//#pragma link off all globals;
//#pragma link off all classes;
//#pragma link off all functions;
#pragma link C++ class HiForest;
#pragma link C++ class Hlts;
#pragma link C++ class Skims;
#pragma link C++ class Evts;
#pragma link C++ class Jets;
#endif

const double pi=acos(-1.);
const double pi2=2*pi -1;


const int knj=12;    //! algo 0: ic5, 1: ak3pu
const char *calgo[knj]={
  "ak3PF","ak4PF","ak5PF",
  "akPu3PF","akPu4PF","akPu5PF",
  "ak3Calo","ak4Calo","ak5Calo",
  "akPu3Calo","akPu4Calo","akPu5Calo"
};

const int maxruns=12;
int Run[maxruns];

void LoadRuns  (const char */*sp*/);
int GetRunIdx  (int /*irun*/);
void GetRunNos(HiForest */*c*/);
bool ascend(int /*run1*/,int /*run2*/);

void GetXsection(const char */*pthat*/,double &/*xsection*/,double &/*maxpthat*/);

//! For jet energy balance
double delphi(double /*phi1*/, double /*phi2*/);

const float ketacut=2.0;
double kptrecocut=30.;
const bool iCountRuns=false;
const double kDelr = 0.3;

TStopwatch timer;
int qa(const char *ctrig="HLT_PAJet4060_NoJetID_v1",
       TString filename="root://eoscms//eos/cms/store/caf/user/yjlee/PP2013_PromptReco_Json_Jet40Jet60_v84_ppTrack/pp_hiForest2_100_1_sAo.root",
       int ifile=0,
       bool isMC=false
)
{

  timer.Start();
  const char *ctracking="pptracking";
  //! Load Lib
  //gSystem->Load("/net/hidsk0001/d00/scratch/pawan/service/CMSSW_6_0_0/src/test/pp/HiForest/V3/hiForest_h.so");
  gSystem->Load("/afs/cern.ch/user/p/pawan/scratch0/CMSSW_6_2_0/src/work/pPb/HiForest/V3/hiForest_h.so");
  const char *ksp="pp";

  //! Load Run numbers 
  for(int ir=0;ir<maxruns;ir++)Run[ir]=-999;
  if(!iCountRuns)LoadRuns(Form("%s",ksp));

  //! Define the input file and HiForest
  //! pp 2013 data w HI Iterative tracking
  TString inname="";
  HiForest *c=0;
  TFile *fout=0;
  

  if(isMC)inname = Form("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt%s/HiForest_v81_merged01/pt%s_pp2013_P01_prod22_v81_merged_forest_0.root",ctrig,ctrig);
  else    inname = filename;
  c = new HiForest(inname,Form("MyForest%s",ksp),cPP);

  if(isMC)fout = new TFile(Form("/net/hidsk0001/d00/scratch/pawan/service/CMSSW_6_0_0/src/test/pp/output/qa/pp_pthat%sGeV_%s_HiForest2_pt%0.0f_2013_prod22_v81_%d.root",ctrig,ctracking,kptrecocut,ifile),"RECREATE");
  else fout = new TFile(Form("output/%s/pp_%s_%s_HiForest2_pt%0.0f_2013_v84_%d.root",ctrig,ctrig,ctracking,kptrecocut,ifile),"RECREATE");        

  Long64_t nb = 0;
  Long64_t nentries = c->GetEntries();

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Running for %s ",ksp)<<std::endl;
  std::cout<<Form("pT  cut for %0.3f ",kptrecocut)<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;

  std::cout<<"Input  file  "<<inname.Data()<<std::endl;
  std::cout<<"HiForest Tree : " <<c->GetName()<<"\t Entries "<<nentries<<std::endl;
  std::cout<<"Output file  "<<fout->GetName()<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  //return 0;


  double wcen=1;
  double wvz=1;

  double xsection=-1;
  double maxpthat=9999;
  double wxs=1;  



  if(iCountRuns){//! switch off all the trees and only keep in evtTree
    GetRunNos(c);
    return 0;
  }else{
    if(isMC){
      c->hasHltTree=0;
      c->hasSkimTree=0;
      c->hasIcPu5JetTree=0;
      c->hasAk1PFJetTree=0;
      c->hasAk2PFJetTree=0;
      c->hasAk6PFJetTree=0;
      c->hasAkPu1PFJetTree=0;
      c->hasAkPu2PFJetTree=0;
      c->hasAkPu6PFJetTree=0;
      c->hasAk1CaloJetTree=0;
      c->hasAk2CaloJetTree=0;
      c->hasAk6CaloJetTree=0;
      c->hasAkPu1CaloJetTree=0;
      c->hasAkPu2CaloJetTree=0;
      c->hasAkPu6CaloJetTree=0;
      c->hasTrackTree=0;

      GetXsection(ctrig,xsection,maxpthat);
      wxs = xsection/(nentries/100000.);

      std::cout<<"MC sample for pT hat "<<ctrig<<"GeV \t maxpthat "<<maxpthat<<"\t xsection : "<<xsection<<"\t weight  : "<<wxs<<std::endl;
      std::cout<<std::endl;
      //! weight  for the merging of the samples for different pT hat bins           

    }else{
      c->hltTree->SetBranchStatus("*",0,0);
      c->hltTree->SetBranchStatus("HLT_PAJet40_NoJetID_v1",1,0);
      c->hltTree->SetBranchStatus("HLT_PAJet60_NoJetID_v1",1,0);
      c->hltTree->SetBranchStatus("HLT_PAJet80_NoJetID_v1",1,0);
      c->hltTree->SetBranchStatus("HLT_PAJet100_NoJetID_v1",1,0);
      
      //! for Skim Tree
      c->skimTree->SetBranchStatus("*",0,0);
      c->skimTree->SetBranchStatus("pcollisionEventSelection",1,0);
      c->skimTree->SetBranchStatus("pPAcollisionEventSelectionPA",1,0);
      c->skimTree->SetBranchStatus("pHBHENoiseFilter",1,0);
      //c->skimTree->SetBranchStatus("phfPosFilter1",1,0);
      //c->skimTree->SetBranchStatus("phfNegFilter1",1,0);
    }
  }
  std::cout<<"Loaded all tree variables "<<std::endl;
  
  
  //! For jets
  Jets *mJets=0;
  int scaFac=5;

  //! Define histograms here
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  TProfile::SetDefaultSumw2();
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  char mrun[6]={0};
  TH1F *hEvt = new TH1F("hEvt","# of events ",4,0,4);
  TH1F *hbin = new TH1F("hbin","centrality bin",40,-0.5,40-0.5);

  //! Data related historgams
  TH2F *hEvtRun = new TH2F("hEvtRun","# of events in that run",maxruns,-0.5,maxruns-0.5,40,-40,40);
  TH1F *hLumiBlock = new TH1F("hLumiBlock","Luminosity block from hlt",scaFac*1500,0,scaFac*1500);
  TH2F *hHFPlusMinus = new TH2F("hHFPlusMinus","HF plus:minus",200,0,400,200,0,400);

  TH1F *hHF     = new TH1F("hHF","HF distribution",200,0,600);
  TH1F *hAvHF     = new TH1F("hAvHF","HF E/hit distribution",100,0,0.1);
  TH1F *hAvHFplus = new TH1F("hAvHFplus","HFplus E/hit distribution",100,0,0.5);
  TH1F *hAvHFminus = new TH1F("hAvHFminus","HFminus E/hit distribution",100,0,0.5);

  
  TH2F *hAvHFRun     = new TH2F("hAvHFRun","HF E/hit distribution Run-by-Run",maxruns,-0.5,maxruns-0.5,100,0,0.5);
  TH2F *hAvHFplusRun = new TH2F("hAvHFplusRun","HFplus E/hit distribution Run-by-Run",maxruns,-0.5,maxruns-0.5,100,0,0.5);
  TH2F *hAvHFminusRun = new TH2F("hAvHFminusRun","HFminus E/hit distribution Run-by-Run",maxruns,-0.5,maxruns-0.5,100,0,0.5);
  TH2F *hRatioHFplusminusRun = new TH2F("hRatioHFplusminusRun","Ratio of HFplus/HFminus vs Run Id",maxruns,-0.5,maxruns-0.5,100,0.,5.);  
  TH2F *hNTracksRun   = new TH2F("hNTracksRun","NTracks vs RunId",maxruns,-0.5,maxruns-0.5,200,0.,400.);    
  TH2F *hNpixelTracksRun = new TH2F("hNpixelTracksRun","hiNpixelTracks vs RunId",maxruns,-0.5,maxruns-0.5,200,0,300);
  TH2F *hLumiRun = new TH2F("hLumiRun","Luminosity Run-by-Run",maxruns,-0.5,maxruns-0.5,scaFac*100,0,scaFac*2000);
  hLumiRun->GetXaxis()->SetLabelSize(0.03);
  for(int ix=1;ix<=hLumiRun->GetNbinsX();ix++){
    if(Run[ix-1]!=-999)sprintf(mrun,"%d",Run[ix-1]);
    else sprintf(mrun,"%s","\0");
    hLumiRun->GetXaxis()->SetBinLabel(ix,mrun);
  }
  TH2F *hJetPtRun[knj];
  TH2F *hJetPURun[knj];
  TH2F *hJetBkRun[knj];
  TH2F *hJetNoRun[knj];

  TH2F *hZDCPlusMinus = new TH2F("hZDCPlusMinus","ZDC plus:minus",500,0,5000,500,0,5000);
  TH1F *hZDC     = new TH1F("hZDC","ZDC distribution",500,0,5000);
  
  //! EE
  TH1F *hET = new TH1F("hET","hiET",200,0,600);
  TH1F *hEB = new TH1F("hEB","hiEB",200,0,1600);
  TH1F *hEE = new TH1F("hEE","hiEE",200,0,600);
  TH2F *hEBEE = new TH2F("hEBEE","hEBEE",200,0,1600,200,0,600);    
  TH2F *hEEPlusMinus = new TH2F("hEEPlusMinus","hEEPlusEEMinus",200,0,400,200,0,400);    
  
  //! vertex
  TH1F *hvz   = new TH1F("hvz","hvz",80,-20,20);
  TH2F *hvxvy = new TH2F("hvxvy","vx vs vy",100,-0.2,0.2,100,-0.2,0.2);
  TH2F *hvzrun   = new TH2F("hvzrun","hvzrun",maxruns,-0.5,maxruns-0.5,80,-20,20);

  //! Pixel
  TH1F *hNpix = new TH1F("hNpix","hiNpix",1000,0,2000);
  TH1F *hNpixelTracks = new TH1F("hNpixelTracks","hiNpixelTracks",200,0,300);
  TH1F *hNTracks = new TH1F("hNTracks","hiNTracks",200,0,400);
  TH2F *hNpixelTracksNTracks = new TH2F("hNpixelTracksNTracks","hNpixelTracksNTracks",200,0,400,200,0,400);

  //! Cross correlations
  TH2F *hHFZDC = new TH2F("hHFZDC","hHFZDC",800,0,600,1000,0,5000);
  TH2F *hNpixelTracksZDC = new TH2F("hNpixelTracksZDC","hNpixelTracksZDC",200,0,300,1000,0,5000);
  TH2F *hHFNpixelTracks  = new TH2F("hHFNpixelTracks","hHFNpixelTracks",800,0,600,200,0,300);
  
  //! jet variables
  TH2F *hjrpt[knj], *hjppt[knj], *hjbpt[knj], *hjeta[knj], *hjphi[knj];
  TH2F *hjrapt[knj];
  TH2F *hetaphi[knj];
  TH1F *hNjets[knj];

  TH3F *hjb[knj];


  TH2F *htrackmax[knj], *hphotonmax[knj];
  TH2F *hchsum[knj], *hnesum[knj], *hphsum[knj];
  TH2F *hhcalsum[knj], *hecalsum[knj];
  
  TH2F *hra_trackmax[knj], *hra_photonmax[knj];
  TH2F *hra_chsum[knj], *hra_nesum[knj], *hra_phsum[knj];
  TH2F *hra_hcalsum[knj], *hra_ecalsum[knj];

  TH1F *hpthat  = new TH1F("hpthat","py-hat distribution",500,0,1000);  
  
  for(int nj=0;nj<knj;nj++){
    //std::cout<<Form("Histograms for algo %s",calgo[nj])<<std::endl;

    //! Define the histograms for jets
    hNjets[nj] = new TH1F(Form("hNjets%d",nj),Form("# of  %s jets",calgo[nj]),3,0.0,3.0);

    hjrpt[nj]  = new TH2F(Form("hjrpt%d",nj),Form("jet reco:raw  p_{T} distribution jet %s",calgo[nj]),500,0,1000,500,0,1000);
    hjppt[nj]  = new TH2F(Form("hjppt%d",nj),Form("jet reco:pilup pt p_{T} distribution jet %s",calgo[nj]),500,0,1000,300,0,600);
    hjbpt[nj]  = new TH2F(Form("hjbpt%d",nj),Form("jet reco:bkgd distribution jet  %s",calgo[nj]),500,0,1000,300,0,600);
    hjeta[nj]  = new TH2F(Form("hjeta%d",nj),Form("jet recopT:eta distribution jet %s",calgo[nj]),500,0,1000,18,-3.0,3.0);
    hjphi[nj]  = new TH2F(Form("hjphi%d",nj),Form("jet recopT:phi distribution jet %s",calgo[nj]),500,0,1000,36,-3.14,3.14);

    hjb[nj]    = new TH3F(Form("hjb%d",nj),Form("jet bkgd distribution jet  %s",calgo[nj]),18,-3.0,3.0,500,0,1000,300,0,600);

    hjrapt[nj]  = new TH2F(Form("hjrapt%d",nj),Form("jet reco:reco/raw  p_{T} distribution jet %s",calgo[nj]),500,0,1000,100,0,5);
    hetaphi[nj]  = new TH2F(Form("hetaphi%d",nj),Form("jet eta:phi distribution jet %s",calgo[nj]),36,-3.0,3.0,36,-3.0,3.0);
    htrackmax[nj]  = new TH2F(Form("htrackmax%d",nj),Form("recopT:trackmax distribution jet %s",calgo[nj]),500,0,1000,200,0,400);
    hphotonmax[nj] = new TH2F(Form("hphotonmax%d",nj),Form("recopT:photonmax distribution jet %s",calgo[nj]),500,0,1000,200,0,400);
    hchsum[nj] = new TH2F(Form("hchsum%d",nj),Form("recopT:chsum distribution jet %s",calgo[nj]),500,0,1000,200,0,400);
    hphsum[nj] = new TH2F(Form("hphsum%d",nj),Form("recopT:phsum distribution jet %s",calgo[nj]),500,0,1000,200,0,400);
    hnesum[nj] = new TH2F(Form("hnesum%d",nj),Form("recopT:nesum distribution jet %s",calgo[nj]),500,0,1000,200,0,400);
    hhcalsum[nj] = new TH2F(Form("hhcalsum%d",nj),Form("hcalsum distribution jet %s",calgo[nj]),500,0,1000,200,0,400);
    hecalsum[nj] = new TH2F(Form("hecalsum%d",nj),Form("ecalsum distribution jet %s",calgo[nj]),500,0,1000,200,0,400);

    hra_trackmax[nj]  = new TH2F(Form("hra_trackmax%d",nj),Form("ratio recopT:trackmax distribution jet %s",calgo[nj]),500,0,1000,100,0,2);
    hra_photonmax[nj] = new TH2F(Form("hra_photonmax%d",nj),Form("ratio recopT:photonmax distribution jet %s",calgo[nj]),500,0,1000,100,0,2);
    hra_chsum[nj] = new TH2F(Form("hra_chsum%d",nj),Form("ratio recopT:chsum distribution jet %s",calgo[nj]),500,0,1000,100,0,2);
    hra_phsum[nj] = new TH2F(Form("hra_phsum%d",nj),Form("ratio recopT:phsum distribution jet %s",calgo[nj]),500,0,1000,100,0,2);
    hra_nesum[nj] = new TH2F(Form("hra_nesum%d",nj),Form("ratio recopT:nesum distribution jet %s",calgo[nj]),500,0,1000,100,0,2);
    hra_hcalsum[nj] = new TH2F(Form("hra_hcalsum%d",nj),Form("ratio hcalsum distribution jet %s",calgo[nj]),500,0,1000,100,0,2);
    hra_ecalsum[nj] = new TH2F(Form("hra_ecalsum%d",nj),Form("ratio ecalsum distribution jet %s",calgo[nj]),500,0,1000,100,0,2);


    //! for jet pt
    hJetPtRun[nj] = new TH2F(Form("hJetPtRun%d",nj),Form("Jet pt run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,80,0,800);
    
    //! for jet pile up
    hJetPURun[nj] = new TH2F(Form("hJetPURun%d",nj),Form("Jet pile up run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,100,0,300);
    
    //! for jet bkgd
    hJetBkRun[nj] = new TH2F(Form("hJetBkRun%d",nj),Form("Jet bkgd run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,100,0,300);
    
    //! # of jets in a run
    hJetNoRun[nj] = new TH2F(Form("hJetNoRun%d",nj),Form("Total jets run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,500,0,10000);
    hJetNoRun[nj]->GetXaxis()->SetLabelSize(0.03); 
    
    
    for(int ix=1;ix<=hJetPtRun[nj]->GetNbinsX();ix++){
      if(Run[ix-1]!=-999)sprintf(mrun,"%d",Run[ix-1]);
      else sprintf(mrun,"%s","\0");
      
      hJetPtRun [nj]->GetXaxis()->SetBinLabel(ix,mrun);
      hJetPURun [nj]->GetXaxis()->SetBinLabel(ix,mrun);
      hJetBkRun [nj]->GetXaxis()->SetBinLabel(ix,mrun);
      hJetNoRun [nj]->GetXaxis()->SetBinLabel(ix,mrun);
      if(nj==0){
	hvzrun->GetXaxis()->SetBinLabel(ix,mrun);
	hAvHFRun->GetXaxis()->SetBinLabel(ix,mrun);
	hAvHFplusRun->GetXaxis()->SetBinLabel(ix,mrun);
	hAvHFminusRun->GetXaxis()->SetBinLabel(ix,mrun);
      }
    }
    std::cout<<Form("Initialized the histograms %s",calgo[nj]) <<std::endl;
  }
  std::cout<<"\t"<<std::endl;
  /////////////////////////////////////////////////////////////////////////////////////////
  //! Vertex z re-weighting function
  TF1 *fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  if(strcmp(ksp,"pp")==0)fVz->SetParameters(1.19756,0.0121871,-0.0048427,-3.3433e-05,4.97766e-06);
  else  fVz->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);

  double nmbj[knj];
  double sumjpt[knj];
  double sumjpu[knj];
  double sumjbk[knj];
  int istat[knj];
  
  double totjetpt[knj][maxruns];
  double totjets [knj][maxruns];
  
  for(int nj=0;nj<knj;nj++){
    for(int ir=0;ir<maxruns;ir++){
      totjetpt[nj][ir]=0;
      totjets[nj][ir]=0;
    }
  }

  hEvt->Fill(2,nentries);  
  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
  //for (Long64_t ievt=0; ievt<500;ievt++) {//! event loop
    //! load the hiForest event
    c->GetEntry(ievt);
    
    int RunNo = c->evt.run;
    int RunId=1;

    if(!isMC){
      RunId = GetRunIdx(RunNo);
      if(RunId<0){
	std::cout<<"Something is wrong with the LoadRuns, Current Run number not listed "<<std::endl;
	exit(0);
      }
    }
    if(ievt%10000==0 && !iCountRuns)std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<RunNo<<"\t HF : "<<c->evt.hiHF<<std::endl;
    //std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<RunNo<<"\t HF : "<<c->evt.hiHF<<std::endl;
    
    //! select good events
    bool goodEvent = false;
    //goodEvent = c->skim.pHBHENoiseFilter && (c->hlt.HLT_PAJet40_NoJetID_v1 || c->hlt.HLT_PAJet60_NoJetID_v1) && c->skim.pPAcollisionEventSelectionPA && fabs(c->evt.vz)<15.;
    if(isMC)goodEvent = fabs(c->evt.vz)<15.;
    else goodEvent = c->skim.pHBHENoiseFilter && (c->hlt.HLT_PAJet40_NoJetID_v1 || c->hlt.HLT_PAJet60_NoJetID_v1) && c->skim.pPAcollisionEventSelectionPA && fabs(c->evt.vz)<15.;

    if(!goodEvent)continue;

    //std::cout<<" ********** Event # " <<ievt<<"\t Vz: "<<c->evt.vz<<"\t HF : "<<c->evt.hiHF<<std::endl;

    //! Event variables
    int /*EvtNo=-999,*/ hiBin=-999, LumiBlock=-999, hiNpix =-999, hiNpixelTracks=-999,hiNtracks=-999;
    float vx=-999,vy=-999,vz=-999, hiHF=-999,hiHFplus=-999,hiHFminus=-999,hiHFhit=-999,hiHFhitplus=-999,
      hiHFhitminus=-999,hiZDC=-999,hiZDCplus=-999,hiZDCminus=-999;
    float hiET=-999,hiEE=-999,hiEB=-999,hiEEplus=-999,hiEEminus=-999;
    
    //EvtNo          = c->evt.evt;
    hiBin          = c->evt.hiBin;
    LumiBlock      = c->evt.lumi;
    hiNpix         = c->evt.hiNpix;
    hiNpixelTracks = c->evt.hiNpixelTracks;
    hiNtracks      = c->evt.hiNtracks;
    
    vx             = c->evt.vx;
    vy             = c->evt.vy;
    vz             = c->evt.vz;
    hiHF           = c->evt.hiHF;
    hiHFplus       = c->evt.hiHFplus;
    hiHFminus      = c->evt.hiHFminus;
    hiHFhit        = c->evt.hiHFhit;
    hiHFhitplus    = c->evt.hiHFhitPlus;
    hiHFhitminus   = c->evt.hiHFhitMinus;
    hiZDC          = c->evt.hiZDC;
    hiZDCplus      = c->evt.hiZDCplus;
    hiZDCminus     = c->evt.hiZDCminus;
    hiET           = c->evt.hiET;
    hiEE           = c->evt.hiEE;  
    hiEB           = c->evt.hiEB;  
    hiEEplus       = c->evt.hiEEplus;  
    hiEEminus      = c->evt.hiEEminus;

    //! vertex re-weighting function
    if(isMC){
      //wvz = 1./fVz->Eval(vz);
      wvz = 1.;
    }

    //if(ievt%10000==0 && !iCountRuns)std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<RunNo<<"\t HF : "<<hiHF<<std::endl;
    //std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<RunNo<<"\t HF : "<<hiHF<<std::endl;            

    //! Jet Algo loop starts here
    //! nj : 0 (icpu5calo) and 1 (akpu3pf)
    for(int nj=0;nj<knj;nj++){

      //! Initialization
      nmbj[nj]=0;
      sumjpt[nj]=0;
      sumjpu[nj]=0;
      sumjbk[nj]=0;
      istat[nj]=0;
      
      //! assign the jets  
      if(nj==0)mJets = &(c->ak3PF);
      else if(nj==1)mJets = &(c->ak4PF);
      else if(nj==2)mJets = &(c->ak5PF);

      else if(nj==3)mJets = &(c->akPu3PF);
      else if(nj==4)mJets = &(c->akPu4PF);
      else if(nj==5)mJets = &(c->akPu5PF);

      else if(nj==6)mJets = &(c->ak3Calo);
      else if(nj==7)mJets = &(c->ak4Calo);
      else if(nj==8)mJets = &(c->ak5Calo);

      else if(nj==9)mJets  = &(c->akPu3Calo);
      else if(nj==10)mJets = &(c->akPu4Calo);
      else if(nj==11)mJets = &(c->akPu5Calo);
      
      if(isMC){
        if(mJets->pthat > maxpthat)continue;
	else hpthat->Fill(mJets->pthat,wxs*wcen*wvz);
      }


      
      //if(nj==0)std::cout<<" ********** pThat cut  " <<mJets->pthat<<"\t nref : "<<mJets->nref<<std::endl;            
      //if(mJets->nref==0)continue;

      //! Loop over # of jets in a given algo
      for(int ij=0; ij<mJets->nref; ij++){
	
 	double recoeta = mJets->jteta[ij];
	bool goodJet=false;
	if(isMC)goodJet = mJets->jtpt[ij] > kptrecocut && mJets->trackMax[ij] > 4.0 && fabs(recoeta) < ketacut && mJets->subid[ij] == 0  && fabs(mJets->refdrjt[ij]) < kDelr && mJets->refpt[ij]>0.;
	else goodJet = mJets->jtpt[ij]>kptrecocut && mJets->trackMax[ij] > 4.0 && fabs(recoeta) < ketacut;
	if(!goodJet)continue;
	//if(fabs(mJets->refdrjt[ij])<kDelr && mJets->subid[ij]!=0)cout<<mJets->jtpt[ij]<<"\t"<<mJets->trackMax[ij]<<"\t"<<recoeta <<"\t"<<mJets->subid[ij]<<"\t"<<mJets->refdrjt[ij]<<"\t"<<mJets->refpt[ij]<<endl;
	
	double recophi = mJets->jtphi[ij];
	double recopt  = mJets->jtpt[ij];
	double rawpt   = mJets->rawpt[ij];
	double pu      = mJets->jtpu[ij];
	double jbkgd   = (mJets->photonSum[ij]+mJets->neutralSum[ij]+mJets->chargedSum[ij]) - mJets->rawpt[ij];
	double trkmax  = mJets->trackMax[ij];
	double phomax  = mJets->photonMax[ij];
	
	double chsum   = mJets->chargedSum[ij];
	double nesum   = mJets->neutralSum[ij];
	double phsum   = mJets->photonSum[ij];
	
	double hcalsum = mJets->hcalSum[ij];
	double ecalsum = mJets->ecalSum[ij];


	//cout<<recopt<<"\t"<<rawpt<<"\t"<<recoeta<<"\t"<<mJets->refdrjt[ij]<<"\t"<<mJets->subid[ij]<<endl;

	istat[nj]=1;
	hNjets [nj]->Fill(1.);
	hetaphi[nj]->Fill(recoeta,recophi,wxs*wcen*wvz);
	
	//! jets within different detector etabins
	hjrpt[nj]->Fill(recopt,rawpt,wxs*wcen*wvz);
	hjppt[nj]->Fill(recopt,pu,wxs*wcen*wvz);
	hjeta[nj]->Fill(recopt,recoeta,wxs*wcen*wvz);
	hjphi[nj]->Fill(recopt,recophi,wxs*wcen*wvz);
	hjbpt[nj]->Fill(recopt,jbkgd,wxs*wcen*wvz);
	hjrapt[nj]->Fill(recopt,rawpt/recopt,wxs*wcen*wvz);

	hjb[nj]->Fill(recoeta,recopt,jbkgd,wxs*wcen*wvz);

	hetaphi[nj]->Fill(recoeta,recophi,wxs*wcen*wvz);

 	htrackmax[nj]->Fill(recopt,trkmax,wxs*wcen*wvz);
	hphotonmax[nj]->Fill(recopt,phomax,wxs*wcen*wvz);
	hchsum[nj]->Fill(recopt,chsum,wxs*wcen*wvz);
	hnesum[nj]->Fill(recopt,nesum,wxs*wcen*wvz);
	hphsum[nj]->Fill(recopt,phsum,wxs*wcen*wvz);
	hhcalsum[nj]->Fill(recopt,hcalsum,wxs*wcen*wvz);
	hecalsum[nj]->Fill(recopt,ecalsum,wxs*wcen*wvz);


 	hra_trackmax[nj]->Fill(recopt,trkmax/recopt,wxs*wcen*wvz);
	hra_photonmax[nj]->Fill(recopt,phomax/recopt,wxs*wcen*wvz);
	hra_chsum[nj]->Fill(recopt,chsum/recopt,wxs*wcen*wvz);
	hra_nesum[nj]->Fill(recopt,nesum/recopt,wxs*wcen*wvz);
	hra_phsum[nj]->Fill(recopt,phsum/recopt,wxs*wcen*wvz);
	hra_hcalsum[nj]->Fill(recopt,hcalsum/recopt,wxs*wcen*wvz);
	hra_ecalsum[nj]->Fill(recopt,ecalsum/recopt,wxs*wcen*wvz);
	
	//! For a given event	  
	sumjpt[nj]+=recopt;
	sumjpu[nj]+=pu;
	sumjbk[nj]+=jbkgd;
	nmbj[nj]++;
	
	//! For a given run
	totjetpt[nj][RunId] += recopt;
	totjets [nj][RunId]++;
    }//! # jet loop
    //if(nj==1 && nmbj[nj]>0)std::cout<<Form(" \t \t \t mbjets %s : ",calgo[nj]) <<nmbj[nj]<<"\t sumjetpt : "<<sumjpt[nj]<<"\t <jtpt> "<<(sumjpt[nj]/nmbj[nj])<<std::endl;
      if(istat[nj]){
	//! # of jets distribution
	hJetPtRun [nj]->Fill(RunId,sumjpt[nj]/nmbj[nj]);
	hJetPURun [nj]->Fill(RunId,sumjpu[nj]/nmbj[nj]);
	hJetBkRun [nj]->Fill(RunId,sumjbk[nj]/nmbj[nj]);
      }
    }//! nj loop

    hvzrun->Fill(RunId,vz);
    //! Luminosity run-by-run
    hLumiBlock->Fill(LumiBlock);
    hEvtRun->Fill(RunId,vz);    
    hLumiRun->Fill(RunId,LumiBlock);
    hNTracksRun->Fill(RunId,hiNtracks);
    hNpixelTracksRun->Fill(RunId,hiNpixelTracks);
    if(hiHFhit>0)hAvHFRun->Fill(RunId,hiHF/hiHFhit);
    if(hiHFhitplus>0)hAvHFplusRun->Fill(RunId,hiHFplus/hiHFhitplus);
    if(hiHFhitminus>0)hAvHFminusRun->Fill(RunId,hiHFplus/hiHFhitminus);
    if(hiHFminus>0)hRatioHFplusminusRun->Fill(RunId,hiHFplus/hiHFminus);
    //! HF
    if(hiHFhit>0){
      hAvHF->Fill(hiHF/hiHFhit,wxs*wvz);
    }
    if(hiHFhitplus>0){
      hAvHFplus->Fill(hiHFplus/hiHFhitplus,wxs*wvz);
    }
    if(hiHFhitminus>0){
      hAvHFminus->Fill(hiHFplus/hiHFhitminus,wxs*wvz);
    }
    hHF->Fill(hiHF);
    hHFPlusMinus->Fill(hiHFplus,hiHFminus,wxs*wvz);
    
    //! ZDC
    hZDC->Fill(hiZDC,wxs*wvz);
    hZDCPlusMinus->Fill(hiZDCplus,hiZDCminus,wxs*wvz);
    
    //! EE
    hET->Fill(hiET,wxs*wvz);
    hEB->Fill(hiEB,wxs*wvz);
    hEE->Fill(hiEE,wxs*wvz);
    hEBEE->Fill(hiEB,hiEE,wxs*wvz);    
    hEEPlusMinus->Fill(hiEEplus,hiEEminus,wxs*wvz);    

    //! centralit bin
    hbin->Fill(hiBin,wxs*wvz);

    //! vertex
    hvz->Fill(vz,wxs*wvz);
    hvxvy->Fill(vx,vy,wxs*wvz);

    //! Pixel
    hNpix->Fill(hiNpix,wxs*wvz);
    hNpixelTracks->Fill(hiNpixelTracks,wxs*wvz);
    hNTracks->Fill(hiNtracks,wxs*wvz);
    hNpixelTracksNTracks->Fill(hiNpixelTracks,hiNtracks,wxs*wvz);

    //! Cross correlations
    hHFZDC->Fill(hiHF,hiZDC,wxs*wvz);
    hNpixelTracksZDC->Fill(hiNpixelTracks,hiZDC,wxs*wvz);
    hHFNpixelTracks->Fill(hiHF,hiNpixelTracks,wxs*wvz);

  }//! event loop ends
  
  if(!iCountRuns){
    
    std::cout<<std::endl;
    std::cout<<std::endl;
    
    for(int nj=0;nj<knj;nj++){
      if(!isMC){
	std::cout<<Form("For %s algo : ",calgo[nj])<<std::endl;
	for(int ir=0;ir<maxruns;ir++){
	  if(Run[ir]<0 || totjets[nj][ir]==0)continue;
	  std::cout<<"\t RunNo : "<<Run[ir]<<"\t Sumjetpt : "<<totjetpt[nj][ir]<<"\t Totjets : "<<totjets[nj][ir]<<"\t  <pT> : "<<totjetpt[nj][ir]/totjets[nj][ir]<<std::endl;
	  hJetNoRun[nj]->Fill(ir+1,totjets[nj][ir]);
	}
      }
    }
    std::cout<<std::endl;
    std::cout<<std::endl;
    for(int nj=0;nj<knj;nj++)std::cout<<"# of Events  : "<<calgo[nj]<<" \t : "<<hvz->Integral()<<"\t # of Jets : "<<hNjets[nj]->Integral()<<std::endl;
  }
  
  //! Write to output file
  fout->cd();
  fout->Write();
  fout->Close();

  //! Check
  timer.Stop();
  float  mbytes = 0.000001*nb;
  double rtime  = timer.RealTime();
  double ctime  = timer.CpuTime();

  std::cout<<"\t"<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<Form("You read %f Mbytes/RealTime seconds",mbytes/rtime)<<std::endl;
  std::cout<<Form("You read %f Mbytes/CpuTime  seconds",mbytes/ctime)<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;

  return 0;
}
int GetRunIdx(int irun)
{
  int runid=-1;
  int istat=1;
  int ir=0;
  while(istat>0){
    if(irun==Run[ir]){
      runid=ir;
      istat=-1;
    }
    ir++;
  }
  return runid;
}
void LoadRuns(const char *sp)
{
  //! 2010
  //Run[0]  = 161366;Run[1]  = 161396;Run[2]  = 161404;Run[3]  = 161439;Run[4]  = 161445;Run[5]  = 161450;Run[6]  = 161454;Run[7]  = 161473;Run[8]  = 161474;
  
  //! 2013
  Run [0]=211739;Run[1]=211740;Run[2]=211752;Run[3]=211760;Run[4]=211765;Run[5]=211792;Run[6]=211797;Run[7]=211812;Run[8]=211821;Run[9]=211822;Run[10]=211823;Run[11]=211831;
  std::cout<<Form("Run # loaded according to the JSON file for %s",sp)<<std::endl;
}
double delphi(double phi1, double phi2)
{
  double dphi = phi2 - phi1;
  dphi = fabs(atan2(sin(dphi),cos(dphi)));
  return dphi;
}

void GetRunNos(HiForest *c){

  /* for getting number of runs and runnumbers */
  c->hasHltTree=0;
  c->hasSkimTree=0;
  
  c->hasAk2CaloJetTree=0;
  c->hasAk3CaloJetTree=0;
  c->hasAk4CaloJetTree=0;
  c->hasAk5CaloJetTree=0;
  
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  c->hasAkPu5CaloJetTree=0;
  
  c->hasAk2PFJetTree=0;
  c->hasAk3PFJetTree=0;
  c->hasAk4PFJetTree=0;
  c->hasAk5PFJetTree=0;
  
  c->hasAkPu2PFJetTree=0;
  c->hasAkPu3PFJetTree=0;
  c->hasAkPu4PFJetTree=0;
  c->hasAkPu5PFJetTree=0;
  
  c->hasTrackTree=0;
  
  c->evtTree->SetBranchStatus("*",0,0);
  c->evtTree->SetBranchStatus("run",1,0);

  int oldrun=-1;
  std::vector <int>runs; 
  runs.clear();

  Long64_t nentries = c->GetEntries();
  std::cout<<"# of entries in TTree : "<<nentries<<std::endl;
  for (Long64_t ievt=0; ievt<nentries;ievt++) {
    c->GetEntry(ievt);
    int newrun=c->evt.run;
    if(newrun!=oldrun){
      int istat=1;
      for (std::vector<int>::iterator it = runs.begin(); it != runs.end(); ++it){
	if((*it)==newrun){
	  istat=0;
	  break;
	}
      }
      if(istat){
	runs.push_back(newrun);
	oldrun=newrun;
	//std::cout<<"iev : "<<ievt<<"\t oldrun : "<<oldrun<<"\t newrun "<<newrun<<std::endl;
      }
    }
  }

  std::sort(runs.begin(),runs.end(),ascend);

  int ic=0;
  for (std::vector<int>::iterator it = runs.begin(); it != runs.end(); ++it){
    std::cout<<"Run["<<ic<<"]="<<(*it)<<";";
    ic++;
  }
  std::cout<<std::endl;
}
bool ascend(int run1,int run2)
{
  return (run1<run2);
}
void GetXsection(const char *pthat,double &xsection,double &maxpthat){
  double xup=0;
  double xsub=0;

  if(strcmp(pthat,"15")==0){
    maxpthat=30;
    xup =2.034e-01;
    xsub=1.079e-02;
  }else if(strcmp(pthat,"30")==0){
    maxpthat=50;
    xup =1.079e-02;
    xsub=1.021e-03;
  }else if(strcmp(pthat,"50")==0){
    maxpthat=80;
    xup =1.021e-03;
    xsub=9.913e-05;
  }else if(strcmp(pthat,"80")==0){
    maxpthat=120;
    xup =9.913e-05;
    xsub=1.128e-05;
  }else if(strcmp(pthat,"120")==0){
    maxpthat=170;
    xup=1.128e-05;
    xsub=1.470e-06;
  }else if(strcmp(pthat,"170")==0){
    maxpthat=220;
    xup=1.470e-06;
    xsub=2.837e-07;
  }else if(strcmp(pthat,"220")==0){
    maxpthat=280;
    xup =2.837e-07;
    xsub=5.323e-08;
  }else if(strcmp(pthat,"280")==0){
    maxpthat=9999;
    xup =5.323e-08;
    xsub=0;
  }
  xsection = xup-xsub;
}
