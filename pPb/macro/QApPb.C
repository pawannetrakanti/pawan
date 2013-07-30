#include "../HiForest/V2/hiForest.h"
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
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
const double pi2=2*pi -1;

const int knj = 16;
const char *cjets[knj] = {"ak2PF","ak3PF","ak4PF","ak5PF",
			  "akPu2PF","akPu3PF","akPu4PF","akPu5PF",
			  "ak2Calo","ak3Calo","ak4Calo","ak5Calo",
			  "akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo"};
const char *calgo[knj] = {"ak2PF","ak3PF","ak4PF","ak5PF",
			  "akPu2PF","akPu3PF","akPu4PF","akPu5PF",
			  "ak2Calo","ak3Calo","ak4Calo","ak5Calo",
			  "akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo"};

const int ncen=5;
const char *ccent[ncen] = {"0-10%","10-30%","30-50%","50-70%","70-100%"};

const int nmult=6;
const int nvtx=4;
const int ketar=4;
const int maxruns=62;

int Run[maxruns];
double Lumi[maxruns];

void LoadRuns  (const char */*sp*/);
int GetRunIdx  (int /*irun*/);

int GetVtxBin  (float /*vz*/);
int GetEtaBin  (float /*eta*/);
int GetCentBin (int /*hiBin*/);
int GetMultBin (int /*ntracks*/);

void FindLeadSubLeadJets(Jets */*mJets*/, int */*ljet*/);

void LoadLib();
TStopwatch timer;
void QApPb(const int iRun=202792,const char *ksp="ppb")
{

  timer.Start();

  LoadLib();

  TString inname="";
  bool ispp=false;
  if(strcmp(ksp,"ppb")==0){
    //!  cmsPfn /store/caf/user/velicanu/PA2013_merged_HiForest/pPb_hiForest2_1_15_test.root //! to get the path from CAF
    
    //! pilot run
    //inname = "root://eoscms//eos/cms/store/caf/user/velicanu/PA2013_merged/PA2013_HiForest_Express_r0_pilot_minbias_v0.root"; // minbias
    //inname = "root://eoscms//eos/cms/store/caf/user/velicanu/PA2013_merged_HiForest/pPb_hiForest2_1_15_test.root"; //! not minbias
    //! Latest
    //inname   = "root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/PA2013_HiForest_Express_r0_pilot_minbias_v0.root"; //! minbias



    //! stable beams
    if(iRun<210760)inname=Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/PA2013_HiForest_Express_r%d_autoforest_v51.root",iRun); 
    else inname=Form("root://eoscms//eos/cms/store/group/phys_heavyions/velicanu/forest/PA2013_HiForest_Express_r%d_autoforest_v63.root",iRun); 

    ispp=false;
  }


  const float ketacut=2.;
  const double kptrecocut=30.;
  const bool iCountRuns=false;

  //! Load Lib
  //gSystem->Load("../HiForest/V2/hiForest_h.so");
  
  char mrun[6]={0};
  //! Load Run numbers 
  for(int ir=0;ir<maxruns;ir++){
    Run[ir]=-999;Lumi[ir]=-999;
  }
  if(!iCountRuns)LoadRuns(Form("%s",ksp));

  //! pt binning
  //!               0    1   2   3   4   5   6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30
  double ptbins[] ={30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,220.,240.,260.,280.,300.,325.,350.,375.,400.,450.,558.,638.,828};
  const int bins = sizeof(ptbins)/sizeof(double) - 1;

  //! Define the input file and HiForest
  //! CMSSW_5_3_3


  //! Define the input file and HiForest
  //! CMSSW_5_3_3
  HiForest *c = new HiForest(inname,Form("Forest%s",ksp),cPPb);

  //! Output file
  //! HIHighPt
  TFile *fout = new TFile(Form("output/test_%s_%d.root",ksp,iRun),"RECREATE");  

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


  //! Define event selection cut

  //! Don't want to loop over trees which is not used in the analysis
  //! event trees
  //c->hasHltTree=0;
  //c->hasSkimTree=0;
  //c->hasEvtTree=0;

  //! jet trees
  //c->hasAk1CaloJetTree=0;
  //c->hasAkPu1CaloJetTree=0;

  //! photon tree
  /*
  c->hasPhotonTree=0;
  c->hasTrackTree=0;
  c->hasPixTrackTree=0;
  c->hasTowerTree=0;
  c->hasHbheTree=0;
  c->hasEbTree=0;
  c->hasGenpTree=0;
  c->hasGenParticleTree=0;
  */

  //! added by pawan
  //! Select only the branches you want to use for the analysis
  //! This increases the speed for running

  //! For Hlt
  //c->hltTree->SetBranchStatus("*",0,0);
  //if(strcmp(ksp,"pp")==0) c->hltTree->SetBranchStatus("HLT_Jet60_v1",1,0);
  //else  c->hltTree->SetBranchStatus("HLT_HIJet80_v1",1,0);

  //! for Skim Tree
  //c->skimTree->SetBranchStatus("*",0,0);
  //c->skimTree->SetBranchStatus("pcollisionEventSelection",1,0);
  //c->skimTree->SetBranchStatus("pHBHENoiseFilter",1,0);

  //! Evt tree
  //c->evtTree->SetBranchStatus("*",0,0);
  //c->evtTree->SetBranchStatus("run",1,0);
  //c->evtTee->SetBranchStatus("lumi",1,0);
  //c->evtTree->SetBranchStatus("hiBin",1,0);//! centrality bin
  //c->evtTree->SetBranchStatus("hiHF",1,0);//! HF
  //c->evtTree->SetBranchStatus("hiHFplus",1,0); //! HF plus
  //c->evtTree->SetBranchStatus("hiHFminus",1,0);//! HF minus
  //c->evtTree->SetBranchStatus("hiZDC",1,0);//!ZDC

  std::cout<<"Loaded all tree variables "<<std::endl;

  //! For jets
  Jets *mJets=0;

  int scaFac=1;
  if(ispp)scaFac=5;

  //! Define histograms here
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TProfile::SetDefaultSumw2();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *hEvt = new TH1F("hEvt","# of events ",4,0,4);
  TH1F *hVz  = new TH1F("hVz","vertex z distribution",80,-20,20);
  TH2F *hEvtRun = new TH2F("hEvtRun","# of events in that run",maxruns,-0.5,maxruns-0.5,40,-40,40);
  TH1F *hLumiBlock = new TH1F("hLumiBlock","Luminosity block from hlt",scaFac*1500,0,scaFac*1500);
  TH2F *hLumiRun = new TH2F("hLumiRun","Luminosity Run-by-Run",maxruns,-0.5,maxruns-0.5,scaFac*100,0,scaFac*2000);
  hLumiRun->GetXaxis()->SetLabelSize(0.03);
  for(int ix=1;ix<=hLumiRun->GetNbinsX();ix++){
    if(Run[ix-1]!=-999)sprintf(mrun,"%d",Run[ix-1]);
    else sprintf(mrun,"%s","\0");
    hLumiRun->GetXaxis()->SetBinLabel(ix,mrun);
  }
  TH2F *hHFPlusMinus = new TH2F("hHFPlusMinus","HF plus:minus",200,0,800,200,0,800);

  TH1F *hAvHF     = new TH1F("hAvHF","HF E/hit distribution",100,0,0.1);
  TH1F *hAvHFplus = new TH1F("hAvHFplus","HFplus E/hit distribution",100,0,0.1);
  TH1F *hAvHFminus = new TH1F("hAvHFminus","HFminus E/hit distribution",100,0,0.2);

  TH2F *hAvHFRun     = new TH2F("hAvHFRun","HF E/hit distribution Run-by-Run",maxruns,-0.5,maxruns-0.5,100,0,0.1);
  TH2F *hAvHFplusRun = new TH2F("hAvHFplusRun","HFplus E/hit distribution Run-by-Run",maxruns,-0.5,maxruns-0.5,100,0,0.1);
  TH2F *hAvHFminusRun = new TH2F("hAvHFminusRun","HFminus E/hit distribution Run-by-Run",maxruns,-0.5,maxruns-0.5,100,0,0.2);

  TH2F *hZDCPlusMinus = new TH2F("hZDCPlusMinus","ZDC plus:minus",1000,0,45000,1000,0,45000);
  TH1F *hZDC     = new TH1F("hZDC","ZDC distribution",5000,0,1e+05);

  TH1F *hHF       = new TH1F("hHF","HF distribution",500,0,500);
  TH1F *hNTracks  = new TH1F("hNTracks","hiNTracks",500,0,500);
  TH2F *hNTracksHF = new TH2F("hNTracksHF","Ntrack vs HF ",500,0,500,500,0,500);
  TH2F *hHFNtrkO  = new TH2F("hHFNtrkO","HFvs NTracks distribution",500,0,500,500,0,500);

  TH1F *hBin        = new TH1F("hBin","Centrality bin",100,-0.5,100-0.5);
  TH2F *hBinHF      = new TH2F("hBinHF","HF vs hiBin",100,-0.5,100-0.5,500,0,500);
  TH2F *hBinNTracks = new TH2F("hBinNTracks","HF vs NTracks",100,-0.5,100-0.5,500,0,500);


  //! EE
  TH1F *hET = new TH1F("hET","hiET",200,0,1600);
  TH1F *hEB = new TH1F("hEB","hiEB",200,0,1200);
  TH1F *hEE = new TH1F("hEE","hiEE",500,0,2500);
  TH2F *hEBEE = new TH2F("hEBEE","hEBEE",200,0,1200,500,0,2500);    
  TH2F *hEEPlusMinus = new TH2F("hEEPlusMinus","hEEPlusEEMinus",200,0,1200,200,0,1200);    

  //! vertex
  TH1F *hvz   = new TH1F("hvz","hvz",40,-40,40);
  hvz->Sumw2();
  TH2F *hvxvy = new TH2F("hvxvy","vx vs vy",100,-0.2,0.2,100,-0.2,0.2);
  TH2F *hvzrun   = new TH2F("hvzrun","hvzrun",maxruns,-0.5,maxruns-0.5,40,-40,40);
  hvzrun->Sumw2();

  //! Pixel
  TH1F *hNpix = new TH1F("hNpix","hiNpix",200,0,3000);
  TH1F *hNpixelTracks = new TH1F("hNpixelTracks","hiNpixelTracks",100,0,300);
  TH2F *hNpixelTracksNTracks = new TH2F("hNpixelTracksNTracks","hNpixelTracksNTracks",500,0,500,100,0,500);

  //! Cross correlations
  TH2F *hHFZDC = new TH2F("hHFZDC","hHFZDC",600,0,6000,1000,0,1e+05);
  TH2F *hNpixelTracksZDC = new TH2F("hNpixelTracksZDC","hNpixelTracksZDC",200,0,3000,1000,0,1e+05);
  TH2F *hHFNpixelTracks  = new TH2F("hHFNpixelTracks","hHFNpixelTracks",1000,0,1000,500,0,500);



  TH1F *hHF_tr[knj];
  TH1F *hTotJets[knj];

  TH2F *hJetPtRun[knj];
  TH2F *hJetPtRun_Mult[knj][nmult];
  TH2F *hJetPURun[knj];
  TH2F *hMBJetsRun[knj];
  TH2F *hTotJetsRun[knj];
  TH2F *hTotJetsCentRun[knj][ncen];
  TH2F *hNevtRun[knj];
  TH1F *hTotEvtRun[knj];

  TH2F *hJetBkgdRun[knj];
  TH2F *hJetBkgdRun_Mult[knj][nmult];

  TH3F *hJetEtaPhiRun[knj];


  //! mb
  TH2F *hjetratiorjptmb[knj];
  TH1F *hjetrptmb[knj], *hjetjptmb[knj], *hjetjpumb[knj];
  TH1F *hjetetamb[knj], *hjetphimb[knj];
  TH2F *hjetptetamb[knj], *hjetptphimb[knj], *hjetetaphimb[knj];
  TH2F *hjetptpumb[knj];
  TH2F *hjetptbkgdmb[knj];
  TH2F *hjetptbkgd_etabmb[knj][ketar];
  TH2F *hJetEnergyScalemb[knj];

  //! vetex bin
  TH1F *hjetrptvz[knj][nvtx], *hjetjptvz[knj][nvtx], *hjetjpuvz[knj][nvtx];
  TH1F *hjetetavz[knj][nvtx], *hjetphivz[knj][nvtx];
  TH2F *hjetptetavz[knj][nvtx], *hjetptphivz[knj][nvtx], *hjetetaphivz[knj][nvtx];
  TH2F *hjetptpuvz[knj][nvtx];

  //! ncen bin
  TH1F *hfrpt[knj][ncen], *hfjpt[knj][ncen], *hfpupt[knj][ncen];
  TH1F *hjeteta  [knj][ncen], *hjetphi[knj][ncen];
  TH2F *hjetpteta[knj][ncen], *hjetptphi[knj][ncen], *hjetetaphi[knj][ncen];
  TH2F *hjetptpu[knj][ncen];
  TH1F *hjetbkgd[knj][ncen];
  TH2F *hjetptbkgd[knj][ncen];


  //! mult bin
  TH1F *hfrpt_mult[knj][nmult], *hfjpt_mult[knj][nmult], *hfpupt_mult[knj][nmult];
  TH1F *hjeteta_mult  [knj][nmult], *hjetphi_mult[knj][nmult];
  TH2F *hjetpteta_mult[knj][nmult], *hjetptphi_mult[knj][nmult], *hjetetaphi_mult[knj][nmult];
  TH2F *hjetptpu_mult[knj][nmult];
  TH1F *hjetbkgd_mult[knj][nmult];
  TH2F *hjetptbkgd_mult[knj][nmult];
  TH2F *hjetptbkgd_etab_mult[knj][nmult][ketar];

  //! eta dependence
  TH1F *hjetrptmb_etab[knj][ketar];
  TH1F *hjetjptmb_etab[knj][ketar];
  TH1F *hjetjpumb_etab[knj][ketar];
  TH2F *hjetptpumb_etab[knj][ketar];

  TH1F *hNevtmb_etab[knj];

  TH1F *hjetrpt_etab[knj][ncen][ketar];
  TH1F *hjetjpt_etab[knj][ncen][ketar];
  TH1F *hjetjpu_etab[knj][ncen][ketar];
  TH2F *hjetptpu_etab[knj][ncen][ketar];
  TH1F *hNevt_etab[knj][ncen];



  TH1F *hNref[knj];
  TH1F *hNjets[knj][ncen], *hNjets_dist[knj][ncen];
  TH1F *hMBJets[knj];

  TH1F *hNevt      [knj][ncen];
  TH1F *hNevtvz    [knj][nvtx];
  TH1F *hNevt_notr [knj][ncen];

  //TH2F *hJetEnergyScale[knj][ncen];

  TH2F *hJetEnergyScale_mult[knj][nmult];

  TH1F *hNtrks_Out;
  TH1F *hNtrks_In;

  TH2F *hNtrksNTracks_In;
  TH2F *hNtrksHF_In;
  TH2F *hNtrksNTracks_Out;
  TH2F *hNtrksHF_Out;

  TH2F *hNtrks_InOut;

  hNtrks_Out = new TH1F("hNtrks_Out","# of tracks outside jet cone",400,0,400);
  hNtrks_In  = new TH1F("hNtrks_In","# of tracks inside jet cone",400,0,400);

  hNtrksNTracks_In = new TH2F("hNtrksNTracks_In","# of tracks inside jet cone",400,0,400,400,0,400);
  hNtrksHF_In = new TH2F("hNtrksHF_In","# of tracks inside jet cone",400,0,400,400,0,400);
  hNtrksNTracks_Out = new TH2F("hNtrksNTracks_Out","# of tracks outside jet cone",400,0,400,400,0,400);
  hNtrksHF_Out = new TH2F("hNtrksHF_Out","# of tracks outside jet cone",400,0,400,400,0,400);
  
  hNtrks_InOut = new TH2F("hNtrks_InOut","# of tracks inside-outside jet cone",400,0,400,400,0,400);

  for(int nj=0;nj<knj;nj++){
    //std::cout<<Form("Histograms for algo %s",calgo[nj])<<std::endl;


    //! Define the histograms for jets
    hTotJets[nj] = new TH1F(Form("hTotJets%d",nj),Form("# of jets w/o cuts with %s",calgo[nj]),4,0,4);
    hMBJets[nj] = new TH1F(Form("hMBJets%d",nj),Form("MB distribution of jets in %s",calgo[nj]),200,-0.5,200-0.5);

    hMBJetsRun[nj] = new TH2F(Form("hMBJetsRun%d",nj),Form("MB distribution of jets run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,100,-0.5,100-0.5);
    hMBJetsRun[nj]->GetXaxis()->SetLabelSize(0.03); 

    hJetEtaPhiRun[nj] = new TH3F(Form("hJetEtaPhiRun%d",nj),Form("Jet eta-phi run-by-run for %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,72,-2.0,2.0,72,-pi,pi);

    //! for jet pt
    hJetPtRun[nj] = new TH2F(Form("hJetPtRun%d",nj),Form("Jet pt run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,500,0,1000);


    //! for jet pile up
    hJetPURun[nj] = new TH2F(Form("hJetPURun%d",nj),Form("Jet pile up run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,100,0,300);
    hJetBkgdRun[nj] = new TH2F(Form("hJetBkgdRun%d",nj),Form("Jet bkgd run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,100,0,300);

    hTotJetsRun[nj] = new TH2F(Form("hTotJetsRun%d",nj),Form("Total jets run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,500,-0.5,10000-0.5);
    hTotJetsRun[nj]->GetXaxis()->SetLabelSize(0.03); 

    hNevtRun[nj] = new TH2F(Form("hNevtRun%d",nj),Form("Total # of events run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,3,-0.5,3-0.5);
    hNevtRun[nj]->GetXaxis()->SetLabelSize(0.03); 

    hTotEvtRun[nj] = new TH1F(Form("hTotEvtRun%d",nj),Form("Total # of events run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5);
    hTotEvtRun[nj]->GetXaxis()->SetLabelSize(0.03); 

    for(int ic=0;ic<ncen;ic++){
      hTotJetsCentRun[nj][ic] = new TH2F(Form("hTotJetsCentRun%d_%d",nj,ic),Form("Total jets run-by-run in %s %s",calgo[nj],ccent[ic]),maxruns,-0.5,maxruns-0.5,500,-0.5,10000-0.5);
      hTotJetsCentRun[nj][ic]->GetXaxis()->SetLabelSize(0.03); 
    }

    for(int im=0;im<nmult;im++){
      hJetPtRun_Mult[nj][im] = new TH2F(Form("hJetPtRun_Mult%d_%d",nj,im),Form("Jet pt run-by-run in  multiplicity %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,500,0,1000);
      hJetBkgdRun_Mult[nj][im] = new TH2F(Form("hJetBkgdRun_Mult%d_%d",nj,im),Form("Jet bkgd run-by-run in multiplicity %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,100,0,300);
    }

    for(int ix=1;ix<=hMBJetsRun[nj]->GetNbinsX();ix++){
      if(Run[ix-1]!=-999)sprintf(mrun,"%d",Run[ix-1]);
      else sprintf(mrun,"%s","\0");

      hMBJetsRun[nj]->GetXaxis()->SetBinLabel(ix,mrun);
      hJetPtRun [nj]->GetXaxis()->SetBinLabel(ix,mrun);

      hJetPURun [nj]->GetXaxis()->SetBinLabel(ix,mrun);
      hTotJetsRun[nj]->GetXaxis()->SetBinLabel(ix,mrun);

      for(int im=0;im<nmult;im++){
	hJetPtRun_Mult [nj][im]->GetXaxis()->SetBinLabel(ix,mrun);
	hJetBkgdRun_Mult [nj][im]->GetXaxis()->SetBinLabel(ix,mrun);
      }
      for(int ic=0;ic<ncen;ic++){
	hTotJetsCentRun[nj][ic]->GetXaxis()->SetBinLabel(ix,mrun);
      }
    }

    hJetEnergyScalemb[nj] = new TH2F(Form("hJetEnergyScalemb%d",nj),Form("hJetEnergyScalemb%d",nj),500,0,2000,100,-1.00,1.00);

    hNref[nj] = new TH1F(Form("hNRef%d",nj),Form("# of jets in %s",calgo[nj]),40,0,100);
    hHF_tr[nj]      = new TH1F(Form("hHF_tr%d",nj),Form("HF triggered events (Jet p_{T} > 85 GeV/c distribution) %s",calgo[nj]),600,0,6000);
    
    hjetetamb[nj] = new TH1F(Form("hjetetamb%d",nj),Form("jet eta distribution jet %s",calgo[nj]),72,-2.0,2.0);
    hjetphimb[nj] = new TH1F(Form("hjetphimb%d",nj),Form("jet phi distribution jet %s",calgo[nj]),72,-pi,pi);
    hjetrptmb[nj] = new TH1F(Form("hjetrptmb%d",nj),Form("jet(raw pt) p_{T} distribution jet %s",calgo[nj]),500,0,1000);
    hjetratiorjptmb[nj] = new TH2F(Form("hjetratiorjptmb%d",nj),Form("corr pT / jet(raw pt) p_{T} distribution jet %s",calgo[nj]),500,0,1000,50,0,10);

    hjetjptmb[nj] = new TH1F(Form("hjetjptmb%d",nj),Form("jet(jt pt) p_{T} distribution jet %s",calgo[nj]),500,0,1000);
    hjetjpumb[nj] = new TH1F(Form("hjetjpumb%d",nj),Form("jet(j pu) p_{T} distribution jet %s",calgo[nj]),100,0,300);
    hjetptpumb[nj] = new TH2F(Form("hjetptpumb%d",nj),Form("jet(pt:pu) distribution jet %s",calgo[nj]),bins,ptbins,100,0,300);
    
    hjetetaphimb[nj] = new TH2F(Form("hjetetaphimb%d",nj),Form("jet eta-phi distribution jet mb %s",calgo[nj]),72,-2.0,2.0,72,-pi,pi);
    hjetptetamb[nj] = new TH2F(Form("hjetptetamb%d",nj),Form("jet pt-eta distribution jet mb %s",calgo[nj]),500,0,1000,72,-2.0,2.0);
    hjetptphimb[nj] = new TH2F(Form("hjetptphimb%d",nj),Form("jet pt-phi distribution jet mb %s",calgo[nj]),500,0,1000,72,-pi,pi);

    hjetptbkgdmb[nj] = new TH2F(Form("hjetptbkgdmb%d",nj),Form("jet(pt:bkgd) distribution jet mb %s",calgo[nj]),500,0,1000,100,0,300);
    hNevtmb_etab[nj] = new TH1F(Form("hNevtmb_etab%d",nj),Form("# of events mb jet %s",calgo[nj]),40,-40,40);

    for(int ie=0;ie<ketar;ie++){
      hjetrptmb_etab[nj][ie] = new TH1F(Form("hjetrptmb_etab%d_%d",nj,ie),Form("jet(raw pt) p_{T} distribution jet %s etabin %d",calgo[nj],ie),bins,ptbins);
      hjetjptmb_etab[nj][ie] = new TH1F(Form("hjetjptmb_etab%d_%d",nj,ie),Form("jet(jt pt) p_{T} distribution jet %s etabin %d",calgo[nj],ie),bins,ptbins);
      hjetjpumb_etab[nj][ie] = new TH1F(Form("hjetjpumb_etab%d_%d",nj,ie),Form("jet(j pu) p_{T} distribution jet %s etabin %d",calgo[nj],ie),100,0,300);
      hjetptpumb_etab[nj][ie] = new TH2F(Form("hjetptpumb_etab%d_%d",nj,ie),Form("jet(pt:pu) distribution jet %s etabin %d",calgo[nj],ie),bins,ptbins,100,0,300);
      
      hjetptbkgd_etabmb[nj][ie] = new TH2F(Form("hjetptbkgd_etabmb%d_%d",nj,ie),Form("jet(pt:bkgd) distribution jet mb %s etabin %d",calgo[nj],ie),500,0,1000,100,0,300);
    }

    for(int iv=0;iv<nvtx;iv++){
      hNevtvz[nj][iv] = new TH1F(Form("hNevtvz%d_%d",nj,iv),Form("# of events cent vz %d %s",iv,calgo[nj]),40,-40,40);
      hjetetavz[nj][iv] = new TH1F(Form("hjetetavz%d_%d",nj,iv),Form("jet eta distribution jet vz %s",calgo[nj]),72,-2.0,2.0);
      hjetphivz[nj][iv] = new TH1F(Form("hjetphivz%d_%d",nj,iv),Form("jet phi distribution jet vz %s",calgo[nj]),72,-pi,pi);
      hjetrptvz[nj][iv] = new TH1F(Form("hjetrptvz%d_%d",nj,iv),Form("jet(raw pt) p_{T} distribution jet vz %s",calgo[nj]),bins,ptbins);
      hjetjptvz[nj][iv] = new TH1F(Form("hjetjptvz%d_%d",nj,iv),Form("jet(jt pt) p_{T} distribution jet vz %s",calgo[nj]),bins,ptbins);
      hjetjpuvz[nj][iv] = new TH1F(Form("hjetjpuvz%d_%d",nj,iv),Form("jet(pu) p_{T} distribution jet vz %s",calgo[nj]),100,0,300);
      hjetptpuvz[nj][iv] = new TH2F(Form("hjetptpuvz%d_%d",nj,iv),Form("jet(pt:pu)  distribution jet vz %s",calgo[nj]),bins,ptbins,100,0,300);
      
      hjetetaphivz[nj][iv] = new TH2F(Form("hjetetaphivz%d_%d",nj,iv),Form("jet eta-phi distribution jet vz %s",calgo[nj]),72,-2.0,2.0,72,-pi,pi);
      hjetptetavz[nj][iv] = new TH2F(Form("hjetptetavz%d_%d",nj,iv),Form("jet pt-eta distribution jet vz %s",calgo[nj]),bins,ptbins,72,-3.0,3.0);
      hjetptphivz[nj][iv] = new TH2F(Form("hjetptphivz%d_%d",nj,iv),Form("jet pt-phi distribution jet vz %s",calgo[nj]),bins,ptbins,72,-pi,pi);
    }

    for(int icen=0;icen<ncen;icen++){
      //hJetEnergyScale[nj][icen] = new TH2F(Form("hJetEnergyScale%d_%d",nj,icen),Form("hJetEnergyScale%d_%d",nj,icen),bins,ptbins,50,-1.00,1.00);

      hjeteta[nj][icen] = new TH1F(Form("hjeteta%d_%d",nj,icen),Form("jet eta distribution jet centb %d %s",icen,calgo[nj]),72,-2.0,2.0);
      hjetphi[nj][icen] = new TH1F(Form("hjetphi%d_%d",nj,icen),Form("jet phi distribution jet centb %d %s",icen,calgo[nj]),72,-pi,pi);

      hfrpt[nj][icen]    = new TH1F(Form("hfrpt%d_%d",nj,icen),Form("jet(raw pt) p_{T} distribution jet centb %d %s",icen,calgo[nj]),500,0.,1000);
      hfjpt[nj][icen]    = new TH1F(Form("hfjpt%d_%d",nj,icen),Form("jet(corrected pt) p_{T} distribution jet centb %d %s",icen,calgo[nj]),500,0.,1000);
      hfpupt[nj][icen]   = new TH1F(Form("hfpupt%d_%d",nj,icen),Form("jet(pilup pt) p_{T} distribution jet centb %d %s",icen,calgo[nj]),250,0.,100);
      hjetptpu[nj][icen] = new TH2F(Form("hjetptpu%d_%d",nj,icen),Form("jet(pt:pu) distribution jet centb %d %s",icen,calgo[nj]),500,0,1000,100,0,300);

      hjetbkgd  [nj][icen] = new TH1F(Form("hjetbkgd%d_%d",nj,icen),Form("jet(pu) p_{T} distribution jet centb %d %s",icen,calgo[nj]),100,0,300);
      hjetptbkgd[nj][icen] = new TH2F(Form("hjetptbkgd%d_%d",nj,icen),Form("jet(pt:bkgd) distribution jet centb %d %s",icen,calgo[nj]),500,0,1000,100,0,300);

      hjetetaphi[nj][icen] = new TH2F(Form("hjetetaphi%d_%d",nj,icen),Form("jet eta-phi distribution jet centb %d %s",icen,calgo[nj]),72,-2.0,2.0,72,-pi,pi);
      hjetpteta[nj][icen] = new TH2F(Form("hjetpteta%d_%d",nj,icen),Form("jet pt-eta distribution jet centb %d %s",icen,calgo[nj]),500,0,1000,72,-2.0,2.0);
      hjetptphi[nj][icen] = new TH2F(Form("hjetptphi%d_%d",nj,icen),Form("jet pt-phi distribution jet centb %d %s",icen,calgo[nj]),500,0,1000,72,-pi,pi);

      hNevt_etab[nj][icen] = new TH1F(Form("hNevt_etab%d_%d",nj,icen),Form("# of events cent %d jet %s",icen,calgo[nj]),15,-15,15);

      for(int ie=0;ie<ketar;ie++){
	hjetrpt_etab[nj][icen][ie] = new TH1F(Form("hjetrpt_etab%d_%d_%d",nj,icen,ie),Form("jet(raw pt) p_{T} distribution centb %d jet %s etabin %d",icen,calgo[nj],ie),500,0,1000);
	hjetjpt_etab[nj][icen][ie] = new TH1F(Form("hjetjpt_etab%d_%d_%d",nj,icen,ie),Form("jet(jt pt) p_{T} distribution centb %d jet %s etabin %d",icen,calgo[nj],ie),500,0,1000);
	hjetjpu_etab[nj][icen][ie] = new TH1F(Form("hjetjpu_etab%d_%d_%d",nj,icen,ie),Form("jet(pu) p_{T} distribution centb %d jet %s etabin %d",icen,calgo[nj],ie),100,0,300);
	hjetptpu_etab[nj][icen][ie] = new TH2F(Form("hjetptpu_etab%d_%d_%d",nj,icen,ie),Form("jet(pt:pu) distribution centb %d jet %s etabin %d",icen,calgo[nj],ie),500,0,1000,100,0,300);

      }

      hNevt[nj][icen] = new TH1F(Form("hNevt%d_%d",nj,icen),Form("# of events cent %d %s",icen,calgo[nj]),15,-15,15);

      hNevt_notr[nj][icen] = new TH1F(Form("hNevt_notr%d_%d",nj,icen),Form("# of events w/o trigger conditioncent %d %d",nj,icen),150,-15,15);

      hNjets[nj][icen] = new TH1F(Form("hNjets%d_%d",nj,icen),Form("# of cent %d jets %s",icen,calgo[nj]),3,0.0,3.0);
      
      hNjets_dist[nj][icen] = new TH1F(Form("hNjets_dist%d_%d",nj,icen),Form("distribution jets in cent %d %s",icen,calgo[nj]),100,0,100);
    }

    //! nmult bin
    for(int imult=0;imult<nmult;imult++){
      hJetEnergyScale_mult[nj][imult] = new TH2F(Form("hJetEnergyScale_mult%d_%d",nj,imult),Form("hJetEnergyScale%d_%d",nj,imult),500,0,2000,100,-1.00,1.00);
      hjeteta_mult[nj][imult] = new TH1F(Form("hjeteta_mult%d_%d",nj,imult),Form("jet eta distribution jet centb %d %s",imult,calgo[nj]),72,-2.0,2.0);
      hjetphi_mult[nj][imult] = new TH1F(Form("hjetphi_mult%d_%d",nj,imult),Form("jet phi distribution jet centb %d %s",imult,calgo[nj]),72,-pi,pi);
      
      hfrpt_mult[nj][imult]    = new TH1F(Form("hfrpt_mult%d_%d",nj,imult),Form("jet(raw pt) p_{T} distribution jet centb %d %s",imult,calgo[nj]),500,0.,1000);
      hfjpt_mult[nj][imult]    = new TH1F(Form("hfjpt_mult%d_%d",nj,imult),Form("jet(corrected pt) p_{T} distribution jet centb %d %s",imult,calgo[nj]),500,0.,1000);
      hfpupt_mult[nj][imult]   = new TH1F(Form("hfpupt_mult%d_%d",nj,imult),Form("jet(pilup pt) p_{T} distribution jet centb %d %s",imult,calgo[nj]),250,0.,100);
      hjetptpu_mult[nj][imult] = new TH2F(Form("hjetptpu_mult%d_%d",nj,imult),Form("jet(pt:pu) distribution jet centb %d %s",imult,calgo[nj]),500,0,1000,100,0,300);

      hjetbkgd_mult[nj][imult]   = new TH1F(Form("hjetbkgd_mult%d_%d",nj,imult),Form("jet(pu) p_{T} distribution jet centb %d %s",imult,calgo[nj]),100,0,300);
      hjetptbkgd_mult[nj][imult] = new TH2F(Form("hjetptbkgd_mult%d_%d",nj,imult),Form("jet(pt:bkgd) distribution jet centb %d %s",imult,calgo[nj]),500,0,1000,100,0,300);

      hjetetaphi_mult[nj][imult] = new TH2F(Form("hjetetaphi_mult%d_%d",nj,imult),Form("jet eta-phi distribution jet centb %d %s",imult,calgo[nj]),72,-2.0,2.0,72,-pi,pi);
      hjetpteta_mult[nj][imult]  = new TH2F(Form("hjetpteta_mult%d_%d",nj,imult),Form("jet pt-eta distribution jet centb %d %s",imult,calgo[nj]),500,0,1000,72,-2.0,2.0);
      hjetptphi_mult[nj][imult]  = new TH2F(Form("hjetptphi_mult%d_%d",nj,imult),Form("jet pt-phi distribution jet centb %d %s",imult,calgo[nj]),500,0,1000,72,-pi,pi);
      for(int ie=0;ie<ketar;ie++){
	hjetptbkgd_etab_mult[nj][imult][ie] = new TH2F(Form("hjetptbkgd_etab_mult%d_%d_%d",nj,imult,ie),Form("jet(pt:bkgd) distribution jet etabin %d centb %d %s",ie,imult,calgo[nj]),500,0,1000,100,0,300);
      }
    }
    std::cout<<Form("Initialized the histograms %s",calgo[nj]) <<std::endl;
  }
  std::cout<<"\t"<<std::endl;
  /////////////////////////////////////////////////////////////////////////////////////////

  /* for getting number of runs and runnumbers */
  int runc=0;
  int newrun=0;
  int oldrun=-1;
  int kstat[maxruns]={0};

  
  float njets[knj][ncen];
  float nmbj[knj];
  float sumjpt[knj];
  float sumjpu[knj];
  float sumjbkg[knj];
  int istat[knj];
  int stat_etab[knj];

  double totjetpt[knj][maxruns];
  double totjets [knj][maxruns];
  double totjetsc[knj][maxruns][ncen];
  int totnevt[knj][maxruns];
  
  for(int nj=0;nj<knj;nj++){
    for(int ir=0;ir<maxruns;ir++){
      totjetpt[nj][ir]=0;
      totjets[nj][ir]=0;
      totnevt[nj][ir]=0;
      for(int ic=0;ic<ncen;ic++){
	totjetsc[nj][ir][ic]=0;
      }
    }
  }

  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;
  hEvt->Fill(2,nentries);
  
  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
  //for (Long64_t ievt=0; ievt<2000;ievt++) {//! event loop
    //! load the hiForest event
    c->GetEntry(ievt);

    //! select good events
    bool goodEvent = false;
    //goodEvent = c->selectEvent() && fabs(c->evt.vz)<15.; 

    //HLT_PAJet100_NoJetID_v1 (105566) 127 
    //HLT_PAJet120_NoJetID_v1 (105567) 127 
    //HLT_PAJet20_NoJetID_v1  (105562) 127 
    //HLT_PAJet40ETM30_v1     (105583) 127 
    //HLT_PAJet40_NoJetID_v1  (105563) 127 
    //HLT_PAJet60ETM30_v1     (105584) 127 
    //HLT_PAJet60_NoJetID_v1  (105564) 127 
    //HLT_PAJet80_NoJetID_v1  (105565) 127 

    
    bool selRun=false;
    if(c->evt.run==210498)selRun = c->evt.lumi>=121 && c->evt.lumi<=351;
    if(c->evt.run==210534)selRun = c->evt.lumi>=23  && c->evt.lumi<=348;
    if(c->evt.run==210614)selRun = c->evt.lumi>=101 && c->evt.lumi<=1652;
    if(c->evt.run==210634)selRun = c->evt.lumi>=97  && c->evt.lumi<=798;    
    if(c->evt.run==210635)selRun = c->evt.lumi>=1   && c->evt.lumi<=935;    
    if(c->evt.run==210638)selRun = c->evt.lumi>=1   && c->evt.lumi<=186;    
    if(c->evt.run==210658)selRun = c->evt.lumi>=100 && c->evt.lumi<=1544;    
    if(c->evt.run==210676)selRun = c->evt.lumi>=690 && c->evt.lumi<=908;    
    if(c->evt.run==210737)selRun = c->evt.lumi>=89  && c->evt.lumi<=284;    
    if(c->evt.run==210738)selRun = c->evt.lumi>=1   && c->evt.lumi<=696;    
    if(c->evt.run==210759)selRun = c->evt.lumi>=12  && c->evt.lumi<=1011;    
    if(c->evt.run==210818)selRun = c->evt.lumi>=99  && c->evt.lumi<=910;    
    if(c->evt.run==210837)selRun = c->evt.lumi>=80  && c->evt.lumi<=1054;    
    if(c->evt.run==210855)selRun = c->evt.lumi>=81  && c->evt.lumi<=1033;    
    if(c->evt.run==210885)selRun = c->evt.lumi>=1   && c->evt.lumi<=906;    
    if(c->evt.run==210906)selRun = c->evt.lumi>=85  && c->evt.lumi<=684;    
    if(c->evt.run==210909)selRun = c->evt.lumi>=1   && c->evt.lumi<=339;    
    if(c->evt.run==210986)selRun = (c->evt.lumi>=133&& c->evt.lumi<=512) || (c->evt.lumi>=806 && c->evt.lumi<=1169);    
    if(c->evt.run==210998)selRun = c->evt.lumi>=129 && c->evt.lumi<=537;    
    if(c->evt.run==211000)selRun = c->evt.lumi>=19  && c->evt.lumi<=434;    
    if(c->evt.run==211001)selRun = c->evt.lumi>=1   && c->evt.lumi<=177;    
    if(c->evt.run==211032)selRun = c->evt.lumi>=175   && c->evt.lumi<=1179;    
    if(c->evt.run==211256)selRun = c->evt.lumi>=80   && c->evt.lumi<=756;    

    //! Reverse beam
    if(c->evt.run==211313)selRun = c->evt.lumi>=55   && c->evt.lumi<=756;    
    if(c->evt.run==211328)selRun = c->evt.lumi>=56   && c->evt.lumi<=734;    
    if(c->evt.run==211347)selRun = c->evt.lumi>=64   && c->evt.lumi<=713;    
    if(c->evt.run==211354)selRun = c->evt.lumi>=99   && c->evt.lumi<=772;    


    //! pilot run
    //bool evSel = c->skim.pHBHENoiseFilter && c->hlt.HLT_PAZeroBiasPixel_SingleTrack_v1 && c->skim.phfPosFilter1 && c->skim.phfNegFilter1 && c->skim.phltPixelClusterShapeFilter && c->skim.pprimaryvertexFilter;
    //selRun=true;
    //! events with Single vertex
    //bool singleVertex = true;
    
    //cout<<"event vz : "<<c->evt.vz<<"\t "<<endl;
    
    bool evSel = c->skim.pHBHENoiseFilter  && c->skim.pPAcollisionEventSelectionPA;
    //! events with Single vertex
    bool singleVertex = c->track.nVtx==1;

    goodEvent  = fabs(c->evt.vz)<15. && evSel && selRun && singleVertex;
    if(!goodEvent)continue;

    //! Event variables
    int RunNo=-999, /*EvtNo=-999,*/ hiBin=-999, LumiBlock=-999, hiNpix =-999, hiNpixelTracks=-999,hiNtracks=-999;
    float vx=-999,vy=-999,vz=-999, hiHF=-999,hiHFplus=-999,hiHFminus=-999,hiHFhit=-999,hiHFhitplus=-999,
      hiHFhitminus=-999,hiZDC=-999,hiZDCplus=-999,hiZDCminus=-999;
    float hiET=-999,hiEE=-999,hiEB=-999,hiEEplus=-999,hiEEminus=-999;
    
    RunNo          = c->evt.run;
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


    int RunId=0;
    /* For getting the number of runs and runnumbers */
    if(iCountRuns){
      newrun=RunNo;
      if(newrun!=oldrun){
	bool mstat=true;
	for(int ir=0;ir<runc;ir++){
	  if(kstat[ir]==newrun)mstat=false;
	}
	if(mstat){
	  oldrun=newrun;
	  kstat[runc]=newrun;
	  runc++;
	  std::cout<<"Run # : " <<newrun<<"\t # of runs : "<<runc<<std::endl;
	}
      }
    }else{
      RunId = GetRunIdx(RunNo);
      if(RunId<0){
	std::cout<<"Something is wrong with the LoadRuns, Current Run number not listed "<<std::endl;
	exit(0);
      }
    }

    //! Centrality bin
    //! 5 or 7 bins
    int centb=-1;
    if(ispp) centb=ncen-1; //pp
    else{ //! PbPb
      centb=GetCentBin(hiBin);
    }

    if(centb<0 || centb>=ncen)continue;
    //centb=0;

    //! Vertex z bin
    //! Only two bins v<0 : 0 and vz>0 : 1
    int vbin = GetVtxBin(vz);
    //if(vbin<0)continue;

    //! Multiplicity bins
    int imult = GetMultBin(hiNtracks);
    
    if(ievt%10000==0 && !iCountRuns)std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<RunNo<<"\t HF : "<<hiHF<<std::endl;




    //! Jet Algo loop starts here
    //! nj : 0 (icpu5calo) and 1 (akpu3pf)
    Int_t iFirst=0;
    Int_t iSecond=0;
    Int_t ntrks_in=0;
    Int_t ntrks_out=0;
    for(int nj=0;nj<knj;nj++){

      //! Initialization
      nmbj[nj]=0;
      sumjpt[nj]=0;
      sumjpu[nj]=0;
      sumjbkg[nj]=0;
      istat[nj]=0;
      stat_etab[nj]=0;
      for(int icen=0;icen<ncen;icen++){
	njets[nj][icen]=0;
      }


      if(nj==0)mJets      = &(c->ak2PF);
      else if(nj==1)mJets = &(c->ak3PF);
      else if(nj==2)mJets = &(c->ak4PF);
      else if(nj==3)mJets = &(c->ak5PF);

      else if(nj==4)mJets = &(c->akPu2PF);
      else if(nj==5)mJets = &(c->akPu3PF);
      else if(nj==6)mJets = &(c->akPu4PF);
      else if(nj==7)mJets = &(c->akPu5PF);

      else if(nj==8)mJets = &(c->ak2Calo);
      else if(nj==9)mJets = &(c->ak3Calo);
      else if(nj==10)mJets = &(c->ak4Calo);
      else if(nj==11)mJets = &(c->ak5Calo);

      else if(nj==12)mJets = &(c->akPu2Calo);
      else if(nj==13)mJets = &(c->akPu3Calo);
      else if(nj==14)mJets = &(c->akPu4Calo);
      else if(nj==15)mJets = &(c->akPu5Calo);



      int *ljet = new int[3];
      FindLeadSubLeadJets(mJets,ljet);
      //! Jet energy scale comparison with data
      if(ljet[0]>=0 && ljet[1]>=0 && mJets->jtpt[ljet[0]]>50. && mJets->jtpt[ljet[1]]>50.){//! atleas a dijet
        int mstat=1;
	double ptdij = (mJets->jtpt[ljet[0]] + mJets->jtpt[ljet[1]])/2.;
        if(ljet[2]>=0){
          //if(iJet->jtpt[ljet[2]]/ptdij > 0.2)mstat=0;
          mstat=0;
        }
        if(mstat){
          double B=-9999;
          double rn1 = gRandom->Rndm();
          double rn2 = gRandom->Rndm();
          if(rn1 > rn2){
            B = (mJets->jtpt[ljet[0]] - mJets->jtpt[ljet[1]])/(mJets->jtpt[ljet[0]] + mJets->jtpt[ljet[1]]);
          }else{
            B = (mJets->jtpt[ljet[1]] - mJets->jtpt[ljet[0]])/(mJets->jtpt[ljet[1]] + mJets->jtpt[ljet[0]]);
          }
          hJetEnergyScale_mult[nj][imult]->Fill(ptdij,B);
          hJetEnergyScalemb[nj]->Fill(ptdij,B);
        }
      }
      delete [] ljet;
      


      //if(nj==1)std::cout<<"  \t \t  "<<Form("\t %s  : ",calgo[nj])<<mJets->nref<<"\t centb : "<<centb<<std::endl;
      //! Loop over # of jets in a given algo
      //if(nj==2)std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<RunNo<<"\t HF : "<<hiHF<<"\t ntracks : "<<hiNtracks<<"\t nref : "<<mJets->nref<<std::endl;

      for(int ij=0; ij<mJets->nref; ij++){

	//! Without any cut how many jets are there
	hTotJets[nj]->Fill(1.);

	//! jet selction cut
	if(mJets->jtpt[ij]<kptrecocut || fabs(mJets->jtphi[ij])>pi)continue;

	//! trackMax cuts
	//if(mJets->trackMax[ij]/mJets->jtpt[ij] < 0.01)continue;

	//! Jet bkgd estimation
	//! (photonSum+neutralSum+chargedSum-rawpt)
	double jbkgd  = (mJets->photonSum[ij]+mJets->neutralSum[ij]+mJets->chargedSum[ij]) - mJets->rawpt[ij];

	//! jets within different etabins
	//! 
	int ebin = GetEtaBin(fabs(mJets->jteta[ij]));
	if(ebin<0)continue;
	stat_etab[nj]=1;
	hjetrptmb_etab [nj][ebin]->Fill(mJets->rawpt[ij]);
	hjetjptmb_etab [nj][ebin]->Fill(mJets->jtpt[ij]);
	hjetjpumb_etab [nj][ebin]->Fill(mJets->jtpu[ij]);
	hjetptpumb_etab[nj][ebin]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);
	hjetptbkgd_etabmb[nj][ebin]->Fill(mJets->jtpt[ij],jbkgd);
	
	hjetrpt_etab [nj][centb][ebin]->Fill(mJets->rawpt[ij]);
	hjetjpt_etab [nj][centb][ebin]->Fill(mJets->jtpt[ij]);
	hjetjpu_etab [nj][centb][ebin]->Fill(mJets->jtpu[ij]);
	hjetptpu_etab[nj][centb][ebin]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);
	

	if(fabs(mJets->jteta[ij]) < ketacut){ //! |eta|<2.0

	  istat[nj]=1;

	  
	  if(mJets->jtpt[ij] > 120 && !iFirst){
	    iFirst = ij;
	  }

	  if(mJets->jtpt[ij] > 30 && !iSecond && !iFirst){
	    iSecond = ij;
	  }
	  

	  /*
	  if(nj==0){
	    ntrks_in=0;
	    ntrks_out=0;

	    double jetphi = mJets->jtphi[ij];
	    double jeteta = mJets->jteta[ij];

	    //! track loop
	    for(int j=0;j<c->track.nTrk;j++) {
	      if(!((c->track.trkPt[j]>0.2)
		   && (fabs(c->track.trkEta[j])>2.0)
		   //&& (c->track.trkEta[j]<0.47)
		   //&& (c->track.trkEta[j]>-1.53)
		   && (c->track.highPurity[j])
		   && (fabs(c->track.trkDz1[j]/c->track.trkDzError1[j])<3)
		   && (fabs(c->track.trkDxy1[j]/c->track.trkDxyError1[j])<3)
		   && (c->track.trkPtError[j]/c->track.trkPt[j]<0.1)
		   ))
		continue;
	      
	      double trkphi = c->track.trkPhi[j];
	      double trketa = c->track.trkEta[j];
	      
	      double dr = sqrt(pow((jetphi - trkphi),2) + pow((jeteta - trketa),2));

	      float dphi  = jetphi - trkphi;
	      if (dphi > pi ) dphi = dphi - 2 * pi;
	      if (dphi < -pi) dphi = dphi + 2 * pi;
	      //if(dphi > 2*pi/3.)continue;

	      //! count the tracks in and out of the jet cone
	      if(dr < 0.5){
		ntrks_in++;
		//std::cout<<"\t \t dr : "<<dr<<"\t dphi : "<< dphi <<std::endl;
	      } else if(dr > 0.5 && dphi < 2*pi/3.){
		ntrks_out++;
	      }
	    }//! track loop ends here
	    //std::cout<<"ntrks_in : "<<ntrks_in<<"\t ntrks_out : "<<ntrks_out<<std::endl;
	  }//! nj==0
	  */
	  
	  //! minbias
	  hjetrptmb [nj]->Fill(mJets->rawpt[ij]);
	  hjetjptmb [nj]->Fill(mJets->jtpt[ij]);
	  hjetratiorjptmb[nj]->Fill(mJets->jtpt[ij],mJets->jtpt[ij]/mJets->rawpt[ij]);
	  hjetjpumb [nj]->Fill(mJets->jtpu[ij]);
	  hjetptpumb[nj]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);
	  hjetptbkgdmb[nj]->Fill(mJets->jtpt[ij],jbkgd);
	  
	  hjetetamb   [nj]->Fill(mJets->jteta[ij]);
	  hjetphimb   [nj]->Fill(mJets->jtphi[ij]);
	  hjetptetamb [nj]->Fill(mJets->jtpt[ij],mJets->jteta[ij]);
	  hjetptphimb [nj]->Fill(mJets->jtpt[ij],mJets->jtphi[ij]);
	  hjetetaphimb[nj]->Fill(mJets->jteta[ij],mJets->jtphi[ij]);
	  
	  if(vbin>=0){
	    //! vertex z dependence
	    hjetrptvz [nj][vbin]->Fill(mJets->rawpt[ij]);
	    hjetjptvz [nj][vbin]->Fill(mJets->jtpt[ij]);
	    hjetjpuvz [nj][vbin]->Fill(mJets->jtpu[ij]);
	    hjetptpuvz[nj][vbin]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);
	    
	    hjetetavz   [nj][vbin]->Fill(mJets->jteta[ij]);
	    hjetphivz   [nj][vbin]->Fill(mJets->jtphi[ij]);
	    hjetptetavz [nj][vbin]->Fill(mJets->jtpt[ij],mJets->jteta[ij]);
	    hjetptphivz [nj][vbin]->Fill(mJets->jtpt[ij],mJets->jtphi[ij]);
	    hjetetaphivz[nj][vbin]->Fill(mJets->jteta[ij],mJets->jtphi[ij]);
	  }
	  
	  //! centrality
	  hNjets [nj][centb]->Fill(1.);

	  hfrpt  [nj][centb]->Fill(mJets->rawpt[ij]);
	  hfjpt  [nj][centb]->Fill(mJets->jtpt[ij]);
	  hfpupt [nj][centb]->Fill(mJets->jtpu[ij]);

	  hjetptpu[nj][centb]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);
	  hjetbkgd  [nj][centb]->Fill(jbkgd);
	  hjetptbkgd[nj][centb]->Fill(mJets->jtpt[ij],jbkgd);

	  hjeteta   [nj][centb]->Fill(mJets->jteta[ij]);
	  hjetphi   [nj][centb]->Fill(mJets->jtphi[ij]);
	  hjetpteta [nj][centb]->Fill(mJets->jtpt[ij],mJets->jteta[ij]);
	  hjetptphi [nj][centb]->Fill(mJets->jtpt[ij],mJets->jtphi[ij]);
	  hjetetaphi[nj][centb]->Fill(mJets->jteta[ij],mJets->jtphi[ij]);
	  


	  if(imult>=0){
	    hfrpt_mult   [nj][imult]->Fill(mJets->rawpt[ij]);
	    hfjpt_mult   [nj][imult]->Fill(mJets->jtpt[ij]);
	    hfpupt_mult  [nj][imult]->Fill(mJets->jtpu[ij]);
	    hjetptpu_mult[nj][imult]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);

	    hjetbkgd_mult  [nj][imult]->Fill(jbkgd);
	    hjetptbkgd_mult[nj][imult]->Fill(mJets->jtpt[ij],jbkgd);
	    
	    hjeteta_mult   [nj][imult]->Fill(mJets->jteta[ij]);
	    hjetphi_mult   [nj][imult]->Fill(mJets->jtphi[ij]);
	    hjetpteta_mult [nj][imult]->Fill(mJets->jtpt[ij],mJets->jteta[ij]);
	    hjetptphi_mult [nj][imult]->Fill(mJets->jtpt[ij],mJets->jtphi[ij]);
	    hjetetaphi_mult[nj][imult]->Fill(mJets->jteta[ij],mJets->jtphi[ij]);

	    int ieta=-1;
	    if(fabs(mJets->jteta[ij])<1.3)ieta=0; //! barrel region
	    else ieta=1; //! HCAL region         
	    hjetptbkgd_etab_mult[nj][imult][ieta]->Fill(mJets->jtpt[ij],jbkgd);
	  }



	  //! For a given event	  
	  njets[nj][centb]++;
	  sumjpt[nj]+=mJets->jtpt[ij];
	  sumjpu[nj]+=mJets->jtpu[ij];
	  sumjbkg[nj]+=jbkgd; 
	  nmbj[nj]++;


	  //! For a given run
	  totjetpt[nj][RunId] += mJets->jtpt[ij];
	  totjets [nj][RunId]++;
	  totjetsc[nj][RunId][centb]++;

	  hJetEtaPhiRun[nj]->Fill(RunId,mJets->jteta[ij],mJets->jtphi[ij]);

	} //! etacut
	//if(centb==4)std::cout<<Form(" \t \t \t mbjets %s : ",calgo[nj]) <<nmbj[nj]<<"\t centb : "<<centb<<"\t jtpt "<<mJets->jtpt[ij]<<std::endl;
      }//! # jet loop

      hNref[nj]->Fill(mJets->nref);
      hNevt_notr[nj][centb]->Fill(vz);

      if(stat_etab[nj]){
	hNevtmb_etab[nj]->Fill(vz);
	hNevt_etab[nj][centb]->Fill(vz);
      }

      if(nj==0){
	if(ntrks_in>0){
	  hNtrks_In->Fill(ntrks_in);
	  hNtrksNTracks_In->Fill(ntrks_in,hiNtracks);
	  hNtrksHF_In->Fill(ntrks_in,hiHF);
	}
	if(ntrks_out>0){
	  hNtrks_Out->Fill(ntrks_out);
	  hNtrksNTracks_Out->Fill(ntrks_out,hiNtracks);
	  hNtrksHF_Out->Fill(ntrks_out,hiHF);
	}
	if(ntrks_in>0 && ntrks_out>0){
	  hNtrks_InOut->Fill(ntrks_in,ntrks_out);
	}
      }

      if(istat[nj]){
	hNevt[nj][centb]->Fill(vz);
	if(vbin>=0)hNevtvz[nj][vbin]->Fill(vz);
	hMBJets[nj]->Fill(nmbj[nj]);

	hHF_tr[nj]->Fill(hiHF);

	//! # of jets distribution
	hNjets_dist[nj][centb]->Fill(njets[nj][centb]);
	hMBJetsRun [nj]->Fill(RunId,nmbj[nj]);
	hJetPtRun  [nj]->Fill(RunId,sumjpt[nj]/nmbj[nj]);
	hJetPtRun_Mult  [nj][imult]->Fill(RunId,sumjpt[nj]/nmbj[nj]);
	hJetBkgdRun [nj]->Fill(RunId,sumjbkg[nj]/nmbj[nj]);
	hJetBkgdRun_Mult  [nj][imult]->Fill(RunId,sumjbkg[nj]/nmbj[nj]);

	hJetPURun  [nj]->Fill(RunId,sumjpu[nj]/nmbj[nj]);
	totnevt[nj][RunId]++;
      }

      if(iFirst && iSecond){
	hHFNtrkO->Fill(hiHF,hiNtracks);
      }

    }//! nj loop

    

    //! Vz distribution
    hVz->Fill(vz);

    //! HF
    if(hiHFhit>0){
      hAvHF->Fill(hiHF/hiHFhit);
      hAvHFRun->Fill(RunId,hiHF/hiHFhit);
    }
    if(hiHFhitplus>0){
      hAvHFplus->Fill(hiHFplus/hiHFhitplus);
      hAvHFplusRun->Fill(RunId,hiHFplus/hiHFhitplus);
    }
    if(hiHFhitminus>0){
      hAvHFminus->Fill(hiHFplus/hiHFhitminus);
      hAvHFminusRun->Fill(RunId,hiHFplus/hiHFhitminus);
    }

    hHFPlusMinus->Fill(hiHFplus,hiHFminus);
    
    //! Luminosity run-by-run
    hLumiRun->Fill(RunId,LumiBlock);
    hLumiBlock->Fill(LumiBlock);
    hBin->Fill(hiBin);
    hEvtRun->Fill(RunId,vz);    
    
    //! ZDC
    hZDC->Fill(hiZDC);
    hZDCPlusMinus->Fill(hiZDCplus,hiZDCminus);
    
    //! EE
    hET->Fill(hiET);
    hEB->Fill(hiEB);
    hEE->Fill(hiEE);
    hEBEE->Fill(hiEB,hiEE);    
    hEEPlusMinus->Fill(hiEEplus,hiEEminus);    

    //! vertex
    hvz->Fill(vz);
    hvzrun->Fill(RunId,vz);
    hvxvy->Fill(vx,vy);

    //! Pixel
    hNpix->Fill(hiNpix);
    hNpixelTracks->Fill(hiNpixelTracks);
    hNpixelTracksNTracks->Fill(hiNpixelTracks,hiNtracks);

    hHF->Fill(hiHF);
    hNTracks->Fill(hiNtracks);
    hNTracksHF->Fill(hiNtracks,hiHF);
    hBinHF->Fill(hiBin,hiHF);
    hBinNTracks->Fill(hiBin,hiNtracks);

    //! Cross correlations
    hHFZDC->Fill(hiHF,hiZDC);
    hNpixelTracksZDC->Fill(hiNpixelTracks,hiZDC);
    hHFNpixelTracks->Fill(hiHF,hiNpixelTracks);

  }//! event loop ends

  if(!iCountRuns){
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    
    if(strcmp(ksp,"pp")==0){
      for(int nj=0;nj<knj;nj++){    
	std::cout<<"# of Events  : "<<calgo[nj]<<" \t : "<<hNevt[nj][ncen-1]->Integral()<<"\t # of Jets : "<<hNjets[nj][ncen-1]->Integral()<<std::endl;
      }
    }else{
      for(int nj=0;nj<knj;nj++){    
	for(int ic=0;ic<ncen;ic++){          
	  std::cout<<"# of Events  : "<<ccent[ic]<<"\t "<<calgo[nj]<<" \t : "<<hNevt[nj][ic]->Integral()<<"\t # of Jets : "<<hNjets[nj][ic]->Integral()<<std::endl;
	}
	std::cout<<"\t"<<std::endl;
      }
    }
    std::cout<<std::endl;
    std::cout<<std::endl;
    
    for(int nj=0;nj<knj;nj++){
      std::cout<<Form("For %s algo : ",calgo[nj])<<std::endl;
      for(int ir=0;ir<maxruns;ir++){
	hTotEvtRun[nj]->Fill(ir+1,nentries);
	if(Run[ir]<0 || totjets[nj][ir]==0)continue;
	std::cout<<"\t RunNo : "<<Run[ir]<<"\t Sumjetpt : "<<totjetpt[nj][ir]<<"\t Totjets : "<<totjets[nj][ir]<<"\t  <pT> : "<<totjetpt[nj][ir]/totjets[nj][ir]<<std::endl;
	hTotJetsRun[nj]->Fill(ir+1,totjets[nj][ir]);
	hNevtRun[nj]->Fill(ir+1,totnevt[nj][ir]);
	for(int ic=0;ic<ncen;ic++){
	  hTotJetsCentRun[nj][ic]->Fill(ir+1,totjetsc[nj][ir][ic]);
	}
      }
    }
  }

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
}
int GetMultBin(int nt)
{
  int ibin=-1;
  if(nt<60)ibin=0;
  else if(nt>=60 && nt<90)ibin=1;
  else if(nt>=90 && nt<110)ibin=2;
  else if(nt>=110&& nt<150)ibin=3;
  else if(nt>=150&& nt<180)ibin=4;
  else if(nt>=180)ibin=5;

  return ibin;
}
int GetCentBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 1% bins of cross section
  //! in 0-100 bins
  if(bin<10)ibin=0; //! 0-10%
  else if(bin>=10 && bin<30)ibin=1;  //! 10-30%
  else if(bin>=30 && bin<50)ibin=2;  //! 30-50%
  else if(bin>=50 && bin<70)ibin=3;  //! 50-70%
  else if(bin>=70 && bin<100)ibin=4;  //! 70-100%
  return ibin;
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
int GetVtxBin(float vz)
{
  int ibin=-1;
  ibin = (int)(nvtx*(vz+15.0)/30.0); //! vertex bin
  if(ibin<0 || ibin>=nvtx)return ibin=-1;
  else return ibin;

}
int GetEtaBin(float eta)
{
  int ibin=-1;

  if(eta>=0 && eta<1.3)ibin=0;       //! barrel
  else if(eta>=1.3 && eta<2.0)ibin=1;//! endcap+tracks
  else if(eta>=2.0 && eta<3.0)ibin=2;//! endcap-notracks
  else if(eta>=3.0 && eta<5.0)ibin=3;//! forward  

  return ibin;
}
void LoadRuns(const char *sp)
{
  if(strcmp(sp,"pbpb")==0){
    /*HiforestDijet_v7.root*/
    Run[0] = 181611; Run[1] = 181683; Run[2] = 181684; Run[3] = 181685; Run[4] = 181686;
    Run[5] = 181687; Run[6] = 181688; Run[7] = 181689; Run[8] = 181690; Run[9] = 181691;

    Run[10] = 181692; Run[11] = 181693; Run[12] = 181695; Run[13] = 181754; Run[14] = 181759;
    Run[15] = 181760; Run[16] = 181777; Run[17] = 181778; Run[18] = 181912; Run[19] = 181913;

    Run[20] = 181914; Run[21] = 181938; Run[22] = 181946; Run[23] = 181950; Run[24] = 181969;
    Run[25] = 181985; Run[26] = 182052; Run[27] = 182065; Run[28] = 182066; Run[29] = 182098;

    Run[30] = 182099; Run[31] = 182124; Run[32] = 182133; Run[33] = 182227; Run[34] = 182228;
    Run[35] = 182239; Run[36] = 182241; Run[37] = 182257; Run[38] = 182296; Run[39] = 182324; 

    Run[40] = 182365; Run[41] = 182382; Run[42] = 182398; Run[43] = 182422; Run[44] = 182536; 
    Run[45] = 182561; Run[46] = 182572; Run[47] = 182591; Run[48] = 182609; Run[49] = 182664; 

    Run[50] = 182785; Run[51] = 182798; Run[52] = 182822; Run[53] = 182838; Run[54] = 182890; 
    Run[55] = 182915; Run[56] = 182916; Run[57] = 182927; Run[58] = 182944; Run[59] = 182960; 

    Run[60] = 182972; Run[61] = 183013;

    //! Recorded Luminosity for these runs 
    //! Luminosity per run (everything is in /mub)
    Lumi[0] = 92.362e-03 ; Lumi[1] = 187.220e-03; Lumi[2] = 219.066e-03; Lumi[3] = 148.002e-03; Lumi[4] = 59.739e-03;
    Lumi[5] = 241.659e-03; Lumi[6] = 109.47e-03 ; Lumi[7] = 35.597e-03 ; Lumi[8] = 128.991e-03; Lumi[9] = 96.059e-03;

    Lumi[10] = 171.437e-03; Lumi[11] = 206.834e-03; Lumi[12] = 83.324e-03 ; Lumi[13] = 38.624e-03 ; Lumi[14] = 100.677e-03; 
    Lumi[15] = 224.282e-03; Lumi[16] = 103.743e-03; Lumi[17] = 169.573e-03; Lumi[18] = 697.353e-03; Lumi[19] =  2.333; 

    Lumi[20] = 47.166e-03; Lumi[21] =  3.184; Lumi[22] = 427.845e-03; Lumi[23] = 709.747e-03; Lumi[24] = 4.541; 
    Lumi[25] = 4.152     ; Lumi[26] =  4.971; Lumi[27] = 189.544e-03; Lumi[28] = 4.286      ; Lumi[29] = 995.874e-03;

    Lumi[30] = 4.351; Lumi[31] = 2.632; Lumi[32] = 1.101      ; Lumi[33] = 198.610e-03; Lumi[34] = 119.643e-03; 
    Lumi[35] = 1.544; Lumi[36] = 3.124; Lumi[37] = 596.649e-03; Lumi[38] = 3.138      ; Lumi[39] = 4.144;

    Lumi[40] = 5.207      ; Lumi[41] = 5.366; Lumi[42] = 105.077e-03; Lumi[43] = 5.123; Lumi[44] = 5.327; 
    Lumi[45] = 816.714e-03; Lumi[46] = 6.188; Lumi[47] = 5.909      ; Lumi[48] = 6.629; Lumi[49] = 1.745; 

    Lumi[50] = 5.539; Lumi[51] = 4.551      ; Lumi[52] = 5.172; Lumi[53] = 3.397; Lumi[54] = 5.002; 
    Lumi[55] = 5.193; Lumi[56] = 264.214e-03; Lumi[57] = 6.613; Lumi[58] = 3.569; Lumi[59] = 5.142; 

    Lumi[60] = 6.264; Lumi[61] = 6.410;
  }
  else if(strcmp(sp,"pp")==0){
    Run[0]  = 161366;Run[1]  = 161396;Run[2]  = 161404;Run[3]  = 161439;Run[4]  = 161445;Run[5]  = 161450;Run[6]  = 161454;Run[7]  = 161473;Run[8]  = 161474;

  }else if(strcmp(sp,"ppb")==0){
    Run[0]   = 202792; Run[1] = 210498; Run[2] = 210534; Run[3]  = 210614; Run[4] = 210634; Run[5] = 210635; Run[6] = 210638;
    Run[7]   = 210658; Run[8] = 210676; Run[9] = 210737; Run[10] = 210738; Run[11]= 210759; Run[12]= 210818; Run[13]= 210837;   
    Run[14]  = 210855; Run[15]= 210885; Run[16]= 210906; Run[17] = 210909; Run[18]= 210986; Run[19]= 210998; Run[20]= 211000; 
    Run[21]  = 211001; Run[22]= 211256; 

    //! Beam revered direction
    Run[23] = 211313; Run[24] = 211328; Run[25] = 211354;
    //! 
    //! Luminosity per run (everything is in /nb)
    Lumi[0]  = 0.00105;Lumi[1] = 0.00896 ; Lumi[2] = 0.148;  Lumi[3] = 1.215; Lumi[4] = 0.780; Lumi[5] = 0.583; Lumi[6] =0.088;
    Lumi[7]  = 1.854;  Lumi[8] = 1.321 ;   Lumi[9] = 0.302;  Lumi[10]= 0.769; Lumi[11]= 1.409; Lumi[12]= 1.000; Lumi[13]=1.177; 
    Lumi[14] = 1.179;  Lumi[15]= 1.138 ;   Lumi[16]= 0.930;  Lumi[17]=0.328;  Lumi[18]=0.958;  Lumi[19]= 0.519; Lumi[20]=0.358;
    Lumi[21] = 0.108;  Lumi[22]= 1.031;

    Lumi[23]= 0.09; Lumi[24] = 0.463; Lumi[25]=0.436;
  }
  std::cout<<"Run # loaded according to increasing order" <<std::endl;
}
void LoadLib()
{
  gSystem->Load("../HiForest/V2/hiForest_h.so");
}
void FindLeadSubLeadJets(Jets *jetc, int *ljet)
{
  ljet[0]=-1; ljet[1]=-2;ljet[2]=-3;

  float tempt=-9;
  //! Get the leading jet                                                                                                                                                                                           
  for(int ij=0; ij<jetc->nref; ij++){
    if(fabs(jetc->jteta[ij])>2.0 || jetc->jtpt[ij]<30)continue;
    float jetpt = jetc->jtpt[ij];
    if(jetpt > tempt){
      tempt = jetpt;
      ljet[0] = ij;
    }
  }

  tempt=-9;
  for(int ij=0; ij<jetc->nref; ij++){
    if(ij==ljet[0])continue;
    if(fabs(jetc->jteta[ij])>2.0 || jetc->jtpt[ij]<30)continue;
    float jetpt = jetc->jtpt[ij];
    float dphi  = jetc->jtphi[ij] - jetc->jtphi[ljet[0]];
    if (dphi > pi ) dphi = dphi - 2 * pi;
    if (dphi < -pi) dphi = dphi + 2 * pi;
    if(dphi < 2*pi/3.)continue;
    if (jetpt > tempt){
      tempt = jetpt;
      ljet[1] = ij;
    }
  }

  tempt=-9;
  for(int ij=0; ij<jetc->nref; ij++){
    if(ij==ljet[0] || ij==ljet[1])continue;
    if(fabs(jetc->jteta[ij])>2.0|| jetc->jtpt[ij]<30)continue;
    float jetpt = jetc->jtpt[ij];
    if (jetpt > tempt){
      tempt = jetpt;
      ljet[2] = ij;
    }
  }
}
