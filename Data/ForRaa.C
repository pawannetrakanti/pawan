#include "/afs/cern.ch/user/p/pawan/scratch0/CMSSW_4_4_2/src/UserCode/HiForest/hiForest.h"

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
#include <cstdlib>
#include <cmath>
#include <iostream>
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

const int ncen=7;
const char *ccent[ncen]={"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","pp"};

//const int knj=7;    //! algo 0: ic5, 1: ak3pu
//const char *calgo[knj]={"icPu5","akPu2PF","akPu3PF","akPu4PF","akPu2Calo","akPu3Calo","akPu4Calo"};

const int knj=3;    //! algo 0: ic5, 1: ak3pu
const char *calgo[knj]={"akPu2PF","akPu3PF","akPu4PF"};

const int nvtx=4;
const int npix=5;
const int ketar=4;
const int maxruns=62;

int Run[maxruns];
double Lumi[maxruns];

void LoadRuns  (const char */*sp*/);
int GetRunIdx  (int /*irun*/);
int GetCentBin (int /*bin*/);
int GetVtxBin  (float /*vz*/);
int GetPixelBin(int /*hiNpixeltracks*/);
int GetEtaBin  (float /*eta*/);


class DuplicateEvents {
public:
  DuplicateEvents(TString infname) {
    inf = TFile::Open(infname);
    t = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
  };
  ~DuplicateEvents() {
    delete inf;
  }
  void MakeList() {
    cout << "Starting Making List to check for duplicate events" << endl;
    evts.clear();
    occurrence.clear();
    int run,evt;
    t->SetBranchAddress("run",&run);
    t->SetBranchAddress("evt",&evt);
    for (int i=0;i<t->GetEntries();i++) {
      t->GetEntry(i);
      if (i%100000==0) cout <<i<<" / "<<t->GetEntries() << " run: " << run << " evt: " << evt << endl;
      int occur = (int)FindOccurrences(run,evt);
      if (occur==0) occurrence.push_back(1);
      else occurrence.push_back(2);
      evts.push_back(std::make_pair(run,evt));
    }
  }
  int FindOccurrences(int run, int evt) {
    int noccur = count(evts.begin(), evts.end(), std::make_pair(run,evt));
    return noccur;
  }
  TFile * inf;
  TTree * t;
  vector<pair<int, int> > evts;
  vector<int> occurrence;
};


//! Smearing factor for only akPu's
/*
//!   Reco pT Cut                                 Centrality
//!   pT>15 GeV                0-5%  , 5-10%  ,  10-30% ,  30-50% ,  50-70% ,   70-90%  ,  pp    
//! with new pp |vz|<15 cm  with full pT hat bins and will be integrated from 80 to 400 GeV in ref pT (in withvtx dir)
double smearf[knj][ncen]    ={{8.883  , 7.446  ,  6.025  ,  2.831  ,  0.956  ,   1.500   ,  0.00},  //! akpu2pf
			      {9.905  , 8.865  ,  7.354  ,  4.664  ,  3.163  ,   0.000   ,  0.00},  //! akpu3pf
			      {11.336 , 10.195 ,  8.456  ,  5.962  ,  4.158  ,   1.585   ,  0.00},  //! akpu4pf
}; 
*/

double ptbins_smear[] = {80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400};
const int bins1 = sizeof(ptbins_smear)/sizeof(double) - 1;
const int smbins = bins1;
double smearf[knj][ncen][smbins] = {
  { {10.6069,9.53263,7.99496,6.47572,6.63517,6.43069,4.79561,5.52368,1.51585,4.73694,5.49985,4.34136,2.69882,3.18559,2.75011,0 }, //! 0-5%  //! ak2PF                
    {9.5939,7.59408,6.41255,2.52901,0,3.03726,3.73099,5.65377,0,3.46223,0,0,2.48996,5.21187,3.39047,4.0801 }, //! 5-10%                                              
    {7.52046,7.09926,4.87256,5.56958,6.40448,5.15708,5.45629,3.48318,3.03069,0,5.56163,4.05211,2.18542,0.828933,4.28571,0.963933 }, //! 10-30%                       
    {4.15823,4.29915,4.08619,0,4.64429,3.1441,0,2.75549,0,3.51214,1.52933,0,2.96407,2.92293,0,4.12055 }, //! 30-50%                                                  
    {3.40076,1.91548,2.24094,4.40189,0,1.06348,3.58938,2.98547,2.11986,0,0,2.39968,0,0,1.17852,0 }, //! 50-70%                                                       
    {3.46271,1.38165,1.47467,0,3.91275,2.09271,0,0,1.40595,0,2.66064,2.75378,1.6508,0,0,3.06773 }, //! 70-90%                                                        
    { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,0.00,0.00,0.00,0.00}}, //! pp                                                           
  
  { {10.0971,10.0556,11.53382,10.09026,11.04904,11.89586,11.82574,10.04593,6.36579,6.71464,6.78906,7.88803,0,7.553219,7.559386,0}, //! 0-5%  //! ak3PF               
    {9.3828,8.2346,8.27719,7.70535,6.25354,7.98441,6.34697,8.95322,8.86133,0,0.097127,0.093448,0,8.376,8.87887,3.3609 }, //! 5-10%                                   
    {7.486,7.78749,7.84899,7.793231,7.62483,7.73366,6.93853,5.0297,3.58286,0.02722,5.53265,6.9743,0.040309,0,6.25054,0.046698 }, //! 10-30%                          
    {5.15368,5.0759,5.3187,4.06908,5.98657,3.83484,3.07333,1.26228,0.50647,3.89258,3.03255,1.04565,1.07938,4.03525,0,4.13655 }, //! 30-50%                           
    {3.99394,3.84754,4.03742,4.56218,4.32862,2.39101,3.83596,4.54976,0.041902,0,0,0.09433,0,2.0601,1.49476,0 }, //! 50-70%                                           
    {2.03394,2.01458,2.79081,0,0.026742,0.0218433,0.021867,0,0,0,0.0281399,0.00403,0.010721,0.048432,0.08341,0.04126 }, //! 70-90%                                   
    { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,0.00,0.00,0.00,0.00}}, //! pp                                                           
  
  { {12.2563,11.5695,10.4316,8.80091,6.58487,7.90287,13.62099,13.63562,0.19621,3.12709,4.27031,6.8802,3.91001,3.8079,3.27438,4.06707 }, //! 0-5%  //! ak4PF          
    {12.0168,9.35782,8.00736,7.29842,7.97529,6.42666,4.56315,4.51409,5.95295,1.24883,4.64867,3.70103,2.31134,6.32349,4.89873,3.1264 }, //! 5-10%                     
    {9.6697,8.94999,13.88295,12.53489,12.92233,15.93666,5.38456,4.00389,4.43362,3.15875,5.23279,4.55506,4.00098,2.07628,4.37486,3.09102 }, //! 10-30%                
    {7.49541,6.63023,5.87943,4.95293,5.67788,4.70276,3.33361,3.45651,2.01514,3.88602,2.04301,2.5694,1.38983,4.30535,2.42315,3.50593 }, //! 30-50%                    
    {6.22191,4.93202,5.07923,4.44174,4.10758,2.37307,3.582,4.4269,3.46693,0,3.01783,2.34337,0,1.90033,2.66977,0 }, //! 50-70%                                        
    {4.56134,3.35514,3.19569,0.1903,0.00174,0.06648,0.045127,0.01317,0.022475,0,0.066526,0.032671,0.082098,0,0.056494,0.098529 }, //! 70-90%                         
    { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,0.00,0.00,0.00,0.00}}, //! pp                                        
};
int GetPtBinSmear(float /*pt*/);




TStopwatch timer;
void ForRaa(const char *ksp="pp")
{

  timer.Start();

  const float ketacut=2.;
  const double kptrecocut=80.;
  const bool iCountRuns=false;

  //! Load Lib
  gSystem->Load("/afs/cern.ch/user/p/pawan/scratch0/CMSSW_4_4_2/src/UserCode/HiForest/hiForest_h.so");
  
  char mrun[6]={0};
  //! Load Run numbers 
  for(int ir=0;ir<maxruns;ir++){
    Run[ir]=-999;Lumi[ir]=-999;
  }
  if(!iCountRuns)LoadRuns(Form("%s",ksp));

  //! pt binning
  //!                   0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   22   23
  //double ptbins[22] ={100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,220.,240.,260.,280.,300.,325.,350.,375.,400.,450.,558.,638.,828};


  //! pt binning
  double ptbins[] ={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
  const int bins = sizeof(ptbins)/sizeof(double) - 1;

  //! Define the input file and HiForest
  //! CMSSW_4_4_2

  TString inname="";
  bool ispp=false;
  if(strcmp(ksp,"pbpb")==0){
    //! merged from prompt reco
    //inname="rfio:/castor/cern.ch/user/f/frankma/forest2/HiForest-promptskim-hihighpt-hltjet80-pt90-v3_part.root";
    inname="/d102/yjlee/hiForest2/HiForest-promptskim-hihighpt-hltjet80-pt90-v3_part.root";

  }else if(strcmp(ksp,"pp")==0){
    //! New production
    //inname="rfio:/castor/cern.ch/user/f/frankma/forest2/HiForest-ppskim-hihighpt-pt90-v1_v3_part.root";
    inname="/d102/yjlee/hiForest2PP/HiForest-ppskim-hihighpt-pt90-v1_v3_part.root";

    ispp=true;
  }
  

  DuplicateEvents dupEvt(inname);
  if(!ispp){
    // Check for duplicate events only for PbPb
    dupEvt.MakeList();
  }


  //! Define the input file and HiForest
  //! CMSSW_4_4_2
  HiForest *c = new HiForest(inname,Form("MyForest%s",ksp),ispp,0);


  //! Output file
  //! HIHighPt
  TFile *fout = 0;
  if(strcmp(ksp,"pp")==0){
    fout = new TFile(Form("Output/HLT_Jet60_v1_ppSmearing_HiForest2_%s_pt%0.0f_2012.root",ksp,kptrecocut),"RECREATE");  
    //fout = new TFile(Form("Output/HLT_Jet60_v1_RandomCone_HiForest2_%s_pt%0.0f_2012.root",ksp,kptrecocut),"RECREATE");  
  }else{
    fout = new TFile(Form("Output/HLT_HIJet80_v1_HiForest2_%s_pt%0.0f_2012.root",ksp,kptrecocut),"RECREATE");  
  }

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
  c->hasIcPu5JetTree=0;
  if(strcmp(ksp,"pp")==0){
    c->hasAkPu2JetTree=0;
    c->hasAkPu3JetTree=0;
    c->hasAkPu4JetTree=0;
  }else{
    c->hasAk2JetTree=0;
    c->hasAk3JetTree=0;
    c->hasAk4JetTree=0;
  }
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;


  //! photon tree
  c->hasPhotonTree=0;

  c->hasTrackTree=0;
  c->hasPixTrackTree=0;
  c->hasTowerTree=0;
  c->hasHbheTree=0;
  c->hasEbTree=0;
  c->hasGenpTree=0;
  c->hasGenParticleTree=0;

  //! added by pawan
  //! Select only the branches you want to use for the analysis
  //! This increases the speed for running

  //! For Hlt
  c->hltTree->SetBranchStatus("*",0,0);
  if(strcmp(ksp,"pp")==0) c->hltTree->SetBranchStatus("HLT_Jet60_v1",1,0);
  else  c->hltTree->SetBranchStatus("HLT_HIJet80_v1",1,0);

  //! for Skim Tree
  c->skimTree->SetBranchStatus("*",0,0);
  c->skimTree->SetBranchStatus("pcollisionEventSelection",1,0);
  c->skimTree->SetBranchStatus("pHBHENoiseFilter",1,0);
  /*
  if(ispp){
  c->skimTree->SetBranchStatus("phiEcalRecHitSpikeFilter",1,0);
  //c->skimTree->SetBranchStatus("phbheReflagNewTimeEnv",1,0);
  //c->skimTree->SetBranchStatus("phcalTimingFilter",1,0);
  //c->skimTree->SetBranchStatus("phfCoincFilter",1,0);
  //c->skimTree->SetBranchStatus("ppurityFractionFilter",1,0);
  }
  */

  //! Evt tree
  //c->evtTree->SetBranchStatus("*",0,0);
  //c->evtTree->SetBranchStatus("run",1,0);
  //c->evtTee->SetBranchStatus("lumi",1,0);
  //c->evtTree->SetBranchStatus("hiBin",1,0);//! centrality bin
  //c->evtTree->SetBranchStatus("hiHF",1,0);//! HF
  //c->evtTree->SetBranchStatus("hiHFplus",1,0); //! HF plus
  //c->evtTree->SetBranchStatus("hiHFminus",1,0);//! HF minus
  //c->evtTree->SetBranchStatus("hiZDC",1,0);//!ZDC

  if(c->hasIcPu5JetTree){
    //! For jets icpu5
    c->icPu5jetTree->SetBranchStatus("*",0,0);
    c->icPu5jetTree->SetBranchStatus("nref",1,0);
    c->icPu5jetTree->SetBranchStatus("rawpt",1,0); 
    c->icPu5jetTree->SetBranchStatus("jtpt",1,0); 
    c->icPu5jetTree->SetBranchStatus("jteta",1,0); 
    c->icPu5jetTree->SetBranchStatus("jtphi",1,0); 
    c->icPu5jetTree->SetBranchStatus("jtpu",1,0); 
    c->icPu5jetTree->SetBranchStatus("trackMax",1,0); 
  }

  //
  //
  //

  if(c->hasAk2JetTree){
    //! For jets ak2pf
    c->ak2jetTree->SetBranchStatus("*",0,0);
    c->ak2jetTree->SetBranchStatus("nref",1,0);
    c->ak2jetTree->SetBranchStatus("rawpt",1,0); 
    c->ak2jetTree->SetBranchStatus("jtpt",1,0); 
    c->ak2jetTree->SetBranchStatus("jteta",1,0); 
    c->ak2jetTree->SetBranchStatus("jtphi",1,0); 
    c->ak2jetTree->SetBranchStatus("jtpu",1,0); 
    c->ak2jetTree->SetBranchStatus("trackMax",1,0); 
  }
  
  if(c->hasAk3JetTree){
    //! For jets ak3pf
    c->ak3jetTree->SetBranchStatus("*",0,0);
    c->ak3jetTree->SetBranchStatus("nref",1,0);
    c->ak3jetTree->SetBranchStatus("rawpt",1,0); 
    c->ak3jetTree->SetBranchStatus("jtpt",1,0); 
    c->ak3jetTree->SetBranchStatus("jteta",1,0); 
    c->ak3jetTree->SetBranchStatus("jtphi",1,0); 
    c->ak3jetTree->SetBranchStatus("jtpu",1,0); 
    c->ak3jetTree->SetBranchStatus("trackMax",1,0); 
  }

  if(c->hasAk4JetTree){
    //! For jets ak4pf
    c->ak4jetTree->SetBranchStatus("*",0,0);
    c->ak4jetTree->SetBranchStatus("nref",1,0);
    c->ak4jetTree->SetBranchStatus("rawpt",1,0); 
    c->ak4jetTree->SetBranchStatus("jtpt",1,0); 
    c->ak4jetTree->SetBranchStatus("jteta",1,0); 
    c->ak4jetTree->SetBranchStatus("jtphi",1,0); 
    c->ak4jetTree->SetBranchStatus("jtpu",1,0); 
    c->ak4jetTree->SetBranchStatus("trackMax",1,0); 
  }


  if(c->hasAkPu2JetTree){
    //! For jets akpu2pf
    c->akPu2jetTree->SetBranchStatus("*",0,0);
    c->akPu2jetTree->SetBranchStatus("nref",1,0);
    c->akPu2jetTree->SetBranchStatus("rawpt",1,0); 
    c->akPu2jetTree->SetBranchStatus("jtpt",1,0); 
    c->akPu2jetTree->SetBranchStatus("jteta",1,0); 
    c->akPu2jetTree->SetBranchStatus("jtphi",1,0); 
    c->akPu2jetTree->SetBranchStatus("jtpu",1,0); 
    c->akPu2jetTree->SetBranchStatus("trackMax",1,0); 
  }

  if(c->hasAkPu3JetTree){
    //! For jets akpu3pf
    c->akPu3jetTree->SetBranchStatus("*",0,0);
    c->akPu3jetTree->SetBranchStatus("nref",1,0);
    c->akPu3jetTree->SetBranchStatus("rawpt",1,0); 
    c->akPu3jetTree->SetBranchStatus("jtpt",1,0); 
    c->akPu3jetTree->SetBranchStatus("jteta",1,0); 
    c->akPu3jetTree->SetBranchStatus("jtphi",1,0); 
    c->akPu3jetTree->SetBranchStatus("jtpu",1,0); 
    c->akPu3jetTree->SetBranchStatus("trackMax",1,0); 
  }

  if(c->hasAkPu4JetTree){
    //! For jets akpu4pf
    c->akPu4jetTree->SetBranchStatus("*",0,0);
    c->akPu4jetTree->SetBranchStatus("nref",1,0);
    c->akPu4jetTree->SetBranchStatus("rawpt",1,0); 
    c->akPu4jetTree->SetBranchStatus("jtpt",1,0); 
    c->akPu4jetTree->SetBranchStatus("jteta",1,0); 
    c->akPu4jetTree->SetBranchStatus("jtphi",1,0); 
    c->akPu4jetTree->SetBranchStatus("jtpu",1,0); 
    c->akPu4jetTree->SetBranchStatus("trackMax",1,0); 
  }



  if(c->hasAkPu2CaloJetTree){
    //! For jets akpu2calo
    c->akPu2CaloJetTree->SetBranchStatus("*",0,0);
    c->akPu2CaloJetTree->SetBranchStatus("nref",1,0);
    c->akPu2CaloJetTree->SetBranchStatus("rawpt",1,0); 
    c->akPu2CaloJetTree->SetBranchStatus("jtpt",1,0); 
    c->akPu2CaloJetTree->SetBranchStatus("jteta",1,0); 
    c->akPu2CaloJetTree->SetBranchStatus("jtphi",1,0); 
    c->akPu2CaloJetTree->SetBranchStatus("jtpu",1,0); 
    c->akPu2CaloJetTree->SetBranchStatus("trackMax",1,0); 
  }
  
  if(c->hasAkPu3CaloJetTree){
    //! For jets akpu3calo
    c->akPu3CaloJetTree->SetBranchStatus("*",0,0);
    c->akPu3CaloJetTree->SetBranchStatus("nref",1,0);
    c->akPu3CaloJetTree->SetBranchStatus("rawpt",1,0); 
    c->akPu3CaloJetTree->SetBranchStatus("jtpt",1,0); 
    c->akPu3CaloJetTree->SetBranchStatus("jteta",1,0); 
    c->akPu3CaloJetTree->SetBranchStatus("jtphi",1,0); 
    c->akPu3CaloJetTree->SetBranchStatus("jtpu",1,0); 
    c->akPu3CaloJetTree->SetBranchStatus("trackMax",1,0); 
  }
  
  if(c->hasAkPu4CaloJetTree){
    //! For jets akpu4calo
    c->akPu4CaloJetTree->SetBranchStatus("*",0,0);
    c->akPu4CaloJetTree->SetBranchStatus("nref",1,0);
    c->akPu4CaloJetTree->SetBranchStatus("rawpt",1,0); 
    c->akPu4CaloJetTree->SetBranchStatus("jtpt",1,0); 
    c->akPu4CaloJetTree->SetBranchStatus("jteta",1,0); 
    c->akPu4CaloJetTree->SetBranchStatus("jtphi",1,0); 
    c->akPu4CaloJetTree->SetBranchStatus("jtpu",1,0); 
    c->akPu4CaloJetTree->SetBranchStatus("trackMax",1,0); 
  }

  std::cout<<"Loaded all tree variables "<<std::endl;


  //! For jets
  Jets *mJets;

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
  TH2F *hHFPlusMinus = new TH2F("hHFPlusMinus","HF plus:minus",600,0,6000,600,0,6000);
  TH1F *hHF     = new TH1F("hHF","HF distribution",600,0,6000);
  TH1F *hAvHF     = new TH1F("hAvHF","HF E/hit distribution",100,0,0.1);
  TH1F *hAvHFplus = new TH1F("hAvHFplus","HFplus E/hit distribution",100,0,0.1);
  TH1F *hAvHFminus = new TH1F("hAvHFminus","HFminus E/hit distribution",100,0,0.2);

  TH2F *hAvHFRun     = new TH2F("hAvHFRun","HF E/hit distribution Run-by-Run",maxruns,-0.5,maxruns-0.5,100,0,0.1);
  TH2F *hAvHFplusRun = new TH2F("hAvHFplusRun","HFplus E/hit distribution Run-by-Run",maxruns,-0.5,maxruns-0.5,100,0,0.1);
  TH2F *hAvHFminusRun = new TH2F("hAvHFminusRun","HFminus E/hit distribution Run-by-Run",maxruns,-0.5,maxruns-0.5,100,0,0.2);

  TH2F *hZDCPlusMinus = new TH2F("hZDCPlusMinus","ZDC plus:minus",1000,0,45000,1000,0,45000);
  TH1F *hZDC     = new TH1F("hZDC","ZDC distribution",5000,0,1e+05);

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
  TH1F *hNpix = new TH1F("hNpix","hiNpix",500,0,80000);
  TH1F *hNpixelTracks = new TH1F("hNpixelTracks","hiNpixelTracks",200,0,3000);
  TH1F *hNTracks = new TH1F("hNTracks","hiNTracks",100,0,1500);
  TH2F *hNpixelTracksNTracks = new TH2F("hNpixelTracksNTracks","hNpixelTracksNTracks",500,0,15000,100,0,1500);

  //! Cross correlations
  TH2F *hHFZDC = new TH2F("hHFZDC","hHFZDC",600,0,6000,1000,0,1e+05);
  TH2F *hNpixelTracksZDC = new TH2F("hNpixelTracksZDC","hNpixelTracksZDC",200,0,3000,1000,0,1e+05);
  TH2F *hHFNpixelTracks  = new TH2F("hHFNpixelTracks","hHFNpixelTracks",600,0,6000,200,0,3000);

  TH1F *hBin  = new TH1F("hBin","Centrality bin",40,-0.5,39.5);

  TH1F *hHF_tr[knj];
  TH1F *hTotJets[knj];

  TH2F *hJetPtRun[knj];
  TH2F *hJetPURun[knj];
  TH2F *hMBJetsRun[knj];
  TH2F *hTotJetsRun[knj];

  //! mb
  TH1F *hjetrptmb[knj], *hjetjptmb[knj], *hjetjpumb[knj];
  TH1F *hjetetamb[knj], *hjetphimb[knj];
  TH2F *hjetptetamb[knj], *hjetptphimb[knj], *hjetetaphimb[knj];
  TH2F *hjetptpumb[knj];

  //! vetex bin
  TH1F *hjetrptvz[knj][nvtx], *hjetjptvz[knj][nvtx], *hjetjpuvz[knj][nvtx];
  TH1F *hjetetavz[knj][nvtx], *hjetphivz[knj][nvtx];
  TH2F *hjetptetavz[knj][nvtx], *hjetptphivz[knj][nvtx], *hjetetaphivz[knj][nvtx];
  TH2F *hjetptpuvz[knj][nvtx];

  //! pixel bin
  TH1F *hjetrpt_pix[knj][npix], *hjetjpt_pix[knj][npix], *hjetjpu_pix[knj][npix];
  TH1F *hjeteta_pix[knj][npix], *hjetphi_pix[knj][npix];
  TH2F *hjetpteta_pix[knj][npix], *hjetptphi_pix[knj][npix], *hjetetaphi_pix[knj][npix];
  TH2F *hjetptpu_pix[knj][npix];

  //! ncen bin
  TH1F *hfrpt[knj][ncen], *hfjpt[knj][ncen], *hfpupt[knj][ncen];

  TH1F *hjetrpt  [knj][ncen], *hjetjpt[knj][ncen], *hjetjpu[knj][ncen], *hjetspt[knj][ncen];
  TH1F *hjeteta  [knj][ncen], *hjetphi[knj][ncen];
  TH2F *hjetpteta[knj][ncen], *hjetptphi[knj][ncen], *hjetetaphi[knj][ncen];
  TH2F *hjetptpu[knj][ncen];

  //! eta dependence
  TH1F *hjetrptmb_etab[knj][ketar];
  TH1F *hjetjptmb_etab[knj][ketar];
  TH1F *hjetjpumb_etab[knj][ketar];
  TH2F *hjetptpumb_etab[knj][ketar];
  TH1F *hNevtmb_etab[knj];

  TH1F *hjetrptpix_etab[knj][npix][ketar];
  TH1F *hjetjptpix_etab[knj][npix][ketar];
  TH1F *hjetjpupix_etab[knj][npix][ketar];
  TH2F *hjetptpupix_etab[knj][npix][ketar];
  TH1F *hNevtpix_etab[knj][npix];

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
  TH1F *hNevtpix   [knj][npix];
  TH1F *hNevt_notr [knj][ncen];


  for(int nj=0;nj<knj;nj++){
    //std::cout<<Form("Histograms for algo %s",calgo[nj])<<std::endl;

    //! Define the histograms for jets
    hTotJets[nj] = new TH1F(Form("hTotJets%d",nj),Form("# of jets w/o cuts with %s",calgo[nj]),4,0,4);
    hMBJets[nj] = new TH1F(Form("hMBJets%d",nj),Form("MB distribution of jets in %s",calgo[nj]),50,0,50);

    hMBJetsRun[nj] = new TH2F(Form("hMBJetsRun%d",nj),Form("MB distribution of jets run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,20,0,20);
    hMBJetsRun[nj]->GetXaxis()->SetLabelSize(0.03); 

    //! for jet pt
    hJetPtRun[nj] = new TH2F(Form("hJetPtRun%d",nj),Form("Jet pt run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,100,0,300);

    //! for jet pile up
    hJetPURun[nj] = new TH2F(Form("hJetPURun%d",nj),Form("Jet pile up run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,100,0,300);

    hTotJetsRun[nj] = new TH2F(Form("hTotJetsRun%d",nj),Form("Total jets run-by-run in %s",calgo[nj]),maxruns,-0.5,maxruns-0.5,500,0,10000);
    hTotJetsRun[nj]->GetXaxis()->SetLabelSize(0.03); 

    for(int ix=1;ix<=hMBJetsRun[nj]->GetNbinsX();ix++){
      if(Run[ix-1]!=-999)sprintf(mrun,"%d",Run[ix-1]);
      else sprintf(mrun,"%s","\0");

      hMBJetsRun[nj]->GetXaxis()->SetBinLabel(ix,mrun);
      hJetPtRun [nj]->GetXaxis()->SetBinLabel(ix,mrun);
      hJetPURun [nj]->GetXaxis()->SetBinLabel(ix,mrun);
      hTotJetsRun[nj]->GetXaxis()->SetBinLabel(ix,mrun);
    }

    hNref[nj] = new TH1F(Form("hNRef%d",nj),Form("# of jets in %s",calgo[nj]),40,0,100);
    hHF_tr[nj]      = new TH1F(Form("hHF_tr%d",nj),Form("HF triggered events (Jet p_{T} > 85 GeV/c distribution) %s",calgo[nj]),600,0,6000);
    
    hjetetamb[nj] = new TH1F(Form("hjetetamb%d",nj),Form("jet eta distribution jet %s",calgo[nj]),18,-3.0,3.0);
    hjetphimb[nj] = new TH1F(Form("hjetphimb%d",nj),Form("jet phi distribution jet %s",calgo[nj]),18,-pi,pi);
    hjetrptmb[nj] = new TH1F(Form("hjetrptmb%d",nj),Form("jet(raw pt) p_{T} distribution jet %s",calgo[nj]),bins,ptbins);
    hjetjptmb[nj] = new TH1F(Form("hjetjptmb%d",nj),Form("jet(jt pt) p_{T} distribution jet %s",calgo[nj]),bins,ptbins);
    hjetjpumb[nj] = new TH1F(Form("hjetjpumb%d",nj),Form("jet(j pu) p_{T} distribution jet %s",calgo[nj]),100,0,300);
    hjetptpumb[nj] = new TH2F(Form("hjetptpumb%d",nj),Form("jet(pt:pu) distribution jet %s",calgo[nj]),bins,ptbins,100,0,300);
    
    hjetetaphimb[nj] = new TH2F(Form("hjetetaphimb%d",nj),Form("jet eta-phi distribution jet mb %s",calgo[nj]),36,-3.0,3.0,36,-pi,pi);
    hjetptetamb[nj] = new TH2F(Form("hjetptetamb%d",nj),Form("jet pt-eta distribution jet mb %s",calgo[nj]),bins,ptbins,36,-3.0,3.0);
    hjetptphimb[nj] = new TH2F(Form("hjetptphimb%d",nj),Form("jet pt-phi distribution jet mb %s",calgo[nj]),bins,ptbins,36,-pi,pi);

    hNevtmb_etab[nj] = new TH1F(Form("hNevtmb_etab%d",nj),Form("# of events mb jet %s",calgo[nj]),40,-40,40);
    for(int ie=0;ie<ketar;ie++){
      hjetrptmb_etab[nj][ie] = new TH1F(Form("hjetrptmb_etab%d_%d",nj,ie),Form("jet(raw pt) p_{T} distribution jet %s etabin %d",calgo[nj],ie),bins,ptbins);
      hjetjptmb_etab[nj][ie] = new TH1F(Form("hjetjptmb_etab%d_%d",nj,ie),Form("jet(jt pt) p_{T} distribution jet %s etabin %d",calgo[nj],ie),bins,ptbins);
      hjetjpumb_etab[nj][ie] = new TH1F(Form("hjetjpumb_etab%d_%d",nj,ie),Form("jet(j pu) p_{T} distribution jet %s etabin %d",calgo[nj],ie),100,0,300);
      hjetptpumb_etab[nj][ie] = new TH2F(Form("hjetptpumb_etab%d_%d",nj,ie),Form("jet(pt:pu) distribution jet %s etabin %d",calgo[nj],ie),bins,ptbins,100,0,300);
    }

    for(int ip=0;ip<npix;ip++){
      hjeteta_pix[nj][ip] = new TH1F(Form("hjeteta_pix%d_%d",nj,ip),Form("jet eta distribution jet pix %d %s",ip,calgo[nj]),18,-3.0,3.0);
      hjetphi_pix[nj][ip] = new TH1F(Form("hjetphi_pix%d_%d",nj,ip),Form("jet phi distribution jet pix %d %s",ip,calgo[nj]),18,-pi,pi);
      hjetrpt_pix[nj][ip] = new TH1F(Form("hjetrpt_pix%d_%d",nj,ip),Form("jet(raw pt) p_{T} distribution jet pix %d %s",ip,calgo[nj]),bins,ptbins);
      hjetjpt_pix[nj][ip] = new TH1F(Form("hjetjpt_pix%d_%d",nj,ip),Form("jet(jt pt) p_{T} distribution jet vz pix %d %s",ip,calgo[nj]),bins,ptbins);
      hjetjpu_pix[nj][ip] = new TH1F(Form("hjetjpu_pix%d_%d",nj,ip),Form("jet(pu) p_{T} distribution jet vz pix %d %s",ip,calgo[nj]),100,0,300);
      hjetptpu_pix[nj][ip] = new TH2F(Form("hjetptpu_pix%d_%d",nj,ip),Form("jet(pt:pu) distribution jet pix %d %s",ip,calgo[nj]),bins,ptbins,100,0,300);

      hjetetaphi_pix[nj][ip] = new TH2F(Form("hjetetaphi_pix%d_%d",nj,ip),Form("jet eta-phi distribution jet pix %d  %s",ip,calgo[nj]),36,-3.0,3.0,36,-pi,pi);
      hjetpteta_pix [nj][ip] = new TH2F(Form("hjetpteta_pix%d_%d",nj,ip),Form("jet pt-eta distribution jet pix %d %s",ip,calgo[nj]),bins,ptbins,18,-3.0,3.0);
      hjetptphi_pix [nj][ip] = new TH2F(Form("hjetptphi_pix%d_%d",nj,ip),Form("jet pt-phi distribution jet pix %d %s",ip,calgo[nj]),bins,ptbins,18,-pi,pi);

      hNevtpix[nj][ip] = new TH1F(Form("hNevtpix%d_%d",nj,ip),Form("# of events cent pix %d %s",ip,calgo[nj]),40,-40,40);

      hNevtpix_etab[nj][ip] = new TH1F(Form("hNevtpix_etab%d_%d",nj,ip),Form("# of events cent pix %d jet %s",ip,calgo[nj]),40,-40,40);
      for(int ie=0;ie<ketar;ie++){
	hjetrptpix_etab[nj][ip][ie] = new TH1F(Form("hjetrptpix_etab%d_%d_%d",nj,ip,ie),Form("jet(raw pt) p_{T} distribution  pix %d jet %s etabin %d",ip,calgo[nj],ie),bins,ptbins);
	hjetjptpix_etab[nj][ip][ie] = new TH1F(Form("hjetjptpix_etab%d_%d_%d",nj,ip,ie),Form("jet(jt pt) p_{T} distribution pix %d jet %s etabin %d",ip,calgo[nj],ie),bins,ptbins);
	hjetjpupix_etab[nj][ip][ie] = new TH1F(Form("hjetjpupix_etab%d_%d_%d",nj,ip,ie),Form("jet(pu) p_{T} distribution pix %d jet %s etabin %d",ip,calgo[nj],ie),100,0,300);
	hjetptpupix_etab[nj][ip][ie] = new TH2F(Form("hjetptpupix_etab%d_%d_%d",nj,ip,ie),Form("jet(pt:pu) distribution pix %d jet %s etabin %d",ip,calgo[nj],ie),bins,ptbins,100,0,300);
      }
    }

    for(int iv=0;iv<nvtx;iv++){
      hNevtvz[nj][iv] = new TH1F(Form("hNevtvz%d_%d",nj,iv),Form("# of events cent vz %d %s",iv,calgo[nj]),40,-40,40);
      hjetetavz[nj][iv] = new TH1F(Form("hjetetavz%d_%d",nj,iv),Form("jet eta distribution jet vz %s",calgo[nj]),18,-3.0,3.0);
      hjetphivz[nj][iv] = new TH1F(Form("hjetphivz%d_%d",nj,iv),Form("jet phi distribution jet vz %s",calgo[nj]),18,-pi,pi);
      hjetrptvz[nj][iv] = new TH1F(Form("hjetrptvz%d_%d",nj,iv),Form("jet(raw pt) p_{T} distribution jet vz %s",calgo[nj]),bins,ptbins);
      hjetjptvz[nj][iv] = new TH1F(Form("hjetjptvz%d_%d",nj,iv),Form("jet(jt pt) p_{T} distribution jet vz %s",calgo[nj]),bins,ptbins);
      hjetjpuvz[nj][iv] = new TH1F(Form("hjetjpuvz%d_%d",nj,iv),Form("jet(pu) p_{T} distribution jet vz %s",calgo[nj]),100,0,300);
      hjetptpuvz[nj][iv] = new TH2F(Form("hjetptpuvz%d_%d",nj,iv),Form("jet(pt:pu)  distribution jet vz %s",calgo[nj]),bins,ptbins,100,0,300);
      
      hjetetaphivz[nj][iv] = new TH2F(Form("hjetetaphivz%d_%d",nj,iv),Form("jet eta-phi distribution jet vz %s",calgo[nj]),36,-3.0,3.0,36,-pi,pi);
      hjetptetavz[nj][iv] = new TH2F(Form("hjetptetavz%d_%d",nj,iv),Form("jet pt-eta distribution jet vz %s",calgo[nj]),bins,ptbins,18,-3.0,3.0);
      hjetptphivz[nj][iv] = new TH2F(Form("hjetptphivz%d_%d",nj,iv),Form("jet pt-phi distribution jet vz %s",calgo[nj]),bins,ptbins,18,-pi,pi);
    }

    for(int icen=0;icen<ncen;icen++){
      hjeteta[nj][icen] = new TH1F(Form("hjeteta%d_%d",nj,icen),Form("jet eta distribution jet centb %d %s",icen,calgo[nj]),18,-3.0,3.0);
      hjetphi[nj][icen] = new TH1F(Form("hjetphi%d_%d",nj,icen),Form("jet phi distribution jet centb %d %s",icen,calgo[nj]),18,-pi,pi);

      hfrpt[nj][icen] = new TH1F(Form("hfrpt%d_%d",nj,icen),Form("jet(raw pt) p_{T} distribution jet centb %d %s",icen,calgo[nj]),500,0.,1000);
      hfjpt[nj][icen] = new TH1F(Form("hfjpt%d_%d",nj,icen),Form("jet(corrected pt) p_{T} distribution jet centb %d %s",icen,calgo[nj]),500,0.,1000);
      hfpupt[nj][icen] = new TH1F(Form("hfpupt%d_%d",nj,icen),Form("jet(pilup pt) p_{T} distribution jet centb %d %s",icen,calgo[nj]),250,0.,100);

      hjetrpt[nj][icen] = new TH1F(Form("hjetrpt%d_%d",nj,icen),Form("jet(raw pt) p_{T} distribution jet centb %d %s",icen,calgo[nj]),bins,ptbins);
      hjetjpt[nj][icen] = new TH1F(Form("hjetjpt%d_%d",nj,icen),Form("jet(jt pt) p_{T} distribution jet centb %d %s",icen,calgo[nj]),bins,ptbins);
      hjetspt[nj][icen] = new TH1F(Form("hjetspt%d_%d",nj,icen),Form("jet(smeared pt) p_{T} distribution jet centb %d %s",icen,calgo[nj]),bins,ptbins);
      hjetjpu[nj][icen] = new TH1F(Form("hjetjpu%d_%d",nj,icen),Form("jet(pu) p_{T} distribution jet centb %d %s",icen,calgo[nj]),100,0,300);
      hjetptpu[nj][icen] = new TH2F(Form("hjetptpu%d_%d",nj,icen),Form("jet(pt:pu) distribution jet centb %d %s",icen,calgo[nj]),bins,ptbins,100,0,300);

      hjetetaphi[nj][icen] = new TH2F(Form("hjetetaphi%d_%d",nj,icen),Form("jet eta-phi distribution jet centb %d %s",icen,calgo[nj]),18,-3.0,3.0,18,-pi,pi);
      hjetpteta[nj][icen] = new TH2F(Form("hjetpteta%d_%d",nj,icen),Form("jet pt-eta distribution jet centb %d %s",icen,calgo[nj]),bins,ptbins,18,-3.0,3.0);
      hjetptphi[nj][icen] = new TH2F(Form("hjetptphi%d_%d",nj,icen),Form("jet pt-phi distribution jet centb %d %s",icen,calgo[nj]),bins,ptbins,18,-pi,pi);

      hNevt_etab[nj][icen] = new TH1F(Form("hNevt_etab%d_%d",nj,icen),Form("# of events cent %d jet %s",icen,calgo[nj]),40,-40,40);
      for(int ie=0;ie<ketar;ie++){
	hjetrpt_etab[nj][icen][ie] = new TH1F(Form("hjetrpt_etab%d_%d_%d",nj,icen,ie),Form("jet(raw pt) p_{T} distribution centb %d jet %s etabin %d",icen,calgo[nj],ie),bins,ptbins);
	hjetjpt_etab[nj][icen][ie] = new TH1F(Form("hjetjpt_etab%d_%d_%d",nj,icen,ie),Form("jet(jt pt) p_{T} distribution centb %d jet %s etabin %d",icen,calgo[nj],ie),bins,ptbins);
	hjetjpu_etab[nj][icen][ie] = new TH1F(Form("hjetjpu_etab%d_%d_%d",nj,icen,ie),Form("jet(pu) p_{T} distribution centb %d jet %s etabin %d",icen,calgo[nj],ie),100,0,300);
	hjetptpu_etab[nj][icen][ie] = new TH2F(Form("hjetptpu_etab%d_%d_%d",nj,icen,ie),Form("jet(pt:pu) distribution centb %d jet %s etabin %d",icen,calgo[nj],ie),bins,ptbins,100,0,300);
      }

      hNevt[nj][icen] = new TH1F(Form("hNevt%d_%d",nj,icen),Form("# of events cent %d %s",icen,calgo[nj]),40,-40,40);

      hNevt_notr[nj][icen] = new TH1F(Form("hNevt_notr%d_%d",nj,icen),Form("# of events w/o trigger conditioncent %d %d",nj,icen),40,-40,40);

      hNjets[nj][icen] = new TH1F(Form("hNjets%d_%d",nj,icen),Form("# of cent %d jets %s",icen,calgo[nj]),3,0.0,3.0);
      
      hNjets_dist[nj][icen] = new TH1F(Form("hNjets_dist%d_%d",nj,icen),Form("distribution jets in cent %d %s",icen,calgo[nj]),100,0,100);
    }
    //std::cout<<Form("Initialized the histograms %s",calgo[nj]) <<std::endl;
    std::cout<<"\t"<<std::endl;
  }
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
  int istat[knj];
  int stat_etab[knj];

  double totjetpt[knj][maxruns];
  double totjets [knj][maxruns];
  
  for(int nj=0;nj<knj;nj++){
    for(int ir=0;ir<maxruns;ir++){
      totjetpt[nj][ir]=0;
      totjets[nj][ir]=0;
    }
  }


  Long64_t nb = 0;
  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;
  hEvt->Fill(2,nentries);
  
  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
    //for (Long64_t ievt=0; ievt<200;ievt++) {//! event loop
    //! load the hiForest event
    nb += c->GetEntry(ievt);

    //! Remove duplicate events
    if(!ispp && dupEvt.occurrence[ievt]==2)continue;



    //! select good events
    bool goodEvent = false;
    if(ispp){
      goodEvent   = c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_Jet60_v1 && fabs(c->evt.vz)<15.;
    }else{
      goodEvent = c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_HIJet80_v1 && fabs(c->evt.vz)<15.; ///! this one to use
    }
    if(!goodEvent)continue;

    //! Event variables
    int RunNo=-999, EvtNo=-999, hiBin=-999, LumiBlock=-999, hiNpix =-999, hiNpixelTracks=-999,hiNtracks=-999;
    float vx=-999,vy=-999,vz=-999, hiHF=-999,hiHFplus=-999,hiHFminus=-999,hiHFhit=-999,hiHFhitplus=-999,
      hiHFhitminus=-999,hiZDC=-999,hiZDCplus=-999,hiZDCminus=-999;
    float hiET=-999,hiEE=-999,hiEB=-999,hiEEplus=-999,hiEEminus=-999;
    
    
    RunNo          = c->evt.run;
    EvtNo          = c->evt.evt;
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

    //! Vertex z bin
    //! Only two bins v<0 : 0 and vz>0 : 1
    int vbin = GetVtxBin(vz);
    //if(vbin<0)continue;

    //! hiPixelTracks bins
    //! 5 bins
    int ipix = GetPixelBin(hiNpixelTracks);
    //if(ipix<0)continue;

    if(ievt%10000==0 && !iCountRuns)std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<RunNo<<"\t HF : "<<hiHF<<std::endl;
    //std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<RunNo<<"\t HF : "<<hiHF<<std::endl;            

    //! Jet Algo loop starts here
    //! nj : 0 (icpu5calo) and 1 (akpu3pf)
    for(int nj=0;nj<knj;nj++){

      //! Initialization
      nmbj[nj]=0;
      sumjpt[nj]=0;
      sumjpu[nj]=0;
      istat[nj]=0;
      stat_etab[nj]=0;
      for(int icen=0;icen<ncen;icen++){
	njets[nj][icen]=0;
      }
      
      //! assign the jets for PbPb
      if(strcmp(ksp,"pp")==0){
	if(nj==0)mJets = &(c->ak2PF);
	else if(nj==1)mJets = &(c->ak3PF);
	else if(nj==2)mJets = &(c->ak4PF);
      }else{
	if(nj==0)mJets = &(c->akPu2PF);
	else if(nj==1)mJets = &(c->akPu3PF);
	else if(nj==2)mJets = &(c->akPu4PF);
      }

      /*if(nj==0)mJets = &(c->icPu5);
	else if(nj==1)mJets = &(c->akPu2PF);
	else if(nj==2)mJets = &(c->akPu3PF);
	else if(nj==3)mJets = &(c->akPu4PF);
	else if(nj==4)mJets = &(c->akPu2Calo);
	else if(nj==5)mJets = &(c->akPu3Calo);
	else if(nj==6)mJets = &(c->akPu4Calo);
      */

      //if(nj==1)std::cout<<"  \t \t  "<<Form("\t %s  : ",calgo[nj])<<mJets->nref<<"\t centb : "<<centb<<std::endl;
      //! Loop over # of jets in a given algo
      for(int ij=0; ij<mJets->nref; ij++){

	//! Without any cut how many jets are there
	hTotJets[nj]->Fill(1.);

	//! jet selction cut
	if(mJets->jtpt[ij]<kptrecocut || fabs(mJets->jtphi[ij])>pi)continue;

	//! trackMax cuts
	if(mJets->trackMax[ij]/mJets->jtpt[ij] < 0.01)continue;

	//! jets within different etabins
	//! 
	int ebin = GetEtaBin(fabs(mJets->jteta[ij]));
	if(ebin<0)continue;
	stat_etab[nj]=1;
	hjetrptmb_etab [nj][ebin]->Fill(mJets->rawpt[ij]);
	hjetjptmb_etab [nj][ebin]->Fill(mJets->jtpt[ij]);
	hjetjpumb_etab [nj][ebin]->Fill(mJets->jtpu[ij]);
	hjetptpumb_etab[nj][ebin]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);

	if(ipix>=0){
	  hjetrptpix_etab [nj][ipix][ebin]->Fill(mJets->rawpt[ij]);
	  hjetjptpix_etab [nj][ipix][ebin]->Fill(mJets->jtpt[ij]);
	  hjetjpupix_etab [nj][ipix][ebin]->Fill(mJets->jtpu[ij]);
	  hjetptpupix_etab[nj][ipix][ebin]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);
	}
	hjetrpt_etab [nj][centb][ebin]->Fill(mJets->rawpt[ij]);
	hjetjpt_etab [nj][centb][ebin]->Fill(mJets->jtpt[ij]);
	hjetjpu_etab [nj][centb][ebin]->Fill(mJets->jtpu[ij]);
	hjetptpu_etab[nj][centb][ebin]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);
	

	if(fabs(mJets->jteta[ij]) < ketacut){
	  istat[nj]=1;
	  
	  //! minbias
	  hjetrptmb [nj]->Fill(mJets->rawpt[ij]);
	  hjetjptmb [nj]->Fill(mJets->jtpt[ij]);
	  hjetjpumb [nj]->Fill(mJets->jtpu[ij]);
	  hjetptpumb[nj]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);
	  
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
	  //! pixel tracks dependence
	  hjetrpt_pix [nj][ipix]->Fill(mJets->rawpt[ij]);
	  hjetjpt_pix [nj][ipix]->Fill(mJets->jtpt[ij]);
	  hjetjpu_pix [nj][ipix]->Fill(mJets->jtpu[ij]);
	  hjetptpu_pix[nj][ipix]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);
	  
	  hjeteta_pix   [nj][ipix]->Fill(mJets->jteta[ij]);
	  hjetphi_pix   [nj][ipix]->Fill(mJets->jtphi[ij]);
	  hjetpteta_pix [nj][ipix]->Fill(mJets->jtpt[ij],mJets->jteta[ij]);
	  hjetptphi_pix [nj][ipix]->Fill(mJets->jtpt[ij],mJets->jtphi[ij]);
	  hjetetaphi_pix[nj][ipix]->Fill(mJets->jteta[ij],mJets->jtphi[ij]);
	  
	  //! centrality
	  hNjets [nj][centb]->Fill(1.);

	  hfrpt  [nj][centb]->Fill(mJets->rawpt[ij]);
	  hfjpt  [nj][centb]->Fill(mJets->jtpt[ij]);
	  hfpupt [nj][centb]->Fill(mJets->jtpu[ij]);


	  hjetrpt [nj][centb]->Fill(mJets->rawpt[ij]);
	  hjetjpt [nj][centb]->Fill(mJets->jtpt[ij]);
	  hjetjpu [nj][centb]->Fill(mJets->jtpu[ij]);
	  hjetptpu[nj][centb]->Fill(mJets->jtpt[ij],mJets->jtpu[ij]);

	  //! pp pt smeared with PbPb resolution
	  if(strcmp(ksp,"pp")==0){
	    int smbin = GetPtBinSmear(mJets->jtpt[ij]);
	    if(smbin>=0){
	      for(int ic=0;ic<ncen-1;ic++){
		float smpt = gRandom->Gaus(mJets->jtpt[ij],smearf[nj][ic][smbin]);
		if(smpt<kptrecocut)continue;
		hjetspt [nj][ic]->Fill(smpt);
	      }
	    }
	  }

	  hjeteta   [nj][centb]->Fill(mJets->jteta[ij]);
	  hjetphi   [nj][centb]->Fill(mJets->jtphi[ij]);
	  hjetpteta [nj][centb]->Fill(mJets->jtpt[ij],mJets->jteta[ij]);
	  hjetptphi [nj][centb]->Fill(mJets->jtpt[ij],mJets->jtphi[ij]);
	  hjetetaphi[nj][centb]->Fill(mJets->jteta[ij],mJets->jtphi[ij]);
	  
	  //! For a given event	  
	  njets[nj][centb]++;
	  sumjpt[nj]+=mJets->jtpt[ij];
	  sumjpu[nj]+=mJets->jtpu[ij];
	  nmbj[nj]++;

	  //! For a given run
	  totjetpt[nj][RunId] += mJets->jtpt[ij];
	  totjets [nj][RunId]++;

	} //! etacut
	//if(centb==4)std::cout<<Form(" \t \t \t mbjets %s : ",calgo[nj]) <<nmbj[nj]<<"\t centb : "<<centb<<"\t jtpt "<<mJets->jtpt[ij]<<std::endl;
      }//! # jet loop

      hNref[nj]->Fill(mJets->nref);
      hNevt_notr[nj][centb]->Fill(vz);

      if(stat_etab[nj]){
	hNevtmb_etab[nj]->Fill(vz);
	hNevtpix_etab[nj][ipix]->Fill(vz);
	hNevt_etab[nj][centb]->Fill(vz);
      }

      if(istat[nj]){
	hNevt[nj][centb]->Fill(vz);
	if(vbin>=0)hNevtvz[nj][vbin]->Fill(vz);
	if(ipix>=0)hNevtpix[nj][ipix]->Fill(vz);
	hMBJets[nj]->Fill(nmbj[nj]);

	hHF_tr[nj]->Fill(hiHF);

	//! # of jets distribution
	hNjets_dist[nj][centb]->Fill(njets[nj][centb]);
	hMBJetsRun [nj]->Fill(RunId,nmbj[nj]);
	hJetPtRun  [nj]->Fill(RunId,sumjpt[nj]/nmbj[nj]);
	hJetPURun  [nj]->Fill(RunId,sumjpu[nj]/nmbj[nj]);

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
    hHF->Fill(hiHF);
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
    hNTracks->Fill(hiNtracks);
    hNpixelTracksNTracks->Fill(hiNpixelTracks,hiNtracks);

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
	for(int ic=0;ic<ncen-1;ic++){          
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
	if(Run[ir]<0 || totjets[nj][ir]==0)continue;
	std::cout<<"\t RunNo : "<<Run[ir]<<"\t Sumjetpt : "<<totjetpt[nj][ir]<<"\t Totjets : "<<totjets[nj][ir]<<"\t  <pT> : "<<totjetpt[nj][ir]/totjets[nj][ir]<<std::endl;
	hTotJetsRun[nj]->Fill(ir+1,totjets[nj][ir]);
      }
    }
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
}
int GetPtBinSmear(float pt)
{
  int ibin=-1;
  for(int ix=0;ix<smbins;ix++){
    if(pt>=ptbins_smear[ix] && pt<ptbins_smear[ix+1]){
      return ix;
    }
  }
  return ibin;
}
int GetCentBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 2.5% bins of cross section
  //! in 0-39 bins
  if(bin<2)ibin=0; //! 0-5%
  else if(bin>=2 && bin<4)ibin=1;   //! 5-10%
  else if(bin>=4 && bin<12)ibin=2;  //! 10-30%
  else if(bin>=12&& bin<20)ibin=3;  //! 30-50%
  else if(bin>=20&& bin<28)ibin=4;  //! 50-70%
  else if(bin>=28&& bin<36)ibin=5;  //! 70-90%
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
int GetPixelBin(int npixtr)
{
  int ibin=-1;
  if(npixtr>=700)ibin=0;
  else if(npixtr>=200 && npixtr<700)ibin=1;
  else if(npixtr>=100 && npixtr<200)ibin=2;
  else if(npixtr>=50  && npixtr<100)ibin=3;
  else if(npixtr<50)ibin=4;

  return ibin;
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
  }
  std::cout<<"Run # loaded according to increasing order" <<std::endl;
}
