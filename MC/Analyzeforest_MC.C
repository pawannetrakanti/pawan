#include "/net/hidsk0001/d00/scratch/pawan/service/CMSSW_6_0_0/src/test/pp/HiForest/V3/hiForest.h"
#include <TF1.h>
#include <TProfile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;


#ifdef _MAKECINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class HiForest;
#pragma link C++ class Evts;
#pragma link C++ class Jets;
#endif


int GetEtaBin(float /*eta*/);
int GetDetEtaBin(float /*eta*/);
int GetPhiBin  (float /*phi*/);
int GetCentBin(int /*hiBin*/);
int GetHFBin   (float /*hiHFsumEta4*/);
int GetPtBin(float /*pt*/);
bool selecJet(Jets */*jetc*/, int /*ij*/);
void FindLeadSubLeadJets(const char */*sp*/,Jets */*jetc*/, int */*ljet*/);
bool isHiPtTrkInJet(const char */*sp*/,Jets */*jetc*/, Tracks */*trks*/,float /*delr*/);
void ShutoffBranches(HiForest */*hi*/,const char */*ksp*/);
double delphi(double /*phi1*/, double /*phi2*/);
void GetXsection(const char */*sp*/,double /*sqrts*/,double /*pthat*/,double &/*xsection*/,double &/*maxpthat*/);

//! pt binning
double ptbins[] ={10,20,30,50,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,280,300,320,340,360,400,450,500,550,600,835};
const int bins = sizeof(ptbins)/sizeof(double) - 1;

//! data pt binning
//double ptbins_data[] ={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
//const int dbins = sizeof(ptbins_data)/sizeof(double) - 1;


//! constants
int iYear=2013;
const double pi=acos(-1.);
const double pi2=2*pi -1;

const int ncent=6; //! for pp fill in 0, for ppb 0-4 and pbpb 0-5
const int ketar=4;
const int maxe =5;
const int maxph=5;
const float ketacut=2.0;
const double kptrecocut=15.;
const double kptgencut =0.;
const double kdRcut=0.3;
const double kVzcut=15.;

const int rbins=50;
double rbinl=0.0;
double rbinh=5.0;


const int knj = 16;
const char *calgo[knj] = {"ak2PF","ak3PF","ak4PF","ak5PF",
			  "akPu2PF","akPu3PF","akPu4PF","akPu5PF",
			  "ak2Calo","ak3Calo","ak4Calo","ak5Calo",
			  "akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo"};
const float kDelR[knj]  = {0.2,0.3,0.4,0.5,
                           0.2,0.3,0.4,0.5,
                           0.2,0.3,0.4,0.5,
                           0.2,0.3,0.4,0.5
};
//! branch selection for jet trees
const char *kbr[] ={"pthat","nref",
		    "refpt","refphi","refeta","refdrjt",
		    "subid","rawpt","jtpt","jteta","jtphi",
		    "trackMax","jtpu","photonSum","chargedSum","neutralSum",
		    "ngen","ngenmatchindex"
};
int nsize = sizeof(kbr)/sizeof(char*);

TStopwatch timer;
int Analyzeforest_MC(const char *ksp="pp",double kPt=80,double kPtCut=30)
{

  timer.Start();

  //! Load Lib
  gSystem->Load("/net/hidsk0001/d00/scratch/pawan/service/CMSSW_6_0_0/src/test/pp/HiForest/V3/hiForest_h.so");
  
  //! Load the hiforest input root file
  HiForest *c = 0;
  TString inname="";
  double sqrts=0;
  if(strcmp(ksp,"ppb")==0){
    c = new HiForest(Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt%0.0f/HiForest_v77_merged01/pt%0.0f_HP04_prod16_v77_merged_forest_0.root",kPt,kPt),"PyHI2012",cPPb);
    sqrts=5.02;
  }else if(strcmp(ksp,"pp")==0){
    //! pp 2.76 TeV  2013
    //! pptracking
    inname=Form("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt%0.0f/HiForest_v81_merged01/pt%0.0f_pp2013_P01_prod22_v81_merged_forest_0.root",kPt,kPt);
    //inname=Form("/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged/pp276Dijet%0.0f_merged.root",kPt);

    //! pp 2.76 TeV 2011
    //if(iYear==2011)inname=Form("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt%0.0f/HiForest_v81_merged01/pt%0.0f_pp2013_P01_prod22_v81_merged_forest_0.root",kPt,kPt);
    //! HI tracking
    //TString inname=Form("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/Signal_Pythia_pt%0.0f/HiForest_v81_merged01/pt%0.0f_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root",kPt,kPt);
    c=new HiForest(inname,Form("Py%0.0f",kPt),cPP);
    sqrts=2.76;
    iYear=2013;
  }else if(strcmp(ksp,"pbpb")==0){
    if(kPt==50 || kPt==80 || kPt==100 || kPt==170)inname=Form("/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet%0.0f_HydjetDrum_v27_mergedV1.root",kPt); 
    else if(kPt==120 || kPt==200 || kPt==250 || kPt==300)inname=Form("/mnt/hadoop/cms/store/user/yenjie/HiForest_v28/Dijet%0.0f_HydjetDrum_v28_mergedV1.root",kPt); 
    c=new HiForest(inname,Form("HI%0.0f",kPt),cPbPb);
    sqrts=2.76;
    iYear=2011;
  }


  c->hasHltTree=0;
  c->hasSkimTree=0;
  c->hasIcPu5JetTree=0;
  c->hasAk1PFJetTree=0;
  c->hasAk6PFJetTree=0;
  c->hasAkPu1PFJetTree=0;
  c->hasAkPu6PFJetTree=0;
  c->hasAk1CaloJetTree=0;
  c->hasAk6CaloJetTree=0;
  c->hasAkPu1CaloJetTree=0;
  c->hasAkPu6CaloJetTree=0;
  c->hasTrackTree=0;

  //! Shut off unwanted branches in jet tree
  //ShutoffBranches(c,ksp);

  //! Get the x-section re-weighting factor
  double xsection=1;
  double maxpthat=9999;
  GetXsection(ksp,sqrts,kPt,xsection,maxpthat);


  std::cout<<std::endl;
  std::cout<<std::endl;
  

  //! To get the jet object from hiforest
  Jets *iJet=0;

  std::cout<<"Loaded all tree variables and # of jet algorithms : "<<knj<<std::endl;
  std::cout<<"\t"<<std::endl;

  //! Open a output file for histos
  //TFile *fout = new TFile(Form("/net/hidsk0001/d00/scratch/pawan/service/CMSSW_6_0_0/src/test/pp/output/jec/MC_%d_%s_pT%0.0fGeV.root",iYear,ksp,kPt),"RECREATE");
  TFile *fout = new TFile(Form("/net/hidsk0001/d00/scratch/pawan/service/CMSSW_6_0_0/src/test/pp/output/jec/MC_%d_v81_prod22_%s_pT%0.0fGeV.root",iYear,ksp,kPt),"RECREATE");

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Running for %s ",ksp)<<std::endl;
  //  std::cout<<Form("pT  cut for %0.3f ",kptrecocut)<<std::endl;
  std::cout<<Form("pT  cut for %0.3f ",kPtCut)<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"My hiForest TTree : " <<c->GetName()<<std::endl;
  std::cout<<"Output file  "<<fout->GetName()<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;


  //! 
  //! Define histograms here
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  TProfile::SetDefaultSumw2();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *hEvt    = new TH1F("hEvt","# of events ",4,0,4);
  TH1F *hVz     = new TH1F("hVz","# of events ",80,-20,20);
  TH1F *hBin    = new TH1F("hBin","Centrality bin",100,-0.5,100-0.5);
  TH1F *hTotEve = new TH1F("hTotEve","# of events in the skimmed files",4,0,4);
  TH1F *hpthat  = new TH1F("hpthat","py-hat distribution",500,0,1000);

  TH1F *hHFPlusEta4 = new TH1F("hHFPlusEta4","HFPlusEta4 distribution",200,0,100);
  TH1F *hHF       = new TH1F("hHF","HF distribution",1000,0,1000);
  TH1F *hNTracks  = new TH1F("hNTracks","hiNTracks",1000,0,1000);
  TH2F *hNTracksHF = new TH2F("hNTracksHF","Ntrack vs HF ",1000,0,1000,1000,0,1000);
  TH2F *hBinHF      = new TH2F("hBinHF","HF vs hiBin",100,-0.5,100-0.5,500,0,500);
  TH2F *hBinNTracks = new TH2F("hBinNTracks","HF vs NTracks",100,-0.5,100-0.5,500,0,500);
  TH2F *hNTracksHFPlusEta4 = new TH2F("hNTracksHFPlusEta4","Ntrack vs HFPlusEta4 ",1000,0,1000,800,0,200);

  TH2F *hrajptrpt [knj][ncent];

  //! Gen matched jets
  TH1F *hgenpt_genm [knj][ncent], *hrecopt_genm[knj][ncent], *hrawpt_genm[knj][ncent];
  TH1F *hjeteta     [knj][ncent], *hjetphi     [knj][ncent];
  TH2F *hjetpteta   [knj][ncent], *hjetptphi   [knj][ncent], *hjetetaphi [knj][ncent];

  //! Ratios of the pt distributions
  TProfile *hrecogen[knj][ncent], *hrecoraw[knj][ncent], *hrawgen[knj][ncent];

  //! Resposnse
  TH2F *hrescrpt_genm [knj][ncent], *hresrrpt_genm [knj][ncent], *hresrcrpt_genm [knj][ncent];
  TH3F *hrescreta_genm[knj][ncent], *hresrreta_genm[knj][ncent], *hresrcreta_genm[knj][ncent];
  TH2F *hratiorawrefpt_eta[knj][ncent][2], *hratiocorrrefpt_eta[knj][ncent][2];

  TH2F *hratiocorrrefpt_genm[knj][ncent];
  TH2F *hratiocorrrefpt_lead[knj][ncent], *hratiocorrrefpt_slead[knj][ncent];

  //! For comparison with data
  TH2F *hJetEnergyScale[knj][ncent];

  TH2F *hpteta[knj][ncent][maxe] ;
  TH2F *hptphi[knj][ncent][maxph] ;


  TH2F *hgenjrecoj[knj][ncent];
  TH1F *hNjets_genm[knj][ncent];
  TH1F *hNevt [knj][ncent];

  //! Background for jets 
  //! (photonSum+neutralSum+chargedSum-rawpt)
  TH2F *hjetptpu_genm[knj][ncent];             //! centrality                                                                       
  TH2F *hjetptpu_etab_genm[knj][ncent][ketar]; //! eta dependence             
  TH1F *hjetbkgd_genm[knj][ncent];
  TH2F *hjetptbkgd_genm[knj][ncent];
  TH2F *hjetptbkgd_etab_genm[knj][ncent][ketar];
  TH2F *hPFFraction_genm[knj][ncent][3];

  //! Efficency histos
  TH1F *hPtAll [knj][ncent], *hPtSel[knj][ncent];
  TH1F *hEtaAll[knj][ncent], *hEtaSel[knj][ncent];
  TH1F *hPhiAll[knj][ncent], *hPhiSel[knj][ncent];

  //! Response vs deltar
  TH1F *hRspVsDeltaR[knj][ncent][25];

  //! DeltaR efficiency                                                                                                                                                              
  TH1F *hDeltaR[knj][ncent];
  TH1F *hDeltaRAll[knj][ncent];
  TH1F *hDeltaRSel[knj][ncent];


  TH3F *hjbkgd[knj][ncent]; 


  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncent;ic++){

      hjbkgd[nj][ic] = new TH3F(Form("hjbkgd%d_%d",nj,ic),Form("hjbkgd%d_%d",nj,ic),12,-ketacut,ketacut,50,0,500,40,0,80);

      hNevt [nj][ic] = new TH1F(Form("hNevt%d_%d",nj,ic),Form("# of events cent %d %s",ic,calgo[nj]),40,-40,40);

      hNjets_genm [nj][ic] = new TH1F(Form("hNjets_genm%d_%d",nj,ic),Form("# of jets cent %d jets %s",ic,calgo[nj]),200,30,630);
      hgenpt_genm [nj][ic] = new TH1F(Form("hgenpt_genm%d_%d",nj,ic),Form("Gen matched gen p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hrecopt_genm[nj][ic] = new TH1F(Form("hrecopt_genm%d_%d",nj,ic),Form("Gen matched reco p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hrawpt_genm [nj][ic] = new TH1F(Form("hrawpt_genm%d_%d",nj,ic),Form("Gen matched raw p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);

      //! Ratios
      hrecogen[nj][ic] = new TProfile(Form("hrecogen%d_%d",nj,ic),Form("reco/gen p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hrecoraw[nj][ic] = new TProfile(Form("hrecoraw%d_%d",nj,ic),Form("reco/raw p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hrawgen[nj][ic]  = new TProfile(Form("hrawgen%d_%d",nj,ic),Form("raw/gen p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      
      //! Gen matched Response and resolution
      hrescrpt_genm[nj][ic]= new TH2F(Form("hrescrpt_genm%d_%d",nj,ic),Form("Gen jet:(Reco/Gen) jet p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,150,rbinl,rbinh);
      hresrrpt_genm[nj][ic]= new TH2F(Form("hresrrpt_genm%d_%d",nj,ic),Form("Gen jet:(Raw/Gen) jet p_{T}  distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,150,rbinl,rbinh);
      hresrcrpt_genm[nj][ic]= new TH2F(Form("hresrcrpt_genm%d_%d",nj,ic),Form("Reco jet:(Reco/Raw) jet p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,150,rbinl,rbinh);

      //! Gen matched response and resolutions in different eta bins
      hrescreta_genm[nj][ic]  = new TH3F(Form("hrescreta_genm%d_%d",nj,ic),Form("Gen jet:(Reco/Gen) jet p_{T} distribution jet centb %d %s",ic,calgo[nj])  ,12,-3.0,3.0,500,0,1000,rbins,rbinl,rbinh);
      hresrreta_genm[nj][ic]  = new TH3F(Form("hresrreta_genm%d_%d",nj,ic),Form("Gen jet:(Raw/Gen) jet p_{T}  distribution jet centb %d %s",ic,calgo[nj])  ,12,-3.0,3.0,500,0,1000,rbins,rbinl,rbinh);
      hresrcreta_genm[nj][ic] = new TH3F(Form("hresrcreta_genm%d_%d",nj,ic),Form("Reco jet:(Reco/Raw) jet p_{T} distribution jet centb %d %s",ic,calgo[nj]),12,-3.0,3.0,500,0,1000,rbins,rbinl,rbinh);

      hjeteta[nj][ic] = new TH1F(Form("hjeteta%d_%d",nj,ic),Form("jet eta distribution jet centb %d %s",ic,calgo[nj]),72,-ketacut,ketacut);
      hjetphi[nj][ic] = new TH1F(Form("hjetphi%d_%d",nj,ic),Form("jet phi distribution jet centb %d %s",ic,calgo[nj]),72,-pi,pi);

      hjetetaphi[nj][ic] = new TH2F(Form("hjetetaphi%d_%d",nj,ic),Form("jet eta-phi distribution jet centb %d %s",ic,calgo[nj]),72,-ketacut,ketacut,72,-pi,pi);
      hjetpteta[nj][ic] = new TH2F(Form("hjetpteta%d_%d",nj,ic),Form("jet pt-eta distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,72,-ketacut,ketacut);
      hjetptphi[nj][ic] = new TH2F(Form("hjetptphi%d_%d",nj,ic),Form("jet pt-phi distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,72,-pi,pi);


      hratiocorrrefpt_lead[nj][ic]= new TH2F(Form("hratiocorrrefpt_lead%d_%d",nj,ic),Form("Leading jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",ic,calgo[nj]),
                                               500,0,1000,rbins,rbinl,rbinh);
      hratiocorrrefpt_slead[nj][ic]= new TH2F(Form("hratiocorrrefpt_slead%d_%d",nj,ic),Form("sub-Leading jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",ic,calgo[nj]),
                                                500,0,1000,rbins,rbinl,rbinh);
      hratiocorrrefpt_genm[nj][ic]= new TH2F(Form("hratiocorrrefpt_genm%d_%d",nj,ic),Form("Gen matched jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",ic,calgo[nj]),
                                               500,0,1000,rbins,rbinl,rbinh);
      hJetEnergyScale[nj][ic] = new TH2F(Form("hJetEnergyScale%d_%d",nj,ic),Form("hJetEnergyScale%d_%d",nj,ic),500,0,1000,50,-1.00,1.00);

      for(int ie=0;ie<2;ie++){
        hratiorawrefpt_eta[nj][ic][ie]= new TH2F(Form("hratiorawrefpt_eta%d_%d_%d",nj,ic,ie),
						   Form("Raw jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",ic,calgo[nj],ie),
						   500,0,1000,rbins,rbinl,rbinh);
        hratiocorrrefpt_eta[nj][ic][ie]= new TH2F(Form("hratiocorrrefpt_eta%d_%d_%d",nj,ic,ie),
						    Form("Reco jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",ic,calgo[nj],ie),
						    500,0,1000,rbins,rbinl,rbinh);
      }
      
      //! Jet background
      hjetptpu_genm[nj][ic] = new TH2F(Form("hjetptpu_genm%d_%d",nj,ic),Form("jet(pt:pu) distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,100,0,300);
      for(int ie=0;ie<ketar;ie++){
        hjetptpu_etab_genm[nj][ic][ie] = new TH2F(Form("hjetptpu_etab_genm%d_%d_%d",nj,ic,ie),Form("jet(pt:pu) distribution jet %s etabin %d cen %d",calgo[nj],ic,ie),
						    500,0,1000,100,0,300);
	hjetptbkgd_etab_genm[nj][ic][ie] = new TH2F(Form("hjetptbkgd_etab_genm%d_%d_%d",nj,ic,ie),Form("jet(pt:bkgd) distribution jet etabin %d centb %d %s",ie,ic,calgo[nj]),500,0,1000,100,0,300);
      }
      hjetbkgd_genm  [nj][ic] = new TH1F(Form("hjetbkgd_genm%d_%d",nj,ic),Form("jet(pu) p_{T} distribution jet centb %d %s",ic,calgo[nj]),100,0,300);
      hjetptbkgd_genm[nj][ic] = new TH2F(Form("hjetptbkgd_genm%d_%d",nj,ic),Form("jet(pt:bkgd) distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,100,0,300);

      //! PF fractions
      for(int ipf=0;ipf<3;ipf++){
	hPFFraction_genm[nj][ic][ipf] = new TH2F(Form("hPFFraction_genm%d_%d_%d",nj,ic,ipf),Form("PF fraction distribution jet centb %d %s %d",ic,calgo[nj],ipf),500,0,1000,15,0,1.5);
      }

      
      for(int m=0;m<maxe;m++){
        hpteta[nj][ic][m] = new TH2F(Form("hpteta%d_%d_%d",nj,ic,m),Form("resolution  pt(eta) distribution cent %d jet %s etabin%d",ic,calgo[nj],m),
				       500,0,1000,rbins,rbinl,rbinh);
      }
      
      for(int m=0;m<maxph;m++){
	hptphi[nj][ic][m] = new TH2F(Form("hptphi%d_%d_%d",nj,ic,m),Form("resolution pt(phi) distribution cent %d jet %s phibin%d",ic,calgo[nj],m),
				       500,0,1000,rbins,rbinl,rbinh);
      }

      hgenjrecoj[nj][ic] =new TH2F(Form("hgenjrecoj%d_%d",nj,ic),Form("gen jet2 : reco jet2 %s cent %d",calgo[nj],ic),500,0.,1000.,500,0.,1000.);


      hrajptrpt[nj][ic] = new TH2F(Form("hrajptrpt%d_%d",nj,ic),Form("corr pT / jet(raw pt) p_{T} distribution %d jet %s",ic,calgo[nj]),500,0,1000,50,0,10);

      //! efficcy histograms
      hPtAll [nj][ic] = new TH1F(Form("hPtAll%d_%d",nj,ic),Form("Denominator pT for algorithm %s cent %d",calgo[nj],ic),40,10,110);
      hEtaAll[nj][ic] = new TH1F(Form("hEtaAll%d_%d",nj,ic),Form("Denominator eta  for algorithm %s cent %d",calgo[nj],ic),20,-ketacut,ketacut);
      hPhiAll[nj][ic] = new TH1F(Form("hPhiAll%d_%d",nj,ic),Form("Denominator  phi  for algorithm %s cent %d",calgo[nj],ic),20,-pi,pi);
      
      hPtSel [nj][ic] = new TH1F(Form("hPtSel%d_%d",nj,ic),Form("Numerator pT for algorithm %s cent %d",calgo[nj],ic),40,10,110);
      hEtaSel[nj][ic] = new TH1F(Form("hEtaSel%d_%d",nj,ic),Form("Numerator eta  for algorithm %s cent %d",calgo[nj],ic),20,-ketacut,ketacut);
      hPhiSel[nj][ic] = new TH1F(Form("hPhiSel%d_%d",nj,ic),Form("Numerator  phi  for algorithm %s cent %d",calgo[nj],ic),20,-pi,pi);

      hDeltaR[nj][ic]    = new TH1F(Form("hDeltaR%d_%d",nj,ic),Form("#DeltaR for algorithm %s cent %d",calgo[nj],ic),100,0,1);
      hDeltaRAll[nj][ic] = new TH1F(Form("hDeltaRAll%d_%d",nj,ic),Form("#DeltaR (all) for algorithm %s cent %d",calgo[nj],ic),100,0,1);
      hDeltaRSel[nj][ic] = new TH1F(Form("hDeltaRSel%d_%d",nj,ic),Form("#DeltaR (sel) for algorithm %s cent %d",calgo[nj],ic),100,0,1);

      for(int ir=0;ir<25;ir++){
	//! Response vs DeltaR
        hRspVsDeltaR[nj][ic][ir] = new TH1F(Form("hRspVsDeltaR%d_%d_%d",nj,ic,ir),Form(" <recopt/refpt> vs. #DeltaR (%d) algorithm %s cent %d",ir,calgo[nj],ic),rbins,rbinl,rbinh);
      }
    }//! ic
  }//! nj
  std::cout<<"Initialized the histograms " <<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////// 


  //! Centrality reweighting function

  //! vertex z reweighting


  Long64_t nentries = c->GetEntries();
  hEvt->Fill(2,nentries);

  //! weight  for the merging of the samples for different pT hat bins
  double wxs = xsection/(nentries/100000.);

  std::cout<<Form("# of entries in TTree for %s  : ",ksp)<<nentries<<"\t and x-section  : "<<xsection<<std::endl;
  std::cout<<std::endl;
  
  Int_t iEvent=0; 
  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
    //for (Long64_t ievt=0; ievt<5;ievt++) {//! event loop
    //! load the hiForest event
    c->GetEntry(ievt);

    int hiBin       = c->evt.hiBin;
    float vz        = c->evt.vz;
    float hiHF      = c->evt.hiHF;
    int   ntracks   = c->evt.hiNtracks;
    
    float hiHFPlusEta4=-999;
    float hiHFMinusEta4=-999;
    float hiHFSum=-999;
    if(strcmp(ksp,"ppb")==0){
      hiHFPlusEta4 = c->evt.hiHFplusEta4; 
      hiHFMinusEta4 = c->evt.hiHFminusEta4; 
      hiHFSum = hiHFPlusEta4 + hiHFMinusEta4 ;
    }

    //! apply vertex cut
    if(fabs(vz)>kVzcut)continue;


    int icent=-1;
    double wcen=1;
    double wvz=1;
    //wxs=1;
    if(strcmp(ksp,"pp")==0){
      icent=0;
    }else if(strcmp(ksp,"pbpb")==0){
      icent = GetCentBin(hiBin);
    }
    else if(strcmp(ksp,"ppb")==0){
      icent = GetHFBin(hiHFSum);
    }


    if(ievt%10000==0)std::cout<<" ********** Event # " <<ievt<<std::endl;
    //std::cout<<" ********** Event # " <<ievt<<"\t vz : "<<vz<<"\t hiBin : "<<hiBin<<"\t wxs : "<<wxs<<std::endl;


    //! Centrality
    if(icent<0 || icent>ncent-1 )continue;

    for(int nj=0;nj<knj;nj++){ //! loop over different jet algorithms
      

      if(nj==0)iJet = &(c->ak2PF);
      else if(nj==1)iJet = &(c->ak3PF);
      else if(nj==2)iJet = &(c->ak4PF);
      else if(nj==3)iJet = &(c->ak5PF);

      else if(nj==4)iJet = &(c->akPu2PF);
      else if(nj==5)iJet = &(c->akPu3PF);
      else if(nj==6)iJet = &(c->akPu4PF);
      else if(nj==7)iJet = &(c->akPu5PF);

      else if(nj==8)iJet = &(c->ak2Calo);
      else if(nj==9)iJet = &(c->ak3Calo);
      else if(nj==10)iJet = &(c->ak4Calo);
      else if(nj==11)iJet = &(c->ak5Calo);

      else if(nj==12)iJet = &(c->akPu2Calo);
      else if(nj==13)iJet = &(c->akPu3Calo);
      else if(nj==14)iJet = &(c->akPu4Calo);
      else if(nj==15)iJet = &(c->akPu5Calo);



      //! xsec-weight
      double pthat = iJet->pthat;
      if(pthat > maxpthat)continue;
      hpthat->Fill(pthat,wxs*wcen*wvz);

      //if(!isHiPtTrkInJet(iJet,&(c->track),kDelR[nj]))continue;
      
      //std::cout<<"\t Jet Algorithm : "<<calgo[nj]<<"\t # of Jets  : "<<iJet->nref<<"\t pthat : "<<pthat<<std::endl;
      if(nj==0)hTotEve->Fill(1); //! akPu3PF      
      

      int *ljet = new int[3];
      FindLeadSubLeadJets(ksp,iJet,ljet);
      if(ljet[0]>=0 && ljet[1]>=0 && iJet->jtpt[ljet[0]]>30. && iJet->jtpt[ljet[1]]>30. && (iJet->trackMax[ljet[0]]>4. || iJet->trackMax[ljet[1]]>4.)){

	//! Leading jet response
	if(ljet[0]>=0)hratiocorrrefpt_lead[nj][icent]->Fill(iJet->refpt[ljet[0]],iJet->jtpt[ljet[0]]/iJet->refpt[ljet[0]],wxs*wcen*wvz);
	//! Sub leading jet response
	if(ljet[1]>=0)hratiocorrrefpt_slead[nj][icent]->Fill(iJet->refpt[ljet[1]],iJet->jtpt[ljet[1]]/iJet->refpt[ljet[1]],wxs*wcen*wvz);
	
	//! Jet energy scale comparison with data
	if(iJet->jtpt[ljet[0]]>80. && iJet->jtpt[ljet[1]]>80.){//! atleas a dijet
	  int mstat=1;
	  double ptdij = (iJet->jtpt[ljet[0]] + iJet->jtpt[ljet[1]])/2.;
	  if(ljet[2]>=0){
	    //if(iJet->jtpt[ljet[2]]/ptdij > 0.2)mstat=0;
	    mstat=0;
	  }
	  if(mstat){
	    double B=-9999;
	    double rn1 = gRandom->Rndm();
	    double rn2 = gRandom->Rndm();
	    if(rn1 > rn2){
	      B = (iJet->jtpt[ljet[0]] - iJet->jtpt[ljet[1]])/(iJet->jtpt[ljet[0]] + iJet->jtpt[ljet[1]]);
	    }else{
	      B = (iJet->jtpt[ljet[1]] - iJet->jtpt[ljet[0]])/(iJet->jtpt[ljet[1]] + iJet->jtpt[ljet[0]]);
	    }
	    hJetEnergyScale[nj][icent]->Fill(ptdij,B,wxs*wcen*wvz);
	  }
	}
      }//! Jet selection

      //! Gen matched jets loop
      int njets  = iJet->nref;
      for(int igen=0; igen<njets; igen++){
	int gj = igen;

	if(iJet->subid[igen] != 0)continue;

        //if(strcmp(ksp,"pbpb")==0 && iJet->subid[igen] != 0)continue;
	//else if(strcmp(ksp,"pp")==0 && iJet->refpt[igen]<0)continue;


	if(iJet->trackMax[gj]<4.0)continue;
	
	float rawpt   = iJet->rawpt[gj];
	float refpt   = iJet->refpt[gj];
	float refeta  = iJet->refeta[gj];
	float refphi  = iJet->refphi[gj];
	float recopt  = iJet->jtpt[gj];
	float recoeta = iJet->jteta[gj];
	float delr    = iJet->refdrjt[gj];

	//if(nj==1)cout<<"subid : "<<iJet->subid[igen]<<endl;

	if(fabs(refeta)<ketacut && refpt>15){
          //! Denominator for matching efficiency
	  hPtAll [nj][icent]->Fill(refpt,wxs*wcen*wvz);
          hEtaAll[nj][icent]->Fill(refeta,wxs*wcen*wvz);
          hPhiAll[nj][icent]->Fill(refphi,wxs*wcen*wvz);
	  
	  //! DeltaR efficiency
	  hDeltaR[nj][icent]->Fill(delr,wxs*wcen*wvz);
	  for (int idrmax=0;idrmax<100;idrmax++) {
	    float drmax = idrmax*0.01+0.005;
	    hDeltaRAll[nj][icent]->Fill(drmax,wxs*wcen*wvz);
	    if (delr<drmax){
	      hDeltaRSel[nj][icent]->Fill(drmax,wxs*wcen*wvz);
	    }
	  }
        }

	if(recopt<kptrecocut || refpt<kptgencut || refpt==0 || fabs(recoeta)>ketacut || fabs(delr)>kdRcut)continue;
	
        if(fabs(refeta)<ketacut && refpt>15){
          //! Numerator for matching efficiency
          hPtSel [nj][icent]->Fill(refpt,wxs*wcen*wvz);
          hEtaSel[nj][icent]->Fill(refeta,wxs*wcen*wvz);
          hPhiSel[nj][icent]->Fill(refphi,wxs*wcen*wvz);
        }
	

	//! Response
	for (int idr=0;idr<25;idr++) {
	  double drcut = 0.0+idr*(0.25-0.00)/(25-1);
	  if (delr>drcut) continue;
	  hRspVsDeltaR[nj][icent][idr]->Fill(recopt/refpt,wxs*wcen*wvz);
	}
	

	hratiocorrrefpt_genm[nj][icent]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hrecogen[nj][icent]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hrecoraw[nj][icent]->Fill(recopt,recopt/rawpt,wxs*wcen*wvz);
        hrawgen [nj][icent]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	hgenjrecoj [nj][icent]->Fill(refpt,recopt,wxs*wcen*wvz);

	int ieta=-1;
	if(fabs(recoeta)<1.3)ieta=0; //! barrel region
	else ieta=1; //! HCAL region

	
	//! Response in eta
	hratiocorrrefpt_eta[nj][icent][ieta]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hratiorawrefpt_eta [nj][icent][ieta]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);

	if(recopt > kPtCut){

	  //! Ratio of recopT / rawpT
	  hrajptrpt[nj][icent] ->Fill(iJet->jtpt[gj],iJet->jtpt[gj]/iJet->rawpt[gj],wxs*wcen*wvz);

	  //! Jet eta, phi, pt, eta-pt, eta-phi and pt-phi
	  hjeteta   [nj][icent]->Fill(iJet->jteta[gj],wxs*wcen*wvz);
	  hjetphi   [nj][icent]->Fill(iJet->jtphi[gj],wxs*wcen*wvz);
	  hjetpteta [nj][icent]->Fill(iJet->jtpt[gj],iJet->jteta[gj],wxs*wcen*wvz);
	  hjetptphi [nj][icent]->Fill(iJet->jtpt[gj],iJet->jtphi[gj],wxs*wcen*wvz);
	  hjetetaphi[nj][icent]->Fill(iJet->jteta[gj],iJet->jtphi[gj],wxs*wcen*wvz);


	  //! Fill the background pT info for the signal matched jets only
	  //! pileup study
	  hjetptpu_genm[nj][icent]->Fill(recopt,iJet->jtpu[gj],wxs*wcen*wvz);
	  hjetptpu_etab_genm[nj][icent][ieta]->Fill(recopt,iJet->jtpu[gj],wxs*wcen*wvz);

	  //! Jet bkgd estimation
	  //! (photonSum+neutralSum+chargedSum-rawpt)
	  double jbkgd  = (iJet->photonSum[gj]+iJet->neutralSum[gj]+iJet->chargedSum[gj]) - rawpt;
	  hjetbkgd_genm  [nj][icent]->Fill(jbkgd,wxs*wcen*wvz);
	  hjetptbkgd_genm[nj][icent]->Fill(recopt,jbkgd,wxs*wcen*wvz);
	  hjetptbkgd_etab_genm[nj][icent][ieta]->Fill(recopt,jbkgd,wxs*wcen*wvz);

	  hjbkgd[nj][icent]->Fill(recoeta,recopt,jbkgd,wxs*wcen*wvz);

	  //! PF fraction	  
	  hPFFraction_genm[nj][icent][0]->Fill(recopt,iJet->photonSum[gj]/rawpt,wxs*wcen*wvz);
	  hPFFraction_genm[nj][icent][1]->Fill(recopt,iJet->neutralSum[gj]/rawpt,wxs*wcen*wvz);
	  hPFFraction_genm[nj][icent][2]->Fill(recopt,iJet->chargedSum[gj]/rawpt,wxs*wcen*wvz);

	}//! recopt cut
	
	hNjets_genm [nj][icent]->Fill(refpt);
	hgenpt_genm [nj][icent]->Fill(refpt,wxs*wcen*wvz);
	hrecopt_genm[nj][icent]->Fill(recopt,wxs*wcen*wvz);	  
	hrawpt_genm [nj][icent]->Fill(rawpt,wxs*wcen*wvz);
	
	//! Very fine bin in ref pt used for response
	hrescrpt_genm [nj][icent]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hresrrpt_genm [nj][icent]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	hresrcrpt_genm[nj][icent]->Fill(recopt,recopt/rawpt,wxs*wcen*wvz);

	//! Very fine bin in ref eta used for response
	hrescreta_genm [nj][icent]->Fill(refeta,refpt,recopt/refpt,wxs*wcen*wvz);
	hresrreta_genm [nj][icent]->Fill(refeta,refpt,rawpt/refpt,wxs*wcen*wvz);
	hresrcreta_genm[nj][icent]->Fill(recoeta,refpt,recopt/rawpt,wxs*wcen*wvz);

        //! Response in different eta and phi bins
        int etabin = GetEtaBin(fabs(refeta));
	int phibin = GetPhiBin(refphi);

        //! Response in eta and phi bins
        if(etabin >= 0 && etabin<maxe){
	  hpteta[nj][icent][etabin]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	}
        if(phibin >= 0 && phibin<maxph){
	  hptphi[nj][icent][phibin]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	}
      }//! igen loop
      
      hNevt[nj][icent]->Fill(vz);
      delete [] ljet;
    }//! nj jet loop
    
    hBin->Fill(hiBin,wxs*wcen*wvz);
    hVz->Fill(vz,wxs*wcen*wvz);
    hHF->Fill(hiHF,wxs*wcen*wvz);
    
    hNTracks->Fill(ntracks,wxs*wcen*wvz);
    hNTracksHF->Fill(ntracks,hiHF,wxs*wcen*wvz);
    hBinHF->Fill(hiBin,hiHF,wxs*wcen*wvz);
    hBinNTracks->Fill(hiBin,ntracks,wxs*wcen*wvz);
    hNTracksHFPlusEta4->Fill(ntracks,hiHFPlusEta4,wxs*wcen*wvz);
    hHFPlusEta4->Fill(hiHFPlusEta4,wxs*wcen*wvz);

    iEvent++;
    //std::cout<<"Completed event #  "<<ievt<<std::endl; 
  }//! event loop ends
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  std::cout<<"Events which passed the pT hat cut : "<<hTotEve->Integral()<<" out of  : "<<hEvt->Integral()
	   <<" efficiency of all the cuts : " <<hTotEve->Integral()/hEvt->Integral()<<std::endl;
  std::cout<<std::endl;

  for(int nj=0;nj<knj;nj++){
    if(strcmp(ksp,"pp")==0)std::cout<<ksp<<"\t"<<calgo[nj]<<"\t # of Jets : "<<hNjets_genm[nj][0]->Integral()<<"\t # of Jets events : "<<hEvt->Integral()<<std::endl;
    else {
      for(int ic=0;ic<ncent;ic++){
	std::cout<<ksp<<"\t"<<calgo[nj]<<"\t # of Jets : "<<hNjets_genm[nj][ic]->Integral()<<"\t # of Jets events : "<<hNevt[nj][ic]->Integral()<<std::endl;
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

  std::cout<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;
  return 1;
}
int GetEtaBin(float eta)
{
  int ibin=-1;
  float min =0.0;
  ibin = (eta - min)*maxe/ketacut;
  return ibin;
}
int GetPhiBin(float phi)
{
  int ibin=-1;
  float min = -pi;
  ibin = (phi - min)*maxph/2*pi;
  return ibin;
}
bool selecJet(Jets *jetc,int it)
{
  bool goodJet = jetc->refpt[it]!=-999 && jetc->refphi[it]!=-999 && jetc->refeta[it]!=-999 && jetc->jtpt[it]!=-999 && jetc->jtphi[it]!=-999 && jetc->jteta[it]!=-999;
  return goodJet;
}
int GetDetEtaBin(float eta)
{
  int ibin=-1;
  if(eta>=0 && eta<1.3)ibin=0;       //! barrel
  else if(eta>=1.3 && eta<2.0)ibin=1;//! endcap+tracks
  else if(eta>=2.0 && eta<3.0)ibin=2;//! endcap-notracks
  else if(eta>=3.0 && eta<5.0)ibin=3;//! forward

  return ibin;
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

void FindLeadSubLeadJets(const char *sp,Jets *jetc, int *ljet)
{
  ljet[0]=-1; ljet[1]=-2;ljet[2]=-3;

  int njs  = jetc->nref;
  float tempt=-9;
  //! Get the leading jet               
  for(int ij=0; ij<njs; ij++){
    if(jetc->subid[ij] != 0)continue;
    //if(strcmp(sp,"pbpb")==0 && jetc->subid[ij] != 0)continue;
    //else if(strcmp(sp,"pp")==0 && jetc->refpt[ij]<0)continue;

    int gj = ij;
    if(fabs(jetc->jteta[gj])>ketacut || jetc->jtpt[gj]<120 || jetc->rawpt[gj]<15 || fabs(jetc->refdrjt[gj])>kdRcut)continue;
    float jetpt = jetc->jtpt[gj];
    if(jetpt > tempt){
      tempt = jetpt;
      ljet[0] = gj;
    }
  }
  if(ljet[0]>=0){
    tempt=-9;
    for(int ij=0; ij<njs; ij++){
      if(jetc->subid[ij] != 0)continue;
      //if(strcmp(sp,"pbpb")==0 && jetc->subid[ij] != 0)continue;
      //else if(strcmp(sp,"pp")==0 && jetc->refpt[ij]<0)continue;
      int gj = ij;
      if(gj==ljet[0])continue;
      if(fabs(jetc->jteta[gj])>ketacut || jetc->jtpt[gj]<30 || jetc->jtpt[gj]>jetc->jtpt[ljet[0]] || jetc->rawpt[gj]<15 || fabs(jetc->refdrjt[gj])>kdRcut)continue;
      float jetpt = jetc->jtpt[gj];
      float dphi  = delphi(jetc->jtphi[gj],jetc->jtphi[ljet[0]]);
      if(dphi < 2*pi/3.)continue;
      if (jetpt > tempt){
        tempt = jetpt;
        ljet[1] = gj;
      }
    }
    if(ljet[1]>=0){
      tempt=-9;
      for(int ij=0; ij<njs; ij++){
	if(jetc->subid[ij] != 0)continue;
	//if(strcmp(sp,"pbpb")==0 && jetc->subid[ij] != 0)continue;
	//else if(strcmp(sp,"pp")==0 && jetc->refpt[ij]<0)continue;
	int gj = ij;
        if(gj==ljet[0] || gj==ljet[1])continue;
        if(fabs(jetc->jteta[gj])>ketacut || jetc->jtpt[gj]<30 || jetc->jtpt[gj]>jetc->jtpt[ljet[0]] || jetc->rawpt[gj]<15 || fabs(jetc->refdrjt[gj])>kdRcut)continue;
        float jetpt = jetc->jtpt[ij];
        if (jetpt > tempt){
          tempt = jetpt;
          ljet[2] = gj;
        }
      }
    }
  }
}

double delphi(double phi1, double phi2)
{
  double dphi = phi2 - phi1;
  dphi = fabs(atan2(sin(dphi),cos(dphi)));
  return dphi;
}
int GetCentBin(int bin)
{
  //! PbPb
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
int GetHFBin(float hfet)
{
  int ibin=-1;
  if(hfet>0 && hfet<20)ibin=0;
  else if(hfet>=20 && hfet<25)ibin=1;
  else if(hfet>=25 && hfet<30)ibin=2;
  else if(hfet>=30 && hfet<40)ibin=3;
  else if(hfet>=40 && hfet<100)ibin=4;

  return ibin;
}

bool isHiPtTrkInJet(const char *sp, Jets *jetc, Tracks *trks,float delr)
{
  bool hiPtTrk=false;
  int njs  = jetc->nref;
  for(int ij=0; ij<njs; ij++){
    if(strcmp(sp,"pbpb")==0 && jetc->subid[ij] != 0)continue;
    else if(strcmp(sp,"pp")==0 && jetc->refpt[ij]<0)continue;
    int gj = ij;
    if(fabs(jetc->jteta[gj])>ketacut || jetc->rawpt[gj]<15)continue;
    double jetphi = jetc->jtphi[gj];
    double jeteta = jetc->jteta[gj];

    for(int j=0;j<trks->nTrk;j++) {
      bool goodTrack =  trks->trkPt[j]>4.0
        && fabs(trks->trkEta[j])<2.0
        && trks->highPurity[j]
	&& fabs(trks->trkDz1[j]/trks->trkDzError1[j])<3
	&& fabs(trks->trkDxy1[j]/trks->trkDxyError1[j])<3
        && trks->trkPtError[j]/trks->trkPt[j]<0.1
        ;
      if(!goodTrack)continue;

      double trkphi = trks->trkPhi[j];
      double trketa = trks->trkEta[j];

      double dr = sqrt(pow(trkphi - jetphi,2) + pow(trketa - jeteta,2));
      if(dr<delr){
        hiPtTrk=true;
        break;
      }
    }
  }
  return hiPtTrk;
}

void ShutoffBranches(HiForest *hi,const char *ksp)
{

  std::cout<<"Shutting off branches for "<<ksp<<std::endl;
  //! added by pawan
  //! Select only the branches you want to use for the analysis
  //! This increases the speed for running

  //! Evt tree
  hi->evtTree->SetBranchStatus("*",0,0);
  //hi->evtTree->SetBranchStatus("run",1,0);
  hi->evtTree->SetBranchStatus("vz",1,0);
  hi->evtTree->SetBranchStatus("hiBin",1,0);
  hi->evtTree->SetBranchStatus("hiHF",1,0);//! HF
  hi->evtTree->SetBranchStatus("hiHFplusEta4",1,0);//! HFplusEta4
  hi->evtTree->SetBranchStatus("hiHFminusEta4",1,0);//! HFminusEta4
  hi->evtTree->SetBranchStatus("hiNtracks",1,0);//! HFplusEta4


  //! Track tree
//  hi->trackTree->SetBranchStatus("*",0,0);
//  hi->trackTree->SetBranchStatus("nTrk",1,0);
//  hi->trackTree->SetBranchStatus("trkPt",1,0);
//  hi->trackTree->SetBranchStatus("trkEta",1,0);
//  hi->trackTree->SetBranchStatus("trkPhi",1,0);
//  hi->trackTree->SetBranchStatus("highPurity",1,0);
//  hi->trackTree->SetBranchStatus("trkDz1",1,0);
//  hi->trackTree->SetBranchStatus("trkDzError1",1,0);
//  hi->trackTree->SetBranchStatus("trkDxy1",1,0);
//  hi->trackTree->SetBranchStatus("trkDxyError1",1,0);
//

  //! PF jet tree
  hi->ak2PFJetTree->SetBranchStatus("*",0,0);
  hi->ak3PFJetTree->SetBranchStatus("*",0,0);
  hi->ak4PFJetTree->SetBranchStatus("*",0,0);
  hi->ak5PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu2PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu3PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu4PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu5PFJetTree->SetBranchStatus("*",0,0);
  hi->ak2CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak3CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak4CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak5CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu2CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu3CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu4CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu5CaloJetTree->SetBranchStatus("*",0,0);

  for(int i=0;i<nsize;i++){
    hi->ak2PFJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->ak3PFJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->ak4PFJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->ak5PFJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->akPu2PFJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->akPu3PFJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->akPu4PFJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->akPu5PFJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);

    hi->ak2CaloJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->ak3CaloJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->ak4CaloJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->ak5CaloJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->akPu2CaloJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->akPu3CaloJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->akPu4CaloJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
    hi->akPu5CaloJetTree->SetBranchStatus(Form("%s",kbr[i]),1,0);
  }

  ///
  std::cout<<"Loaded all tree variables "<<std::endl;

}
void GetXsection(const char *sp,double sqrts,double pthat,double &xsection,double &maxpthat){

  double xup=0;
  double xsub=0;

  if(sqrts==2.76){

    if(strcmp(sp,"pbpb")==0){
      if(pthat==15){
	maxpthat=30;
	xup=2.034e-01;
	xsub=1.079e-02;
      }
      else if(pthat==30){
	maxpthat=50;
	xup=1.079e-02;
	xsub=1.021e-03;
      }
      else if(pthat==50){
	maxpthat=80;
	xup=1.021e-03;
	xsub=9.913e-05;
      }
      else if(pthat==80){
	maxpthat=100;
	xup=9.913e-05;
	xsub=3.069e-05;
      }
      else if(pthat==100){
	maxpthat=120;	
	xup=3.069e-05;
	xsub=1.129e-05;
      }
      else if(pthat==120){
	maxpthat=170;
	xup=1.129e-05;
	xsub=1.470e-06;
      }
      else if(pthat==170){
	maxpthat=200;
	xup=1.470e-06;
	xsub=5.310e-07;
      }
      else if(pthat==200){
	maxpthat=250;
	xup=5.310e-07;
	xsub=1.192e-07;
      }
      else if(pthat==250){
	maxpthat=300;
	xup=1.192e-07;
	xsub=3.176e-08;
      }
      else if(pthat==300){
	maxpthat=9999;
	xup=3.176e-08;
	xsub=0;
      }
    }//! pbpb
    else  if(strcmp(sp,"pp")==0){
      if(pthat==15){
	maxpthat=30;
	xup=2.034e-01;
	xsub=1.079e-02;
      }
      else if(pthat==30){
	maxpthat=50;
	xup=1.079e-02;
	xsub=1.021e-03;
      }
      else if(pthat==50){
	maxpthat=80;
	xup=1.021e-03;
	xsub=9.913e-05;
      }
      else if(pthat==80){
	maxpthat=120;
	xup=9.913e-05;
	xsub=1.128e-05;
      }
      else if(pthat==120){
	maxpthat=170;
	xup =1.128e-05;
	xup=1.470e-06;
      }
      else if(pthat==170){
	maxpthat=220;
	xup=1.470e-06;
	xsub=2.837e-07;
      }
      else if(pthat==220){
	maxpthat=280;
	xup =2.837e-07;
	xsub=5.323e-08;
      }
      else if(pthat==280){
	maxpthat=9999;
	xup =5.323e-08;
	xsub=0;
      }
    }//! pp
  }else if (sqrts==5.02){ //! pPb
    if(pthat==15){
      maxpthat=30;
      xup =5.335e-01;
      xsub=3.378e-02;
    }else if(pthat==30){
      maxpthat=50;
      xup =3.378e-02;
      xsub=3.778e-03;
    }else if(pthat==50){
      maxpthat=80;
      xup =3.778e-03;
      xsub=4.412e-04;
    }else if(pthat==80){
      maxpthat=120;
      xup =4.412e-04;
      xsub=6.147e-05;
    }else if(pthat==120){
      maxpthat=170;
      xup=6.147e-05;
      xsub=1.018e-05;
    }else if(pthat==170){
      maxpthat=220;
      xup=1.018e-05;
      xsub=2.477e-06;
    }else if(pthat==220){
      maxpthat=280;
      xup =2.477e-06;
      xsub=6.160e-07;
    }else if(pthat==280){
      if(strcmp(sp,"pp")==0){
	maxpthat=9999;
	xup =6.160e-07;
	xsub=0;
      }else{
	xup =6.160e-07;
	maxpthat=370;
	xsub = 1.088e-07;
      }
    }else if(pthat==370){
      maxpthat=9999;
      xup=1.088e-07;
      xsub=0;
    }
  }
  xsection = xup-xsub;
}
