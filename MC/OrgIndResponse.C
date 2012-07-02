#include "hiForest.h"
#include <TF1.h>
#include <TProfile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
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


int GetCentBin(int /*bin*/);
int GetEtaBin(float /*eta*/);
int GetDetEtaBin(float /*eta*/);
int GetPhiBin(float /*phi*/);
int GetPtBin(float /*pt*/);
bool selecJet(Jets */*jetc*/, int /*ij*/);
void FindLeadSubLeadJets(Jets */*jetc*/, int */*ljet*/);

//! pt binning
double ptbins[] ={30.,40.,50.,60.,70.,80.,90.,
		  100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,
		  200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,
		  300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,
		  400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,
		  500.};
const int bins = sizeof(ptbins)/sizeof(double) - 1;

//! constants
const int iYear=2012;
const double pi=acos(-1.);
const double pi2=2*pi -1;
const int ncen=6;
const int ketar=4;
const int maxe=5;
const int maxph=5;
const float ketacut=2.;
const double kptrecocut=15.;
const double kptgencut =0.;
const double kdRcut=0.3;
const double kVzcut=15.;
const int rbins=50;

/* Full list 
const int knj=25;
const char *calgo[knj]= {
  "icPu5Calo",
  "ak1PF","ak2PF","ak3PF","ak4PF","ak5PF","ak6PF",
  "ak1Calo","ak2Calo","ak3Calo","ak4Calo","ak5Calo","ak6Calo",
  "akPu1PF","akPu2PF","akPu3PF","akPu4PF","akPu5PF","akPu6PF",
  "akPu1Calo","akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo","akPu6Calo"
};
*/

TStopwatch timer;
int IndResponse(double kPt=80,const char *ksp="pbpb")
{

  timer.Start();

  //! Load Lib
  gSystem->Load("hiForest_h.so");
  
  //! Load the hiforest input root file
  HiForest *c = 0;
  if(strcmp("pp",ksp)==0)c = new HiForest(Form("/net/hisrv0001/home/icali/hadoop/Pythia/Z2/ppDijet_merged/pp276Dijet%0.0f_merged.root",kPt),"pp2012",1,1);
  else {
    if(kPt==30)c = new HiForest("/net/hisrv0001/home/icali/hadoop/Hydjet1.8/Z2/Dijet_merged/Pythia30_HydjetDrum_mix01_HiForest2_v19.root","pbpb2012",0,1);
    else if(kPt==50)c = new HiForest("/net/hisrv0001/home/icali/hadoop/Hydjet1.8/Z2/Dijet_merged/Pythia50_HydjetDrum_mix01_HiForest2_v19.root","pbpb2012",0,1);
    else if(kPt==80)c = new HiForest("/net/hisrv0001/home/icali/hadoop/Hydjet1.8/Z2/Dijet_merged/Pythia80_HydjetDrum_mix01_HiForest2_v20.root","pbpb2012",0,1);
    else if(kPt==120)c = new HiForest("/net/hisrv0001/home/icali/hadoop/Hydjet1.8/Z2/Dijet_merged/Pythia120_HydjetDrum_mix01_HiForest2_v21_ivan.root","pbpb2012",0,1);
    else if(kPt==170)c = new HiForest("/net/hisrv0001/home/icali/hadoop/Hydjet1.8/Z2/Dijet_merged/Pythia170_HydjetDrum_mix01_HiForest2_v19.root","pbpb2012",0,1);
    else if(kPt==200)c = new HiForest("/net/hisrv0001/home/icali/hadoop/Hydjet1.8/Z2/Dijet_merged/Pythia200_HydjetDrum_mix01_HiForest2_v21_ivan.root","pbpb2012",0,1);
    else if(kPt==250)c = new HiForest("/net/hisrv0001/home/icali/hadoop/Hydjet1.8/Z2/Dijet_merged/Pythia250_HydjetDrum_mix01_HiForest2_v21_ivan.root","pbpb2012",0,1);
  }


  double xsection=0;
  double maxpthat=9999;
  if(kPt==15){
    maxpthat=30;
    xsection=1.566e-01;
  }
  else if(kPt==30){
    maxpthat=50;
    xsection=1.079e-02;
  }
  else if(kPt==50){
    maxpthat=80;
    xsection=1.021e-03;
  }
  else if(kPt==80){
    maxpthat=120;
    xsection=9.913e-05;
  }
  else if(kPt==120){
    maxpthat=170;
    xsection=1.128e-05;
  }
  else if(kPt==170){
    maxpthat=200;
    xsection=1.470e-06;
  }
  else if(kPt==200){
    maxpthat=250;
    xsection=5.310e-07;
  }
  else if(kPt==250){
    maxpthat=300;
    xsection=9.351e-08;
  }
  else if(kPt==300){
    maxpthat=9999;
    xsection=2.447e-08;
  }


  std::cout<<std::endl;
  std::cout<<std::endl;


  //! Don't want to loop over trees which is not used in the analysis
  //! event trees
  //c->hasEvtTree=0;

  //! jet trees
  //! Switch on only the jet trees which you  require
  /*const char *cjets[5] = {"icPu5","akPu1PF","akPu2PF","akPu1Calo","akPu2Calo"};
  c->SelectJetAlgo(cjets,5);
  */

  //! photon tree
  //c->hasPhotonTree=0;


  //! Select only the jet branches which you are going to use
  //! This increases the speed for running over trees with lot of branches.
  //! This is currently for only jet algos and evtTree  uniformly applied to all the 

  //! Event Tree 
  const char *evlist[]={"hiBin","vz","hiHF"};
  const int kevbr = sizeof(evlist)/sizeof(const char *);
  c->SelectBranches("evtTree",evlist,kevbr);

  //! jet Tree algorithms
  const char *jtlist[]={"nref","pthat","rawpt","jtpt","jteta","jtphi","jtpu","refpt","refeta","refphi","refdrjt","refparton_flavor"};
  const int kjtbr = sizeof(jtlist)/sizeof(const char *);
  c->SelectBranches("JetTree",jtlist,kjtbr);

  std::cout<<"Selected the branches of need from evtTree and JetTrees : "<<std::endl;

  //! To get the jet object from hiforest
  Jets *iJet=0;
  const int knj = c->GetNAlgo(); //! # of jet algorithms in this hiforest 
  std::cout<<"Loaded all tree variables and # of jet algorithms : "<<knj<<std::endl;
  std::cout<<"\t"<<std::endl;

  //! Away-side jet definition
  const char *cdphi="2pi3";
  double kdphicut = 7.*pi/8.;
  if(strcmp(cdphi,"2pi3")==0)kdphicut=2.*pi/3.;
  else if(strcmp(cdphi,"1pi4")==0)kdphicut=1.*pi/4.;

  //! Open a output file for histos
  TFile *fout = new TFile(Form("Output/%s/Response_newHiForest_DJ_%0.0fGeV_%s_%d.root",ksp,kPt,ksp,iYear),"RECREATE");
  //TFile *fout = new TFile(Form("Output/%s/novtxcut/Response_newHiForest_DJ_%0.0fGeV_%s_%d.root",ksp,kPt,ksp,iYear),"RECREATE");
  //TFile *fout = new TFile(Form("test_newHiForest_DJ_%0.0fGeV_%s_%d.root",kPt,ksp,iYear),"RECREATE");

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Running for %s ",ksp)<<std::endl;
  std::cout<<Form("pT  cut for %0.3f ",kptrecocut)<<std::endl;
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
  TProfile::SetDefaultSumw2();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *hEvt    = new TH1F("hEvt","# of events ",4,0,4);
  TH1F *hVz     = new TH1F("hVz","# of events ",80,-20,20);
  TH1F *hBin    = new TH1F("hBin","Centrality bin",40,0,40);
  TH1F *hHF     = new TH1F("hHF","Centrality variable from HF",600,0,6000);
  TH1F *hTotEve = new TH1F("hTotEve","# of events in the skimmed files",4,0,4);

  TH1F *hgenpt [knj][ncen], *hrecopt[knj][ncen], *hrawpt[knj][ncen];
  TH1F *hgenptC [knj][ncen], *hrecoptC[knj][ncen], *hrawptC[knj][ncen];
  TH1F *hgeneta[knj][ncen], *hrecoeta[knj][ncen];
  TH1F *hgenphi[knj][ncen], *hrecophi[knj][ncen];

  //! Ratios of the pt distributions
  TProfile *hrecogen[knj][ncen], *hrecoraw[knj][ncen], *hrawgen[knj][ncen], *hrecoraw_ref[knj][ncen];

  //! Resposnse
  TH2F *hcorrptrefpt[knj][ncen], *hrawptrefpt[knj][ncen], *hcorrptrawpt[knj][ncen];
  TH2F *hrescrpt[knj][ncen], *hresrrpt[knj][ncen], *hresrcrpt[knj][ncen];
  TH2F *hratiorawrefpt[knj][ncen], *hratiocorrrefpt[knj][ncen], *hratiocorrrawpt[knj][ncen];
  TH2F *hratiorawrefpt_eta[knj][ncen][2], *hratiocorrrefpt_eta[knj][ncen][2];
  TH2F *hratiorawrefpt_mb_eta[knj][2], *hratiocorrrefpt_mb_eta[knj][2];

  TH2F *hdiffpteta[knj][ncen][maxe], *hdiffptphi[knj][ncen][maxph];
  TH2F *hpteta[knj][ncen][maxe], *hptphi[knj][ncen][maxph];
  TH2F *hdiffpt[knj][ncen], *hdiffeta[knj][ncen], *hdiffphi[knj][ncen];


  TH2F *hgenj1genj2  [knj][ncen];
  TH2F *hrecoj1recoj2[knj][ncen];

  TH2F *hgenj1recoj1 [knj][ncen];
  TH2F *hgenj2recoj2 [knj][ncen];

  TH2F *hgenj1recoj2 [knj][ncen];
  TH2F *hgenj2recoj1 [knj][ncen];

  TH2F *hgenjrecoj[knj][ncen];
  TH2F *hgenjrecoj_pflavor[knj][ncen];
  TH2F *hgenjrecoj_awayside[knj][ncen];


  //! Away side jet energy scale centrality dependence
  TH2F *hjesas    [knj][ncen];
  TH2F *hjesasr   [knj][ncen];
  TH2F *hjesas_eta[knj][ncen][maxe];

  TH2F *hFracjets[knj][ncen];
  TH2F *hAvjets[knj][ncen];
  TH2F *hMeanPtjets[knj][ncen];
  TH1F *hNjets[knj][ncen];
  TH1F *hNevt [knj][ncen];

  //! Pileup effect study
  TProfile *hjetptpumb[knj];             //! mb 
  TH2F *hjetptpu[knj][ncen];             //! centrality                                                                       
  TH2F *hjetptpumb_etab[knj][ketar];     //! eta dependence                                                                   
  TH2F *hjetptpu_etab[knj][ncen][ketar]; //! eta dependence             

  //! xj histograms
  TH2F *hxgj [knj][ncen];
  TH2F *hdpgj[knj][ncen];
  TH2F *hxdp [knj][ncen][bins];

  //! Efficency histos
  TH1F *hPtAll [knj][ncen], *hPtSel[knj][ncen];
  TH1F *hEtaAll[knj][ncen], *hEtaSel[knj][ncen];
  TH1F *hPhiAll[knj][ncen], *hPhiSel[knj][ncen];

  //! Response vs deltar
  TH1F *hRspVsDeltaR[knj][ncen][25];

  //! Fake jet study
  TH1F *hFakePtAll [knj][ncen], *hFakePtSel[knj][ncen];
  TH1F *hFakeEtaAll[knj][ncen], *hFakeEtaSel[knj][ncen];
  TH1F *hFakePhiAll[knj][ncen], *hFakePhiSel[knj][ncen];

  //! DeltaR efficiency                                                                                                                                                              
  TH1F *hDeltaR[knj][ncen];
  TH1F *hDeltaRAll[knj][ncen];
  TH1F *hDeltaRSel[knj][ncen];

  for(int nj=0;nj<knj;nj++){
    const char *algoname = c->GetAlgoName(nj);

    //! pileup study
    hjetptpumb[nj] = new TProfile(Form("hjetptpumb%d",nj),Form("jet(pt:pu) distribution jet %s",algoname),bins,ptbins,0,300);
    for(int ie=0;ie<ketar;ie++){
      hjetptpumb_etab[nj][ie] = new TH2F(Form("hjetptpumb_etab%d_%d",nj,ie),Form("jet(pt:pu) distribution jet %s etabin %d",algoname,ie),bins,ptbins,100,0,300);
    }

    for(int iea=0;iea<2;iea++){
      hratiorawrefpt_mb_eta[nj][iea]= new TH2F(Form("hratiorawrefpt_mb_eta%d_%d",nj,iea),Form("Raw jet / Gen jet p_{T} (raw) distribution jet mb %s etabin%d",algoname,iea),
					       bins,ptbins,rbins,0.0,5.0);
      hratiocorrrefpt_mb_eta[nj][iea]= new TH2F(Form("hratiocorrrefpt_mb_eta%d_%d",nj,iea),Form("Reco jet / Gen jet p_{T} (raw) distribution jet mb %s etabin%d",algoname,iea),
						bins,ptbins,rbins,0.0,5.0);
    }


    for(int icen=0;icen<ncen;icen++){
      hFracjets[nj][icen]   = new TH2F(Form("hFracjets%d_%d",nj,icen),Form("Fraction of jets in given pt hat cent %d %s",icen,algoname),500,0,1000,500,0.,1000.);
      hAvjets[nj][icen]     = new TH2F(Form("hAvjets%d_%d",nj,icen),Form("<#> of jets cent %d %s",icen,algoname),500,0,1000,30,0,30);
      hMeanPtjets[nj][icen] = new TH2F(Form("hMeanPtjets%d_%d",nj,icen),Form("<pT> of jets cent %d %s",icen,algoname),500,0,1000,500,0,1000);

      hNjets[nj][icen] = new TH1F(Form("hNjets%d_%d",nj,icen),Form("# of jets cent %d jets %s",icen,algoname),3,0.0,3.0);
      hNevt [nj][icen] = new TH1F(Form("hNevt%d_%d",nj,icen),Form("# of events cent %d %s",icen,algoname),40,-40,40);
      hgeneta[nj][icen] = new TH1F(Form("hgeneta%d_%d",nj,icen),Form("gen eta distribution jet centb %d %s",icen,algoname),60,-3.0,3.0);
      hrecoeta[nj][icen] = new TH1F(Form("hrecoeta%d_%d",nj,icen),Form("reco eta distribution jet centb %d %s",icen,algoname),60,-3.0,3.0);

      hgenphi[nj][icen] = new TH1F(Form("hgenphi%d_%d",nj,icen),Form("gen phi distribution jet centb %d %s",icen,algoname),36,-pi,pi);
      hrecophi[nj][icen] = new TH1F(Form("hrecophi%d_%d",nj,icen),Form("reco phi distribution jet centb %d %s",icen,algoname),36,-pi,pi);

      hgenpt[nj][icen]  = new TH1F(Form("hgenpt%d_%d",nj,icen),Form("gen p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000);
      hrecopt[nj][icen] = new TH1F(Form("hrecopt%d_%d",nj,icen),Form("reco p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000);
      hrawpt[nj][icen]  = new TH1F(Form("hrawpt%d_%d",nj,icen),Form("raw p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000);

      //! with pt bins
      hgenptC[nj][icen]  = new TH1F(Form("hgenptC%d_%d",nj,icen),Form("gen p_{T} distribution jet centb %d %s",icen,algoname),bins,ptbins);
      hrecoptC[nj][icen] = new TH1F(Form("hrecoptC%d_%d",nj,icen),Form("reco p_{T} distribution jet centb %d %s",icen,algoname),bins,ptbins);
      hrawptC[nj][icen]  = new TH1F(Form("hrawptC%d_%d",nj,icen),Form("raw p_{T} distribution jet centb %d %s",icen,algoname),bins,ptbins);

      //! Ratios
      hrecogen[nj][icen] = new TProfile(Form("hrecogen%d_%d",nj,icen),Form("reco/gen p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000);
      hrecoraw[nj][icen] = new TProfile(Form("hrecoraw%d_%d",nj,icen),Form("reco/raw p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000);
      hrecoraw_ref[nj][icen]  = new TProfile(Form("hrecoraw_ref%d_%d",nj,icen),Form("reco/raw : ref p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000);
      hrawgen[nj][icen]  = new TProfile(Form("hrawgen%d_%d",nj,icen),Form("raw/gen p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000);


      hcorrptrefpt[nj][icen]= new TH2F(Form("hcorrptrefpt%d_%d",nj,icen),Form("Gen jet:Reco jet p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000,500,0,1000);
      hrawptrefpt[nj][icen]= new TH2F(Form("hrawptrefpt%d_%d",nj,icen),Form("Gen jet:Raw jet p_{T}  distribution jet centb %d %s",icen,algoname),500,0,1000,500,0,1000);
      hcorrptrawpt[nj][icen]= new TH2F(Form("hcorrptrawpt%d_%d",nj,icen),Form("Reco jet Corr:Raw jet p_{T}  distribution jet centb %d %s",icen,algoname),500,0,1000,500,0,1000);


      hrescrpt[nj][icen]= new TH2F(Form("hrescrpt%d_%d",nj,icen),Form("Gen jet:(Reco/Gen) jet p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000,150,0,5);
      hresrrpt[nj][icen]= new TH2F(Form("hresrrpt%d_%d",nj,icen),Form("Gen jet:(Raw/Gen) jet p_{T}  distribution jet centb %d %s",icen,algoname),500,0,1000,150,0,5);
      hresrcrpt[nj][icen]= new TH2F(Form("hresrcrpt%d_%d",nj,icen),Form("Reco jet:(Reco/Raw) jet p_{T} distribution jet centb %d %s",icen,algoname),500,0,1000,150,0,5);

      hratiorawrefpt[nj][icen]= new TH2F(Form("hratiorawrefpt%d_%d",nj,icen),Form("Raw jet / Gen jet p_{T} (raw) distribution jet centb %d %s",icen,algoname),
					 bins,ptbins,rbins,0.0,5.0);
      hratiocorrrefpt[nj][icen]= new TH2F(Form("hratiocorrrefpt%d_%d",nj,icen),Form("Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,algoname),
					  bins,ptbins,rbins,0.0,5.0);
      hratiocorrrawpt[nj][icen]= new TH2F(Form("hratiocorrrawpt%d_%d",nj,icen),Form("Correc. jet / Raw jet jet p_{T} distribution jet centb %d %s",icen,algoname),
                                          bins,ptbins,rbins,0.0,5.0);


      hjesas[nj][icen]= new TH2F(Form("hjesas%d_%d",nj,icen),Form("Reco jet / Gen jet p_{T} (corr.) distribution away-side centb %d %s",icen,algoname),
				 bins,ptbins,rbins,0.0,5.0);
      hjesasr[nj][icen]= new TH2F(Form("hjesasr%d_%d",nj,icen),Form("Rawjet / Gen jet p_{T} (corr.) distribution away-side centb %d %s",icen,algoname),
				  bins,ptbins,rbins,0.0,5.0);

      for(int ie=0;ie<2;ie++){
        hratiorawrefpt_eta[nj][icen][ie]= new TH2F(Form("hratiorawrefpt_eta%d_%d_%d",nj,icen,ie),
						   Form("Raw jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",icen,algoname,ie),
						   bins,ptbins,rbins,0.0,5.0);
        hratiocorrrefpt_eta[nj][icen][ie]= new TH2F(Form("hratiocorrrefpt_eta%d_%d_%d",nj,icen,ie),
						    Form("Reco jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",icen,algoname,ie),
						    bins,ptbins,rbins,0.0,5.0);
      }
      
      
      hdiffpt[nj][icen]= new TH2F(Form("hdiffpt%d_%d",nj,icen),Form(" diff p_{T} distribution jet centb %d %s",icen,algoname),bins,ptbins,500,-5.0,5.0);
      hdiffeta[nj][icen] = new TH2F(Form("hdiffeta%d_%d",nj,icen),Form("diff eta(pt) distribution jet centb %d %s",icen,algoname),bins,ptbins,100,-0.5,0.5);
      hdiffphi[nj][icen] = new TH2F(Form("hdiffphi%d_%d",nj,icen),Form("diff phi(pt) distribution jet centb %d %s",icen,algoname),bins,ptbins,100,-0.5,0.5);
      hjetptpu[nj][icen] = new TH2F(Form("hjetptpu%d_%d",nj,icen),Form("jet(pt:pu) distribution jet centb %d %s",icen,algoname),bins,ptbins,100,0,300);

      for(int ie=0;ie<ketar;ie++){
        hjetptpu_etab[nj][icen][ie] = new TH2F(Form("hjetptpu_etab%d_%d_%d",nj,icen,ie),Form("jet(pt:pu) distribution jet %s etabin %d cen %d",algoname,icen,ie),
					       bins,ptbins,100,0,300);
      }
      
      for(int m=0;m<maxe;m++){
        hdiffpteta[nj][icen][m] = new TH2F(Form("hdiffpteta%d_%d_%d",nj,icen,m),Form("diff pt(eta) distribution cent %d jet %s etabin%d",icen,algoname,m),
					   bins,ptbins,500,-5.0,5.0);
        hpteta[nj][icen][m] = new TH2F(Form("hpteta%d_%d_%d",nj,icen,m),Form("resolution  pt(eta) distribution cent %d jet %s etabin%d",icen,algoname,m),
				       bins,ptbins,rbins,0.0,5.0);
        hjesas_eta[nj][icen][m]= new TH2F(Form("hjesas_eta%d_%d_%d",nj,icen,m),Form("Reco jet / Gen jet p_{T} (corr.) distribution away-side centb %d %s etabin%d",
										    icen,algoname,m),bins,ptbins,rbins,0.0,5.0);
      }
      
      for(int m=0;m<maxph;m++){
	hdiffptphi[nj][icen][m] = new TH2F(Form("hdiffptphi%d_%d_%d",nj,icen,m),Form("diff pt(phi) distribution cent %d jet %s phibin%d",icen,algoname,m),
					   bins,ptbins,500,-5.0,5.0);
	hptphi[nj][icen][m] = new TH2F(Form("hptphi%d_%d_%d",nj,icen,m),Form("resolution pt(phi) distribution cent %d jet %s phibin%d",icen,algoname,m),
				       bins,ptbins,rbins,0.0,5.0);
      }

      hxgj [nj][icen] = new TH2F(Form("hxgj%d_%d",nj,icen),Form("x=jet/gamma %s cent %d",algoname,icen),bins,ptbins,48,0.,1.2);
      hdpgj[nj][icen] = new TH2F(Form("hdpgj%d_%d",nj,icen),Form("dphi %s cent %d ",algoname,icen),bins,ptbins,36,0.,pi);

      for(int ib=0;ib<bins;ib++){
        hxdp [nj][icen][ib]  = new TH2F(Form("hxdp%d_%d_%d",nj,icen,ib),Form("x:dphi %s cent %d ptbin %d",algoname,icen,ib),48,0,1.2,36,0.,pi);
      }

      hgenj1genj2[nj][icen] =new TH2F(Form("hgenj1genj2%d_%d",nj,icen),Form("gen jet1(lead) : gen jet2(sublead) %s cent %d",algoname,icen),600,0.,300.,600,0.,300.);
      hrecoj1recoj2[nj][icen] =new TH2F(Form("hrecoj1recoj2%d_%d",nj,icen),Form("reco jet1(lead) : reco jet2(sublead) %s cent %d",algoname,icen),600,0.,300.,600,0.,300.);
      hgenj1recoj1[nj][icen] =new TH2F(Form("hgenj1recoj1%d_%d",nj,icen),Form("gen jet1(lead) : reco jet1(lead) %s cent %d",algoname,icen),600,0.,300.,400,0.,300.);
      hgenj2recoj2[nj][icen] =new TH2F(Form("hgenj2recoj2%d_%d",nj,icen),Form("gen jet2(subleadlead) : reco jet2(sublead) %s cent %d",algoname,icen),600,0.,300.,400,0.,300.);
      hgenj1recoj2[nj][icen] =new TH2F(Form("hgenj1recoj2%d_%d",nj,icen),Form("gen jet1(lead) : reco jet2(sublead) %s cent %d",algoname,icen),600,0.,300.,600,0.,300.);
      hgenj2recoj1[nj][icen] =new TH2F(Form("hgenj2recoj1%d_%d",nj,icen),Form("gen jet2(sublead) : reco jet1(lead) %s cent %d",algoname,icen),600,0.,300.,600,0.,300.);
      
      hgenjrecoj[nj][icen] =new TH2F(Form("hgenjrecoj%d_%d",nj,icen),Form("gen jet2 : reco jet2 %s cent %d",algoname,icen),500,0.,1000.,500,0.,1000.);
      hgenjrecoj_pflavor[nj][icen] =new TH2F(Form("hgenjrecoj_pflavor%d_%d",nj,icen),Form("gen jet : reco jet with parton flavor %s cent %d",algoname,icen),500,0.,1000.,500,0.,1000.);
      hgenjrecoj_awayside[nj][icen] =new TH2F(Form("hgenjrecoj_awayside%d_%d",nj,icen),Form("gen jet : reco jet away-side (%s) %s cent %d",cdphi,algoname,icen),500,0.,1000.,500,0.,1000.);

      //! efficency histograms
      hPtAll [nj][icen] = new TH1F(Form("hPtAll%d_%d",nj,icen),Form("Denominator pT for algorithm %s cent %d",algoname,icen),bins,ptbins);
      hEtaAll[nj][icen] = new TH1F(Form("hEtaAll%d_%d",nj,icen),Form("Denominator eta  for algorithm %s cent %d",algoname,icen),20,-2.0,2.0);
      hPhiAll[nj][icen] = new TH1F(Form("hPhiAll%d_%d",nj,icen),Form("Denominator  phi  for algorithm %s cent %d",algoname,icen),20,-pi,pi);
      
      hPtSel [nj][icen] = new TH1F(Form("hPtSel%d_%d",nj,icen),Form("Numerator pT for algorithm %s cent %d",algoname,icen),bins,ptbins);
      hEtaSel[nj][icen] = new TH1F(Form("hEtaSel%d_%d",nj,icen),Form("Numerator eta  for algorithm %s cent %d",algoname,icen),20,-2.0,2.0);
      hPhiSel[nj][icen] = new TH1F(Form("hPhiSel%d_%d",nj,icen),Form("Numerator  phi  for algorithm %s cent %d",algoname,icen),20,-pi,pi);

      //! Fake jets
      hFakePtAll [nj][icen] = new TH1F(Form("hFakePtAll%d_%d",nj,icen),Form("Denominator Fake jets pT for algorithm %s cent %d",algoname,icen),bins,ptbins);
      hFakeEtaAll[nj][icen] = new TH1F(Form("hFakeEtaAll%d_%d",nj,icen),Form("Denominator Fake jets eta  for algorithm %s cent %d",algoname,icen),20,-2.0,2.0);
      hFakePhiAll[nj][icen] = new TH1F(Form("hFakePhiAll%d_%d",nj,icen),Form("Denominator Fake  phi  for algorithm %s cent %d",algoname,icen),20,-pi,pi);

      hFakePtSel [nj][icen] = new TH1F(Form("hFakePtSel%d_%d",nj,icen),Form("Numerator Fake jet pT for algorithm %s cent %d",algoname,icen),bins,ptbins);
      hFakeEtaSel[nj][icen] = new TH1F(Form("hFakeEtaSel%d_%d",nj,icen),Form("Numerator Fake jet eta  for algorithm %s cent %d",algoname,icen),20,-2.0,2.0);
      hFakePhiSel[nj][icen] = new TH1F(Form("hFakePhiSel%d_%d",nj,icen),Form("Numerator Fake jet  phi  for algorithm %s cent %d",algoname,icen),20,-pi,pi);

      hDeltaR[nj][icen]    = new TH1F(Form("hDeltaR%d_%d",nj,icen),Form("#DeltaR for algorithm %s cent %d",algoname,icen),100,0,1);
      hDeltaRAll[nj][icen] = new TH1F(Form("hDeltaRAll%d_%d",nj,icen),Form("#DeltaR (all) for algorithm %s cent %d",algoname,icen),100,0,1);
      hDeltaRSel[nj][icen] = new TH1F(Form("hDeltaRSel%d_%d",nj,icen),Form("#DeltaR (sel) for algorithm %s cent %d",algoname,icen),100,0,1);

      for(int ir=0;ir<25;ir++){
	//! Response vs DeltaR
        hRspVsDeltaR[nj][icen][ir] = new TH1F(Form("hRspVsDeltaR%d_%d_%d",nj,icen,ir),Form(" <recopt/refpt> vs. #DeltaR (%d) algorithm %s cent %d",ir,algoname,icen),rbins,0.0,5.0);
      }
    }//! icen
  }//! nj
  std::cout<<"Initialized the histograms " <<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////// 


  //! Centrality reweighting function
  TF1* fcen = new TF1("fcen","exp(-1.0*pow(x+1.11957e+01,2)/pow(1.34120e+01,2)/2)",0,40);

  Long64_t nb = 0;
  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;
  std::cout<<std::endl;
  hEvt->Fill(2,nentries);

  //! weight  for the merging of the samples for different pT hat bins
  Float_t wxs = xsection/nentries;
  
  Int_t iEvent=0; 
  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
  //for (Long64_t ievt=0; ievt<1;ievt++) {//! event loop
    //! load the hiForest event
    nb += c->GetEntry(ievt);

    int hiBin       = c->evt.hiBin;
    float vz        = c->evt.vz;
    float hiHF      = c->evt.hiHF;

    //! testing
    //if(hiBin>4 && strcmp(ksp,"pbpb")==0)continue;

    //! apply vertex cut
    if(fabs(vz)>kVzcut)continue;


    //! Centrality bin
    if(hiBin<0 || hiBin>39)continue;
    int centb=-1;
    double wcen=1;
    if(strcmp(ksp,"pbpb")==0){
      centb=GetCentBin(hiBin);
      wcen = fcen->Eval(hiBin);
    }else{
      centb=0; //! pp
      wcen=1;
    }

    if(ievt%10000==0)std::cout<<" ********** Event # " <<ievt<<std::endl;
    //std::cout<<" ********** Event # " <<ievt<<"\t vz : "<<vz<<"\t hiBin : "<<hiBin<<"\t wxs : "<<wxs<<std::endl;


    //! Centrality from 0-90% 
    if(centb==-1 || centb==ncen)continue;


    //! Methods to get the jet objects of a given algorithm
    /* Get the Jet object of a given type of algorithm according to the number from Full list
       iJet = c->GetJet(15);
       std::cout<<"\t GetJet()       : "<<c->GetAlgoName(15)<<"\t # of Jets  : "<<iJet->nref<<std::endl;
    */
    /* Get the Jet object of a given type of algorithm according to the NAME from Full list
    iJet = c->GetJetByAlgo("akPu3PF");
    std::cout<<"\t GetJetByAlgo() : "<<c->GetAlgoName(15)<<"\t # of Jets  : "<<iJet->nref<<std::endl;
    */

    int istat=0;
    for(int nj=0;nj<knj;nj++){ //! loop over different jet algorithms
      //for(int nj=15;nj<16;nj++){ //! loop over different jet algorithms
      //! Get the jet object
      iJet = c->GetJet(nj);

      //! xsec-weight
      double pthat = iJet->pthat;
      if(pthat > maxpthat)continue;
      istat=1;
      
      //std::cout<<"\t Jet Algorithm : "<<c->GetAlgoName(nj)<<"\t # of Jets  : "<<iJet->nref<<"\t pthat : "<<pthat<<std::endl;
      if(nj==0)hTotEve->Fill(1); //! akPu3PF      

      float njets=0.;
      float meanpt=0.;
      int *ljet = new int[2];
      FindLeadSubLeadJets(iJet,ljet);

      //if(ljet[0]>=0 && ljet[1]>=0 && iJet->refpt[ljet[0]]< (maxpthat+10) && iJet->refpt[ljet[1]]<(maxpthat+10)){
      if(ljet[0]>=0 && ljet[1]>=0){
	
        hgenj1genj2  [nj][centb]->Fill(iJet->refpt[ljet[0]],iJet->refpt[ljet[1]],wxs*wcen);
        hrecoj1recoj2[nj][centb]->Fill(iJet->jtpt [ljet[0]],iJet->jtpt [ljet[1]],wxs*wcen);
	
        hgenj1recoj1 [nj][centb]->Fill(iJet->refpt[ljet[0]],iJet->jtpt [ljet[0]],wxs*wcen);
        hgenj2recoj2 [nj][centb]->Fill(iJet->refpt[ljet[1]],iJet->jtpt [ljet[1]],wxs*wcen);
	
        hgenj1recoj2 [nj][centb]->Fill(iJet->refpt [ljet[0]],iJet->jtpt[ljet[1]],wxs*wcen);
        hgenj2recoj1 [nj][centb]->Fill(iJet->refpt [ljet[1]],iJet->jtpt[ljet[0]],wxs*wcen);
      }


      //int istat=0;
      for(int ij=0; ij<iJet->nref; ij++){
	//if(!selecJet(iJet,ij))continue;

	float refpt   = iJet->refpt[ij];
        float recopt  = iJet->jtpt[ij];
        float rawpt   = iJet->rawpt[ij];
	float refeta  = iJet->refeta[ij];
        float recoeta = iJet->jteta[ij];
        float refphi  = iJet->refphi[ij];
        float recophi = iJet->jtphi[ij];
	int refparton_flavor = iJet->refparton_flavor[ij];

	//! Reference parton flavor
	//if(abs(refparton_flavor)>21 || refparton_flavor==-999 || refpt==0)continue;

	//! Matching criteria 
        float delr = iJet->refdrjt[ij];

	//! extra condition to remove the stray ref pt which
	//! gets extra weight factor
	//if(refpt > (maxpthat+10))continue;


        //! For fake jets
        if(fabs(recoeta)<ketacut){
          hFakePtAll[nj][centb]->Fill(recopt,wxs*wcen);
          hFakeEtaAll[nj][centb]->Fill(recoeta,wxs*wcen);
          hFakePhiAll[nj][centb]->Fill(recophi,wxs*wcen);
          if((refparton_flavor==-999 || refparton_flavor==0) && recopt>0){
            hFakePtSel[nj][centb]->Fill(recopt,wxs*wcen);
            hFakeEtaSel[nj][centb]->Fill(recoeta,wxs*wcen);
            hFakePhiSel[nj][centb]->Fill(recophi,wxs*wcen);
          }
        }

	//if(refparton_flavor==-999 || refparton_flavor==0)continue;
	if(recopt<kptrecocut || refpt<kptgencut || refpt==0)continue;
	//std::cout<<"\t \t "<<ij<<"\t refpt : "<<refpt<<"\t recopt : "<<recopt<<"\t raw pT : "<<rawpt<<std::endl;

	if(fabs(refeta)<ketacut && refpt>80){
          //! Denominator for matching efficiency
	  hPtAll [nj][centb]->Fill(refpt,wxs*wcen);
          hEtaAll[nj][centb]->Fill(refeta,wxs*wcen);
          hPhiAll[nj][centb]->Fill(refphi,wxs*wcen);
	  
	  
	  //! Response
	  for (int idr=0;idr<25;idr++) {
	    double drcut = 0.0+idr*(0.25-0.00)/(25-1);
	    if (delr>drcut) continue;
	    hRspVsDeltaR[nj][centb][idr]->Fill(recopt/refpt,wxs*wcen);
	  }
	  
	  //! DeltaR efficiency
	  hDeltaR[nj][centb]->Fill(delr,wxs*wcen);
	  for (int idrmax=0;idrmax<100;idrmax++) {
	    float drmax = idrmax*0.01+0.005;
	    hDeltaRAll[nj][centb]->Fill(drmax,wxs*wcen);
	    if (delr<drmax) hDeltaRSel[nj][centb]->Fill(drmax,wxs*wcen);
	  }
        }

        //! Matching cut for gen-jet and reco-jet
	if(delr > kdRcut)continue;


        if(fabs(iJet->refeta[ij])<ketacut && refpt>80){
          //! Numerator for matching efficiency
          hPtSel [nj][centb]->Fill(refpt,wxs*wcen);
          hEtaSel[nj][centb]->Fill(refeta,wxs*wcen);
          hPhiSel[nj][centb]->Fill(refphi,wxs*wcen);
        }

        //! pile up eta dependence
        int ebin = GetDetEtaBin(fabs(recoeta));
        hjetptpumb_etab[nj][ebin]->Fill(recopt,iJet->jtpu[ij],wxs*wcen);
	hjetptpu_etab[nj][centb][ebin]->Fill(recopt,iJet->jtpu[ij],wxs*wcen);

	//! jet selction cut
        if(fabs(recoeta)>ketacut)continue;

	/*
	if(refpt>140 && refpt<180 && nj==15){
	  std::cout<<"\t \t pthat : "<<iJet->pthat<<"\t"<<ij<<" refpt : "<<refpt<<"\t rawpt : "<<rawpt<<"\t recopt : "<<recopt<<"\t scale : "<<recopt/refpt<<std::endl;



	}
	*/


	
	//istat=1;
	//! 1D distributions
	njets++;
	meanpt += refpt;
        hNjets [nj][centb]->Fill(1.);
        hgenpt [nj][centb]->Fill(refpt,wxs*wcen);
        hgeneta[nj][centb]->Fill(refeta,wxs*wcen);
	hgenphi[nj][centb]->Fill(refphi,wxs*wcen);


	hrecopt [nj][centb]->Fill(recopt,wxs*wcen);
        hrecoeta[nj][centb]->Fill(recoeta,wxs*wcen);
        hrecophi[nj][centb]->Fill(recophi,wxs*wcen);
        hrawpt [nj][centb]->Fill(rawpt,wxs*wcen);


        hgenptC [nj][centb]->Fill(refpt,wxs*wcen);
	hrecoptC [nj][centb]->Fill(recopt,wxs*wcen);
        hrawptC [nj][centb]->Fill(rawpt,wxs*wcen);

	hrecogen[nj][centb]->Fill(refpt,recopt/refpt,wxs*wcen);
        hrawgen [nj][centb]->Fill(refpt,rawpt/refpt,wxs*wcen);
        hrecoraw[nj][centb]->Fill(rawpt,recopt/rawpt,wxs*wcen);
        hrecoraw_ref [nj][centb]->Fill(refpt,recopt/rawpt,wxs*wcen);

	hFracjets[nj][centb]->Fill(pthat,refpt,wxs*wcen);

	/*
	if(nj==2 && centb==0){
	  std::cout<<"ievt : "<<ievt<<"\t refpt : "<<refpt<<"\t recopt : "<<recopt<<"\t rawpt : "<<rawpt
		   <<"\t scale : "<<(recopt/refpt)<<"\t wxs : "<<wxs<<"\t wcen : "<<wcen<<std::endl;
	}
	*/

	//! Response  (ratio of recopt/refpt)
        hratiocorrrefpt[nj][centb]->Fill(refpt,recopt/refpt,wxs*wcen);
	hratiorawrefpt [nj][centb]->Fill(refpt,rawpt/refpt,wxs*wcen);
	hratiocorrrawpt [nj][centb]->Fill(recopt,recopt/rawpt,wxs*wcen);

	//! 2D correlation between refpt and recopt
        hcorrptrefpt[nj][centb]->Fill(refpt,recopt,wxs*wcen);
        hrawptrefpt [nj][centb]->Fill(refpt,rawpt,wxs*wcen);
	hcorrptrawpt[nj][centb]->Fill(rawpt,recopt,wxs*wcen);

        //! Very fine bin in ref pt
	hrescrpt[nj][centb]->Fill(refpt,recopt/refpt,wxs*wcen);
	hresrrpt[nj][centb]->Fill(refpt,rawpt/refpt,wxs*wcen);
	hresrcrpt[nj][centb]->Fill(recopt,recopt/rawpt,wxs*wcen);

	int ieta=-1;
	if(fabs(recoeta)<1.3)ieta=0; //! barrel region
        else ieta=1; //! HCAL region
        if(ieta>=0){
          hratiocorrrefpt_mb_eta[nj][ieta]->Fill(refpt,recopt/refpt,wxs*wcen);
          hratiorawrefpt_mb_eta [nj][ieta]->Fill(refpt,rawpt/refpt,wxs*wcen);
          hratiocorrrefpt_eta[nj][centb][ieta]->Fill(refpt,recopt/refpt,wxs*wcen);
          hratiorawrefpt_eta [nj][centb][ieta]->Fill(refpt,rawpt/refpt,wxs*wcen);
	}

	//! pileup study
        hjetptpumb[nj]->Fill(recopt,iJet->jtpu[ij],wxs*wcen);
        hjetptpu[nj][centb]->Fill(recopt,iJet->jtpu[ij],wxs*wcen);

        //! Response in different eta and phi bins
        int etabin = GetEtaBin(fabs(refeta));
	int phibin = GetPhiBin(refphi);

	if(etabin < 0 || etabin>=maxe || phibin < 0 || phibin>=maxph)continue;

        //! Response in eta and phi bins
        hpteta[nj][centb][etabin]->Fill(refpt,recopt/refpt,wxs*wcen);
        hptphi[nj][centb][phibin]->Fill(refpt,recopt/refpt,wxs*wcen);

        //! Relative resolution
	hdiffpt[nj][centb]        ->Fill(refpt,(recopt - refpt)/refpt,wxs*wcen);
        hdiffpteta[nj][centb][etabin]->Fill(refpt,(recopt - refpt)/refpt,wxs*wcen);
        hdiffptphi[nj][centb][phibin]->Fill(refpt,(recopt - refpt)/refpt,wxs*wcen);

        //! eta and phi resolutions
        hdiffeta[nj][centb]->Fill(refpt,(recoeta - refeta),wxs*wcen);
	hdiffphi[nj][centb]->Fill(refpt,(recophi - refphi),wxs*wcen);

        //! Gen:Reco correlation plots
	hgenjrecoj [nj][centb]->Fill(refpt,recopt,wxs*wcen);
	if(fabs(iJet->refparton_flavor[ij])<=21)hgenjrecoj_pflavor [nj][centb]->Fill(refpt,recopt,wxs*wcen);

	//! Away-side jet energy scale
        //! leading jet selection 
        //! for away-side jet selection 
        if(ljet[0]>=0 && ljet[0]!=ij){
          float dphi = recophi - iJet->jtphi[ljet[0]];
          if (dphi > pi ) dphi = dphi - 2 * pi;
          if (dphi < -pi) dphi = dphi + 2 * pi;
          //if (dphi < 0.5)continue;

          if(fabs(dphi)>kdphicut){ //! for dijet-jet
            hgenjrecoj_awayside [nj][centb]->Fill(refpt,recopt,wxs*wcen);
	    
            hjesasr[nj][centb]->Fill(iJet->refpt[ljet[0]],rawpt/refpt,wxs*wcen);
            hjesas [nj][centb]->Fill(iJet->refpt[ljet[0]],recopt/refpt,wxs*wcen);
            hjesas_eta[nj][centb][etabin]->Fill(iJet->refpt[ljet[0]],recopt/refpt,wxs*wcen);
	    
            float x = recopt/iJet->jtpt[ljet[0]];
            hxgj [nj][centb]->Fill(x,wxs*wcen);
            hdpgj[nj][centb]->Fill(dphi,wxs*wcen);

            int ibin = GetPtBin(iJet->jtpt[ljet[0]]);
	    if(ibin>=0)hxdp [nj][centb][ibin]->Fill(x,dphi,wxs*wcen);
          }
	}
      }//! ij loop

      hNevt[nj][centb]->Fill(vz);
      if(njets>0){
	hAvjets[nj][centb]->Fill(pthat,njets,wxs*wcen);
	hMeanPtjets[nj][centb]->Fill(pthat,meanpt/njets,wxs*wcen);
      }
      delete [] ljet;


    }//! nj jet loop

    if(istat){
      hBin->Fill(hiBin);
      hVz->Fill(vz);
      hHF->Fill(hiHF,wxs*wcen);
      iEvent++;
    }
    //std::cout<<"Completed event #  "<<ievt<<std::endl; 
  }//! event loop ends
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  std::cout<<"Events which passed the pT hat cut : "<<hTotEve->Integral()<<" out of  : "<<hEvt->Integral()
	   <<" efficiency of all the cuts : " <<hTotEve->Integral()/hEvt->Integral()<<std::endl;
  std::cout<<std::endl;

  if(strcmp(ksp,"pp")==0){
    for(int nj=0;nj<knj;nj++){
      std::cout<<"# of Events for : "<<c->GetAlgoName(nj)<<"\t"<<hNevt[nj][0]->Integral()<<"\t # of Jets : "<<hNjets[nj][0]->Integral()<<std::endl;
    }
  }else{
    for(int nj=0;nj<knj;nj++){
      for(int icen=0;icen<ncen;icen++){
	std::cout<<"# of Events for : "<<c->GetAlgoName(nj)<<"\t icen : "<<icen<<"\t"<<hNevt[nj][icen]->Integral()<<"\t # of Jets : "<<hNjets[nj][icen]->Integral()<<std::endl;
      }
      std::cout<<std::endl;
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

  std::cout<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<Form("You read %f Mbytes/RealTime seconds",mbytes/rtime)<<std::endl;
  std::cout<<Form("You read %f Mbytes/CpuTime  seconds",mbytes/ctime)<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;
  return 1;
}
int GetCentBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 2.5% bins of cross section
  //! in 0-39 bins

  if(ncen==6){
    /*
    if(bin<4)ibin=0; //! 0-10%
    else if(bin>=4 && bin<8)ibin=1;    //! 10-20%
    else if(bin>=8 && bin<12)ibin=2;   //! 20-30%
    else if(bin>=12 && bin<20)ibin=3;  //! 30-50%
    else if(bin>=20 && bin<28)ibin=4;  //! 50-70%
    else if(bin>=28 && bin<40)ibin=5;  //! 70-100%
    */

    if(bin<2)ibin=0; //! 0-5%
    else if(bin>=2 && bin<4)ibin=1;   //! 5-10%
    else if(bin>=4 && bin<12)ibin=2;  //! 10-30%
    else if(bin>=12&& bin<20)ibin=3;  //! 30-50%
    else if(bin>=20&& bin<28)ibin=4;  //! 50-70%
    else if(bin>=28&& bin<36)ibin=5;  //! 70-90%

  } else if(ncen==5){
    if(bin<4)ibin=0; //! 0-10%
    else if(bin>=4  && bin<8)ibin=1;  //! 10-20%
    else if(bin>=8  && bin<12)ibin=2; //! 20-30%
    else if(bin>=12 && bin<20)ibin=3; //! 30-50%
    else if(bin>=20 && bin<40)ibin=4; //! 50-100%

  }else if(ncen==7){
    if(bin<4)ibin=0; //! 0-10% 
    else if(bin>=4 && bin<8)ibin=1;    //! 10-20%
    else if(bin>=8 && bin<12)ibin=2;   //! 20-30%
    else if(bin>=12 && bin<16)ibin=3;  //! 30-40%
    else if(bin>=16 && bin<20)ibin=4;  //! 40-50%
    else if(bin>=20 && bin<24)ibin=5;  //! 50-60%
    else if(bin>=24 && bin<32)ibin=6;  //! 60-80%
    else if(bin>=32)ibin=7;            //! 80-100% 
  }
  return ibin;
}
int GetEtaBin(float eta)
{
  int ibin=-1;
  float min =0.0;
  ibin = (eta - min)*maxe/2.;
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
void FindLeadSubLeadJets(Jets *jetc, int *ljet)
{
  ljet[0]=-1; ljet[1]=-2;

  float tempt1=-9;
  float tempt2=-9;

  //! Get the leading jet                                                                                                                                                                                      
  for(int ij=0; ij<jetc->nref; ij++){
    //if(!selecJet(jetc,ij))continue;
    if(fabs(jetc->refdrjt[ij])>kdRcut || jetc->refparton_flavor[ij]==-999 || jetc->refparton_flavor[ij]==0)continue;
    if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<kptrecocut)continue;

    float jetpt = jetc->jtpt[ij];
    if(jetpt > tempt1){
      tempt1 = jetpt;
      ljet[0] = ij;
    }
  }

  for(int ij=0; ij<jetc->nref; ij++){
    //if(!selecJet(jetc,ij) || ij==ljet[0])continue;
    if(fabs(jetc->refdrjt[ij])>kdRcut || jetc->refparton_flavor[ij]==-999 || jetc->refparton_flavor[ij]==0 || ij==ljet[0])continue;
    if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<kptrecocut)continue;

    float jetpt = jetc->jtpt[ij];
    float dphi  = jetc->jtphi[ij] - jetc->jtphi[ljet[0]];
    if (dphi > pi ) dphi = dphi - 2 * pi;
    if (dphi < -pi) dphi = dphi + 2 * pi;
    //if(dphi < 2*pi/3.)continue;                                                                                                                                                                              

    if (jetpt > tempt2){
      tempt2 = jetpt;
      ljet[1] = ij;
    }
  }
}

