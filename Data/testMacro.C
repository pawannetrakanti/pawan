#include "HiForest/hiForest.h"
#include "SmearingFactors_New.h"
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

void testMacro(const char *ksp="pp")
{

  //! Load Lib
  gSystem->Load("HiForest/hiForest_h.so");

  //! Load Smearing factors
  LoadParameters();

  float ketacut=2.;
  double kptrecocut=60.;


  //! Define the input file and HiForest
  TString inname="";
  bool ispp=false;
  if(strcmp(ksp,"pbpb")==0){
    //! merged from prompt reco
    inname="/net/hidsk0001/d00/scratch/yjlee/merge/pbpbDijet_v20/promptskim-hihighpt-hltjet80-pt90-v20.root";
    kptrecocut=80;
  }else if(strcmp(ksp,"pp")==0){
    //! New production
    inname="/net/hisrv0001/home/yenjie/scratch/hiForest/prod/productionPP/CMSSW_4_4_2_patch6/test/ppForest2/pp_merged_full.root";
    ispp=true;
  }

  HiForest *c = new HiForest(inname,Form("MyForest%s",ksp),ispp,0);

  //! shutting off trees which we do not use now
  c->hasIcPu5JetTree=0;
  c->hasAkPu2JetTree=0;
  c->hasAkPu4JetTree=0;
  c->hasAk2JetTree=0;
  c->hasAk4JetTree=0;
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  c->hasAk2CaloJetTree=0;
  c->hasAk3CaloJetTree=0;
  c->hasAk4CaloJetTree=0;
  c->hasPhotonTree=0;
  c->hasTrackTree=0;
  c->hasPixTrackTree=0;
  c->hasTowerTree=0;
  c->hasHbheTree=0;
  c->hasEbTree=0;
  c->hasGenpTree=0;
  c->hasGenParticleTree=0;
  c->hasPFTree=0;
  c->hasTrackTree=0;
  c->hasMetTree=0;

  if(strcmp(ksp,"pp")==0)c->hasAkPu3JetTree=0;
  else c->hasAk3JetTree=0;

  Jets *mJets=0; //! to grab jet objects

  Long64_t nb = 0;
  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;
  
  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
    //! load the hiForest event
    nb += c->GetEntry(ievt);

    int hiBin          = c->evt.hiBin;
    if(ievt%10000==0)std::cout<<" ********** Event # " <<ievt<<"\t hiBin : "<<hiBin<<std::endl;            

    //! select good events
    bool goodEvent = false;
    if(ispp){
      goodEvent = c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_Jet60_v1 && fabs(c->evt.vz)<15.;
    }else{
      goodEvent = c->skim.pcollisionEventSelection && c->skim.pHBHENoiseFilter && c->hlt.HLT_HIJet80_v1 && fabs(c->evt.vz)<15.; 
    }
    if(!goodEvent)continue;

    
    if(strcmp(ksp,"pp")==0){      //! assign the jets for pp
      mJets = &(c->ak3PF);
    }else{      //! assign the jets for PbPb
      mJets = &(c->akPu3PF);
    }

    //! Loop over # of jets in a given algo
    for(int ij=0; ij<mJets->nref; ij++){

      //! jet selction cut
      if(mJets->jtpt[ij]<kptrecocut || fabs(mJets->jteta[ij])>ketacut)continue;
      
	//! trackMax cuts
      if(mJets->trackMax[ij]/mJets->jtpt[ij] < 0.01)continue;

      //! For MC we have used the refpt for the
      //! Matched gen-jet only: gensubid==0 && genamatchedindex>0; in ngen loop for PbPb

      float recopt = mJets->jtpt[ij];

      //! pp pt smeared with PbPb resolution
      if(strcmp(ksp,"pp")==0){
	for(int ic=0;ic<6;ic++){

	  //! The first argument (2,..,..) is for ak3PF jets in pp
	  //! Currently implemented for only ak3PF jets 

	  //! Shifts the jet energy scale and smears the reco jet pt
	  float smpt0  = GetSmearedPtData(2,ic,recopt,0,"no");  //! Data
	  //float smpt0  = GetSmearedPtMC(2,ic,recopt,refpt);   //! MC

	  //! All no jet energy scale shift only smearing;
	  float smpt1  = GetSmearedPtData_NoMeanShift(2,ic,recopt,0,"no"); //! Data
	  //float smpt1  = GetSmearedPtMC_NoMeanShift(2,ic,recopt,refpt);  //! MC

	  //! All only jet energy scale shift no smearing;
	  float smpt2  = GetSmearedPtData_OnlyMeanShift(2,ic,recopt,0,"no");  //! Data
	  //float smpt2  = GetSmearedPtMC_OnlyMeanShift(2,ic,recopt,refpt);   //! MC

	  //! All the after burners are shut off (both for resolution and scale);
	  float smpt3  = GetSmearedPtData_woAfBurn(2,ic,recopt,0,"no");  //! Data
	  //float smpt3  = GetSmearedPtMC_woAfBurn(2,ic,recopt,refpt);   //! MC 
	  
	  //! If you are using 4 centrality bins (0-10%, 10-30%, 30-50% and 50-100%) then use these factors
	  //! Re-weighting factor for data
	  //! If you want you can use your own re-weighting factors and ignore them.
	  float rewe0 = GetReWeight(2,ic,smpt0);             //! JFF & JS
	  float rewe1 = GetReWeight_NoMeanShift(2,ic,smpt1); //! JFF & JS

	  //! Add these piece of code in your smearing loop
	  int nic=-1;
	  if(ic==0 || ic==1)nic=0; //! 0-10%
	  else if(ic==2 || ic==3)nic=ic-1; //! 10-30% and 30-50%
	  else if(ic==4 || ic==5)nic=3; //! 50-100%

	  //! To use them
	  //! hMyHist[nic]->Fill(smpt0,rewe0);
	}
      }else  if(strcmp(ksp,"pbpb")==0){
	//! If you want to correct the jet energy scale for PbPb:
	//! Note it takes input as hiBin from the hiForest and not the centrality bins which
	//! we use.
	recopt = GetPbPbCorrectedScaleData(2,hiBin,recopt);       //! Data
	//recopt = GetPbPbCorrectedScaleMC(2,hiBin,recopt,refpt); //! MC
      }
    }//! # jet loop
  }//! event loop ends
}
