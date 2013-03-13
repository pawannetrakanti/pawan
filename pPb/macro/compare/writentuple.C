// From the prompt reco
#include "../../HiForest/V2/hiForest.h"
#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
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

const float ketacut=3.;
const double kptrecocut=30.;

void LoadLib();
void ShutoffBranches(HiForest */*hi*/,const char */*ksp*/);
void FindLeadSubLeadJets(Jets */*mJets*/, int */*ljet*/);

TStopwatch timer;
int writentuple(char *ksp="ppb")
{

  timer.Start();

  LoadLib();

  TString inname="";
  if(strcmp(ksp,"ppb")==0){
    //! Prompt reco for Jet80
    //!  cmsPfn /store/caf/user/velicanu/PA2013_merged_HiForest/pPb_hiForest2_1_15_test.root //! to get the path from CAF
    //inname="root://eoscms//eos/cms/store/caf/user/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_JSonPPb_forestv72.root";
    
    inname="root://eoscms//eos/cms/store/caf/user/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_JSonPPb_forestv77.root";

    //! Minbias
    //inname="root://eoscms//eos/cms/store/group/phys_heavyions/kjung/MinBiasUPCForest_v71/hiForest_MinBiasUPC_v71_Merged1.root";
    //inname="root://eoscms//eos/cms/store/group/phys_heavyions/kjung/MinBiasUPCForest_v71/hiForest_MinBiasUPC_v71_Merged2.root";
  }else  if(strcmp(ksp,"pbpb")==0){
    inname="root://eoscms//eos/cms/store/caf/user/yjlee/PbPHiForest2_PbPbPAHighPtJet80_cent50-100_pprereco.root";
  }else if(strcmp(ksp,"pp")==0){
    inname="root://eoscms//eos/cms/store/caf/user/yjlee/pp/hiForest2_pp_ppreco_415_90percent.root";
  }

  //! Load Lib
  //gSystem->Load("../HiForest/V2/hiForest_h.so");


  //! Define the input file and HiForest
  //! CMSSW_5_3_3
  HiForest *c = new HiForest(inname,Form("Forest%s",ksp),cPPb);
  cout<<"Loaded the hiforest tree : "<<c->GetName()<<endl;
  ShutoffBranches(c,ksp);



  //! Output file
  //! HIHighPt
  TFile *fout = new TFile("ntuple_pawan.root","RECREATE");  

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
  c->hasAk3CaloJetTree=0;
  c->hasAk5CaloJetTree=0;
  c->hasAk6CaloJetTree=0;

  c->hasAkPu1CaloJetTree=0;
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  c->hasAkPu5CaloJetTree=0;
  c->hasAkPu6CaloJetTree=0;

  c->hasAk1PFJetTree=0;
  c->hasAk2PFJetTree=0;
  c->hasAk4PFJetTree=0;
  c->hasAk5PFJetTree=0;
  c->hasAk6PFJetTree=0;

  c->hasAkPu1PFJetTree=0;
  c->hasAkPu2PFJetTree=0;
  c->hasAkPu4PFJetTree=0;
  c->hasAkPu5PFJetTree=0;
  c->hasAkPu6PFJetTree=0;

  c->hasTrackTree=0;

  //! For jets
  Jets *mJets=0;
  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;

  string jetVars = "";
  jetVars += "evt:run:vz:jet80:noise:ntrk:npix:npixtrk:hf:hfp:hfm:hfp4:hfm4:zdc:zdcp:zdcm:jtpt1:raw1:jteta1:jtphi1:trkMax1:jtpt2:raw2:jteta2:jtphi2:trkMax2:jtpt3:raw3:jteta3:jtphi3:trkMax3";
  TNtuple *ntjet=0;
  ntjet = new TNtuple("ntjet","",jetVars.data());

  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
    //! load the hiForest event
    c->GetEntry(ievt);
    
    //! only selecting pPb runs
    if(c->evt.run>211256)continue;

    //! events with Single vertex
    bool evSel = false;
    evSel = fabs(c->evt.vz)<15. && c->skim.pHBHENoiseFilter  && c->skim.pPAcollisionEventSelectionPA && c->skim.pVertexFilterCutGplus && c->hlt.HLT_PAJet80_NoJetID_v1;
    if(!evSel)continue;


    double noise=-1;
    double pt1 = -9, pt2 = -9, pt3 = -9,
      raw1 = -9, raw2 = -9, raw3 = -9,
      eta1 = -9,eta2 = -9, eta3  = -9,
      phi1 = -9,phi2 = -9, phi3  = -9,
      trkMax1 = -9,trkMax2=-9,trkMax3=-9;
    


    double run   = c->evt.run;
    double evt   = c->evt.evt;
    double vz    = c->evt.vz;
    double jet80 = c->hlt.HLT_PAJet80_NoJetID_v1;

    double ntrk    = c->evt.hiNtracks;
    double npix    = c->evt.hiNpix;
    double npixtrk = c->evt.hiNpixelTracks;

    double hf   = c->evt.hiHF;
    double hfp  = c->evt.hiHFplus;
    double hfm  = c->evt.hiHFminus;
    double hfp4 = c->evt.hiHFplusEta4;
    double hfm4 = c->evt.hiHFminusEta4;

    double zdc  = c->evt.hiZDC;
    double zdcp = c->evt.hiZDCplus;
    double zdcm = c->evt.hiZDCminus;

    if(ievt%10000==0)std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<run<<std::endl;
    
    mJets = &(c->akPu3PF);
    
    int *ljet = new int[3];
    FindLeadSubLeadJets(mJets,ljet);

    int jtLead = -1, jtSubLead = -1, jtThird = -1;

    if(ljet[0] >=0 ) jtLead    = ljet[0];
    if(ljet[1] >=0 ) jtSubLead = ljet[1];    
    if(ljet[2] >=0 ) jtThird   = ljet[2];    
    
    if(jtLead<0)continue;
    if(mJets->jtpt[jtLead] < 80)continue;

    if(jtLead > -1){
      pt1     = mJets->jtpt[jtLead];
      eta1    = mJets->jteta[jtLead];
      phi1    = mJets->jtphi[jtLead];
      raw1    = mJets->rawpt[jtLead];
      trkMax1 = mJets->trackMax[jtLead];
    }

    if(jtSubLead > -1){
      pt2     = mJets->jtpt[jtSubLead];
      eta2    = mJets->jteta[jtSubLead];
      phi2    = mJets->jtphi[jtSubLead];
      raw2    = mJets->rawpt[jtSubLead];
      trkMax2 = mJets->trackMax[jtSubLead];
    }

    if(jtThird > -1){
      pt3     = mJets->jtpt[jtThird];
      eta3    = mJets->jteta[jtThird];
      phi3    = mJets->jtphi[jtThird];
      raw3    = mJets->rawpt[jtThird];
      trkMax3 = mJets->trackMax[jtThird];
    }

    float jentry[] = {evt,run,vz,jet80,noise,ntrk,npix,npixtrk,hf,hfp,hfm,hfp4,hfm4,zdc,zdcp,zdcm,
		      pt1,raw1,eta1,phi1,trkMax1,
		      pt2,raw2,eta2,phi2,trkMax2,
		      pt3,raw3,eta3,phi3,trkMax3};
    
    ntjet->Fill(jentry);
    delete [] ljet;
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
void FindLeadSubLeadJets(Jets *jetc, int *ljet)
{
  ljet[0]=-1; ljet[1]=-2; ljet[2]=-3;
  
  double tempt=-9;
  //! Get the leading jet
  for(int ij=0; ij<jetc->nref; ij++){
    //if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<100 || jetc->rawpt[ij]<15 || jetc->trackMax[ij]<4.)continue;
    if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<100 || jetc->rawpt[ij]<15)continue;
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
      //if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15 || jetc->trackMax[ij]<4.)continue;
      if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15)continue;
      double jetpt = jetc->jtpt[ij];
      if (jetpt > tempt){
	tempt = jetpt;
	ljet[1] = ij;
      }
    }

    if(ljet[1]>=0){
      // third jet
      tempt=-9;
      for(int ij=0; ij<jetc->nref; ij++){
	if(ij==ljet[0] || ij==ljet[1])continue;
	//if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15 || jetc->trackMax[ij]<4.)continue;
	if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15)continue;
	double jetpt = jetc->jtpt[ij];
	if (jetpt > tempt){
	  tempt = jetpt;
	  ljet[2] = ij;
	}
      }
    }//! third jet
  }
}
void LoadLib()
{
  gSystem->Load("../../HiForest/V2/hiForest_h.so");
}
void ShutoffBranches(HiForest *hi,const char *ksp)
{

  //! added by pawan
  //! Select only the branches you want to use for the analysis
  //! This increases the speed for running

  //! For Hlt
  hi->hltTree->SetBranchStatus("*",0,0);
  if(strcmp(ksp,"ppb")==0){//! pPb
    hi->hltTree->SetBranchStatus("HLT_PAJet40_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet60_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet80_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet100_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAZeroBiasPixel_SingleTrack_v1",1,0);
  }else if(strcmp(ksp,"pbpb")==0){  //! PbPb
    hi->hltTree->SetBranchStatus("HLT_HIJet80_v1",1,0);
  }
  else if(strcmp(ksp,"pp")==0){  //! pp
    hi->hltTree->SetBranchStatus("HLT_Jet60_v1",1,0);
  }

  //! for Skim Tree
  hi->skimTree->SetBranchStatus("*",0,0);
  hi->skimTree->SetBranchStatus("pHBHENoiseFilter",1,0);
  if(strcmp(ksp,"ppb")==0){
    hi->skimTree->SetBranchStatus("pVertexFilterCutGplus",1,0);
    hi->skimTree->SetBranchStatus("pPAcollisionEventSelectionPA",1,0);
  }else if(strcmp(ksp,"pbpb")==0 || strcmp(ksp,"pp")==0){
    hi->skimTree->SetBranchStatus("pcollisionEventSelection",1,0);
  }




  //! Evt tree
  hi->evtTree->SetBranchStatus("*",0,0);
  hi->evtTree->SetBranchStatus("run",1,0);
  hi->evtTree->SetBranchStatus("evt",1,0);
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
  hi->ak3PFJetTree->SetBranchStatus("*",0,0);
  hi->ak3PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak3PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak3PFJetTree->SetBranchStatus("trackMax",1,0);
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

  hi->akPu3PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu3PFJetTree->SetBranchStatus("nref",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("trackMax",1,0);

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

  hi->ak3CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak3CaloJetTree->SetBranchStatus("nref",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("trackMax",1,0);

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

  hi->akPu3CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu3CaloJetTree->SetBranchStatus("nref",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("trackMax",1,0);

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
