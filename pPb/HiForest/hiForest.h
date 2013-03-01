#ifndef HIFOREST_H
#define HIFOREST_H

#include <iostream>
#include <vector>
#include <algorithm>

#include "commonTool.h"
#include "SetupPhotonTree.h"
#include "SetupJetTree.h"
#include "SetupHltTree.h"
#include "SetupSkimTree.h"
#include "SetupTrackTree.h"
#include "SetupHitTree.h"
#include "SetupEvtTree.h"
#include "SetupMetTree.h"
#include "SetupGenpTree.h"
#include "SetupPFTree.h"
#include "SetupGenParticleTree.h"
#include "TrackingCorrectionsv6.h"
#include "TrackingParam.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TF1.h>
#include <TCut.h>

#include "DummyJetCorrector.h"

// ==========================================================
// Main class which can be used to read the hiForest trees
//
// Author: Yen-Jie Lee
//
// ==========================================================

namespace names{
  enum Algo{
  };
  string AlgoRename[100] = {
  };
  string AlgoAnalyzer[100] = {
  };
}

enum collisionType { cPbPb, cPP, cPPb };

class HiForest : public TNamed
{

  using TObject::Draw;

  public: 
   HiForest(const char *file, const char *name="forest", collisionType cMode = cPbPb, bool ismc = 0, bool isrecorrected = 0);
  ~HiForest();

  //==================================================================================================================================
  // Utility functions
  //==================================================================================================================================
  void CheckTree(TTree *t,const char *title);			// Check the status of a tree
  //void CheckArraySizes();					// Check if array size is large enough
  void GetEntry(int i);
  int  GetEntries();  						// Get the number of entries 

  //void GetEnergyScaleTable(char *fileNameTable);                // Get photon energy scale table
  //void InitTree();						// Initialize track correction

  void PrintStatus();						// Print the status of the hiForest
  void SetOutputFile(const char *name);               		// Set output file name for skim
  void AddCloneTree(TTree* t, const char *dirName, const char *treeName);   // Add a clone tree to the clone forest
  void FillOutput();						// Fill output forest  
  
  Long64_t Draw(const char* varexp, const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0){
     return tree->Draw(varexp,selection,option,nentries,firstentry);
  }

  //==================================================================================================================================
  // Event filtering utility functions
  //==================================================================================================================================
  bool selectEvent();
  TCut eventSelection();

  //==================================================================================================================================
  // Photon utility functions
  //==================================================================================================================================
  /*
  bool isSpike(int i);                          		// return true if it is considered as a spike candidate
  bool isLooseEGamma(int i);
  bool isLoosePhoton(int i);
  bool isGoodPhoton(int i);                   			// return true if it is considered as a hiGoodPhoton candidate
  bool isIsolatedPhoton(int i);                     		// return true if it is considered as a hiGoodPhoton candidate
  bool isMCSignal(int i);
  bool isDirectPhoton(int i);
  bool isFragPhoton(int i);
  float getCorrEt(int j);                                       // Get Corrected photon Et
  */

  //==================================================================================================================================
  // Jet utility functions
  //==================================================================================================================================
  /*
  void sortJets(TTree* jetTree, Jets& jets, double etaMax = 2, double ptMin = 40, bool allEvents = 1, int smearType = -1);
  int leadingJet();
  int subleadingJet();
  int thirdJet();
  double deltaPhiDijet(Jets& jets);
  bool hasDiJet(Jets& jets, double pt1 = 100, double pt2 = 40, double dphiMin = 2.*3.1415926/3.);
  void fakeRejection(TTree *jetTree, Jets &jets, bool allEvents);
  */
  //==================================================================================================================================
  // Track utility functions
  //==================================================================================================================================
  /*
  int getMatchedCaloTowerAllowReuse(int j);
  int getMatchedHBHEAllowReuse(int j);
  void matchTrackCalo(bool allEvents = 1);
  double getTrackCorrectionPara(int j);
  double getTrackCorrection(int j);
  bool selectTrack(int j);
  */

  //==================================================================================================================================
  // Get track-jet correlated variables. Not needed if correlatePF is run.
  //==================================================================================================================================
  //void correlateTracks(TTree* jetTree, Jets& jets, bool allEvents = 1, bool smeared = 0);

  // TFile
  TFile *inf; 					// Input file 
  TFile *outf;                                  // Output file if we want to export the forest

  // Trees

  TTree *evtTree;                               // Event Tree
  TTree *icPu5jetTree;				// Jet Tree with icPu5 algorithm, see branches in SetupJetTree.h

  TTree *ak1PFJetTree;			// Jet Tree with ak1PF algorithm, see branches in SetupJetTree.h
  TTree *ak2PFJetTree;			// Jet Tree with ak2PF algorithm, see branches in SetupJetTree.h
  TTree *ak3PFJetTree;			// Jet Tree with ak3PF algorithm, see branches in SetupJetTree.h
  TTree *ak4PFJetTree;			// Jet Tree with ak4PF algorithm, see branches in SetupJetTree.h
  TTree *ak5PFJetTree;			// Jet Tree with ak5PF algorithm, see branches in SetupJetTree.h
  TTree *ak6PFJetTree;			// Jet Tree with ak6PF algorithm, see branches in SetupJetTree.h

  TTree *akPu1PFJetTree;			// Jet Tree with akPu1PF algorithm, see branches in SetupJetTree.h
  TTree *akPu2PFJetTree;			// Jet Tree with akPu2PF algorithm, see branches in SetupJetTree.h
  TTree *akPu3PFJetTree;			// Jet Tree with akPu3PF algorithm, see branches in SetupJetTree.h
  TTree *akPu4PFJetTree;			// Jet Tree with akPu4PF algorithm, see branches in SetupJetTree.h
  TTree *akPu5PFJetTree;			// Jet Tree with akPu5PF algorithm, see branches in SetupJetTree.h
  TTree *akPu6PFJetTree;			// Jet Tree with akPu6PF algorithm, see branches in SetupJetTree.h


  TTree *akPu1CaloJetTree;			// Jet Tree with akPu1Calo algorithm, see branches in SetupJetTree.h
  TTree *akPu2CaloJetTree;			// Jet Tree with akPu2Calo algorithm, see branches in SetupJetTree.h
  TTree *akPu3CaloJetTree;			// Jet Tree with akPu3Calo algorithm, see branches in SetupJetTree.h
  TTree *akPu4CaloJetTree;			// Jet Tree with akPu4Calo algorithm, see branches in SetupJetTree.h
  TTree *akPu5CaloJetTree;			// Jet Tree with akPu5Calo algorithm, see branches in SetupJetTree.h
  TTree *akPu6CaloJetTree;			// Jet Tree with akPu6Calo algorithm, see branches in SetupJetTree.h

  TTree *ak1CaloJetTree;			// Jet Tree with ak1Calo algorithm, see branches in SetupJetTree.h
  TTree *ak2CaloJetTree;			// Jet Tree with ak2Calo algorithm, see branches in SetupJetTree.h
  TTree *ak3CaloJetTree;			// Jet Tree with ak3Calo algorithm, see branches in SetupJetTree.h
  TTree *ak4CaloJetTree;			// Jet Tree with ak4Calo algorithm, see branches in SetupJetTree.h
  TTree *ak5CaloJetTree;			// Jet Tree with ak5Calo algorithm, see branches in SetupJetTree.h
  TTree *ak6CaloJetTree;			// Jet Tree with ak5Calo algorithm, see branches in SetupJetTree.h

  TTree *hltTree;				// OpenHLT Tree, see branches in SetupHltTree.h
  TTree *skimTree;				// Skim Tree, contains event selection info, see branches in SetupSkimTree.h
  TTree *trackTree;				// Track Tree, see branches in SetupTrackTree.h

  /*
  TTree *photonTree;				// Photon Tree, see branches in SetupPhotonTree.h
  TTree *pixtrackTree;				// Track Tree, see branches in SetupTrackTree.h
  TTree *towerTree;                             // Tower Tree
  TTree *hbheTree;                              // HCAL HBHE Tree
  TTree *ebTree;                                // ECAL eb Tree

  TTree *metTree;                               // MET Tree
  TTree *pfTree;                                // PF candidate Tree, see branches in SetupPFTree.h
  TTree *genpTree;                              // Gen photon of the signal event Tree
  TTree *genParticleTree;                       // Stable Gen particles
  */

  TTree *tree;					// Pointer to the available tree, all trees in the forest are friended to each other

  //vector<TTree*> jetTrees;
  vector<TTree*> cloneForest;                   // Vector of clones for skim

  TF1* fGauss;

  // Branches
  Evts evt;
  Hlts hlt;
  Skims skim;

  //vector<Jets> alljets;

  Jets icPu5;

  Jets akPu1PF;
  Jets akPu2PF;
  Jets akPu3PF;
  Jets akPu4PF;
  Jets akPu5PF;
  Jets akPu6PF;

  Jets ak1PF;
  Jets ak2PF;
  Jets ak3PF;
  Jets ak4PF;
  Jets ak5PF;
  Jets ak6PF;


  Jets akPu1Calo;
  Jets akPu2Calo;
  Jets akPu3Calo;
  Jets akPu4Calo;
  Jets akPu5Calo;
  Jets akPu6Calo;

  Jets ak1Calo;
  Jets ak2Calo;
  Jets ak3Calo;
  Jets ak4Calo;
  Jets ak5Calo;
  Jets ak6Calo;

  Tracks track;

  /*
  Photons photon;
  Tracks pixtrack;
  Hits tower;
  Hits hbhe;
  Hits eb;

  Mets met;
  PFs pf;
  Genps genp;
  GenParticles genparticle;
  */

  // Booleans
  bool hasEvtTree;
  bool hasHltTree;
  bool hasSkimTree;
  bool hasTrackTree;

  bool hasIcPu5JetTree;

  bool hasAk1PFJetTree;
  bool hasAk2PFJetTree;
  bool hasAk3PFJetTree;
  bool hasAk4PFJetTree;
  bool hasAk5PFJetTree;
  bool hasAk6PFJetTree;

  bool hasAkPu1PFJetTree;
  bool hasAkPu2PFJetTree;
  bool hasAkPu3PFJetTree;
  bool hasAkPu4PFJetTree;
  bool hasAkPu5PFJetTree;
  bool hasAkPu6PFJetTree;


  bool hasAk1CaloJetTree;
  bool hasAk2CaloJetTree;
  bool hasAk3CaloJetTree;
  bool hasAk4CaloJetTree;
  bool hasAk5CaloJetTree;
  bool hasAk6CaloJetTree;


  bool hasAkPu1CaloJetTree;
  bool hasAkPu2CaloJetTree;
  bool hasAkPu3CaloJetTree;
  bool hasAkPu4CaloJetTree;
  bool hasAkPu5CaloJetTree;
  bool hasAkPu6CaloJetTree;

  /*
  bool hasPhotonTree;
  bool hasMetTree;
  bool hasPFTree;
  bool hasPixTrackTree;
  bool hasTowerTree;
  bool hasHbheTree;
  bool hasEbTree;
  bool hasGenpTree;
  bool hasGenParticleTree;
  */

  bool setupOutput;
  bool verbose;
  collisionType collisionMode;
  bool mc;


  bool doJetCorrection;
  bool doTrackCorrections;
  bool doTrackingSeparateLeadingSubleading;

  /*
  // Extra variables
  Float_t* towerEt;
  Float_t* towerdR;
  Float_t* hbheEt;
  Float_t* hbhedR;

  Float_t* jtChg;
  Float_t* jtNeut;
  Float_t* jtEM;

  Float_t* jtChgGen;
  Float_t* jtNeutGen;
  Float_t* jtEMGen;

  Float_t* jtPtMax;
  Float_t* jtPtMean;
  Float_t* jtPtMeanW;

  Int_t* jtLeadType;

  Int_t jtLead;
  Int_t jtSubLead;
  bool jtHasDijet;
  bool jtHasLeadingJet;

  Float_t* tjDeltaEtaLead;
  Float_t* tjDeltaPhiLead;
  Float_t* zLead;

  Float_t* tjDeltaEtaSubLead;
  Float_t* tjDeltaPhiSubLead;
  Float_t* zSubLead;

  Float_t* zOldLead;
  Float_t* zOldSubLead;

  Float_t* zSingleLead;
  Float_t* zLabLead;
  Float_t* zSingleSubLead;
  Float_t* zLabSubLead;
  Float_t* tjDeltaThetaLead;
  Float_t* tjDeltaThetaLabLead;
  Float_t* tjDeltaThetaSingleLead;

  Float_t* tjDeltaThetaSubLead;
  Float_t* tjDeltaThetaLabSubLead;
  Float_t* tjDeltaThetaSingleSubLead;

  Float_t* corrLead;
  Float_t* corrSubLead;
  */

  Float_t minJetPtForTrkCor;
  Float_t leadingJetPtForTrkCor;
  Float_t subleadingJetPtForTrkCor;
  Float_t leadingJetEtaForTrkCor;
  Float_t subleadingJetEtaForTrkCor;
  Float_t leadingJetPhiForTrkCor;
  Float_t subleadingJetPhiForTrkCor;

  int nEntries;
  int currentEvent;
  double cone;
  bool initialized;

  //vector<JetCorrectorParameters> vpar_HI310x;
  //  FactorizedJetCorrector *_JEC_HI310X;

//  TrackingParam trackCorrFromParam("TrkCorr_3DHistos_pThat80_Inclusive.root","CentralityWeights.root","PtResidualWeights.root");
//  TrackingParam *trackCorrFromParam;

//  vector<TrackingCorrections*> trackCorrections;
  //TF1* fEnergyScale[2][10];  // [a][b],  a =0 for unconverted,  a=1 for converted.   b: 1,2,3 is centrality bin. b=0 is empty                                                                                       
  
 private:
  
  ClassDef(HiForest,8)  
};


HiForest::HiForest(const char *infName, const char* name, collisionType cMode, bool ismc, bool recjec):
   tree(0),
   fGauss(0),
   verbose(0),
   collisionMode(cMode),
   mc(ismc),
   nEntries(0),
   currentEvent(0)
{
   tree = new TTree("tree","");
  SetName(name);
  // Input file
  inf = TFile::Open(infName);

  cone = 0.3;
  doTrackCorrections = 0;

  // Track correction initialized?
  initialized = 0;

  // Print out collision mode:
  cout <<"Collision Mode:";
  if (collisionMode == cPP) cout <<" P+P"<<endl;  
  if (collisionMode == cPPb) cout <<" P+Pb"<<endl;  
  if (collisionMode == cPbPb) cout <<" Pb+Pb"<<endl;  

  // Load trees. Hard coded for the moment
  hltTree          = (TTree*) inf->Get("hltanalysis/HltTree");
  skimTree         = (TTree*) inf->Get("skimanalysis/HltTree");
  //photonTree       = (TTree*) inf->Get("multiPhotonAnalyzer/photon");
  if (collisionMode == cPbPb || collisionMode == cPP) trackTree        = (TTree*) inf->Get("anaTrack/trackTree");
  if (collisionMode == cPPb) trackTree        = (TTree*) inf->Get("ppTrack/trackTree");
  //towerTree        = (TTree*) inf->Get("rechitanalyzer/tower");

  icPu5jetTree     = (TTree*) inf->Get("icPu5JetAnalyzer/t");

  ak1PFJetTree = (TTree*) inf->Get("ak1PFJetAnalyzer/t");
  ak2PFJetTree = (TTree*) inf->Get("ak2PFJetAnalyzer/t");
  ak3PFJetTree = (TTree*) inf->Get("ak3PFJetAnalyzer/t");
  ak4PFJetTree = (TTree*) inf->Get("ak4PFJetAnalyzer/t");
  ak5PFJetTree = (TTree*) inf->Get("ak5PFJetAnalyzer/t");
  ak6PFJetTree = (TTree*) inf->Get("ak6PFJetAnalyzer/t");
  akPu1PFJetTree = (TTree*) inf->Get("akPu1PFJetAnalyzer/t");
  akPu2PFJetTree = (TTree*) inf->Get("akPu2PFJetAnalyzer/t");
  akPu3PFJetTree = (TTree*) inf->Get("akPu3PFJetAnalyzer/t");
  akPu4PFJetTree = (TTree*) inf->Get("akPu4PFJetAnalyzer/t");
  akPu5PFJetTree = (TTree*) inf->Get("akPu5PFJetAnalyzer/t");
  akPu6PFJetTree = (TTree*) inf->Get("akPu6PFJetAnalyzer/t");


  ak1CaloJetTree = (TTree*) inf->Get("ak1CaloJetAnalyzer/t");
  ak2CaloJetTree = (TTree*) inf->Get("ak2CaloJetAnalyzer/t");
  ak3CaloJetTree = (TTree*) inf->Get("ak3CaloJetAnalyzer/t");
  ak4CaloJetTree = (TTree*) inf->Get("ak4CaloJetAnalyzer/t");
  ak5CaloJetTree = (TTree*) inf->Get("ak5CaloJetAnalyzer/t");
  ak6CaloJetTree = (TTree*) inf->Get("ak6CaloJetAnalyzer/t");

  akPu1CaloJetTree = (TTree*) inf->Get("akPu1CaloJetAnalyzer/t");
  akPu2CaloJetTree = (TTree*) inf->Get("akPu2CaloJetAnalyzer/t");
  akPu3CaloJetTree = (TTree*) inf->Get("akPu3CaloJetAnalyzer/t");
  akPu4CaloJetTree = (TTree*) inf->Get("akPu4CaloJetAnalyzer/t");
  akPu5CaloJetTree = (TTree*) inf->Get("akPu5CaloJetAnalyzer/t");
  akPu6CaloJetTree = (TTree*) inf->Get("akPu6CaloJetAnalyzer/t");


  //hbheTree         = (TTree*) inf->Get("rechitanalyzer/hbhe");
  //ebTree           = (TTree*) inf->Get("rechitanalyzer/eb");
  evtTree          = (TTree*) inf->Get("hiEvtAnalyzer/HiTree");
  //metTree          = (TTree*) inf->Get("anaMET/metTree");
  //pfTree           = (TTree*) inf->Get("pfcandAnalyzer/pfTree");
  //genpTree         = (TTree*) inf->Get("genpana/photon");
  
  // doesn't load genParticle by default
  //genParticleTree  = (TTree*) inf->Get("HiGenParticleAna/hi");

  // Check the validity of the trees.
  //hasPhotonTree        = (photonTree       	!= 0);
  //hasPFTree            = (pfTree     		!= 0);
  hasEvtTree           = (evtTree      		!= 0);
  //hasMetTree           = (metTree      		!= 0);
  hasIcPu5JetTree      = (icPu5jetTree 		!= 0);


  hasAk1PFJetTree  = (ak1PFJetTree 	!= 0);
  hasAk2PFJetTree  = (ak2PFJetTree 	!= 0);
  hasAk3PFJetTree  = (ak3PFJetTree 	!= 0);
  hasAk4PFJetTree  = (ak4PFJetTree 	!= 0);
  hasAk5PFJetTree  = (ak5PFJetTree 	!= 0);
  hasAk6PFJetTree  = (ak6PFJetTree 	!= 0);
  

  hasAkPu1PFJetTree  = (akPu1PFJetTree 	!= 0);
  hasAkPu2PFJetTree  = (akPu2PFJetTree 	!= 0);
  hasAkPu3PFJetTree  = (akPu3PFJetTree 	!= 0);
  hasAkPu4PFJetTree  = (akPu4PFJetTree 	!= 0);
  hasAkPu5PFJetTree  = (akPu5PFJetTree 	!= 0);
  hasAkPu6PFJetTree  = (akPu6PFJetTree 	!= 0);


  hasAk1CaloJetTree  = (ak1CaloJetTree 	!= 0);
  hasAk2CaloJetTree  = (ak2CaloJetTree 	!= 0);
  hasAk3CaloJetTree  = (ak3CaloJetTree 	!= 0);
  hasAk4CaloJetTree  = (ak4CaloJetTree 	!= 0);
  hasAk5CaloJetTree  = (ak5CaloJetTree 	!= 0);
  hasAk6CaloJetTree  = (ak6CaloJetTree 	!= 0);
  

  hasAkPu1CaloJetTree  = (akPu1CaloJetTree 	!= 0);
  hasAkPu2CaloJetTree  = (akPu2CaloJetTree 	!= 0);
  hasAkPu3CaloJetTree  = (akPu3CaloJetTree 	!= 0);
  hasAkPu4CaloJetTree  = (akPu4CaloJetTree 	!= 0);
  hasAkPu5CaloJetTree  = (akPu5CaloJetTree 	!= 0);
  hasAkPu6CaloJetTree  = (akPu6CaloJetTree 	!= 0);

  hasTrackTree     = (trackTree    		!= 0);
  hasHltTree       = (hltTree    		!= 0);
  hasSkimTree      = (skimTree   		!= 0);
  //hasTowerTree     = (towerTree    		!= 0);
  //hasHbheTree      = (hbheTree     		!= 0);
  //hasEbTree        = (ebTree       		!= 0);
  //hasGenpTree	   = (genpTree     		!=0);
  //hasGenParticleTree = (genParticleTree   	!=0);
  setupOutput = false;
  
  // Setup branches. See also Setup*.h

  if (hasEvtTree) {
    evtTree->SetName("evtTree");
    if (tree == 0) tree = evtTree; else tree->AddFriend(evtTree);
    setupEvtTree(evtTree,evt);
  }


  if (hasHltTree) {
    hltTree->SetName("hltTree");
    if (tree == 0) tree = hltTree; else tree->AddFriend(hltTree);
    setupHltTree(hltTree,hlt);
  }

  if (hasSkimTree) {
    skimTree->SetName("skim");
    if (tree == 0) tree = skimTree; else tree->AddFriend(skimTree);
    setupSkimTree(skimTree,skim);
  }

  if (hasIcPu5JetTree) {
    icPu5jetTree->SetName("icPu5");
    if (tree == 0) tree = icPu5jetTree; else tree->AddFriend(icPu5jetTree);
    setupJetTree(icPu5jetTree,icPu5);
  }
 
  //////// PF jets ////////
  if (hasAkPu1PFJetTree) {
    akPu1PFJetTree->SetName("akPu1PF");
    if (tree == 0) tree = akPu1PFJetTree; else tree->AddFriend(akPu1PFJetTree);
    setupJetTree(akPu1PFJetTree,akPu1PF);
  }
  if (hasAkPu2PFJetTree) {
    akPu2PFJetTree->SetName("akPu2PF");
    if (tree == 0) tree = akPu2PFJetTree; else tree->AddFriend(akPu2PFJetTree);
    setupJetTree(akPu2PFJetTree,akPu2PF);
  }

  if (hasAkPu3PFJetTree) {
    akPu3PFJetTree->SetName("akPu3PF");
    if (tree == 0) tree = akPu3PFJetTree; else tree->AddFriend(akPu3PFJetTree);
    setupJetTree(akPu3PFJetTree,akPu3PF);
  }

  if (hasAkPu4PFJetTree) {
    akPu4PFJetTree->SetName("akPu4PF");
    if (tree == 0) tree = akPu4PFJetTree; else tree->AddFriend(akPu4PFJetTree);
    setupJetTree(akPu4PFJetTree,akPu4PF);
  }

  if (hasAkPu5PFJetTree) {
    akPu5PFJetTree->SetName("akPu5PF");
    if (tree == 0) tree = akPu5PFJetTree; else tree->AddFriend(akPu5PFJetTree);
    setupJetTree(akPu5PFJetTree,akPu5PF);
  }
  if (hasAkPu6PFJetTree) {
    akPu6PFJetTree->SetName("akPu6PF");
    if (tree == 0) tree = akPu6PFJetTree; else tree->AddFriend(akPu6PFJetTree);
    setupJetTree(akPu6PFJetTree,akPu6PF);
  }
  if (hasAk1PFJetTree) {
    ak1PFJetTree->SetName("ak1PF");
    if (tree == 0) tree = ak1PFJetTree; else tree->AddFriend(ak1PFJetTree);
    setupJetTree(ak1PFJetTree,ak1PF);
  }

  if (hasAk2PFJetTree) {
    ak2PFJetTree->SetName("ak2PF");
    if (tree == 0) tree = ak2PFJetTree; else tree->AddFriend(ak2PFJetTree);
    setupJetTree(ak2PFJetTree,ak2PF);
  }

  if (hasAk3PFJetTree) {
    ak3PFJetTree->SetName("ak3PF");
    if (tree == 0) tree = ak3PFJetTree; else tree->AddFriend(ak3PFJetTree);
    setupJetTree(ak3PFJetTree,ak3PF);
  }

  if (hasAk4PFJetTree) {
    ak4PFJetTree->SetName("ak4PF");
    if (tree == 0) tree = ak4PFJetTree; else tree->AddFriend(ak4PFJetTree);
    setupJetTree(ak4PFJetTree,ak4PF);
  }

  if (hasAk5PFJetTree) {
    ak5PFJetTree->SetName("ak5PF");
    if (tree == 0) tree = ak5PFJetTree; else tree->AddFriend(ak5PFJetTree);
    setupJetTree(ak5PFJetTree,ak5PF);
  }
  if (hasAk6PFJetTree) {
    ak6PFJetTree->SetName("ak6PF");
    if (tree == 0) tree = ak6PFJetTree; else tree->AddFriend(ak6PFJetTree);
    setupJetTree(ak6PFJetTree,ak6PF);
  }
  ///////////////


  if (hasAkPu1CaloJetTree) {
    akPu1CaloJetTree->SetName("akPu1Calo");
    if (tree == 0) tree = akPu1CaloJetTree; else tree->AddFriend(akPu1CaloJetTree);
    setupJetTree(akPu1CaloJetTree,akPu1Calo);
  }

  if (hasAkPu2CaloJetTree) {
    akPu2CaloJetTree->SetName("akPu2Calo");
    if (tree == 0) tree = akPu2CaloJetTree; else tree->AddFriend(akPu2CaloJetTree);
    setupJetTree(akPu2CaloJetTree,akPu2Calo);
  }

  if (hasAkPu3CaloJetTree) {
    akPu3CaloJetTree->SetName("akPu3Calo");
    if (tree == 0) tree = akPu3CaloJetTree; else tree->AddFriend(akPu3CaloJetTree);
    setupJetTree(akPu3CaloJetTree,akPu3Calo);
  }

  if (hasAkPu4CaloJetTree) {
    akPu4CaloJetTree->SetName("akPu4Calo");
    if (tree == 0) tree = akPu4CaloJetTree; else tree->AddFriend(akPu4CaloJetTree);
    setupJetTree(akPu4CaloJetTree,akPu4Calo);
  }

  if (hasAkPu5CaloJetTree) {
    akPu5CaloJetTree->SetName("akPu5Calo");
    if (tree == 0) tree = akPu5CaloJetTree; else tree->AddFriend(akPu5CaloJetTree);
    setupJetTree(akPu5CaloJetTree,akPu5Calo);
  }
  if (hasAkPu6CaloJetTree) {
    akPu6CaloJetTree->SetName("akPu6Calo");
    if (tree == 0) tree = akPu6CaloJetTree; else tree->AddFriend(akPu6CaloJetTree);
    setupJetTree(akPu6CaloJetTree,akPu6Calo);
  }


  if (hasAk1CaloJetTree) {
    ak1CaloJetTree->SetName("ak1Calo");
    if (tree == 0) tree = ak1CaloJetTree; else tree->AddFriend(ak1CaloJetTree);
    setupJetTree(ak1CaloJetTree,ak1Calo);
  }

  if (hasAk2CaloJetTree) {
    ak2CaloJetTree->SetName("ak2Calo");
    if (tree == 0) tree = ak2CaloJetTree; else tree->AddFriend(ak2CaloJetTree);
    setupJetTree(ak2CaloJetTree,ak2Calo);
  }

  if (hasAk3CaloJetTree) {
    ak3CaloJetTree->SetName("ak3Calo");
    if (tree == 0) tree = ak3CaloJetTree; else tree->AddFriend(ak3CaloJetTree);
    setupJetTree(ak3CaloJetTree,ak3Calo);
  }

  if (hasAk4CaloJetTree) {
    ak4CaloJetTree->SetName("ak4Calo");
    if (tree == 0) tree = ak4CaloJetTree; else tree->AddFriend(ak4CaloJetTree);
    setupJetTree(ak4CaloJetTree,ak4Calo);
  }

  if (hasAk5CaloJetTree) {
    ak5CaloJetTree->SetName("ak5Calo");
    if (tree == 0) tree = ak5CaloJetTree; else tree->AddFriend(ak5CaloJetTree);
    setupJetTree(ak5CaloJetTree,ak5Calo);
  }
  if (hasAk6CaloJetTree) {
    ak6CaloJetTree->SetName("ak6Calo");
    if (tree == 0) tree = ak6CaloJetTree; else tree->AddFriend(ak6CaloJetTree);
    setupJetTree(ak6CaloJetTree,ak6Calo);
  }


  if (hasTrackTree) {
    trackTree->SetName("track");
    if (tree == 0) tree = trackTree; else tree->AddFriend(trackTree);
    trackTree->SetAlias("mergedGeneral","(trkAlgo<4||(highPurity))");
    trackTree->SetAlias("mergedSelected","(trkAlgo<4||(highPurity&&trkAlgo==4)))");
    setupTrackTree(trackTree,track);
  }
   
   

  /*
  if (hasMetTree) {
    metTree->SetName("met");
    if (tree == 0) tree = metTree; else tree->AddFriend(metTree);
    setupMetTree(metTree,met);
  }

  if (hasPixTrackTree) {
    pixtrackTree->SetName("pixtrack");
    if (tree == 0) tree = pixtrackTree; else tree->AddFriend(pixtrackTree);
    setupTrackTree(pixtrackTree,pixtrack);
  }

  if (hasPhotonTree) {
    photonTree->SetName("photon");
    photonTree->SetAlias("swiss","1-(eRight+eLeft+eTop+eBottom)/eMax");
    if (tree == 0) tree = photonTree; else tree->AddFriend(photonTree);
    setupPhotonTree(photonTree,photon);
  }

  if (hasPFTree) {
    pfTree->SetName("pf");
    if (tree == 0) tree = pfTree; else tree->AddFriend(pfTree);
    setupPFTree(pfTree,pf);
  }

  if (hasTowerTree) {
    towerTree->SetName("tower");
    if (tree == 0) tree = towerTree; else tree->AddFriend(towerTree);
    setupHitTree(towerTree,tower);
  }

  if (hasHbheTree) {
    hbheTree->SetName("hbhe");
    if (tree == 0) tree = hbheTree; else tree->AddFriend(hbheTree);
    setupHitTree(hbheTree,hbhe);
  }

  if (hasEbTree) {
    ebTree->SetName("eb");
    if (tree == 0) tree = ebTree; else tree->AddFriend(ebTree);
    setupHitTree(ebTree,eb);
  }
  
  if (hasGenpTree) {
    genpTree->SetName("genp");
    if (tree == 0) tree = genpTree; else tree->AddFriend(genpTree);
    setupGenpTree(genpTree,genp);
  }

  if (hasGenParticleTree) {
    genParticleTree->SetName("genParticle");
    if (tree == 0) tree = genParticleTree; else tree->AddFriend(genParticleTree);
    setupGenParticleTree(genParticleTree,genparticle);
  }
  */
  tree->SetMarkerStyle(20);

  // Print the status of thre forest
  PrintStatus();

}

HiForest::~HiForest()
{
  if (setupOutput) {
    for (unsigned int i=0; i<cloneForest.size(); i++)
    {
      cloneForest[i]->AutoSave();
    }
    outf->Close();
  }
}

void HiForest::GetEntry(int i)
{

   currentEvent = i;
  // get the entry of the available trees
  if (hasEvtTree)      evtTree   ->GetEntry(i);
  if (hasHltTree)      hltTree      ->GetEntry(i);
  if (hasSkimTree)     skimTree     ->GetEntry(i);
  if (hasTrackTree)    trackTree    ->GetEntry(i);

  if (hasIcPu5JetTree) icPu5jetTree ->GetEntry(i);

  if (hasAk1PFJetTree) ak1PFJetTree ->GetEntry(i);
  if (hasAk2PFJetTree) ak2PFJetTree ->GetEntry(i);
  if (hasAk3PFJetTree) ak3PFJetTree ->GetEntry(i);
  if (hasAk4PFJetTree) ak4PFJetTree ->GetEntry(i);
  if (hasAk5PFJetTree) ak5PFJetTree ->GetEntry(i);
  if (hasAk6PFJetTree) ak6PFJetTree ->GetEntry(i);


  if (hasAkPu1PFJetTree) akPu1PFJetTree ->GetEntry(i);
  if (hasAkPu2PFJetTree) akPu2PFJetTree ->GetEntry(i);
  if (hasAkPu3PFJetTree) akPu3PFJetTree ->GetEntry(i);
  if (hasAkPu4PFJetTree) akPu4PFJetTree ->GetEntry(i);
  if (hasAkPu5PFJetTree) akPu5PFJetTree ->GetEntry(i);
  if (hasAkPu6PFJetTree) akPu6PFJetTree ->GetEntry(i);



  if (hasAk1CaloJetTree) ak1CaloJetTree ->GetEntry(i);
  if (hasAk2CaloJetTree) ak2CaloJetTree ->GetEntry(i);
  if (hasAk3CaloJetTree) ak3CaloJetTree ->GetEntry(i);
  if (hasAk4CaloJetTree) ak4CaloJetTree ->GetEntry(i);
  if (hasAk5CaloJetTree) ak5CaloJetTree ->GetEntry(i);
  if (hasAk6CaloJetTree) ak6CaloJetTree ->GetEntry(i);


  if (hasAkPu1CaloJetTree) akPu1CaloJetTree ->GetEntry(i);
  if (hasAkPu2CaloJetTree) akPu2CaloJetTree ->GetEntry(i);
  if (hasAkPu3CaloJetTree) akPu3CaloJetTree ->GetEntry(i);
  if (hasAkPu4CaloJetTree) akPu4CaloJetTree ->GetEntry(i);
  if (hasAkPu5CaloJetTree) akPu5CaloJetTree ->GetEntry(i);
  if (hasAkPu6CaloJetTree) akPu6CaloJetTree ->GetEntry(i);


  /*
  if (hasPhotonTree)   photonTree   ->GetEntry(i);
  if (hasPFTree)       pfTree   ->GetEntry(i);
  if (hasMetTree)      metTree   ->GetEntry(i);
  if (hasPixTrackTree) pixtrackTree ->GetEntry(i);
  if (hasTowerTree)    towerTree    ->GetEntry(i);
  if (hasHbheTree)     hbheTree     ->GetEntry(i);
  if (hasEbTree)       ebTree     ->GetEntry(i);
  if (hasGenpTree)     genpTree   ->GetEntry(i);
  if (hasGenParticleTree) genParticleTree   ->GetEntry(i);
  */

  minJetPtForTrkCor = 40;
  leadingJetPtForTrkCor = -100;
  subleadingJetPtForTrkCor = -100;
  leadingJetEtaForTrkCor = -100;
  subleadingJetEtaForTrkCor = -100;
  leadingJetPhiForTrkCor = -100;
  subleadingJetPhiForTrkCor = -100;
}

int HiForest::GetEntries()
{
  // get the entries of the available trees
  return nEntries;
}

/*
void HiForest::InitTree()
{
   // Setup Track Corrections 	 
   if(doTrackCorrections){

      trackCorrFromParam = new TrackingParam();

      if (collisionMode==cPbPb) {
         trackCorrections.push_back(new TrackingCorrections("Forest2STAv14","Forest2_MergedGeneral_jetfine"));
      } else {
         trackCorrections.push_back(new TrackingCorrections("Forest2STApp","Forest2_MergedGeneral_jetfine"));
      }
//       trackCorrections.push_back(new TrackingCorrections("Forest2STAv12","Forest2_MergedGeneral_j1"));
//       trackCorrections.push_back(new TrackingCorrections("Forest2STAv12","Forest2_MergedGeneral_j2"));

      for(int i = 0; i < trackCorrections.size(); ++i){
         if (collisionMode==cPbPb) {
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/IterTrkCorrv14XSec_hy18dj80to100_akPu3PF_100_-1_-1000_genJetMode0.root",80);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/IterTrkCorrv14XSec_hy18dj100to170_akPu3PF_100_-1_-1000_genJetMode0.root",100);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/IterTrkCorrv14XSec_hy18dj170to200_akPu3PF_100_-1_-1000_genJetMode0.root",170);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/IterTrkCorrv14XSec_hy18dj200to250_akPu3PF_100_-1_-1000_genJetMode0.root",200);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/IterTrkCorrv14XSec_hy18dj250to300_akPu3PF_100_-1_-1000_genJetMode0.root",250);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/IterTrkCorrv14XSec_hy18dj300to9999_akPu3PF_100_-1_-1000_genJetMode0.root",300);
           // v12: frozen for preapproval
//            trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv12XSec/IterTrkCorrv12XSec_hy18dj100_akPu3PF_100_40_2749_genJetMode0.root",100);
//            trackCorrections[i]->AddNormFile("trkcorr/IterTrkCorrv12XSec/IterTrkCorrv12XSec_hy18dj100_akPu3PF_-1_-1_-1000_genJetMode0.root");
//            trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv12XSec/IterTrkCorrv12XSec_hy18dj200_akPu3PF_100_40_2749_genJetMode0.root",200);
//            trackCorrections[i]->AddNormFile("trkcorr/IterTrkCorrv12XSec/IterTrkCorrv12XSec_hy18dj200_akPu3PF_-1_-1_-1000_genJetMode0.root");
         } else {
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj80to120_akPu3PF_100_-1_-1000_genJetMode0.root",80);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj120to170_akPu3PF_100_-1_-1000_genJetMode0.root",120);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj170to200_akPu3PF_100_-1_-1000_genJetMode0.root",170);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj200to250_akPu3PF_100_-1_-1000_genJetMode0.root",200);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj250to9999_akPu3PF_100_-1_-1000_genJetMode0.root",250);
         }
         trackCorrections[i]->weightSamples_ = true;
         trackCorrections[i]->smoothLevel_ = 0;
         trackCorrections[i]->trkPhiMode_ = false;
         trackCorrections[i]->ppMode_ = (collisionMode==cPP);
         trackCorrections[i]->Init();
      }
      minJetPtForTrkCor = 40;
      initialized = 1;
      if (doTrackingSeparateLeadingSubleading&&trackCorrections.size()<3) {
         cout << "Fatal Error: Not enough correction tables to do separate leading/subleading jet tracking correction" << endl;
         exit(1);
      }
   }
}
*/

void HiForest::CheckTree(TTree *t,const char *title)
{
  int entries = t->GetEntries();
  if (nEntries==0) nEntries = entries;
  cout <<title<<": "<<entries<<" entries loaded.";
  if (entries != nEntries) {
    cout <<" Inconsistent number of entries!!"<<endl;
  } else {
    cout <<endl;
  }
}

/*
void HiForest::CheckArraySizes(){

  vector<int> trackOverflow(0);
  vector<int> objectOverflow(0);

  for(int ie = 0; ie < nEntries; ++ie){
    GetEntry(ie);
    if(track.nTrk > maxEntryTrack) trackOverflow.push_back(track.nTrk);
    if(akPu3PF.nref > maxEntry) objectOverflow.push_back(akPu3PF.nref);
    if(icPu5.nref > maxEntry) objectOverflow.push_back(icPu5.nref);
  }

  for(int i = 0; i < trackOverflow.size(); ++i){
    cout<<trackOverflow[i]<<endl;
  }

  if(trackOverflow.size() == 0) cout<<"Track sizes OK"<<endl;
  else cout<<"tracks crash"<<endl; // TODO : really crash

  for(int i = 0; i < objectOverflow.size(); ++i){
    cout<<objectOverflow[i]<<endl;
  }

  if(objectOverflow.size() == 0) cout<<"Object sizes OK"<<endl;
  else cout<<"objects crash"<<endl; // TODO : really crash
}
*/

void HiForest::PrintStatus()
{
  if (hasHltTree)      CheckTree(hltTree,      "HltTree");
  if (hasSkimTree)     CheckTree(skimTree,     "SkimTree");
  if (hasEvtTree)      CheckTree(evtTree,      "EvtTree");

  if (hasIcPu5JetTree) CheckTree(icPu5jetTree, "IcPu5jetTree");

  if (hasAk1PFJetTree) CheckTree(ak1PFJetTree, "Ak1PFJetTree");
  if (hasAk2PFJetTree) CheckTree(ak2PFJetTree, "Ak2PFJetTree");
  if (hasAk3PFJetTree) CheckTree(ak3PFJetTree, "Ak3PFJetTree");
  if (hasAk4PFJetTree) CheckTree(ak4PFJetTree, "Ak4PFJetTree");
  if (hasAk5PFJetTree) CheckTree(ak5PFJetTree, "Ak5PFJetTree");
  if (hasAk6PFJetTree) CheckTree(ak6PFJetTree, "Ak6PFJetTree");

  if (hasAkPu1PFJetTree) CheckTree(akPu1PFJetTree, "AkPu1PFJetTree");
  if (hasAkPu2PFJetTree) CheckTree(akPu2PFJetTree, "AkPu2PFJetTree");
  if (hasAkPu3PFJetTree) CheckTree(akPu3PFJetTree, "AkPu3PFJetTree");
  if (hasAkPu4PFJetTree) CheckTree(akPu4PFJetTree, "AkPu4PFJetTree");
  if (hasAkPu5PFJetTree) CheckTree(akPu5PFJetTree, "AkPu5PFJetTree");
  if (hasAkPu6PFJetTree) CheckTree(akPu6PFJetTree, "AkPu6PFJetTree");



  if (hasAk1CaloJetTree) CheckTree(ak1CaloJetTree, "Ak1CaloJetTree");
  if (hasAk2CaloJetTree) CheckTree(ak2CaloJetTree, "Ak2CaloJetTree");
  if (hasAk3CaloJetTree) CheckTree(ak3CaloJetTree, "Ak3CaloJetTree");
  if (hasAk4CaloJetTree) CheckTree(ak4CaloJetTree, "Ak4CaloJetTree");
  if (hasAk5CaloJetTree) CheckTree(ak5CaloJetTree, "Ak5CaloJetTree");
  if (hasAk6CaloJetTree) CheckTree(ak6CaloJetTree, "Ak6CaloJetTree");

  if (hasAkPu1CaloJetTree) CheckTree(akPu1CaloJetTree, "AkPu1CaloJetTree");
  if (hasAkPu2CaloJetTree) CheckTree(akPu2CaloJetTree, "AkPu2CaloJetTree");
  if (hasAkPu3CaloJetTree) CheckTree(akPu3CaloJetTree, "AkPu3CaloJetTree");
  if (hasAkPu4CaloJetTree) CheckTree(akPu4CaloJetTree, "AkPu4CaloJetTree");
  if (hasAkPu5CaloJetTree) CheckTree(akPu5CaloJetTree, "AkPu5CaloJetTree");
  if (hasAkPu6CaloJetTree) CheckTree(akPu6CaloJetTree, "AkPu6CaloJetTree");


  if (hasTrackTree)    CheckTree(trackTree,    "TrackTree");

  /*
  if (hasPixTrackTree) CheckTree(pixtrackTree, "PixTrackTree");
  if (hasPhotonTree)   CheckTree(photonTree,   "PhotonTree");
  if (hasPFTree)       CheckTree(pfTree,   "PFTree");
  if (hasMetTree)      CheckTree(metTree,   "MetTree");
  if (hasTowerTree)    CheckTree(towerTree,    "TowerTree");
  if (hasHbheTree)     CheckTree(hbheTree,     "HbheTree");
  if (hasEbTree)       CheckTree(ebTree,     "EbTree");
  if (hasGenpTree)      CheckTree(genpTree,   "GenpTree");
  if (hasGenParticleTree)      CheckTree(genParticleTree,   "GenParticleTree");
  */
}

void HiForest::SetOutputFile(const char *name)
{
  outf = new TFile(name,"recreate");
  if (hasHltTree)      AddCloneTree(hltTree,      "hltanalysis",        "HltTree");
  if (hasSkimTree)     AddCloneTree(skimTree,     "skimanalysis",       "HltTree");
  if (hasIcPu5JetTree) AddCloneTree(icPu5jetTree, "icPu5JetAnalyzer",   "t");


  if (hasAk1PFJetTree) AddCloneTree(ak1PFJetTree, "ak1PFJetAnalyzer", "t");
  if (hasAk2PFJetTree) AddCloneTree(ak2PFJetTree, "ak2PFJetAnalyzer", "t");
  if (hasAk3PFJetTree) AddCloneTree(ak3PFJetTree, "ak3PFJetAnalyzer", "t");
  if (hasAk4PFJetTree) AddCloneTree(ak4PFJetTree, "ak4PFJetAnalyzer", "t");
  if (hasAk5PFJetTree) AddCloneTree(ak5PFJetTree, "ak5PFJetAnalyzer", "t");
  if (hasAk6PFJetTree) AddCloneTree(ak6PFJetTree, "ak6PFJetAnalyzer", "t");

  if (hasAkPu1PFJetTree) AddCloneTree(akPu1PFJetTree, "akPu1PFJetAnalyzer", "t");
  if (hasAkPu2PFJetTree) AddCloneTree(akPu2PFJetTree, "akPu2PFJetAnalyzer", "t");
  if (hasAkPu3PFJetTree) AddCloneTree(akPu3PFJetTree, "akPu3PFJetAnalyzer", "t");
  if (hasAkPu4PFJetTree) AddCloneTree(akPu4PFJetTree, "akPu4PFJetAnalyzer", "t");
  if (hasAkPu5PFJetTree) AddCloneTree(akPu5PFJetTree, "akPu5PFJetAnalyzer", "t");
  if (hasAkPu6PFJetTree) AddCloneTree(akPu6PFJetTree, "akPu6PFJetAnalyzer", "t");


  if (hasAk1CaloJetTree) AddCloneTree(ak1CaloJetTree, "ak1CaloJetAnalyzer", "t");
  if (hasAk2CaloJetTree) AddCloneTree(ak2CaloJetTree, "ak2CaloJetAnalyzer", "t");
  if (hasAk3CaloJetTree) AddCloneTree(ak3CaloJetTree, "ak3CaloJetAnalyzer", "t");
  if (hasAk4CaloJetTree) AddCloneTree(ak4CaloJetTree, "ak4CaloJetAnalyzer", "t");
  if (hasAk5CaloJetTree) AddCloneTree(ak5CaloJetTree, "ak5CaloJetAnalyzer", "t");
  if (hasAk6CaloJetTree) AddCloneTree(ak6CaloJetTree, "ak6CaloJetAnalyzer", "t");

  if (hasAkPu1CaloJetTree) AddCloneTree(akPu1CaloJetTree, "akPu1CaloJetAnalyzer", "t");
  if (hasAkPu2CaloJetTree) AddCloneTree(akPu2CaloJetTree, "akPu2CaloJetAnalyzer", "t");
  if (hasAkPu3CaloJetTree) AddCloneTree(akPu3CaloJetTree, "akPu3CaloJetAnalyzer", "t");
  if (hasAkPu4CaloJetTree) AddCloneTree(akPu4CaloJetTree, "akPu4CaloJetAnalyzer", "t");
  if (hasAkPu5CaloJetTree) AddCloneTree(akPu5CaloJetTree, "akPu5CaloJetAnalyzer", "t");
  if (hasAkPu6CaloJetTree) AddCloneTree(akPu6CaloJetTree, "akPu6CaloJetAnalyzer", "t");


  if (hasTrackTree)    AddCloneTree(trackTree,    "ppTrack",           "trackTree");

  /*
  if (hasPixTrackTree) AddCloneTree(pixtrackTree, "anaPixTrack",        "trackTree");
  if (hasPhotonTree)   AddCloneTree(photonTree,   "multiPhotonAnalyzer",            "photon");
  if (hasPFTree)   AddCloneTree(pfTree,   "pfcandAnalyzer",            "pfTree");
  if (hasEvtTree)      AddCloneTree(evtTree,      "hiEvtAnalyzer",            "HiTree");
  if (hasMetTree)      AddCloneTree(metTree,      "anaMET",            "metTree");
  if (hasTowerTree)    AddCloneTree(towerTree,    "rechitanalyzer",              "tower");
  if (hasHbheTree)     AddCloneTree(hbheTree,     "rechitanalyzer",               "hbhe");
  if (hasEbTree)       AddCloneTree(ebTree,       "rechitanalyzer",               "eb");
  */

  setupOutput = true;
}

void HiForest::AddCloneTree(TTree* t, const char *dirName, const char *treeName)
{
  // Make directory
  outf->cd();
  outf->mkdir(dirName);
  outf->cd(dirName);

  // Add a clone tree to the clone forest
  TTree *tClone = t->CloneTree(0);
  tClone->SetMaxTreeSize(4000000000);
  tClone->SetName(treeName);
  
  cloneForest.push_back(tClone);
}

void HiForest::FillOutput()
{
  if (setupOutput) {
     for (unsigned int i=0; i<cloneForest.size(); i++)
     {
       cloneForest[i]->Fill();
     } 
  } else {
       cout <<"ERROR: Specify an output file by hiForest.SetOutputFile(filename)!"<<endl;
  }
}

// ====================== Event Utilities ========================

bool HiForest::selectEvent(){
  /*
   bool select = skim.phbheReflagNewTimeEnv 
      && 
      skim.phcalTimingFilter 
      && 
      skim.pHBHENoiseFilter 
      && 
      skim.phiEcalRecHitSpikeFilter;
  */

   bool select = skim.pHBHENoiseFilter || mc;
   if(collisionMode==cPbPb){
      select = select && skim.pcollisionEventSelection;
   }else if(collisionMode==cPPb){
      select = select &&
	 hlt.HLT_PAZeroBiasPixel_SingleTrack_v1 &&
	 //	 skim.pPAcollisionEventSelection &&
	 skim.phfPosFilter1&&
	 skim.phfNegFilter1&&
	skim.phltPixelClusterShapeFilter&&
	skim.pprimaryvertexFilter;
   }else if(collisionMode==cPP){
     //select = select && skim.phfCoincFilter && skim.ppurityFractionFilter;
   }
   return select;
}

TCut HiForest::eventSelection(){
  //   TCut select("skim.phbheReflagNewTimeEnv && skim.phcalTimingFilter && skim.pHBHENoiseFilter && skim.phiEcalRecHitSpikeFilter");
  TCut select("skim.pHBHENoiseFilter");
  if(collisionMode==cPbPb){
    select = select && "skim.pcollisionEventSelection";
  }else{
    //      select = select && "skim.phfCoincFilter && skim.ppurityFractionFilter";
  }
   return select;
}

/*
void HiForest::GetEnergyScaleTable(char *fileNameTable) {
   TFile* f = new TFile(fileNameTable);
   fEnergyScale[0][1] = (TF1*)f->Get("fit_hscale_r9gt94_1");
   fEnergyScale[0][2] = (TF1*)f->Get("fit_hscale_r9gt94_2");
   fEnergyScale[0][3] = (TF1*)f->Get("fit_hscale_r9gt94_3");
   fEnergyScale[1][1] = (TF1*)f->Get("fit_hscale_r9lt94_1");
   fEnergyScale[1][2] = (TF1*)f->Get("fit_hscale_r9lt94_2");
   fEnergyScale[1][3] = (TF1*)f->Get("fit_hscale_r9lt94_3");
}
*/

// ====================== Track Utilities ========================
//#include "TrackUtilities.C"

// ======================= Jet Utilities =========================
//#include "JetUtilities.C"

// ====================== Photon Utilities ========================
//#include "PhotonUtilities.C"


#endif
