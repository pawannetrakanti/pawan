#ifndef HIFOREST_H
#define HIFOREST_H

#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>

#include "SetupEvtTree.h"
#include "SetupJetTree.h"


#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TF1.h>
#include <TCut.h>

using namespace std;

typedef std::map<const char  *,int > JetMap;
typedef std::pair<const char *,int > JetPair;

// ==========================================================
// Main class which can be used to read the hiForest trees
//
// Author: Yen-Jie Lee
//
// ==========================================================

const int ndir=25;
const char *calgo[ndir]= {
  "icPu5",
  "ak1PF","ak2PF","ak3PF","ak4PF","ak5PF","ak6PF",
  "ak1Calo","ak2Calo","ak3Calo","ak4Calo","ak5Calo","ak6Calo",
  "akPu1PF","akPu2PF","akPu3PF","akPu4PF","akPu5PF","akPu6PF",
  "akPu1Calo","akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo","akPu6Calo"
};


class HiForest : public TNamed
{

  using TObject::Draw;
  using TObject::GetName;

  public: 
  HiForest(const char *file, const char *name="forest", bool ispp = 0, bool ismc = 0);
  ~HiForest();

  // Utility functions
  int  GetEntry(int i);
  int  GetEntries();  						// Get the number of entries 
  void CheckTree(TTree *t,const char *title);			// Check the status of a tree
  void PrintStatus();						// Print the status of the hiForest

  Long64_t Draw(const char* varexp, const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0){
    return tree->Draw(varexp,selection,option,nentries,firstentry);
  }
  const char *GetName(){
    return fname;
  }

  //! Jet Algorithms
  int GetNAlgo(); 
  const char *GetAlgoName(int i);
  Jets *GetJet(int i);
  Jets *GetJetByAlgo(const char *algo);
  void FillMap();
  void EmptyMap();
  void SelectJetAlgo (const char * algo[],const int len);
  void SelectBranches(const char *tname,const char * blst[],const int len);
  
  // TFile
  TFile *inf; 					// Input file 

  // Trees
  TTree *evtTree;                               // Event Tree
  TTree *JetTree[ndir];                         // Jet Trees

  TTree *tree;					// Pointer to the available tree, all trees in the forest are friended to each other

  // Branches

  // event tree
  Evts evt;
  
  //! Jet Trees
  Jets *mJets;
  JetMap JetContainer;
  JetMap::iterator it;

  // Booleans
  bool hasEvtTree;
  bool hasJetTree[ndir];

  bool pp;
  bool mc;

  int nEntries;
  int currentEvent;

 private:
  const char *fname;

  ClassDef(HiForest,7)
};
#endif

HiForest::HiForest(const char *infName, const char* name, bool ispp, bool ismc):
   tree(0),
   pp(ispp),
   mc(ismc),
   nEntries(0),
   currentEvent(0),
   fname(name)
{
  tree = new TTree("tree",name);

  //! Initialize the jet collection
  mJets = new Jets[ndir];
  
  // Input file
  inf = TFile::Open(infName);

  // Load trees. Hard coded for the moment
  evtTree      = (TTree*) inf->Get("hiEvtAnalyzer/HiTree");

  // Check the validity of the trees.
  hasEvtTree       = (evtTree   != 0);

  // Setup branches. See also Setup*.h
  if (hasEvtTree) {
    evtTree->SetName("event");
    if (tree == 0) tree = evtTree; else tree->AddFriend(evtTree);
    setupEvtTree(evtTree,evt);
    //std::cout<<"SetupEvt Tree : " <<std::endl;
  }

  //! Jet trees
  for(int idir=0;idir<ndir;idir++){
    JetTree[idir] = (TTree*)inf->Get(Form("%sJetAnalyzer/t",calgo[idir]));
    hasJetTree[idir] = (JetTree[idir]!=0);
    if(hasJetTree[idir]){
      JetTree[idir]->SetName(Form("%s",calgo[idir]));
      if (tree == 0) tree = JetTree[idir]; else tree->AddFriend(JetTree[idir]);
      setupJetTree(JetTree[idir],mJets[idir]);      
      //std::cout<<Form("SetupJet Tree : %s",calgo[idir])<<"\t # of entries : "<<JetTree[idir]->GetEntries()<<std::endl;
    }
  }

  tree->SetMarkerStyle(20);
  std::cout<<"All Trees set up " <<std::endl;

  // Print the status of thre forest
  PrintStatus();

  //! Fill map for association of jet algo and the trees
  FillMap();
}

HiForest::~HiForest()
{
  delete [] mJets;
  EmptyMap();
}

int HiForest::GetEntry(int i)
{
  int nb=0;
  currentEvent = i;
  // get the entry of the available trees
  if (hasEvtTree)      nb+=evtTree      ->GetEntry(i);
  
  for(int idir=0;idir<ndir;idir++){
    if (hasJetTree[idir]) nb+= JetTree[idir]->GetEntry(i);
  }
  return nb;
}

int HiForest::GetEntries()
{
  // get the entries of the available trees
  return nEntries;
}

void HiForest::CheckTree(TTree *t,const char *title)
{
  int entries = t->GetEntries();
  if (nEntries==0) nEntries = entries;
  std::cout <<title<<"\t : "<<entries<<" entries loaded.";
  if (entries != nEntries) {
    std::cout <<" Inconsistent number of entries!!"<<std::endl;
  } else {
    std::cout <<std::endl;
  }
}
void HiForest::PrintStatus()
{
  if (hasEvtTree)      CheckTree(evtTree,      "EvtTree");
  for(int i=0;i<ndir;i++){
    if(hasJetTree[i])CheckTree(JetTree[i],Form("%sJetTree",calgo[i]));
    //std::cout<<"Jet Tree : "<<calgo[i]<<"\t "<<algo[j]<<std::endl;
  }
}
// ====================== Jet related functions ==================
int HiForest::GetNAlgo(){
  //return ndir;
  return JetContainer.size();
}
const char *HiForest::GetAlgoName(int i){
  return calgo[i];
}
Jets *HiForest::GetJet(int i){
  return &(mJets[i]);
}

void HiForest::FillMap(){
  std::cout <<std::endl;
  std::cout <<std::endl;
  for(int idir=0;idir<ndir;idir++){
    if(!hasJetTree[idir])continue;
    JetContainer.insert(JetPair(calgo[idir],idir));
  }
  std::cout<<" FillMap() : Size of the Jet Container : " <<JetContainer.size()<<std::endl;
}
Jets *HiForest::GetJetByAlgo(const char *algo){
  int iFound=-1;
  for (it = JetContainer.begin(); it != JetContainer.end(); ++it) {
    if(strcmp((*it).first,algo)==0){
      iFound=(*it).second;
      break;
    }
  }
  return &(mJets[iFound]);
}
void HiForest::SelectJetAlgo(const char *algo[],const int len){

  for(int i=0;i<ndir;i++){
    hasJetTree[i]=0;
    for(int j=0;j<len;j++){
      if(strcmp(algo[j],calgo[i])==0){
	hasJetTree[i]=1;
	CheckTree(JetTree[i],Form("%sJetTree",calgo[i]));
	//std::cout<<"Jet Tree : "<<calgo[i]<<"\t "<<algo[j]<<std::endl;
      }
    }
  }
}
void HiForest::SelectBranches(const char *tname,const char *blst[],const int len){
  
  if(strcmp(tname,"JetTree")==0){
    for(int i=0;i<ndir;i++){
      //std::cout<<"Jet Tree status : "<<calgo[i]<<"\t"<<hasJetTree[i]<<std::endl;
      if(!hasJetTree[i])continue;
      JetTree[i]->SetBranchStatus("*",0);
      //std::cout<<calgo[i]<<"\t # of entries : "<<JetTree[i]->GetEntries()<<std::endl;
      for(int j=0;j<len;j++){  
	JetTree[i]->SetBranchStatus(blst[j],1);
      //std::cout<<"\t branches : "<<blst[j]<<"\t status : "<<JetTree[i]->GetBranchStatus(blst[j])<<std::endl;
      }  
    }
  }

  if(strcmp(tname,"evtTree")==0){
    if (hasEvtTree){
      evtTree->SetBranchStatus("*",0);
      for(int j=0;j<len;j++){        
	evtTree->SetBranchStatus(blst[j],1);
      }
    }else std::cout<<"eveTree not available : "<<std::endl;
  }
}
void HiForest::EmptyMap(){
  if (JetContainer.empty()) return;
  JetContainer.clear();
}
// ====================== Event Utilities ========================
