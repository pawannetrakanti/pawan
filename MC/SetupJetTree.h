#include "commonSetup.h"
#include <iostream>
#include <string>
#include <TTree.h>

using namespace std;

class Jets
{
  public:

  Jets(){};
  ~Jets(){};

  // Declaration of leaf types
  int           evt;
  float         b;
  int           nref;
  float         rawpt[MaxNumber];   //[nref]
  float         jtpt[MaxNumber];   //[nref]
  float         jteta[MaxNumber];   //[nref]
  float         jty[MaxNumber];   //[nref]
  float         jtphi[MaxNumber];   //[nref]
  float         jtpu[MaxNumber];   //[nref]

  float         discr_fr01[MaxNumber];   //[nref]
  float         trackMax[MaxNumber];   //[nref]
  float         trackSum[MaxNumber];   //[nref]
  int           trackN[MaxNumber];   //[nref]
  float         trackHardSum[MaxNumber];   //[nref]
  int           trackHardN[MaxNumber];   //[nref]
  float         chargedMax[MaxNumber];   //[nref]
  float         chargedSum[MaxNumber];   //[nref]
  int           chargedN[MaxNumber];   //[nref]
  float         chargedHardSum[MaxNumber];   //[nref]
  int           chargedHardN[MaxNumber];   //[nref]
  float         photonMax[MaxNumber];   //[nref]
  float         photonSum[MaxNumber];   //[nref]
  int           photonN[MaxNumber];   //[nref]
  float         photonHardSum[MaxNumber];   //[nref]
  int           photonHardN[MaxNumber];   //[nref]
  float         neutralMax[MaxNumber];   //[nref]
  float         neutralSum[MaxNumber];   //[nref]
  int           neutralN[MaxNumber];   //[nref]
  float         eMax[MaxNumber];   //[nref]
  float         eSum[MaxNumber];   //[nref]
  int           eN[MaxNumber];   //[nref]
  float         muMax[MaxNumber];   //[nref]
  float         muSum[MaxNumber];   //[nref]
  int           muN[MaxNumber];   //[nref]
  float         matchedPt[MaxNumber];   //[nref]
  float         matchedR[MaxNumber];   //[nref]
  int           beamId1;
  int           beamId2;

  float         pthat;
  float         refpt[MaxNumber];   //[nref]
  float         refeta[MaxNumber];   //[nref]
  float         refy[MaxNumber];   //[nref]
  float         refphi[MaxNumber];   //[nref]
  float         refdphijt[MaxNumber];   //[nref]
  float         refdrjt[MaxNumber];   //[nref]
  float         refparton_pt[MaxNumber];   //[nref]
  int           refparton_flavor[MaxNumber];   //[nref]
  int           refparton_flavorForB[MaxNumber];   //[nref]

  int           ngen;
  int           genmatchindex[MaxNumber];   //[ngen]
  float         genpt[MaxNumber];   //[ngen]
  float         geneta[MaxNumber];   //[ngen]
  float         geny[MaxNumber];   //[ngen]
  float         genphi[MaxNumber];   //[ngen]
  float         gendphijt[MaxNumber];   //[ngen]
  float         gendrjt[MaxNumber];   //[ngen]
  int           gensubid[MaxNumber];   //[ngen]

  // List of branches
  TBranch        *b_evt;   //!
  TBranch        *b_b;   //!
  TBranch        *b_nref;   //!
  TBranch        *b_rawpt;   //!
  TBranch        *b_jtpt;   //!
  TBranch        *b_jteta;   //!
  TBranch        *b_jty;   //!
  TBranch        *b_jtphi;   //!
  TBranch        *b_jtpu;   //!

  TBranch        *b_discr_fr01;   //!
  TBranch        *b_trackMax;   //!
  TBranch        *b_trackSum;   //!
  TBranch        *b_trackN;   //!
  TBranch        *b_trackHardSum;   //!
  TBranch        *b_trackHardN;   //!
  TBranch        *b_chargedMax;   //!
  TBranch        *b_chargedSum;   //!
  TBranch        *b_chargedN;   //!
  TBranch        *b_chargedHardSum;   //!
  TBranch        *b_chargedHardN;   //!
  TBranch        *b_photonMax;   //!
  TBranch        *b_photonSum;   //!
  TBranch        *b_photonN;   //!
  TBranch        *b_photonHardSum;   //!
  TBranch        *b_photonHardN;   //!
  TBranch        *b_neutralMax;   //!
  TBranch        *b_neutralSum;   //!
  TBranch        *b_neutralN;   //!
  TBranch        *b_eMax;   //!
  TBranch        *b_eSum;   //!
  TBranch        *b_eN;   //!
  TBranch        *b_muMax;   //!
  TBranch        *b_muSum;   //!
  TBranch        *b_muN;   //!
  TBranch        *b_matchedPt;   //!
  TBranch        *b_matchedR;   //!
  TBranch        *b_beamId1;   //!
  TBranch        *b_beamId2;   //!

  TBranch        *b_pthat;   //!
  TBranch        *b_refpt;   //!
  TBranch        *b_refeta;   //!
  TBranch        *b_refy;   //!
  TBranch        *b_refphi;   //!
  TBranch        *b_refdphijt;   //!
  TBranch        *b_refdrjt;   //!
  TBranch        *b_refparton_pt;   //!
  TBranch        *b_refparton_flavor;   //!
  TBranch        *b_refparton_flavorForB;   //!

  TBranch        *b_ngen;   //!
  TBranch        *b_genmatchindex;   //!
  TBranch        *b_genpt;   //!
  TBranch        *b_geneta;   //!
  TBranch        *b_geny;   //!
  TBranch        *b_genphi;   //!
  TBranch        *b_gendphijt;   //!
  TBranch        *b_gendrjt;   //!
  TBranch        *b_gensubid;   //!

};

void setupJetTree(TTree *t, Jets &tjets,bool doCheck=1)
{

  std::string name = t->GetName();

  t->SetBranchAddress("evt", &tjets.evt, &tjets.b_evt);
  t->SetBranchAddress("b", &tjets.b, &tjets.b_b);
  t->SetBranchAddress("nref", &tjets.nref, &tjets.b_nref);
  t->SetBranchAddress("rawpt", tjets.rawpt, &tjets.b_rawpt);
  t->SetBranchAddress("jtpt", tjets.jtpt, &tjets.b_jtpt);
  t->SetBranchAddress("jteta", tjets.jteta, &tjets.b_jteta);
  t->SetBranchAddress("jty", tjets.jty, &tjets.b_jty);
  t->SetBranchAddress("jtphi", tjets.jtphi, &tjets.b_jtphi);
  t->SetBranchAddress("jtpu", tjets.jtpu, &tjets.b_jtpu);

  t->SetBranchAddress("discr_fr01", tjets.discr_fr01, &tjets.b_discr_fr01);
  t->SetBranchAddress("trackMax", tjets.trackMax, &tjets.b_trackMax);
  t->SetBranchAddress("trackSum", tjets.trackSum, &tjets.b_trackSum);
  t->SetBranchAddress("trackN", tjets.trackN, &tjets.b_trackN);
  t->SetBranchAddress("trackHardSum", tjets.trackHardSum, &tjets.b_trackHardSum);
  t->SetBranchAddress("trackHardN", tjets.trackHardN, &tjets.b_trackHardN);
  t->SetBranchAddress("chargedMax", tjets.chargedMax, &tjets.b_chargedMax);
  t->SetBranchAddress("chargedSum", tjets.chargedSum, &tjets.b_chargedSum);
  t->SetBranchAddress("chargedN", tjets.chargedN, &tjets.b_chargedN);
  t->SetBranchAddress("chargedHardSum", tjets.chargedHardSum, &tjets.b_chargedHardSum);
  t->SetBranchAddress("chargedHardN", tjets.chargedHardN, &tjets.b_chargedHardN);
  t->SetBranchAddress("photonMax", tjets.photonMax, &tjets.b_photonMax);
  t->SetBranchAddress("photonSum", tjets.photonSum, &tjets.b_photonSum);
  t->SetBranchAddress("photonN", tjets.photonN, &tjets.b_photonN);
  t->SetBranchAddress("photonHardSum", tjets.photonHardSum, &tjets.b_photonHardSum);
  t->SetBranchAddress("photonHardN", tjets.photonHardN, &tjets.b_photonHardN);
  t->SetBranchAddress("neutralMax", tjets.neutralMax, &tjets.b_neutralMax);
  t->SetBranchAddress("neutralSum", tjets.neutralSum, &tjets.b_neutralSum);
  t->SetBranchAddress("neutralN", tjets.neutralN, &tjets.b_neutralN);
  t->SetBranchAddress("eMax", tjets.eMax, &tjets.b_eMax);
  t->SetBranchAddress("eSum", tjets.eSum, &tjets.b_eSum);
  t->SetBranchAddress("eN", tjets.eN, &tjets.b_eN);
  t->SetBranchAddress("muMax", tjets.muMax, &tjets.b_muMax);
  t->SetBranchAddress("muSum", tjets.muSum, &tjets.b_muSum);
  t->SetBranchAddress("muN", tjets.muN, &tjets.b_muN);
  t->SetBranchAddress("matchedPt", tjets.matchedPt, &tjets.b_matchedPt);
  t->SetBranchAddress("matchedR", tjets.matchedR, &tjets.b_matchedR);
  t->SetBranchAddress("beamId1", &tjets.beamId1, &tjets.b_beamId1);
  t->SetBranchAddress("beamId2", &tjets.beamId2, &tjets.b_beamId2);

  t->SetBranchAddress("pthat", &tjets.pthat, &tjets.b_pthat);
  t->SetBranchAddress("refpt", tjets.refpt, &tjets.b_refpt);
  t->SetBranchAddress("refeta", tjets.refeta, &tjets.b_refeta);
  t->SetBranchAddress("refy", tjets.refy, &tjets.b_refy);
  t->SetBranchAddress("refphi", tjets.refphi, &tjets.b_refphi);
  t->SetBranchAddress("refdphijt", tjets.refdphijt, &tjets.b_refdphijt);
  t->SetBranchAddress("refdrjt", tjets.refdrjt, &tjets.b_refdrjt);
  t->SetBranchAddress("refparton_pt", tjets.refparton_pt, &tjets.b_refparton_pt);
  t->SetBranchAddress("refparton_flavor", tjets.refparton_flavor, &tjets.b_refparton_flavor);
  t->SetBranchAddress("refparton_flavorForB", tjets.refparton_flavorForB, &tjets.b_refparton_flavorForB);


  if(name.find("Pu") != string::npos && name.find("PF") != string::npos){
    t->SetBranchAddress("ngen", &tjets.ngen, &tjets.b_ngen);
    t->SetBranchAddress("genmatchindex", tjets.genmatchindex, &tjets.b_genmatchindex);
    t->SetBranchAddress("genpt", tjets.genpt, &tjets.b_genpt);
    t->SetBranchAddress("geneta", tjets.geneta, &tjets.b_geneta);
    t->SetBranchAddress("geny", tjets.geny, &tjets.b_geny);
    t->SetBranchAddress("genphi", tjets.genphi, &tjets.b_genphi);
    t->SetBranchAddress("gendphijt", tjets.gendphijt, &tjets.b_gendphijt);
    t->SetBranchAddress("gendrjt", tjets.gendrjt, &tjets.b_gendrjt);
    t->SetBranchAddress("gensubid", tjets.gensubid, &tjets.b_gensubid);
  }
  
  if (doCheck) {
    if (t->GetMaximum("nref")>MaxNumber || t->GetMaximum("ngen")>MaxNumber) std::cout <<"FATAL ERROR: Arrary size of nref || ngen too small!!!  "<<t->GetMaximum("nref")<<std::endl;
  }

}
