

#include "/afs/cern.ch/user/p/pawan/scratch0/CMSSW_5_3_3/src/pPb/HiForest/V3/hiForest.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "TMath.h"

#include "TCut.h"
#include <string>

using namespace std;

void check(const char* infname = "root://eoscms//eos/cms/store/caf/user/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_JSonPPb_forestv77.root",
	   const char* outname = "ntuple_doga.root",
	   bool MC = 0,
	   bool PbPb = 0,
	   double jetEtaMax = 3.0){

  cout<<"Begin"<<endl;
  
  int Nevents = 50000;
  Nevents = -1;
  
  bool usePF = 1;

  double leadPtMin = 100;
  double jetTrackMin = 4;

  //double ptSubLeadMin = 30.;
  //double ptLeadMin    = 120.;


  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();


  
  HiForest * t;
  if(PbPb){
    t = new HiForest(infname,"",cPbPb,MC);
  }else{
    t = new HiForest(infname,"",cPPb,MC);
  }

  t->hasPhotonTree *= 0;
  t->hasMetTree *= 0;
  t->hasPFTree *= 0;
  t->hasPFTree *= 0;
  
  
  t->hasAk2JetTree *= 0;
  t->hasAk3JetTree *= 0;
  t->hasAk4JetTree *= 0;
  t->hasAk5JetTree *= 0;
  
  t->hasAkPu2JetTree *= 0;
  t->hasAkPu3JetTree *= 1;
  t->hasAkPu4JetTree *= 0;
  t->hasAkPu5JetTree *= 0;
  
  t->hasAk2CaloJetTree *= 0;
  t->hasAk3CaloJetTree *= 0;
  t->hasAk4CaloJetTree *= 0;
  t->hasAk5CaloJetTree *= 0;
  
  t->hasAkPu2CaloJetTree *= 0;
  t->hasAkPu3CaloJetTree *= 0;
  t->hasAkPu4CaloJetTree *= 0;
  t->hasAkPu5CaloJetTree *= 0;
  
  t->hasTrackTree *= 0;
  t->hasPixTrackTree *= 0;
  t->hasTowerTree *= 0;
  t->hasHbheTree *= 0;
  t->hasEbTree *= 0;
  t->hasGenpTree *= 0;
  t->hasGenParticleTree *= MC;

  t->InitTree();
  
  cout<<"a"<<endl;
  
  
  if(Nevents > 0){
    t->nEntries = Nevents;
  }else{
    Nevents = t->nEntries;
  }

  TFile* outf = new TFile(outname,"recreate");
  string jetVars = "";
  jetVars += "evt:run:hfp:hfm:zdc:vz:jet80:ntrk:npix:noise:jtpt1:raw1:jteta1:jtphi1:trkMax1:jtpt2:raw2:jteta2:jtphi2:trkMax2:jtpt3:raw3:jteta3:jtphi3:trkMax3";
  cout<<"Filling jet variables   : "<<jetVars.data()<<endl;
  TNtuple *ntjet=0;
  ntjet = new TNtuple("ntjet","",jetVars.data());
  
  outf->cd();
  
  vector<JetIndex> vecs;
  vecs.reserve(maxEntry);
  
  Jets *jets1 = 0;
  
  if(usePF){
    jets1 = &(t->akPu3PF);
  }else{
    jets1 = &(t->akPu3Calo);
  }
  cout<<"a"<<endl;

  for(int iev = 0; iev < Nevents; ++iev){
    
    if(iev%1000==0){ 
      cout<<"Processing entry : "<<iev<<endl;
    }
    t->GetEntry(iev);
    int evt = t->hlt.Event;
    int run = t->hlt.Run;
    if(run>211256) continue;
    
    if(!(t->skim.pPAcollisionEventSelectionPA && t->skim.pHBHENoiseFilter && t->hlt.HLT_PAJet80_NoJetID_v1 && t->skim.pVertexFilterCutGplus && fabs(t->evt.vz)<15.)) continue;
    bool jet80 = t->hlt.HLT_PAJet80_NoJetID_v1;

    // add other selection
    double noise=-1;

    double pt1 = -9, pt2 = -9, pt3 = -9,
      raw1 = -9, raw2 = -9, raw3 = -9,
      eta1 = -9,eta2 = -9, eta3  = -9,
      phi1 = -9,phi2 = -9, phi3  = -9,
      trkMax1 = -9,trkMax2=-9,trkMax3=-9;
    
     double hfp  = t->evt.hiHFplusEta4;
     double hfm  = t->evt.hiHFminusEta4;
     double zdc  = t->evt.hiZDCminus;
     double vz   = t->evt.vz;
     //double vz   = t->track.zVtx[t->track.maxVtx];
     double ntrk = t->evt.hiNtracks;
     double npix = t->evt.hiNpix;

     vecs.clear();

     for(int j = 0; j < jets1->nref; ++j){
       if(jets1->rawpt[j] < 15) continue;
       
       if( fabs(jets1->jteta[j]) > jetEtaMax ) continue;
       if(jets1->trackMax[j]<jetTrackMin) continue;

       JetIndex entry;
       entry.pt = jets1->jtpt[j];
       entry.index = j;
       vecs.push_back(entry);
     }

     sort(vecs.begin(),vecs.end(),comparePt);

     int jtLead = -1, jtSubLead = -1, jtThird = -1;

     if(vecs.size() > 0) jtLead    = vecs[0].index;
     if(vecs.size() > 1) jtSubLead = vecs[1].index;
     if(vecs.size() > 2) jtThird   = vecs[2].index;

     if((vecs.size() < 1 || jets1->jtpt[jtLead] < leadPtMin)) continue;
          


     if(jtLead > -1){
       pt1     = jets1->jtpt[jtLead];
       eta1    = jets1->jteta[jtLead];
       phi1    = jets1->jtphi[jtLead];
       raw1    = jets1->rawpt[jtLead];
       trkMax1 = jets1->trackMax[jtLead];
     }

     if(jtSubLead > -1){
       pt2   =jets1->jtpt[jtSubLead];
       eta2  =jets1->jteta[jtSubLead];
       phi2  =jets1->jtphi[jtSubLead];
       raw2  =jets1->rawpt[jtSubLead];
       trkMax2 = jets1->trackMax[jtSubLead];
     }

     if(jtThird > -1){
       pt3     = jets1->jtpt[jtThird];
       eta3    = jets1->jteta[jtThird];
       phi3    = jets1->jtphi[jtThird];
       raw3    = jets1->rawpt[jtThird];
       trkMax3 = jets1->trackMax[jtThird];
     }

     
     float jentry[] = {evt,run,hfp,hfm,zdc,vz,jet80,ntrk,npix,noise,
		       pt1,raw1,eta1,phi1,trkMax1,
		       pt2,raw2,eta2,phi2,trkMax2,
		       pt3,raw3,eta3,phi3,trkMax3};
     
     ntjet->Fill(jentry);
  }
  
  outf->Write();
  outf->Close();

  cout<<"Congrats!!!"<<endl;
}




