#ifndef analyzeTrackingCorrection_h
#define analyzeTrackingCorrection_h
#include "TChain.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include <iostream>
#include <cmath>
#include "analyzeDiJetMPT.h"

typedef struct
{
   Int_t ids;
   Float_t etas;
   Float_t phis;
   Float_t pts;
   Int_t hits;
   Int_t status;
   Int_t acc;
   Int_t nrec;
   Float_t ptr;
   Float_t dz;
   Float_t d0;
   Float_t pterr;
   Float_t d0err;
   Float_t dzerr;
   Int_t hitr;
   Int_t algo;
   Float_t jet;
   Float_t jeta;
   Float_t jdr;
} SimTrack_t;

typedef struct
{
   Int_t charge;
   Float_t etar;
   Float_t ptr;
   Float_t phir;
   Float_t dz;
   Float_t d0;
   Float_t pterr;
   Float_t d0err;
   Float_t dzerr;
   Int_t hitr; 
   Int_t algo;
   Int_t nsim;
   Int_t status;
   Int_t ids;
   Int_t parids;
   Float_t etas;
   Float_t pts;
   Float_t jet;
   Float_t jeta;
   Float_t jdr;
} RecTrack_t;

class TrkCorrHisAna
{
  public:
    vector<Double_t> ptBins;
    vector<Double_t> etaBins;
    vector<Double_t> phiBins;
    vector<Double_t> jetBins;
    vector<Double_t> zBins;
    vector<Int_t> centBins;

    TString name_;
    bool trkPhiMode_;
    bool isPP_;
    TFile * outFile_;
    float jetPtMin;

    // SimTrack
    TH2F* hsim;
    TH2F* hacc;
    TH2F* heff;
    TH2F* hmul;

    TH3F* hsim3D;
    TH3F* heff3D;
    TH3F* hmul3D;
    TH3F* hresStoR3D;

    // RecTrack
    TH2F* hrec;
    TH2F* hfak;
    TH2F* hsec;

    TH3F* hrec3D;
    TH3F* hfak3D;
    TH3F* hsec3D;

    // vector of histograms
    std::vector<TH3F*> vhsim3D;
    std::vector<TH3F*> vheff3D;
    std::vector<TH3F*> vhmul3D;

    std::vector<TH3F*> vhrec3D;
    std::vector<TH3F*> vhfak3D;
    std::vector<TH3F*> vhsec3D;

    std::vector<TH3F*> vhresStoR3D;

    // monitors
    std::vector<TH2F*> vhTrkJetPtDr;
    std::vector<TH2F*> vhSimJetPtDr;

    // methods
    TrkCorrHisAna(TString name, TFile * outf, float jetPtMin=40, bool pp=false);
    void DeclareHistograms();
    void FillRecHistograms(const EvtSel & evt, const DiJet & gj, const RecTrack_t & r);
    void FillSimHistograms(const EvtSel & evt, const DiJet & gj, const SimTrack_t & s);
};


TrkCorrHisAna::TrkCorrHisAna(TString name, TFile * outf, float jPtMin, bool pp) :
  name_(name),
  isPP_(pp),
  trkPhiMode_(false),
  jetPtMin(jPtMin)
{
   outFile_ = outf;
   
   // pt bins
   const double small = 1e-3;
   double pt;
   for(pt =   0.2  ; pt <   1.2-small; pt +=  0.05) ptBins.push_back(pt); // 20 bins
   for(pt =   1.2; pt <   2.0-small; pt +=  0.1 ) ptBins.push_back(pt); // 8 bins
   const Int_t numHigtPtBins=18;
   Float_t highPtBins[numHigtPtBins+1] = {2,2.5,3,4,5,7.5,10,12,15,20,25,30,45,60,90,120,180,300,500};
   ptBins.insert(ptBins.end(),highPtBins,highPtBins+numHigtPtBins+1);   // eta bins

   // eta bins
   double etaMin   = -2.4;
   double etaMax   =  2.4;
   double etaWidth =  0.4;
   for(double eta = etaMin; eta < etaMax + etaWidth/2; eta += etaWidth)
    etaBins.push_back(eta);

   // phi bins
   double phiMin   = -TMath::Pi();
   double phiMax   =  TMath::Pi();
   double phiWidth =  (phiMax - phiMin)/10;
   for(double phi = phiMin; phi < phiMax + phiWidth/2; phi += phiWidth)
      phiBins.push_back(phi);
       
   //jet bins
//    const Int_t numJetBins=10;
//    Float_t jBins[numJetBins+1] = {0,20,40,60,80,120,160,200,250,500,1000};
   const Int_t numJetBins=5;
   Float_t jBins[numJetBins+1] = {0,jetPtMin,80,120,200,1000};
   jetBins.insert(jetBins.end(),jBins,jBins+numJetBins+1);

   //centrality bins
   if (!isPP_) {
      centBins.push_back(0);
      centBins.push_back(12);
      centBins.push_back(40);
   } else {
      centBins.push_back(0);
      centBins.push_back(40);
   }
}

void TrkCorrHisAna::DeclareHistograms()
{
   ////////////////////////////////////////////
   // Set 3rd Dimension
   ////////////////////////////////////////////
   zBins = jetBins;
   TString zTitle = "jet E_{T} (GeV/c)";

   if (trkPhiMode_) {
      zBins = phiBins;   
      zTitle = "#phi";
   }
   
   cout << "===== " << name_ << " =====" << endl;
	cout << endl << "Pt " << ptBins.size()-1 << " bins:";
	for (Int_t i=0; i<ptBins.size(); ++i) cout << ptBins[i] << " ";
	cout << endl << "Eta " << etaBins.size()-1 << " bins:";
	for (Int_t i=0; i<etaBins.size(); ++i) cout << etaBins[i] << " ";
	cout << endl << "zaxis: " << zBins.size()-1 << " bins:";
	for (Int_t i=0; i<zBins.size(); ++i) cout << zBins[i] << " ";
	cout << endl;
	
   // Setup output dir
   if (!outFile_) {
      cout << "No outfile defined" << endl;
      exit(1);
   }
   outFile_->mkdir(name_);
   outFile_->cd(name_);
   
   // simulated
   hsim = new TH2F("hsim","Sim Tracks;#eta;p_{T} (GeV/c)",
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0]);
   
   // accepted
   hacc = new TH2F("hacc","Accepted Tracks;#eta;p_{T} (GeV/c)",
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0]);
   
   // efficiency
   heff = new TH2F("heff","Effic Rec Tracks;#eta;p_{T} (GeV/c)",
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0]);
   
   // multiply reconstructed
   hmul = new TH2F("hmul","Mult Rec Tracks;#eta;p_{T} (GeV/c)",
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0]);
   
   // reconstructed
   hrec = new TH2F("hrec","Rec Tracks;#eta;p_{T} (GeV/c)",
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0]);
   
   // fakes
   hfak = new TH2F("hfak","Fake Tracks;#eta;p_{T} (GeV/c)",
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0]);
   
   // secondary
   hsec = new TH2F("hsec","Secondary Tracks;#eta;p_{T} (GeV/c)",
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0]);
   
   // simulated 3D 
   hsim3D = new TH3F("hsim3D","Sim Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0],
      zBins.size()-1, &zBins[0]);
   
   // efficiency  3D 
   heff3D = new TH3F("heff3D","Effic Rec Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0],
      zBins.size()-1, &zBins[0]);
   
   // multiply reconstructed 3D 
   hmul3D = new TH3F("hmul3D","Mult Rec Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0],
      zBins.size()-1, &zBins[0]);
   
   
   // reconstructed 3D 
   hrec3D = new TH3F("hrec3D","Rec Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0],
      zBins.size()-1, &zBins[0]);
   
   // fakes 3D 
   hfak3D = new TH3F("hfak3D","Fake Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0],
      zBins.size()-1, &zBins[0]);
   
   // secondary
   hsec3D = new TH3F("hsec3D","Secondary Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0],
      zBins.size()-1, &zBins[0]);
   
   // mom resolution (Sim to Rec) 
   hresStoR3D = new TH3F("hresStoR3D","Momentum resolution (sim to rec);#eta;sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
      etaBins.size()-1, &etaBins[0],
      ptBins.size()-1, &ptBins[0],
      ptBins.size()-1, &ptBins[0]);

   for(unsigned i=0;i<centBins.size()-1;i++){
      vhsim3D.push_back(new TH3F("","Sim Tracks;#eta;p_{T} (GeV/c);"+zTitle, etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], zBins.size()-1, &zBins[0]) );
      vheff3D.push_back(new TH3F("","Effic Rec Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], zBins.size()-1, &zBins[0]) );
      vhmul3D.push_back(new TH3F("","Mult Rec Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], zBins.size()-1, &zBins[0]) );
      vhrec3D.push_back(new TH3F("","Rec Rec Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], zBins.size()-1, &zBins[0]) );
      vhfak3D.push_back(new TH3F("","Fake Rec Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], zBins.size()-1, &zBins[0]) );
      vhsec3D.push_back(new TH3F("","Secondary Tracks;#eta;p_{T} (GeV/c);"+zTitle,
      etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], zBins.size()-1, &zBins[0]) );
      vhresStoR3D.push_back(new TH3F("","Momentum resolution (sim to rec);#eta;sim p_{T} (GeV/c);rec p_{T} (GeV/c)", etaBins.size()-1, &etaBins[0], ptBins.size()-1, &ptBins[0], ptBins.size()-1, &ptBins[0]) );

      vhsim3D[i]->SetName(Form("hsim3D_cbin%dto%d",centBins[i],centBins[i+1]));
      vheff3D[i]->SetName(Form("heff3D_cbin%dto%d",centBins[i],centBins[i+1]));
      vhrec3D[i]->SetName(Form("hrec3D_cbin%dto%d",centBins[i],centBins[i+1]));
      vhfak3D[i]->SetName(Form("hfak3D_cbin%dto%d",centBins[i],centBins[i+1]));
      vhmul3D[i]->SetName(Form("hmul3D_cbin%dto%d",centBins[i],centBins[i+1]));
      vhsec3D[i]->SetName(Form("hsec3D_cbin%dto%d",centBins[i],centBins[i+1]));
      vhresStoR3D[i]->SetName(Form("hresStoR3D_cbin%dto%d",centBins[i],centBins[i+1]));
   }

   // monitors
for (int j=0; j<2; ++j) vhTrkJetPtDr.push_back(new TH2F(Form("hTrkJet%dDr",j+1),Form(";#DeltaR(trk,jet%d);p_{T} (GeV/c);",j+1),50,0,8,50,0,100));
for (int j=0; j<2; ++j) vhSimJetPtDr.push_back(new TH2F(Form("hSimJet%dDr",j+1),Form(";#DeltaR(simtrk,jet%d);p_{T} (GeV/c);",j+1),50,0,8,50,0,100));
}

void TrkCorrHisAna::FillRecHistograms(const EvtSel & evt, const DiJet & gj, const RecTrack_t & r)
{
   float zvar = r.jet;
   if (trkPhiMode_) zvar = r.phir;

   // monitor
   vhTrkJetPtDr[0]->Fill(deltaR(r.etar,r.phir,gj.eta1,gj.phi1),r.ptr);
   vhTrkJetPtDr[1]->Fill(deltaR(r.etar,r.phir,gj.eta2,gj.phi2),r.ptr);

   // corrections
   hrec->Fill(r.etar, r.ptr);
   hrec3D->Fill(r.etar, r.ptr, zvar);
   if(!r.nsim) hfak->Fill(r.etar, r.ptr), hfak3D->Fill(r.etar, r.ptr, zvar);
   if(r.nsim>0 && r.status<0) hsec->Fill(r.etar, r.ptr), hsec3D->Fill(r.etar, r.ptr, zvar); // nsim>0 redudant?
   
   // filling histogram in vector
   for(unsigned i=0;i<centBins.size()-1;i++){
      if(evt.cBin>=centBins[i] && evt.cBin<centBins[i+1]){
         vhrec3D[i]->Fill(r.etar, r.ptr, zvar);
         if(!r.nsim) vhfak3D[i]->Fill(r.etar, r.ptr, zvar);
         if(r.nsim>0 && r.status<0) vhsec3D[i]->Fill(r.etar, r.ptr, zvar);
      }
   } // end of vector loop
}

void TrkCorrHisAna::FillSimHistograms(const EvtSel & evt, const DiJet & gj, const SimTrack_t & s) {
   float zvar = s.jet;
   if (trkPhiMode_) zvar = s.phis;

   if(s.status>0) {
      // monitor
      vhSimJetPtDr[0]->Fill(deltaR(s.etas,s.phis,gj.eta1,gj.phi1),s.pts);
      vhSimJetPtDr[1]->Fill(deltaR(s.etas,s.phis,gj.eta2,gj.phi2),s.pts);
      
      // corrections
      hsim->Fill(s.etas, s.pts);
      hsim3D->Fill(s.etas, s.pts, zvar);
      if(s.acc)    hacc->Fill(s.etas, s.pts);
      if(s.nrec>0) heff->Fill(s.etas, s.pts), heff3D->Fill(s.etas, s.pts, zvar);
      if(s.nrec==1) hresStoR3D->Fill(s.etas, s.pts, s.ptr);
      if(s.nrec>1) hmul->Fill(s.etas, s.pts), hmul3D->Fill(s.etas, s.pts, zvar);
      
      // filling histogram in vector 
      for(unsigned i=0;i<centBins.size()-1;i++){
         if(evt.cBin>=centBins[i] && evt.cBin<centBins[i+1]){
            vhsim3D[i]->Fill(s.etas, s.pts, zvar);
            if(s.nrec>0) vheff3D[i]->Fill(s.etas, s.pts, zvar);
            if(s.nrec==1) vhresStoR3D[i]->Fill(s.etas, s.pts, s.ptr);
            if(s.nrec>1) vhmul3D[i]->Fill(s.etas, s.pts, zvar);
         }
      } // end of vector loop 
   } // end of (s.status) loop 
}


class TrkReso {
public:
   // Tracking Resolution
   TF1 * fGaus;
   TF1 * fReso;
   TrkReso(float pt=-1) {
      fGaus = new TF1("fGaus","gaus",-1,3);
      fGaus->SetParameter(0,1); // normalization
      fGaus->SetParameter(1,1); // mean
      fReso = new TF1("fReso","[0]*pow(x,[3])/(1+exp([1]*(x+[2]))) + [4]*pow(x,[5])",0.8,120);
      fReso->FixParameter(0,-0.351);
      fReso->FixParameter(1,0.0177);
      fReso->FixParameter(2,9.084);
      fReso->FixParameter(3,-0.174);
      fReso->FixParameter(4,0.221);
      fReso->FixParameter(5,-0.299);
      if (pt>0.5) {
         float reso = fReso->Eval(pt);
         fGaus->SetParameter(2,reso); // resolution
      }
     }
   float GetSmear(float pt) {
      float reso = fReso->Eval(pt);
      fGaus->SetParameter(2,reso);
      float sm = fGaus->GetRandom();
      return sm;
   }
   float GetSmear() {
      return fGaus->GetRandom();
   }
};
#endif //analyzeTrackingCorrection_h
