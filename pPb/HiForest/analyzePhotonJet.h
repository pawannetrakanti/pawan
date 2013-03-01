#ifndef analyzePhotonJet_h
#define analyzePhotonJet_h

// Convinient Output Classes
class EvtSel {
public:
   int run,evt,nOccur,cBin;
   int nG,nJ,nT;
   bool trig,offlSel,noiseFilt,anaEvtSel;
   float vz,weight,npart,ncoll,sampleWeight,evtPlane,pthat,samplePtHat;
   TString leaves;
   EvtSel(){
      leaves="run/I:evt:nOccur:cBin:nG:nJ:nT:trig/O:offlSel:noiseFilt:anaEvtSel:vz/F:weight:npart:ncoll:sampleWeight:evtPlane:pthat:samplePtHat";
   }
};

static const int MAXTRK = 10000;

class GammaJet{
public:
   GammaJet() :
   photonEt(-99),photonEta(0),photonPhi(0),
   jetEt(-99),jetEta(0),jetPhi(0),
   deta(-99),dphi(-99), Aj(-99),
   sigmaIetaIeta(-99),
   isEle(false),
   nTrk(0),nJet(0) {
      leaves = "photonEt/F:photonRawEt:photonEta:photonPhi:jetEt:jetRawEt:jetEta:jetPhi:deta:dphi:Agj:hovere:sigmaIetaIeta:sumIsol"
      ":phoMatJetEt:phoMatJetEta:phoMatJetPhi"
      ":ltrkPt:ltrkEta:ltrkPhi:ltrkJetDr:jltrkPt:jltrkEta:jltrkPhi:jltrkJetDr"
      ":jlpfPt:jlpfEta:jlpfPhi:jlpfJetDr:jlpfId:pfPhoPt"
      ":refPhoPt:refPhoEta:refPhoPhi:refPhoFlavor:genCalIsoDR04:refJetEt:refJetEta:refJetPhi:refPartonPt:refPartonFlavor"
      ":genPhoPt:genPhoEta:genPhoPhi:genJetPt:genJetEta:genJetPhi"
      ":isEle/O";
   }
   float photonEt,photonRawEt,photonEta,photonPhi;
   float jetEt,jetRawEt,jetEta,jetPhi;
   float deta,dphi,Aj;
   float hovere,sigmaIetaIeta,sumIsol;
   float phoMatJetEt,phoMatJetEta,phoMatJetPhi;
   float ltrkPt,ltrkEta,ltrkPhi,ltrkJetDr,jltrkPt,jltrkEta,jltrkPhi,jltrkJetDr;
   float jlpfPt,jlpfEta,jlpfPhi,jlpfJetDr,jlpfId,pfPhoPt;
   float refPhoPt,refPhoEta,refPhoPhi,refPhoFlavor,genCalIsoDR04,refJetEt,refJetEta,refJetPhi,refPartonPt,refPartonFlavor;
   float genPhoPt,genPhoEta,genPhoPhi,genJetPt,genJetEta,genJetPhi;
   bool isEle;
   int nTrk;
   float trkPt[MAXTRK];
   float trkEta[MAXTRK];
   float trkPhi[MAXTRK];   
   float trkJetDr[MAXTRK];
   int nJet;
   float inclJetPt[MAXTRK];
   float inclJetEta[MAXTRK];
   float inclJetPhi[MAXTRK];   
   float inclJetRefPt[MAXTRK];
   float inclJetRefPartonPt[MAXTRK];
   float inclJetRefResp[MAXTRK];
   int nGenJet;
   float inclGenJetPt[MAXTRK];
   float inclGenJetEta[MAXTRK];
   float inclGenJetPhi[MAXTRK];   
   float inclGenJetResp[MAXTRK];
   int nPf;
   float pfPt[MAXTRK];
   float pfEta[MAXTRK];
   float pfPhi[MAXTRK];
   int pfId[MAXTRK];
   int nGenp;
   float genpPt[MAXTRK];
   float genpEta[MAXTRK];
   float genpPhi[MAXTRK];
   int genpCh[MAXTRK];
   TString leaves;
   void clear() {
      photonEt=-99; photonEta=-99; photonPhi=-99;
      jetEt=-99; jetRawEt=-99; jetEta=-99; jetPhi=-99;
      deta=-99; dphi=-99; Aj=-99;
      sigmaIetaIeta=-99;
      phoMatJetEt=-99; phoMatJetEta=-99; phoMatJetPhi=-99;
      ltrkPt=-99; ltrkEta=-99; ltrkPhi=-99; ltrkJetDr=-99;
      jltrkPt=-99; jltrkEta=-99; jltrkPhi=-99; jltrkJetDr=-99;
      jlpfPt=-99; jlpfEta=-99; jlpfPhi=-99; jlpfJetDr=-99; jlpfId=-99; pfPhoPt=0;
      refPhoPt=-99; refPhoFlavor=-99; refJetEt=-99; refJetEta=-99; refJetPhi=-99; refPartonPt=-99; refPartonFlavor=-99;
      genPhoPt=-99; genPhoEta=-99; genPhoPhi=-99; genJetPt=-99; genJetEta=-99; genJetPhi=-99;
      isEle=false;
      nTrk=0; nJet=0; nGenJet=0; nPf=0; nGenp=0;
   }
};

class Isolation{
public:
   float cc1,cc2,cc3,cc4,cc5;
   float cr1,cr2,cr3,cr4,cr5;
   float ct1PtCut20,ct2PtCut20,ct3PtCut20,ct4PtCut20,ct5PtCut20;
   TString leaves;
   Isolation() {
      leaves = "cc1:cc2:cc3:cc4:cc5:cr1:cr2:cr3:cr4:cr5:ct1PtCut20:ct2PtCut20:ct3PtCut20:ct4PtCut20:ct5PtCut20";
   }
   void Set(HiForest * c, int j) {
      cc1=c->photon.cc1[j];
      cc2=c->photon.cc2[j];
      cc3=c->photon.cc3[j];
      cc4=c->photon.cc4[j];
      cc5=c->photon.cc5[j];
      cr1=c->photon.cr1[j];
      cr2=c->photon.cr2[j];
      cr3=c->photon.cr3[j];
      cr4=c->photon.cr4[j];
      cr5=c->photon.cr5[j];
      ct1PtCut20=c->photon.ct1PtCut20[j];
      ct2PtCut20=c->photon.ct2PtCut20[j];
      ct3PtCut20=c->photon.ct3PtCut20[j];
      ct4PtCut20=c->photon.ct4PtCut20[j];
      ct5PtCut20=c->photon.ct5PtCut20[j];
   }
};


void BookGJBranches(TTree * tgj, EvtSel & evt, GammaJet & gj, Isolation & isol) {
   tgj->Branch("evt",&evt.run,evt.leaves);
   tgj->Branch("jet",&gj.photonEt,gj.leaves);
   tgj->Branch("isolation",&isol.cc1,isol.leaves);
   tgj->Branch("nTrk",&gj.nTrk,"nTrk/I");
   tgj->Branch("trkPt",gj.trkPt,"trkPt[nTrk]/F");
   tgj->Branch("trkEta",gj.trkEta,"trkEta[nTrk]/F");
   tgj->Branch("trkPhi",gj.trkPhi,"trkPhi[nTrk]/F");
   tgj->Branch("trkJetDr",gj.trkJetDr,"trkJetDr[nTrk]/F");
   tgj->Branch("nJet",&gj.nJet,"nJet/I");
   tgj->Branch("inclJetPt",gj.inclJetPt,"inclJetPt[nJet]/F");
   tgj->Branch("inclJetEta",gj.inclJetEta,"inclJetEta[nJet]/F");
   tgj->Branch("inclJetPhi",gj.inclJetPhi,"inclJetPhi[nJet]/F");
   tgj->Branch("inclJetRefPt",gj.inclJetRefPt,"inclJetRefPt[nJet]/F");
   tgj->Branch("inclJetRefPartonPt",gj.inclJetRefPartonPt,"inclJetRefPartonPt[nJet]/F");
   tgj->Branch("inclJetRefResp",gj.inclJetRefResp,"inclJetRefResp[nJet]/F");
   tgj->Branch("nGenJet",&gj.nGenJet,"nGenJet/I");
   tgj->Branch("inclGenJetPt",gj.inclGenJetPt,"inclGenJetPt[nGenJet]/F");
   tgj->Branch("inclGenJetEta",gj.inclGenJetEta,"inclGenJetEta[nGenJet]/F");
   tgj->Branch("inclGenJetPhi",gj.inclGenJetPhi,"inclGenJetPhi[nGenJet]/F");
   tgj->Branch("inclGenJetResp",gj.inclGenJetResp,"inclGenJetResp[nGenJet]/F");
   tgj->Branch("nPf",&gj.nPf,"nPf/I");
   tgj->Branch("pfPt",gj.pfPt,"pfPt[nPf]/F");
   tgj->Branch("pfEta",gj.pfEta,"pfEta[nPf]/F");
   tgj->Branch("pfPhi",gj.pfPhi,"pfPhi[nPf]/F");
   tgj->Branch("pfId",gj.pfId,"pfId[nPf]/I");
   tgj->Branch("nGenp",&gj.nGenp,"nGenp/I");
   tgj->Branch("genpPt",gj.genpPt,"genpPt[nGenp]/F");
   tgj->Branch("genpEta",gj.genpEta,"genpEta[nGenp]/F");
   tgj->Branch("genpPhi",gj.genpPhi,"genpPhi[nGenp]/F");
   tgj->Branch("genpCh",gj.genpCh,"genpCh[nGenp]/I");
}

class CentralityReWeight {
public:
   CentralityReWeight(TString data, TString mc,TCut s) : 
   datafname(data),mcfname(mc),sel(s) {}
   void Init()
   {
      cout << "Reweight Centrality: " << endl;
      cout << "Data file: " << datafname << endl;
      cout << "MC file:   " << mcfname << endl;
      TChain * tdata = new TChain("tgj");
      TChain * tmc = new TChain("tgj");
      tdata->Add(datafname);
      tmc->Add(mcfname);
      
      hCentData = new TH1D("hCentData","",40,0,40);
      hCentMc = new TH1D("hCentMc","",40,0,40);
      hReWt = new TH1D("hReWt","",40,0,40);
      
      //cout << "data: " << tdata->GetName() << " " << tdata->GetEntries() << endl;
      //cout << "mc: " << tmc->GetName() << " " << tmc->GetEntries() << endl;
      tdata->Project("hCentData","cBin",sel&&"anaEvtSel");
      tmc->Project("hCentMc","cBin",sel);
      hCentData->Scale(1./hCentData->Integral());
      hCentMc->Scale(1./hCentMc->Integral());
      hReWt->Divide(hCentData,hCentMc);
   }
   float GetWeight(int cBin) {
      int bin=cBin+1;
      if (hCentData->GetBinContent(bin)==0 || hCentMc->GetBinContent(bin)==0) {
         return 0; 
      }
      return hCentData->GetBinContent(bin)/hCentMc->GetBinContent(bin);
   }
   TString datafname,mcfname;
   TCut sel;
   TH1D * hCentData;
   TH1D * hCentMc;
   TH1D * hReWt;
};

float getNpart(int cBin) { 
   if (cBin ==0) return 393.633;
   if (cBin ==1) return 368.819;
   if (cBin ==2) return 343.073;
   if (cBin ==3) return 317.625;
   if (cBin ==4) return 292.932;
   if (cBin ==5) return 271.917;
   if (cBin ==6) return 249.851;
   if (cBin ==7) return 230.72;
   if (cBin ==8) return 212.465;
   if (cBin ==9) return 194.752;
   if (cBin ==10) return 178.571;
   if (cBin ==11) return 163.23;
   if (cBin ==12) return 149.187;
   if (cBin ==13) return 136.011;
   if (cBin ==14) return 123.414;
   if (cBin ==15) return 111.7;
   if (cBin ==16) return 100.831;
   if (cBin ==17) return 90.7831;
   if (cBin ==18) return 80.9823;
   if (cBin ==19) return 72.6236;
   if (cBin ==20) return 64.1508;
   if (cBin ==21) return 56.6284;
   if (cBin ==22) return 49.9984;
   if (cBin ==23) return 43.3034;
   if (cBin ==24) return 37.8437;
   if (cBin ==25) return 32.6659;
   if (cBin ==26) return 27.83;
   if (cBin ==27) return 23.7892;
   if (cBin ==28) return 20.1745;
   if (cBin ==29) return 16.8453;
   if (cBin ==30) return 14.0322;
   if (cBin ==31) return 11.602;
   if (cBin ==32) return 9.52528;
   if (cBin ==33) return 7.6984;
   if (cBin ==34) return 6.446;
   if (cBin ==35) return 4.96683;
   if (cBin ==36) return 4.23649;
   if (cBin ==37) return 3.50147;
   if (cBin ==38) return 3.16107;
   if (cBin ==39) return 2.7877;
   return -1;
}

float getNcoll(int cBin) { 
   if (cBin == 0) return  1747.86 ;
   if (cBin == 1) return  1567.53 ;
   if (cBin == 2) return  1388.39 ;
   if (cBin == 3) return  1231.77 ;
   if (cBin == 4) return  1098.2 ;
   if (cBin == 5) return  980.439 ;
   if (cBin == 6) return  861.609 ;
   if (cBin == 7) return  766.042 ;
   if (cBin == 8) return  676.515 ;
   if (cBin == 9) return  593.473 ;
   if (cBin == 10) return  521.912 ;
   if (cBin == 11) return  456.542 ;
   if (cBin == 12) return  398.546 ;
   if (cBin == 13) return  346.647 ;
   if (cBin == 14) return  299.305 ;
   if (cBin == 15) return  258.344 ;
   if (cBin == 16) return  221.216 ;
   if (cBin == 17) return  188.677 ;
   if (cBin == 18) return  158.986 ;
   if (cBin == 19) return  134.7 ;
   if (cBin == 20) return  112.547 ;
   if (cBin == 21) return  93.4537 ;
   if (cBin == 22) return  77.9314 ;
   if (cBin == 23) return  63.5031 ;
   if (cBin == 24) return  52.0469 ;
   if (cBin == 25) return  42.3542 ;
   if (cBin == 26) return  33.9204 ;
   if (cBin == 27) return  27.3163 ;
   if (cBin == 28) return  21.8028 ;
   if (cBin == 29) return  17.2037 ;
   if (cBin == 30) return  13.5881 ;
   if (cBin == 31) return  10.6538 ;
   if (cBin == 32) return  8.35553 ;
   if (cBin == 33) return  6.40891 ;
   if (cBin == 34) return  5.13343 ;
   if (cBin == 35) return  3.73215 ;
   if (cBin == 36) return  3.06627 ;
   if (cBin == 37) return  2.41926 ;
   if (cBin == 38) return  2.11898 ;
   if (cBin == 39) return  1.76953 ;
   return -1;
}

class Response {
public:
   // jet energy studies
   // Default and user constructor.
   TF1 * fGaus;
   TF1 * vfReso[4];
   Response() {
      fGaus = new TF1("fGaus","gaus",-1,3);
      fGaus->SetParameter(0,1); // normalization
      fGaus->SetParameter(1,1); // mean
      vfReso[0] = new TF1("fres0","sqrt(pow(0.0246,2)+pow(1.213/sqrt(x),2)+pow(5.23/x,2))",0,300);
      vfReso[1] = new TF1("fres1","sqrt(pow(0.0246,2)+pow(1.213/sqrt(x),2)+pow(5.10/x,2))",0,300);
      vfReso[2] = new TF1("fres2","sqrt(pow(0.0246,2)+pow(1.213/sqrt(x),2)+pow(3.88/x,2))",0,300);
      vfReso[3] = new TF1("fres3","sqrt(pow(0.0246,2)+pow(1.213/sqrt(x),2)+pow(0.001/x,2))",0,300);
   }
   float GetSmear(int cBin, float pt) {
      int mybin;
      if (cBin<4) mybin=0;
      else if (cBin<12) mybin=1;
      else if (cBin<20) mybin=2;
      else mybin=3;

      float reso = vfReso[mybin]->Eval(pt);
      fGaus->SetParameter(2,reso);
      float sm = fGaus->GetRandom();
//      cout << "GetSmear for: " << mybin << " " << pt << " reso: " << reso << " gaus: " << smpt/pt << " smeared: " << smpt << endl;
      return sm;
   }
};
#endif