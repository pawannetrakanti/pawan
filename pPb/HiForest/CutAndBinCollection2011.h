#ifndef CutAndBinCollection_C
#define CutAndBinCollection_C
double a1, a2, a3, a4, a5, a6;

#include <TCut.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <TMath.h>
#include <math.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "commonUtility.h"
#include "multiTreeUtilPhoton2011.h"
#include "histFunctionD.C"
#define PI 3.141592653589

// DON'T FORGET TO APPLY HOE CUT SEPARATELY    
// Convinient Output Classes                                                                                                                 

struct valPair {
  double val;
  double err;
};



class fitResult {
 public:
   double nSig; 
   double nSigErr;
   double purity010;
};

TCut FisherCut = "((4.6774452168946995e-01 +(cc5-cc1) *9.6796013515455164e-04 +(cc4-cc1) *-1.2788583705647016e-02 +(cc3-cc1) *1.3674667235554151e-02 +(cc2-cc1) *-7.4576842527504350e-02 +(cr5) *-1.2105656031270820e-02 +(cr4) *-1.8158418924831903e-03 +(cr3) *-1.4267772594659891e-02 +(cr2) *-3.0555858981100050e-02 +(ct5PtCut20) *6.4351309460660500e-03 +(ct4PtCut20) *-7.3308097775112357e-03 +(ct3PtCut20) *-2.2250052480332189e-03 +(ct2PtCut20) *-2.4645948244900417e-02 +(ct1PtCut20) *2.1959860851889978e-03 > 0.3))";

TString fisherVar = "4.6774452168946995e-01 +(cc5-cc1) *9.6796013515455164e-04 +(cc4-cc1) *-1.2788583705647016e-02 +(cc3-cc1) *1.3674667235554151e-02 +(cc2-cc1) *-7.4576842527504350e-02 +(cr5) *-1.2105656031270820e-02 +(cr4) *-1.8158418924831903e-03 +(cr3) *-1.4267772594659891e-02 +(cr2) *-3.0555858981100050e-02 +(ct5PtCut20) *6.4351309460660500e-03 +(ct4PtCut20) *-7.3308097775112357e-03 +(ct3PtCut20) *-2.2250052480332189e-03 +(ct2PtCut20) *-2.4645948244900417e-02 +(ct1PtCut20) *2.1959860851889978e-03";

TCut iso3dCut  = "cc4 < 6.9 && ct4PtCut20 < 3.00 && cr4<5"; //cc4 <2 && cr4 < 2.2 && ct4PtCut20 < 2.4";
TCut isoSumCut  = "(cc4+cr4+ct4PtCut20)/0.9 <1";

TCut sbIsoCut =" (cc4+cr4+ct4PtCut20)/0.9>10 && (cc4+cr4+ct4PtCut20)/0.9 < 20 ";

int nBinsExt = 2500;
const double lumiPP = 231.;
const double nMBHI  = 55.69e6 ; //52.851e6  ; //695692e6 ; 
double nEventsHI[4] = {0,   nMBHI* 0.10,  nMBHI * 0.20,  nMBHI * 0.70 };
double taa[5]       = {0,          23.2e-9,          11.6e-9,           3.92e-09, 4.69e-10  };  // pico barn^1         



TString fnameDATA_2010   = "/d100/yjlee/hiForest/merged_HI2010_SD_Photon40_prod02.root";
TString fnameDATA_2011   = "/d100/velicanu/forest/HiForestPhoton35_Skim.root";


//const int nPtBin = 4;
//double ptBin[nPtBin+1] = {50,60,80,120,200};

//const int nPtBin = 1;
//double ptBin[nPtBin+1] = {60,300};

const int nPtBin = 2;    
double ptBin[nPtBin+1] = {50,60,80};


const int nCent_std = 4;
double centBin_std[nCent_std+1] = {0,4,12,20,40};

TCut univHoECutPP = "(hadronicOverEm<0.05)";
TCut univHoECutHI = "(hadronicOverEm<0.1)";

float isolationCut = 5.0;

TCut isFragment = "abs(genMomId)<22";
TCut isPrompt = "abs(genMomId)==22";

TCut genMatchCut0      = "isGenMatched && abs(genMomId)<=22";
TCut genMatchCut1      = Form("isGenMatched && genMomId==22 && genCalIsoDR04 < %.1f",isolationCut);
TCut genMatchCut      = Form("(isGenMatched && abs(genMatchedEta)<1.44 && abs(etCorrected/genMatchedPt-1)<.3 && abs(genMomId) <= 22 && genCalIsoDR04 < %.1f)",isolationCut);
TCut genMatchCutBkg      = "(isGenMatched && abs(genMatchedEta)<1.44 && abs(etCorrected/genMatchedPt-1)<.6)  &&  ( (abs(genMomId) > 22) || (genCalIsoDR04 > 5.0) ) ";


TCut genPhotonCut     = Form("( abs(gpEta) < 1.44 && abs(gpId)==22 && abs(gpMomId) <= 22 && gpCollId ==0  && gpIsoDR04 < %.3f)",isolationCut);


TString swissCrx      = "(1 - (eRight+eLeft+eTop+eBottom)/eMax)";
TCut hiSpikeCutMC     = Form("(  %s < 0.90 && sigmaIetaIeta>0.002 && sigmaIphiIphi>0.002)",swissCrx.Data());
TCut ppSpikeCutMC     = Form("(  %s < 0.95 && sigmaIetaIeta>0.002 && sigmaIphiIphi>0.002)",swissCrx.Data());

//TCut hiSpikeCutNoPhi  = Form("( ( %s < 0.90 && sigmaIetaIeta>0.002) ",swissCrx.Data());      

//TCut hiSpikeCutNoPhi  = Form("( ( %s < 0.90 && sigmaIetaIeta>0.002) ",swissCrx.Data());


TCut seedTimeCut      = "abs(seedTime)<3";
TCut hiSpikeCutData   = hiSpikeCutMC && seedTimeCut;
TCut ppSpikeCutData   = ppSpikeCutMC && seedTimeCut;

TCut etaCut       = "abs(eta)<1.44 && abs(scEta)<1.479";

//TCut etaCut           = " (abs(scEta) < 1.479 && abs(eta)<1.44) && (rawEnergy/energy > 0.5)";// && (!isEBGap&&!isEEGap&&!isEBEEGap)";

TCut genEtaCut  =       "                      (abs(eta) < 1.44)";
TCut vtxCut     = "abs(vtxZ)<15";
TCut finalCutSigHI  = genMatchCut     &&  hiSpikeCutMC && etaCut && vtxCut ;
TCut finalCutBkgHI  = !genMatchCut  &&  hiSpikeCutMC && etaCut && vtxCut ;
TCut finalCutDataHI =                     hiSpikeCutData && etaCut  && vtxCut;

TCut finalCutSigPP  = genMatchCut     &&  ppSpikeCutMC && etaCut && vtxCut ;
TCut finalCutBkgPP  = !genMatchCut  &&  ppSpikeCutMC && etaCut && vtxCut;
TCut finalCutDataPP =                     ppSpikeCutData && etaCut  && vtxCut;



TCut finalCutGen  =  genPhotonCut && vtxCut ;


/*      
	photon15 :   5.414e-05                                                        
	photon30 :   4.971e-06                                                        
	photon50 :   6.650e-07                                                        
	
	dijet 15 :   2.034e-01                                                        
	dijet 30 :   1.078e-02                                                        
	dijet 50 :   1.024e-03                                                        
	dijet 80 :   9.968e-05                                                        
*/

//photons
const double csPho15 =   3.761e-05 ; 
const double csPho30 =   3.844e-06 ;
const double csPho50 =   5.987e-07 ;
const double csPho80 =   8.566e-08 ;  
//dijets

const double csDij15 =   1.928e-01;
const double csDij20 =   5.111e-02;
const double csDij30 =   9.715e-03;
const double csDij50 =   9.198e-04;
const double csDij80 =   9.968e-05;
const double csDij130=   7.167e-06;

const double emFilter15 =  0.00930641;
const double emFilter20 =  0.00777401;
const double emFilter30 =  0.0091598;
const double emFilter50 =  0.0533272;
const double emFilter80 =  0.203695;
const double emFilter130=  0.465538;

double csEmj20  =  csDij20  * emFilter20 ;
double csEmj30  =  csDij30  * emFilter30 ;
double csEmj50  =  csDij50  * emFilter50 ;
double csEmj80  =  csDij80  * emFilter80 ;
double csEmj130 =  csDij130 * emFilter130;


// good photon efficiency
const double emFilter15pp =  3.233e-3;
const double emFilter30pp =  3.617e-2;
const double emFilter50pp =  9.5397e-2;
const double emFilter80pp =  1.16679e-1;

const double photonFilter15pp = 0.3644;
const double photonFilter30pp = 0.9404;
const double photonFilter50pp = 0.9563;


TCut ptHatCutPho15 = "ptHat>15 && ptHat <= 30";
TCut ptHatCutPho30 = "ptHat>30 && ptHat <= 50";
TCut ptHatCutPho50 = "ptHat>50";
//TCut ptHatCutPho80 = "ptHat>80";

TCut ptHatCutEmj15 = "ptHat>15 && ptHat <= 30";
//TCut ptHatCutEmj18 = "ptHat>18 && ptHat <= 30";
TCut ptHatCutEmj20  = "ptHat>20 && ptHat <= 30";
TCut ptHatCutEmj30  = "ptHat>30 && ptHat <= 50";
TCut ptHatCutEmj50  = "ptHat>50 && ptHat <= 80";
TCut ptHatCutEmj80  = "ptHat>80 && ptHat <= 130";
TCut ptHatCutEmj130 = "ptHat>130";


//float scaleDij = weightDij15 / weightPho15;

//TCut hoeOnlyCut         = Form("(hadronicOverEm<%f)",hoeCut);
//TCut recoIsoCut         = Form("(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<%f",caloCut);
//TCut recoSideBandCut    = Form("(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)>%f",  sbCut);


//int nIsoBin = 30;
//float isoMin(-1), isoMax(100), cIsoMin(-20), cIsoMax(30);
//int nShowerBin = 30;
//float ssMin(0), ssMax(0.1);


int color[9] =  {0,1,2,4,8,1,1,1,1};  // the for centrality bins...
int colorEt[9]= {0,1,2,4,8,1,1,1,1};


void mcStyle(TH1* h=0) {
   h->SetLineColor(kPink);
   h->SetFillColor(kOrange+7);
   h->SetFillStyle(3004);
}

void sbStyle(TH1* h=0) {
   h->SetLineColor(kGreen+4);
   h->SetFillColor(kGreen+1);
   h->SetFillStyle(3001);
}

void dijetStyle(TH1* h=0) {
   h->SetLineColor(kBlue);
   h->SetFillColor(kAzure-8);
   h->SetFillStyle(3001);
}

double getNoEmc (TString theFname="", TCut theCut="") {
   TFile *fd = new TFile(theFname.Data());
   TTree *ana = (TTree*) fd->FindObjectAny("Analysis");
   cout << "number of events of " << theCut.GetTitle() << "    : " << ana->GetEntries( theCut ) << endl;
   return ana->GetEntries( theCut ) ;
}



float  getHnaprt(int ibin) {
  if (ibin ==0) return  393.633;
  if (ibin ==1) return  368.819;
  if (ibin ==2) return  343.073;
  if (ibin ==3) return  317.625;
  if (ibin ==4) return  292.932;
  if (ibin ==5) return  271.917;
  if (ibin ==6) return  249.851;
  if (ibin ==7) return  230.72;
  if (ibin ==8) return  212.465;
  if (ibin ==9) return  194.752;
  if (ibin ==10) return  178.571;
  if (ibin ==11) return  163.23;
  if (ibin ==12) return  149.187;
  if (ibin ==13) return  136.011;
  if (ibin ==14) return  123.414;
  if (ibin ==15) return  111.7;
  if (ibin ==16) return  100.831;
  if (ibin ==17) return  90.7831;
  if (ibin ==18) return  80.9823;
  if (ibin ==19) return  72.6236;
  if (ibin ==20) return  64.1508;
  if (ibin ==21) return  56.6284;
  if (ibin ==22) return  49.9984;
  if (ibin ==23) return  43.3034;
  if (ibin ==24) return  37.8437;
  if (ibin ==25) return  32.6659;
  if (ibin ==26) return  27.83;
  if (ibin ==27) return  23.7892;
  if (ibin ==28) return  20.1745;
  if (ibin ==29) return  16.8453;
  if (ibin ==30) return  14.0322;
  if (ibin ==31) return  11.602;
  if (ibin ==32) return  9.52528;
  if (ibin ==33) return  7.6984;
  if (ibin ==34) return  6.446;
  if (ibin ==35) return  4.96683;
  if (ibin ==36) return  4.23649;
  if (ibin ==37) return  3.50147;
  if (ibin ==38) return  3.16107;
  if (ibin ==39) return  2.7877;
  return -100000;
}



class CutAndBinCollection
{
 public:
   CutAndBinCollection() {
   }
   ~CutAndBinCollection() {
      
   }
   //   float getHoECut( float percentBin=1, float et=30 ); 
   TString getTmplName()  { return tmplName_;};
   TString getTmplVar()   { return tmplVar_;};
   TString getTmplVarMC() { return tmplVarMC_;};
   TString getTmplVarBkg() { return tmplVarBkg_;};


   TString getTmplLeg()   { return tmplLeg_;};
   TString getTmplXTitle(){ return tmplXTitle_;};
   TCut getIBCut() {        return IBCut_;};
   TCut getSBCut() {        return SBCut_;};
   int  getNTmplBin() {     return nTmplBin_;};
   float getTmplMin() {     return tmplMin_;};
   float getTmplMax() {     return tmplMax_;};
   float getFitMin()  {     return fitMin_;};
   float getFitMax()  {     return fitMax_;};
   bool isHI()      {       return hiOrPp_;};
   TCut getIBCutWithPreWOEleVeto(); 
   TCut getIBCutWithPre();
   TCut getSBCutWithPre();
   
   void setTmpl(int tmplOpt);
 private:
   
   TString tmplName_;
   TString tmplVar_;
   TString tmplVarMC_;
   TString tmplVarBkg_;
   
   TString tmplLeg_;
   TString tmplXTitle_;
   TCut IBCutWOEleVeto_;
   TCut IBCut_;
   TCut SBCut_;
   int nTmplBin_;
   float tmplMin_;
   float tmplMax_;
   float fitMin_;
   float fitMax_;
   bool hiOrPp_;
};

void CutAndBinCollection::setTmpl( int tmplOpt ) { 
   
   if ( tmplOpt==100) {  // sigma Ieta Ieta                                                                            
      hiOrPp_        = true;
      tmplName_        = "showerShape";
      // tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";                        
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      
      IBCutWOEleVeto_ =" ( cc4 + cr4+ ct4j20 )/0.9 < 5";
      IBCut_  =  "!isEle" && IBCutWOEleVeto_;
      SBCut_  = " !isEle && ( cc4 + cr4+ ct4j20 )/0.9 > 6  && ( cc4 + cr4 + ct4j20 )/0.9 < 11";
      
      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                         
      tmplMax_             = 0.025;
      fitMin_              = 0.005;
      fitMax_              = 0.025;
   }

   if ( tmplOpt==1001) {  // sigma Ieta Ieta                                            
      hiOrPp_        = true;
      tmplName_        = "showerShape";
      // tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";            
      
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      
      IBCutWOEleVeto_ =" ( cc4 + cr4+ ct4j20 ) < 5*0.9";
      IBCut_  = "!isEle" && IBCutWOEleVeto_;
      SBCut_  = "!isEle &&  ( cc4 + cr4+ ct4j20 ) > 6  && ( cc4 + cr4 + ct4j20 ) < 11";
      
      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARA
      tmplMax_             = 0.025;
      fitMin_              = 0.008;
      fitMax_              = 0.020;
   }


   if ( tmplOpt==101) {  // sigma Ieta Ieta                                                                                          
      hiOrPp_        = true;
      tmplName_        = "showerShape";
      // tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";                                    
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";

      IBCutWOEleVeto_ =" (cc4 + cr4 + ct4j20)/0.9 < 5";
      IBCut_  = "!isEle" && IBCutWOEleVeto_;
      SBCut_  = "!isEle &&  (cc4 + cr4 + ct4j20)/0.9 > 6  && (cc4 + cr4 + ct4j20)/0.9  < 11";

      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                                     
      tmplMax_             = 0.025;
      fitMin_              = 0.003;
      fitMax_              = 0.025;
   }
   if ( tmplOpt==102) {  // sigma Ieta Ieta                                                                                                                                                          
      hiOrPp_        = true;
      tmplName_        = "showerShape";

      tmplVar_        = "sieie47";
      tmplVarMC_      = "sieie47";
      tmplVarBkg_     = "sieie47";
      tmplLeg_        = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_     = "Shower Shape (#sigma_{I#etaI#eta})";

      IBCutWOEleVeto_ = " (cc4 + cr4 + ct4PtCut)/0.9 < 5";
      IBCut_  = "!isEle" &&  IBCutWOEleVeto_;
      SBCut_  = "!isEle && (cc4 + cr4 + ct4PtCut)/0.9 > 6  && (cc4 + cr4 + ct4PtCut)/0.9  < 11";

      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SE
      tmplMax_             = 0.025;
      fitMin_              = 0.003;
      fitMax_              = 0.025;
   }

   if ( tmplOpt==200) {  // sigma Ieta Ieta                                                         
      hiOrPp_        = false;
      tmplName_        = "showerShape";
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      
      IBCutWOEleVeto_ =" (cc4j + cr4j + ct4j20 ) < 1.5";
      IBCut_  = "!hasPixelSeed" && IBCutWOEleVeto_;
      SBCut_  = "!hasPixelSeed &&  (cc4j + cr4j + ct4j20 ) > 2.0 && (cc4j + cr4j + ct4j20 ) < 5.0";
      
      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT                                                                                                                             
      tmplMax_             = 0.025;
      fitMin_              = 0.008;
      fitMax_              = 0.017;
   }
   
   if ( tmplOpt==201) {  // sigma Ieta Ieta                                                                                                                                                       
      hiOrPp_        = false;
      tmplName_        = "showerShape";
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      
      IBCut_  =    " (!hasPixelSeed) && (ecalRecHitSumEtConeDR04 < 4.2 + 0.003*pt) && ( hcalTowerSumEtConeDR04  < 2.2 + 0.001*pt) && ( trkSumPtHollowConeDR04  < 2.0 + 0.001*pt)";
      SBCut_  =    " (!hasPixelSeed) &&  (ecalRecHitSumEtConeDR04 < 4.2 + 0.003*pt) && (hcalTowerSumEtConeDR04   < 2.2 + 0.001*pt) && (trkSumPtHollowConeDR04 > 2.0 && trackIso <5.0 )";
      
      nTmplBin_            = 25;
      tmplMin_             = 0.00;  
      tmplMax_             = 0.025;
      fitMin_              = 0.003;
      fitMax_              = 0.025;
   }
   
   
   

   if ( tmplOpt==107) {  // sigma Ieta Ieta                                                                                                                                                                       
      hiOrPp_        = true;
      tmplName_        = "showerShape";
      // tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";                                                                                                                 
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";

      IBCut_  = " !hasPixelSeed && cc4j<3 && cr4j<2 && ct4jPtCut < 1.5";
      SBCut_  = " !hasPixelSeed && cc4j<3 && cr4j<2 && ct4jPtCut >1.5 && ct4jPtCut <3.5";

      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                                                                                                                  
      tmplMax_             = 0.025;
      fitMin_              = 0.004;
      fitMax_              = 0.014;
   }
   
   
   if ( tmplOpt== 207) {  // sigma Ieta Ieta              
      
      hiOrPp_        = false;
      tmplName_        = "showerShape";
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";

      IBCut_  = " !hasPixelSeed && cc4j<3 &&  cr4j<2 && ct4jPtCut<2.0";
      SBCut_  = " !hasPixelSeed && cc4j<3 &&  cr4j<2  && ct4jPtCut>2.0 && ct4jPtCut<5.0 ";
      
      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE C                                                      
      tmplMax_             = 0.025;
      fitMin_              = 0.002;
      fitMax_              = 0.020;
   }
   
   
   // old ones
   if ( tmplOpt==1) {
      tmplName_        = "showerShape";
      // tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";                        
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      //      IBCut_  = " (t1PtCut/et/r9<00.5) && (cc4 + cr4) < 5  &&  (ct4PtCut < 3 )";                                 
      //   SBCut_  = " (t1PtCut/et/r9<00.5) && (cc4 + cr4) < 5  &&  (ct4PtCut > 6 ) ";                                   
      
      // isEle is removed for closure test                                                                               
      IBCut_  = "  (cc4 + cr4) < 5";
      SBCut_  = "  (cc4 + cr4) > 6  &&  (cc4+cr4 ) < 10";
      
      nTmplBin_            = 40;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                         
      tmplMax_             = 0.025;
      fitMin_              = 0.006;
      fitMax_              = 0.014;
      
   }
   if ( tmplOpt==146) {
      tmplName_        = "showerShape";
      tmplVar_       = "sieie46";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      IBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<5";
      SBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)>6";
      nTmplBin_            = 40;
      tmplMin_             = 0;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                             
      tmplMax_             = 0.025;
      fitMin_              = 0.005;
      fitMax_              = 0.022;
   }
   if ( tmplOpt==147) {
      tmplName_        = "showerShape";
      tmplVar_       = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      IBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<5";
      SBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)>6";
      //  IBCut_        = "( (ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<5 )";
      //  SBCut_        = "dr41 < 0.1  && ( (ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)>6 )";

      nTmplBin_            = 40;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                             
      tmplMax_             = 0.025;
      fitMin_              = 0.005;
      fitMax_              = 0.012;
   }
   if ( tmplOpt==777) {
      tmplName_        = "showerShape";
      tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      IBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<5";
      SBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)>6";
      //  IBCut_        = "( (ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<5 )";                                                                                                 
      //  SBCut_        = "dr41 < 0.1  && ( (ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)>6 )";                                                                                  

      nTmplBin_            = 40;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                                                                                                                  
      tmplMax_             = 0.025;
      fitMin_              = 0.005;
      fitMax_              = 0.012;
   }
   
   
   if ( tmplOpt==12047) {  // trackIsohi20 used
      tmplName_        = "showerShape";
      tmplVar_       = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      IBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<15 && trackIsohi20<6";
      SBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<15 && trackIsohi20>10";
      
      nTmplBin_            = 40;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                                                                                                                  
      tmplMax_             = 0.025;
      fitMin_              = 0.005;
      fitMax_              = 0.012;
   }
      
   

   
   

   
   if ( tmplOpt==1277780) {  // sigma Ieta Ieta                                                                           
      tmplName_        = "showerShape";
      // tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";                       
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      
      IBCut_  = " !isEle &&  (cc4 + cr4+ ct4PtCut80)/.9 < 5";
      SBCut_  = " !isEle &&  (cc4 + cr4+ ct4PtCut80)/.9 > 6  && (cc4 + cr4+ ct4PtCut95)/.9 < 11";
      
      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                        
      tmplMax_             = 0.025;
      fitMin_              = 0.004;
      fitMax_              = 0.014;
   }

   if ( tmplOpt==1277790) {  // sigma Ieta Ieta                                                                                                 
      tmplName_        = "showerShape";
      // tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";                                               
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      
      IBCut_  = " !isEle &&  (cc4 + cr4+ ct4PtCut90)/.9 < 5";
      SBCut_  = " !isEle &&  (cc4 + cr4+ ct4PtCut90)/.9 > 6  && (cc4 + cr4+ ct4PtCut90)/.9 < 11";
      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                  
      tmplMax_             = 0.025;
      fitMin_              = 0.004;
      fitMax_              = 0.014;
   }




   if ( tmplOpt==127770) {  // sigma Ieta Ieta                                                                           
      tmplName_        = "showerShape";
      // tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";                       
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = "sieie47";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      //      IBCut_  = " (t1PtCut/et/r9<00.5) && (cc4 + cr4) < 5  &&  (ct4PtCut < 3 )";                                
      //   SBCut_  = " (t1PtCut/et/r9<00.5) && (cc4 + cr4) < 5  &&  (ct4PtCut > 6 ) ";                                  
      
      // isEle is removed fr closure test                                                                               
      IBCut_  = " !isEle &&  (cc4 + cr4+ ct4PtCut)/.9 < 5";
      //  IBCut_  = " !isEle &&  (ecalIso-compEcalIso + hcalIso-compHcalIso + ct4jPtCut) < 5";                          
      
      SBCut_  = " !isEle &&  (cc4 + cr4+ ct4PtCut)/.9 > 6  && (cc4 + cr4+ ct4PtCut)/.9 < 11";
      
      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                        
      tmplMax_             = 0.025;
      fitMin_              = 0.004;
      fitMax_              = 0.014;
   }
  


   if ( tmplOpt==127772) {  // sigma Ieta Ieta                                                                                               
      tmplName_        = "showerShape";
      
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = " ((sieie47 - 0.0105)*0.97 + 0.0105 + -0.000332)*( 0 <= cBin && cBin <= 3) +  ((sieie47 - 0.0095)*0.98 + 0.0095 + -0.000278)*( 4 <= cBin && cBin <= 11)+  ((sieie47 - 0.0095)*0.94 + 0.0095 + -0.000243)*( 12 <= cBin && cBin <= 39) ";



      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      IBCut_  = " !hasPixelSeed &&  (cc4 + cr4+ ct4PtCut)/.9 < 5";
      SBCut_  = " !hasPixelSeed &&  (cc4 + cr4+ ct4PtCut)/.9 > 6  && (cc4 + cr4+ ct4PtCut)/.9 < 11";
      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CU
      tmplMax_             = 0.025;
      fitMin_              = 0.004;
      fitMax_              = 0.014;
   }
   
   
   if ( tmplOpt==127773) {  // sigma Ieta Ieta                                                                         \
                                                                                                                        
      tmplName_        = "showerShape";
      
      tmplVar_       = "sieie47";
      tmplVarMC_     = "sieie47";
      tmplVarBkg_    = " ((sieie47 - 0.0105)*0.87 + 0.009 + -0.000837)*( 0 <= cBin && cBin <= 3) +   ((sieie47 - 0.0095)*0.98 + 0.0095 + -0.000429)*( 4 <= cBin && cBin <= 11) +  ((sieie47 - 0.0095)*0.87 + 0.0095 + -0.000578)*( 12 <= cBin && cBin <= 39)";



      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      IBCut_  = " !isEle &&  (cc4 + cr4+ ct4PtCut)/.9 < 5";
      SBCut_  = " !isEle &&  (cc4 + cr4+ ct4PtCut)/.9 > 6  && (cc4 + cr4+ ct4PtCut)/.9 < 11";
      nTmplBin_            = 25;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CU                                                    
      tmplMax_             = 0.025;
      fitMin_              = 0.004;
      fitMax_              = 0.014;
   }


   
   if ( tmplOpt==12778) {  // Isolation 
      tmplName_        = "Calo Iso";
      tmplVar_       = "cc4-cc1+cr4-cr1";
      tmplLeg_       = "Calo Iso(GeV)";
      tmplXTitle_    = "Bkg. subtracted Calo Iso (GeV)";
      IBCut_  = " sieie47<0.013  &&  (ct4PtCut < 5) ";
      SBCut_  = " (ct4PtCut < 5)  &&  sieie47>0.013  &&  sieie47<0.014";

      nTmplBin_            =  60;
      tmplMin_             = -40; 
      tmplMax_             =  80;
      fitMin_              = -20;
      fitMax_              =  40;
   }
   
   
   if ( tmplOpt==20001) { 
     tmplName_        = "showerShape";
     tmplVar_       = "sieie46*( (cBin<8) && (pt<30) ) + sieie47*( (cBin>=8) || (pt>=30) ) ";
     tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
     tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
     IBCut_  = "(t1PtCut/et/r9<00.5&&abs(eta)<1.44&&rawEnergy/energy>0.5) && (cr4+cc4)<5 && (ct4PtCut - ct1PtCut < 3)";
     SBCut_  = "(t1PtCut/et/r9<00.5&&abs(eta)<1.44&&rawEnergy/energy>0.5) && (cr4+cc4)<5 && (ct4PtCut - ct1PtCut > 6)";
     
     
     nTmplBin_            = 40;
     tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY \ \
     tmplMax_             = 0.025;
     fitMin_              = 0.005;
     fitMax_              = 0.012;

   }

   
   
   if ( tmplOpt==145) {
      tmplName_        = "showerShape";
      tmplVar_       = "sieie45";
      tmplLeg_       = "Shower Shape (#sigma_{I#etaI#eta})";
      tmplXTitle_    = "Shower Shape (#sigma_{I#etaI#eta})";
      //IBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<5";
      // SBCut_  = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)>6";
      IBCut_        = "( (ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<5 )";                        
      SBCut_        = "( (ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)>6 && (ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)<10)";         

      nTmplBin_            = 40;
      tmplMin_             = 0.00;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                                         
      tmplMax_             = 0.025;
      fitMin_              = 0.005;
      fitMax_              = 0.020;
   }


   else if ( tmplOpt==2) { 
      tmplName_        = "rechitIso";
      tmplVar_      = "(ecalRecHitSumEtConeDR04 - compEcalIso + hcalTowerSumEtConeDR04 - compHcalIso)";
      tmplLeg_       = "Calo Isolation";
      tmplXTitle_    = "Calo Iso. (GeV)";
      IBCut_  = Form("sigmaIetaIeta > 0.003 && sigmaIetaIeta < 0.0105");
      SBCut_  = Form("sigmaIetaIeta > 0.011");
      nTmplBin_            = 40;
      tmplMin_             = -20;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                        
      tmplMax_             = 30;
      fitMin_              = -15;
      fitMax_              = 20;
   }    			       
   else if ( tmplOpt==3) {
      tmplName_        = "basicClusIso";
      tmplVar_         = "(cc4 + cr4)";
      tmplLeg_         = "Calo Isolation (cx)";
      tmplXTitle_      = "Calo Iso. (GeV)";
      IBCut_  = Form("sigmaIetaIeta > 0.003 && sigmaIetaIeta < 0.0105");
      SBCut_  = Form("sigmaIetaIeta > 0.011");
      nTmplBin_            = 40;
      tmplMin_             = -20;   // DON'T FORGET TO APPLY HOE CUT SEPARATELY                                                        
      tmplMax_             = 30;
      fitMin_              = -15;
      fitMax_              = 20;
   }
   
   
   
   else 
      cout << "  we dont have such template option " << endl;
}



TCut CutAndBinCollection::getIBCutWithPre() {
   //   float thehoe = getHoECut( percentBin, et) ;                                                                                      
   //  TCut thehoeCut = Form("hadronicOverEm< %f", (float)thehoe);                                                                       
   if ( hiOrPp_) 
      return univHoECutHI && IBCut_;
   else 
      return univHoECutPP && IBCut_;

}
TCut CutAndBinCollection::getSBCutWithPre(){
   //  float thehoe = getHoECut( percentBin, et) ;                                                                                       
   //   TCut thehoeCut = Form("hadronicOverEm< %f", (float)thehoe);                                                                      
   if ( hiOrPp_)
      return univHoECutHI && SBCut_;
   else
      return univHoECutPP && SBCut_;
}
TCut CutAndBinCollection::getIBCutWithPreWOEleVeto() { 
   if ( hiOrPp_)
      return univHoECutHI && IBCutWOEleVeto_;
   else
      return univHoECutPP && IBCutWOEleVeto_;

}


/*
float CutAndBinCollection::getHoECut(   float percentBin, float et ) {  // 100% centrality used   
float thehoe = 0.1;
if ( percentBin <= 20 ) {   // should change univHoECut  accordingly
if         ( et <=30 )                thehoe = 0.15;
else                                  thehoe = 0.1;
}
if ( percentBin > 20 ) {
if         ( et <=30 )                thehoe = 0.1;
else                                  thehoe = 0.05;
}

cout << endl << " centrality = " << percentBin << "%" << "    et = " << et << "     hoeCut = " << thehoe << endl;
return thehoe;
}

*/


void addCentralityFriend(TTree *tSig, TTree *tData,TCut selectionCut)
{
   //copied from /d101/yjlee/HIPhoton/ana/photonAna2011/common.h
   static int counter=0;
   TH1D *hSigCent = new TH1D("hSigCent","",40,-0.5,39.5);
   TH1D *hDataCent = new TH1D("hDataCent","",40,-0.5,39.5);

   hSigCent->SetLineColor(2);
   tSig->Project("hSigCent","cBin",selectionCut);
   tData->Project("hDataCent","cBin",selectionCut);
   hDataCent->Scale(1./hDataCent->GetEntries());
   hSigCent->Scale(1./hSigCent->GetEntries());
   hDataCent->Divide(hSigCent);
   TNtuple *nt = new TNtuple(Form("ntCentFriend%d",counter),"","cBinWeight");

   Int_t cBin;
   tSig->SetBranchAddress("cBin",&cBin);

   for (int i=0;i<tSig->GetEntries();i++)
      {
	 tSig->GetEntry(i);
	 int bin = hDataCent->FindBin(cBin);
	 //cout <<cBin<<" "<<hDataCent->GetBinContent(bin)<<endl;                                                                                                                                                    
	 nt->Fill(hDataCent->GetBinContent(bin));
      }
   counter++;
   delete hSigCent;
   delete hDataCent;
   tSig->AddFriend(nt);
}


void getNColl( float* ncoll) {

ncoll[0] = 1747.86 ; 
ncoll[1] = 1567.53 ; 
ncoll[2] = 1388.39 ; 
ncoll[3] = 1231.77 ; 
ncoll[4] = 1098.2 ; 
ncoll[5] = 980.439 ; 
ncoll[6] = 861.609 ; 
ncoll[7] = 766.042 ; 
ncoll[8] = 676.515 ; 
ncoll[9] = 593.473 ; 
ncoll[10] = 521.912 ; 
ncoll[11] = 456.542 ; 
ncoll[12] = 398.546 ; 
ncoll[13] = 346.647 ; 
ncoll[14] = 299.305 ; 
ncoll[15] = 258.344 ; 
ncoll[16] = 221.216 ; 
ncoll[17] = 188.677 ; 
ncoll[18] = 158.986 ; 
ncoll[19] = 134.7 ; 
ncoll[20] = 112.547 ; 
ncoll[21] = 93.4537 ; 
ncoll[22] = 77.9314 ; 
ncoll[23] = 63.5031 ; 
ncoll[24] = 52.0469 ; 
ncoll[25] = 42.3542 ; 
ncoll[26] = 33.9204 ; 
ncoll[27] = 27.3163 ; 
ncoll[28] = 21.8028 ; 
ncoll[29] = 17.2037 ; 
ncoll[30] = 13.5881 ; 
ncoll[31] = 10.6538 ; 
ncoll[32] = 8.35553 ; 
ncoll[33] = 6.40891 ; 
ncoll[34] = 5.13343 ; 
ncoll[35] = 3.73215 ; 
ncoll[36] = 3.06627 ; 
ncoll[37] = 2.41926 ; 
ncoll[38] = 2.11898 ; 
 ncoll[39] = 1.76953 ; 
 
}



fitResult doFit(TH1D* hSig=0, TH1D* hBkg=0, TH1D* hData1=0, double &nSig=a1, double &nSigErr=a2, float varLow=0.001, float varHigh=0.028, bool drawLeg=true, bool drawHist=false,double &chisq=a5,double &purity011=a6) {
   
   TH1D* hDatatmp = (TH1D*)hData1->Clone(Form("%s_datatmp",hData1->GetName()));
   double realNev = hDatatmp->GetEntries();
   int nBins = hDatatmp->GetNbinsX();
   histFunction2 *myFits = new histFunction2(hSig,hBkg);
   TF1 *f = new TF1("f",myFits,&histFunction2::evaluate,varLow,varHigh,2);
   f->SetParameters( hDatatmp->Integral(1,nBins+1), 0.3);
   f->SetParLimits(1,0,1);
   hDatatmp->Fit("f","LL M O Q","",varLow,varHigh);
   hDatatmp->Fit("f","LL M O Q","",varLow,varHigh);
   chisq = (double)f->GetChisquare()/ f->GetNDF()  ;
   
   //   cout <<" cs = " << chisq << endl;                                                                                                                                                                           
   fitResult res;
   res.nSig =0;
   double nev = f->GetParameter(0);
   double ratio = f->GetParameter(1);
   double ratioErr = f->GetParError(1);
   res.nSig    = nev * ratio;
   res.nSigErr = nev * ratioErr;

   TH1F *hSigPdf = (TH1F*)hSig->Clone(Form("%s_tmp",hSig->GetName()));
   hSigPdf->Scale(res.nSig/hSigPdf->Integral(1,nBins+1));

   //   cout << " all events(real) =   " <<  realNev << endl;                                                                                                                                                       
   //  cout << " all events(fitted) =   " <<  nev << endl;                                                                                                                                                          
   //  cout << " signals    =   " << nSig << endl;                                                                                                                                                                  

   TH1F *hBckPdf = (TH1F*)hBkg->Clone(Form("%s_tmp",hBkg->GetName()));
   hBckPdf->Scale((nev-res.nSig)/hBckPdf->Integral(1,nBins+1));

   double ss1 = hSigPdf->Integral(1, hSigPdf->FindBin(0.00999),"width");
   double bb1 = hBckPdf->Integral(1, hBckPdf->FindBin(0.00999),"width");
   //   cout <<"  hte bin = " <<hSigPdf->FindBin(0.00999) << endl;
   res.purity010 = ss1/(ss1+bb1);
   cout << "purity = " << res.purity010 << endl;
   hSigPdf->Add(hBckPdf);
   handsomeTH1(hSigPdf);
   mcStyle(hSigPdf);
   sbStyle(hBckPdf);
   cleverRange(hSigPdf,1.5);
   hSigPdf->SetNdivisions(510);

   hSigPdf->SetYTitle("Entries");
   hSigPdf->DrawCopy("hist");
   hBckPdf->DrawCopy("same hist");
   hData1->DrawCopy("same");
   TH1D* temphSigPdf = (TH1D*)hSigPdf->Clone("temp1");
   TH1D* temphBckPdf = (TH1D*)hBckPdf->Clone("temp2");
   if(drawLeg){
      TLegend *t3=new TLegend(0.5402006,0.5963235,0.9186019,0.7853466,NULL,"brNDC");
      t3->AddEntry(hData1,"Pb+Pb  #sqrt{s}_{_{NN}}=2.76 TeV","pl");
      t3->AddEntry(temphSigPdf,"Signal","lf");
      t3->AddEntry(temphBckPdf,"Background","lf");
      t3->SetFillColor(0);
      t3->SetBorderSize(0);
      t3->SetFillStyle(0);
      t3->SetTextFont(63);
      t3->SetTextSize(15);
      t3->Draw();
      drawCMS2011(0.53,0.9,150,16);
   }
   

   //   delete hSigPdf;                                                                                                                                                                                             
   //   delete hBckPdf;                                                                                                                                                                                             

   return res;

}



#endif
