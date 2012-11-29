#ifndef SMEARINGFACTORS_H
#define SMEARINGFACTORS_H

#include <iostream>
#include <utility>
#include <TROOT.h>
#include <TStyle.h>

#include <TString.h>
#include <TRandom3.h>
#include <TF1.h>



void  LoadParameters();
int   GetCBin(int /*bin*/);

float GetPtBin(float /*smpt*/);

float GetSmFactor(int /*nj*/,int /*ic*/,float /*recopt*/);
float GetMeanShift(int /*nj*/,int /*ic*/,float /*recopt*/);
float AfterBurnMean(int /*nj*/,int /*ic*/,float /*smpt*/,float /*refpt*/);

float GetSmearedPtMC(int /*nj*/,int /*ic*/,float /*recopt*/,float /*refpt*/);
float GetSmearedPtMC_woAfBurn(int /*nj*/,int /*ic*/,float /*recopt*/,float /*refp*/);
float GetSmearedPtMC_NoMeanShift(int /*nj*/,int /*ic*/,float /*recopt*/,float /*refpt*/); //! JFF & JS
float GetSmearedPtMC_OnlyMeanShift(int /*nj*/,int /*ic*/,float /*recopt*/,float /*refpt*/);

float GetSmearedPtData(int /*nj*/,int /*ic*/,float /*recopt*/,float /*fpercent*/,const char */*csys*/);
float GetSmearedPtData_woAfBurn(int /*nj*/,int /*ic*/,float /*recopt*/,float /*fpercent*/,const char */*csys*/);
float GetSmearedPtData_NoMeanShift(int /*nj*/,int /*ic*/,float /*recopt*/,float /*fpercent*/,const char */*csys*/); //! JFF & JS
float GetSmearedPtData_OnlyMeanShift(int /*nj*/,int /*ic*/,float /*recopt*/,float /*fpercent*/,const char */*csys*/);

//! For PbPb jet energy scale correction for Pu algorithms
float GetPbPbCorrectedScaleMC  (int /*nj*/,int /*hiBin*/,float /*recopt*/,float /*refpt*/); //! JFF & JS
float GetPbPbCorrectedScaleData(int /*nj*/,int /*hiBin*/,float /*recopt*/);  //! JFF & JS

//! Re-weighting factor for data
float GetReWeight(int /*nj*/,int /*ic*/, float /*smpt*/); //! JFF & JS
float GetReWeight_NoMeanShift(int /*nj*/,int /*ic*/, float /*smpt*/); //! JFF & JS

const int NCEN=7;
const int KNJ =7;

//! Smearing function
TF1 *fresol[KNJ][NCEN], *fscale[KNJ][NCEN];
TF1 *fasmf [KNJ][NCEN];
TF1 *fReWe; //! only for 0-10% and 50-100% and ak3PF pp jets

//!                             0-5%    5-10%  10-30%  30-50%   50-70%  70-90%   pp 
double resol[KNJ][NCEN][3]={ {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo
			     {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2PF
			     {
			       {0.0296102 ,1.51239 ,6.37142},//! 0-5% ak3PF
			       {0.021402  ,1.51753 ,5.13989}, //! 5-10%  
			       {0.042533  ,1.38629 ,4.81684}, //! 10-30%  
			       {0.0395597 ,1.33222 ,3.03127}, //! 30-50%  
			       {0.0520218 ,1.17195 ,4.46752}, //! 50-70%  
			       {0.0602054 ,1.03093 ,5.09772}, //! 70-90%	  
			       {0.055736  ,1.18204 ,2.26452}},//!  pp ak3PF

			     {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak4PF
			     {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu2PF
			     {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu3PF
			     {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}  //! akPu4PF
};


double scale[KNJ][NCEN][3]={ {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2PF
                             {
			       {1.009   ,-5.03813 ,191.723}, //! 0-5% 
			       {1.00859 ,-5.26335 ,193.504}, //! 5-10%  
			       {1.00378 ,-3.84755 ,102.459}, //! 10-30%  
			       {1.00489 ,-3.87293 ,106.204}, //! 30-50%  
			       {1.0032  ,-3.05924 ,74.5900}, //! 50-70%  
			       {1.00091 ,-1.46321 ,24.0954}, //! 70-90%   
			       {1.00175 ,-0.60521 ,26.2456}}, //! pp ak3PF	  
			      
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak4PF
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu2PF
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu3PF
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}  //! akPu4PF
};



double afsmf[KNJ][NCEN][3] = {{{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2PF
                              {{8.58808e+00,-1.03265e+03,6.01265e+03}, //! 0-5%
			       {6.88475e+00,-5.91181e+02,-1.05710e+04},//! 5-10%
			       {6.58010e+00,-6.82048e+02,-1.45276e+03},//! 10-30%
			       {5.58927e+00,-6.70886e+02, 8.61490e+03},//! 30-50%
			       {5.05257e+00,-4.52675e+02,-3.14922e+03},//! 50-70% 
			       {2.77350e-02,5.17173e+02 ,-3.51948e+04},//! 70-90%
			       {0,0,0}},//! pp //! ak3PF                                                                                  
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak4PF               
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu2PF             
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu3PF             
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu4PF             
};

//! after burner for mean diff                       
double amdiff[KNJ][NCEN][33]={ {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! icpu5calo
                               {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! ak2PF
                               {
				 {-0.0957503,-0.0366432,-0.00239193,-0.00111985,0.00406867,0.0129364,-0.000904024,0.00290388,0.000446916,0.014991,0.00441539,0.0102446,7.05719e-05,-0.0119217,-0.00174958,-0.000980139,0.00811815,0.00433797,0.00240582,-0.00321567,-0.00202787,0.00674027,0.000503898,3.97563e-05,0.00317174,0.000580072,0.00688094,0.00317574,-0.000688612,-0.00319052,0.0101621,0.00435334,0.00570643},//! 0-5%
				 {-0.105525,-0.0316963,-0.00745112,-0.00189686,0.00352544,0.00309801,-0.00600392,-0.0121201,0.00319189,0.0059008,0.000226438,0.00754851,0.00883555,-0.00923926,0.0188677,-0.00176173,0.00513905,0.00692612,0.00295156,-0.00486153,-0.00638604,-0.0018059,0.00472921,-0.00532085,-0.000431478,0.00263381,0.00789589,0.000686705,0.00072521,0.00738513,0.00688565,-0.00357831,0.00438887},//! 5-10%
				 {-0.0516163,-0.00960845,-0.00638938,0.00395107,0.00338864,0.00748146,0.00172335,0.000775218,-0.0170699,-0.00944746,-0.00334263,-0.000624597,0.00544232,-0.000833631,0.0170863,0.00166899,0.0073728,0.0096156,0.00508732,-0.0119,-0.00246751,0.0054419,0.00258797,-0.00117266,-0.000197947,-0.00679743,0.00264871,0.00251925,-0.000521064,-0.00134331,0.0041877,-0.000174522,0.004668},//! 10-30%
				 {-0.0296276,-0.0102345,-0.00182575,0.00346065,-0.000382781,0.00336576,0.00171405,0.000446558,-0.00463557,-0.00347996,-0.000578165,0.000472009,0.00233442,-0.00186092,0.01033,0.00327784,0.00852287,0.00558299,0.00221169,-0.00364447,-0.00356013,0.000103354,0.00445223,0.000342906,-0.00261116,0.00263554,0.0114818,0.00180966,-0.00305378,-0.0103016,0.00147146,0.00506431,0.00155097},//! 30-50%
				 {-0.031548,-0.00918972,-0.00248939,-0.000692368,0.000444174,0.00123191,0.0023036,0.000940561,-0.00794965,-0.00277925,-0.00262636,0.000455439,0.00198525,-0.00150663,0.0165374,0.000724912,0.00977331,0.00610918,-0.00144237,-0.00246561,-0.00128436,0.00214124,0.00414962,0.0038625,1.8537e-05,-0.000538766,0.00695884,0.00556296,-0.00397176,-0.00839287,0.00705308,0.000589132,0.00280684},//! 50-70%
				 {-0.018268,-0.00804049,0.000407219,0.00383049,0.00274009,-0.00042206,-0.000943899,-0.00426054,-0.00859624,-0.00543714,-0.00393707,0.000441194,0.00764209,-0.00200003,0.0171912,0.00492692,0.00745159,0.00339669,0.00219131,-0.0038476,0.00133491,0.00189805,0.00218612,0.00452137,0.00087893,-0.00327235,0.00359207,0.00202942,0.00302291,-0.00834006,0.00341874,-0.00160706,-0.00174284},//! 70-90%
				 {0,0,0,0,0,0,0,0,0,0}}, //!pp //! ak3PF
                               {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! ak4PF
                               {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! akPu2PF
                               {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! akPu3PF
                               {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}} //! akPu4PF
};


//! Reweighting factor for 0-10% and 50-100% only for data
double fwe010[2][3] = {{2.32745e+00, -9.63227e+01, 5.88044e+03}, //! 0-5%                                                                  
		       {1.13736e+00,  2.20014e+02,-1.28424e+04}  //! 5-10%                                                                 
};
double fwe50100[2][3] = {{2.30339e+00,  -8.18276e+01,  5.35947e+03}, //! 50-70%                                                            
			 {1.50452e+00,   1.20980e+02, -7.33869e+03}  //! 70-100%                                                           
};
//! Reweighting factor for 0-10% and 50-100% only for data - No Mean Shift for pp
double fwe010_noms[2][3] = {{1.94440e+00,  9.70262e+00, -3.85305e+02}, //! 0-5%
			    {1.72225e+00,  6.63995e+01, -3.86893e+03}  //! 5-10%
};
double fwe50100_noms[2][3] = {{2.05212e+00, -1.58706e+01,  9.43296e+02}, //! 50-70%
			      {1.94988e+00,  1.51774e+01, -8.94539e+02}  //! 70-100%
};



double PT_BINS[34] = {30,40,50,60,70,80,90,100,
		      110,120,130,140,150,160,170,180,190,200,
		      210,220,230,240,250,260,270,280,290,300,
		      310,320,330,340,350,400};
const int NBINS=34;

void LoadParameters()
{
  for(int nj=0;nj<KNJ;nj++){
    for(int i=0;i<NCEN;i++){
      fresol[nj][i]  = new TF1(Form("fresol%d_%d",nj,i),
			      "sqrt(pow([0],2)+pow([1]/sqrt(x),2)+pow([2]/x,2))",30,400);
      fscale[nj][i] = new TF1(Form("fscale%d_%d",nj,i),"[0] + [1]/x + [2]/x/x"  ,30,400);
      fasmf [nj][i] = new TF1(Form("fasmfe%d_%d",nj,i),"[0] + [1]/x + [2]/x/x"  ,30,400);      
      for(int im=0;im<3;im++){
	fresol[nj][i]->SetParameter(im,resol[nj][i][im]);
	fscale[nj][i]->SetParameter(im,scale[nj][i][im]);
	fasmf[nj][i]->SetParameter (im,afsmf[nj][i][im]);
      }
    }
  }

  //! Reweight factor for 0-10% and 50-100%
  fReWe = new TF1("fReWe","[0] + [1]/x + [2]/x/x",30,400);
}



//////////////////////////////////////////////////MC related functions///////////////////////////////////////////////////////////////
float GetSmearedPtMC(int nj,int ic,float recopt,float refpt)
{
  int icen = ic;

  //! Get the jet energy scale
  float mpp   = fscale[nj][NCEN-1]->Eval(refpt);
  float mpbpb = fscale[nj][icen]->Eval(refpt);

  //! Calculate the shift in scale
  float mdf = mpp - mpbpb;

  //! Now shift scale first
  float smpt = recopt - mdf*recopt;

  //! afterburn to adjust the remaining residual
  smpt = AfterBurnMean(nj,icen,smpt,refpt);
  
  //! Get resolutions
  float rpp   = fresol[nj][NCEN-1]->Eval(refpt); //! pp
  float rpbpb = fresol[nj][icen]->Eval(refpt);     //! PbPb
  
  //! Calculate the smearing factor
  float smf = 0;
  if(rpp>rpbpb){
    smf=0;
  }else{
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*100;
    float af = fasmf[nj][icen]->Eval(refpt); 
    smf += af;
  }
  //! Smearing the recopt
  smpt = gRandom->Gaus(smpt,smf);
  return smpt;
}
float GetSmearedPtMC_woAfBurn(int nj,int ic,float recopt,float refpt)
{
  int icen = ic;

  //! Get the jet energy scale
  float mpp   = fscale[nj][NCEN-1]->Eval(refpt);
  float mpbpb = fscale[nj][icen]->Eval(refpt);

  //! Calculate the shift in scale
  float mdf = mpp - mpbpb;

  //! Now shift scale first
  float smpt = recopt - mdf*recopt;

  //! afterburn to adjust the remaining residual
  //smpt = AfterBurnMean(nj,icen,smpt,refpt);
  
  //! Get resolutions
  float rpp   = fresol[nj][NCEN-1]->Eval(refpt); //! pp
  float rpbpb = fresol[nj][icen]->Eval(refpt);     //! PbPb
  
  //! Calculate the smearing factor
  float smf = 0;
  if(rpp>rpbpb){
    smf=0;
  }else{
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*100;
    //float af = fasmf[nj][icen]->Eval(refpt); 
    //smf += af;
  }
  //! Smearing the recopt
  smpt = gRandom->Gaus(smpt,smf);
  return smpt;
}
float GetPbPbCorrectedScaleMC(int nj,int hibin,float recopt,float refpt)
{
  if(nj!=2)return recopt; //! currently only for akPu3PF jets

  int icen = GetCBin(hibin);
  //! Get the jet energy scale
  float mpbpb = fscale[nj][icen]->Eval(refpt);

  //! Now correct scale first
  float corr_recopt = recopt/mpbpb;

  return corr_recopt;
}
float GetSmearedPtMC_NoMeanShift(int nj,int ic,float recopt,float refpt)
{
  int icen = ic;
  float smpt = recopt;
  //! Get resolutions
  float rpp   = fresol[nj][NCEN-1]->Eval(refpt); //! pp
  float rpbpb = fresol[nj][icen]->Eval(refpt);   //! PbPb
  
  //! Calculate the smearing factor
  float smf = 0;
  if(rpp>rpbpb){
    smf=0;
  }else{
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*100;
    float af = fasmf[nj][icen]->Eval(refpt); 
    smf += af;
  }
  //! Smearing the recopt
  smpt = gRandom->Gaus(smpt,smf);
  return smpt;
}
float GetSmearedPtMC_OnlyMeanShift(int nj,int ic,float recopt,float refpt)
{
  int icen = ic;
  //! Get the jet energy scale
  float mpp   = fscale[nj][NCEN-1]->Eval(refpt);
  float mpbpb = fscale[nj][icen]->Eval(refpt);

  //! Calculate the shift in scale
  float mdf = mpp - mpbpb;

  //! Now shift scale first
  float smpt = recopt - mdf*recopt;
  return smpt;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////// Data related functions ////////////////////////////////////////////
float GetSmearedPtData(int nj,int ic,float recopt,float fpercent,const char *csys)
{
  int icen = ic;
  //! Mean shift
  float mpp   = fscale[nj][NCEN-1]->Eval(recopt);
  float mpbpb = fscale[nj][icen]->Eval(recopt);
  float mdf   = mpp - mpbpb;

  //! Get resolutions
  float rpp   = fresol[nj][NCEN-1]->Eval(recopt); //! pp
  float rpbpb = fresol[nj][icen]->Eval(recopt);   //! PbPb

  float smf=0;
  float smpt=0;

  //! Calculate the smearing factor
  if(rpp>rpbpb){
    smf=0;
  }else{
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*100;
    float af = fasmf[nj][icen]->Eval(recopt); 
    smf += af;
  }

  if(strcmp(csys,"low")==0){

    smpt = recopt - ((mdf - (2.*fpercent/100.)*mdf)*recopt);
    smpt = AfterBurnMean(nj,icen,smpt,recopt);
    smf  = smf - (fpercent/100.)*smf;

  }else if(strcmp(csys,"up")==0){

    smpt = recopt - ((mdf + (2.*fpercent/100.)*mdf)*recopt);
    smpt = AfterBurnMean(nj,icen,smpt,recopt);
    smf  = smf + (fpercent/100.)*smf;

  }else{

    //! Now shift scale first
    smpt = recopt - mdf*recopt;
    
    //! afterburn to adjust the remaining residual in mean shift
    smpt = AfterBurnMean(nj,icen,smpt,recopt);
  }

  //! Smearing the recopt
  smpt = gRandom->Gaus(smpt,smf);  
  return smpt;
}
float GetSmearedPtData_woAfBurn(int nj,int ic,float recopt,float fpercent,const char *csys)
{
  int icen = ic;
  //! Mean shift
  float mpp   = fscale[nj][NCEN-1]->Eval(recopt);
  float mpbpb = fscale[nj][icen]->Eval(recopt);
  float mdf   = mpp - mpbpb;

  //! Get resolutions
  float rpp   = fresol[nj][NCEN-1]->Eval(recopt); //! pp
  float rpbpb = fresol[nj][icen]->Eval(recopt);   //! PbPb

  float smf=0;
  float smpt=0;

  //! Calculate the smearing factor
  if(rpp>rpbpb){
    smf=0;
  }else{
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*100;
    //float af = fasmf[nj][icen]->Eval(recopt); 
    //smf += af;
  }

  if(strcmp(csys,"low")==0){

    smpt = recopt - ((mdf - (2.*fpercent/100.)*mdf)*recopt);
    //smpt = AfterBurnMean(nj,icen,smpt,recopt);
    smf  = smf - (fpercent/100.)*smf;

  }else if(strcmp(csys,"up")==0){

    smpt = recopt - ((mdf + (2.*fpercent/100.)*mdf)*recopt);
    //smpt = AfterBurnMean(nj,icen,smpt,recopt);
    smf  = smf + (fpercent/100.)*smf;

  }else{

    //! Now shift scale first
    smpt = recopt - mdf*recopt;
    
    //! afterburn to adjust the remaining residual
    //smpt = AfterBurnMean(nj,icen,smpt,recopt);
  }

  //! Smearing the recopt
  smpt = gRandom->Gaus(smpt,smf);  
  return smpt;
}
float GetSmearedPtData_NoMeanShift(int nj,int ic,float recopt,float fpercent,const char *csys)
{
  int icen = ic;
  float smpt = recopt;
  //! Get resolutions
  float rpp   = fresol[nj][NCEN-1]->Eval(recopt); //! pp
  float rpbpb = fresol[nj][icen]->Eval(recopt);     //! PbPb
  
  //! Calculate the smearing factor
  float smf = 0;
  if(rpp>rpbpb){
    smf=0;
  }else{
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*100;
    float af = fasmf[nj][icen]->Eval(recopt); 
    smf += af;
  }

  if(strcmp(csys,"low")==0){
    smf  = smf - (fpercent/100.)*smf;
  }else if(strcmp(csys,"up")==0){
    smf  = smf + (fpercent/100.)*smf;
  }

  //! Smearing the recopt
  smpt = gRandom->Gaus(smpt,smf);
  return smpt;
}
float GetSmearedPtData_OnlyMeanShift(int nj,int ic,float recopt,float fpercent,const char *csys)
{
  int icen = ic;
  //! Get the jet energy scale
  float mpp   = fscale[nj][NCEN-1]->Eval(recopt);
  float mpbpb = fscale[nj][icen]->Eval(recopt);

  //! Calculate the shift in scale
  float mdf = mpp - mpbpb;
  float smpt=recopt;

  if(strcmp(csys,"low")==0){
    smpt = recopt - ((mdf - (2.*fpercent/100.)*mdf)*recopt);
    smpt = AfterBurnMean(nj,icen,smpt,recopt);

  }else if(strcmp(csys,"up")==0){
    smpt = recopt - ((mdf + (2.*fpercent/100.)*mdf)*recopt);
    smpt = AfterBurnMean(nj,icen,smpt,recopt);

  }else{
    smpt = recopt - mdf*recopt;
    smpt = AfterBurnMean(nj,icen,smpt,recopt);
  }
  return smpt;
}
float GetPbPbCorrectedScaleData(int nj,int hibin,float recopt)
{
  if(nj!=2)return recopt;

  int icen = GetCBin(hibin);
  //! Get the jet energy scale
  float mpbpb = fscale[nj][icen]->Eval(recopt);

  //! Now correct scale first
  float corr_recopt = recopt/mpbpb;
  //cout<<"\t \t GetPbPbCScale : "<<recopt<<"\t corr_recopt : "<<corr_recopt<<"\t mpbpb : "<<mpbpb<<endl;
  return corr_recopt;
}
float GetReWeight(int nj,int ic,float smpt)
{
  float rewe=1;
  if(nj!=2)return rewe;

  if(ic==0 || ic==1){ //! 0-10%
    fReWe->SetParameters(fwe010[ic][0],fwe010[ic][1],fwe010[ic][2]);
    rewe = 1./fReWe->Eval(smpt);
  }
  else if(ic==4 || ic==5){//! 50-100%
    fReWe->SetParameters(fwe50100[ic-4][0],fwe50100[ic-4][1],fwe50100[ic-4][2]);
    rewe = 1./fReWe->Eval(smpt);
  }
  return rewe;
}
float GetReWeight_NoMeanShift(int nj,int ic,float smpt)
{
  float rewe=1;
  if(nj!=2)return rewe;

  if(ic==0 || ic==1){ //! 0-10%
    fReWe->SetParameters(fwe010_noms[ic][0],fwe010_noms[ic][1],fwe010_noms[ic][2]);
    rewe = 1./fReWe->Eval(smpt);
  }
  else if(ic==4 || ic==5){//! 50-100%
    fReWe->SetParameters(fwe50100_noms[ic-4][0],fwe50100_noms[ic-4][1],fwe50100_noms[ic-4][2]);
    rewe = 1./fReWe->Eval(smpt);
  }
  return rewe;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







////////////////////////////////////////////////Common functions////////////////////////////////////////////////////////////
float GetSmFactor(int nj,int ic,float recopt)
{
  int icen = ic;
  //! Get resolutions
  float rpp   = fresol[nj][NCEN-1]->Eval(recopt); //! pp
  float rpbpb = fresol[nj][icen]->Eval(recopt);     //! PbPb

  float smf=0;
  //! Calculate the smearing factor
  if(rpp>rpbpb){
    smf=0;
  }else{
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*100;
    float af = fasmf[nj][ic]->Eval(recopt); 
    smf += af;
  }
  return smf;
}
float GetMeanShift(int nj,int ic,float recopt)
{
  int icen = ic;
  //! Mean shift
  float mpp   = fscale[nj][NCEN-1]->Eval(recopt);
  float mpbpb = fscale[nj][icen]->Eval(recopt);
  float mdf   = mpp - mpbpb;

  return mdf;
}
int GetCBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 2.5% bins of cross section
  //! in 0-39 bins
  if(bin<2)ibin=0;                  //! 0-5%
  else if(bin>=2 && bin<4)ibin=1;   //! 5-10%
  else if(bin>=4 && bin<12)ibin=2;  //! 10-30%
  else if(bin>=12&& bin<20)ibin=3;  //! 30-50%
  else if(bin>=20&& bin<28)ibin=4;  //! 50-70%
  else if(bin>=28&& bin<36)ibin=5;  //! 70-90%
  else ibin=5; //! 90-100% use the same 70-90% factors
  return ibin;
}
float AfterBurnMean(int nj,int ic,float smpt,float refpt)
{
  int icen = ic;
  int ib = GetPtBin(refpt);
  if(smpt<50)return smpt; //! do not shift this
  else{
    smpt += smpt*amdiff[nj][icen][ib];
  }
  return smpt;
}
float GetPtBin(float pt)
{
  for(int i=0;i<NBINS-1;i++){
    if(pt>=PT_BINS[i] && pt<PT_BINS[i+1])return i;
  }
  return -1;
}
#endif
