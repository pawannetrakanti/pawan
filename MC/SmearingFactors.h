#ifndef SMEARINGFACTORS_H
#define SMEARINGFACTORS_H

#include <iostream>
#include <utility>
#include <TROOT.h>
#include <TStyle.h>

#include <TString.h>
#include <TRandom3.h>
#include <TF1.h>

float GetSmearedPtMC(int /*nj*/,int /*ic*/,float /*recopt*/,float /*refpt*/);
//float GetSmearedPtData(int /*nj*/,int /*ic*/,float /*recopt*/);
float GetSmearedPtData(int /*nj*/,int /*ic*/,float /*recopt*/,float /*fpercent*/,const char */*csys*/);
void  LoadParameters();
int   GetCBin(int /*bin*/);
float AfterBurnMean(int /*nj*/,int /*ic*/,float /*smpt*/,float /*refpt*/);
float GetPtBin(float /*smpt*/);
const int NCEN=7;
const int KNJ =7;

//! Smearing function
TF1 *fsm[KNJ][NCEN], *fmd[KNJ][NCEN];
//!                             0-5%    5-10%  10-30%  30-50%   50-70%  70-90%   pp 
double smearf[KNJ][NCEN][3]={ {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2PF
                              {{16.3179,-0.0862868,0.000154046}, //! 0-5%
                               {15.1583,-0.073777,8.97108e-05},  //! 5-10%
                               {13.6676,-0.0717466,0.000118877}, //! 10-30%
                               {9.33318,-0.0541488,0.000102715}, //! 30-50%
                               {10.0646,-0.0982037,0.000269721}, //! 50-70%
                               {5.6704,-0.0423724,9.28418e-05},  //! 70-90%
                               {0,0,0}}, //! ak3PF 
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak4PF
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu2PF
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu3PF
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}  //! akPu4PF
};

//! for low pT <90 Gev/c 30-50,50-70 and 70-90                                                                                                                   
double lptsmf[KNJ][NCEN][3]={ {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo                                                           
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2PF                                                               
			      {{6.72286,8.58926,6.89564}, //! 0-5%                                                                                                
                               {6.59894,8.81884,8.41884}, //! 5-10%                                                                                               
                               {6.09014,7.65841,7.17648}, //! 10-30%                                                                                              
                               {5.14246,5.28939,5.316377},//! 30-50%                                                                                             
                               {4.55183,4.29964,3.00092}, //! 50-70%                                                                                             
                               {4.41342,4.02342,1.11234}, //! 70-90%                                                                                              
                               {0,0,0}}, //! ak3PF                                                                                                               
	
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak4PF                                                               
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu2PF                                                             
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu3PF                                                             
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}} //! akPu4PF                                                              
};

double afsmf[KNJ][NCEN][4] = {{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}, //! icpu5calo 
			      {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}, //! ak2PF 
			      {{4.3297,7.94448,5.63229,0.977756}, //! 0-5%
			       {3.08093,3.77259,3.80761,1.92774}, //! 5-10%
			       {3.17783,3.04095,3.75831,0}, //! 10-30%
			       {4.54152,4.49648,2.05938,0}, //! 30-50%
			       {4.06814,3.18181,3.29295,0}, //! 50-70%
			       {1.85248,2.71036,0,0}, //! 70-90%
			       {0,0,0,0}}, //! pp //! ak3PF 
			      {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}, //! ak4PF 
			      {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}, //! akPu2PF 
			      {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}, //! akPu3PF 
			      {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}, //! akPu4PF 
};

double mdiff[KNJ][NCEN][3]={ {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo                                                            
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2PF                                                                
                             {{0.00268504,0.00016407,-4.40974e-07},
                              {0.0159945,4.53507e-05,-1.95167e-07},
                              {0.0351176,-9.87183e-05,7.58988e-08},
                              {0.0430083,-0.000215562,3.35325e-07},
                              {0.0499618,-0.000235293,3.07907e-07},
                              {0.0387144,-0.000228153,4.00751e-07},
                              {0,0,0}}, //! ak3PF                                                                                                                
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak4PF                                                                
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu2PF                                                              
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu3PF                                                              
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}  //! akPu4PF                                                              
};

//! after burner for mean diff
double amdiff[KNJ][NCEN][10]={ {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! icpu5calo
			       {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! ak2PF
			       {{0.0190468,-0.0109088,-0.00455594,-0.00919259,-0.00619966,-0.00380635,0.00996578,0.00677246,0.00889736,0.00211281},
				{0.0164925,-0.0207178,-0.00572103,-0.000149071,0.00611407,0.0017544,-0.000367522,0.00357312,0.00718009,0.00362355},
				{0.0115435,-0.00132334,0.0056141,0.00313061,-0.0107589,0.00898504,0.00247216,0.00880492,0.00123292,0.00262469},
				{0.00472796,-0.00598991,0.00342047,0.00623858,0.00498122,0.00545985,-0.000427663,0.00338912,0.00128466,0.00184447},
				{-0.00560087,0.00213766,0.0134992,0.00633717,0.00128788,0.00414634,0.00403112,0.00647879,0.00192851,0.00288677},
				{-0.00115108,0.00171649,0.0130042,0.00296229,0.00481486,0.00578582,0.00068301,0.000454724,-0.00379604,0.00160635},
				{0,0,0,0,0,0,0,0,0,0}}, //! ak3PF                                                                                                                
			       {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! ak4PF
			       {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! akPu2PF
			       {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! akPu3PF
			       {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}  //! akPu4PF
};

double PT_BINS[11] = {30,50,70,90,120,160,200,240,280,340,400};


void LoadParameters()
{
  for(int nj=0;nj<KNJ;nj++){
    for(int i=0;i<NCEN;i++){
      fsm[nj][i] = new TF1(Form("fsm%d_%d",nj,i),"pol2",20,800);
      fmd[nj][i] = new TF1(Form("fmd%d_%d",nj,i),"pol2",20,800);
      
      for(int im=0;im<3;im++){
	fsm[nj][i]->SetParameter(im,smearf[nj][i][im]);
	fmd[nj][i]->SetParameter(im,mdiff [nj][i][im]);
      }
    }
  }
}

float GetSmearedPtMC(int nj,int ic,float recopt,float refpt)
{
  //! Get Centrality bin
  //int icen = GetCBin(ic);
  int icen=ic;

  //! Smearing currently doing from 30 GeV onwards
  float smpt = recopt - fmd[nj][icen]->Eval(refpt)*recopt;
  smpt = AfterBurnMean(nj,ic,smpt,refpt);
  int ib=-1;
  if(smpt>=90){
    float tsmf = fsm[nj][icen]->Eval(refpt);
    if(refpt>120 && refpt<=240){
      ib=-1;
      ib = (int)(refpt - 120.)/40.;
      if(ib<0)tsmf+=0;
      else tsmf += afsmf[nj][icen][ib];
    }
    smpt = gRandom->Gaus(smpt,tsmf);
  }else{
    ib=-1;
    ib = (int)(recopt - 30.)/20.;
    if(ib<0)return recopt;
    else smpt = gRandom->Gaus(smpt,lptsmf[nj][icen][ib]);
  }
  return smpt;
}

float GetSmearedPtData(int nj,int ic,float recopt,float fpercent,const char *csys)
{
  //! Get Centrality bin
  //int icen = GetCBin(ic);
  int icen=ic;
  float smpt = recopt;
  int ib=-1;
  float tsmf = fsm[nj][icen]->Eval(recopt);

  if(strcmp(csys,"low")==0){
    tsmf  = fsm[nj][icen]->Eval(recopt) - (fpercent/100.)*fsm[nj][icen]->Eval(recopt);
    if(recopt<=90){
      ib = (int)(recopt - 30.)/20.;
      if(ib<0)return smpt;
      tsmf  = lptsmf[nj][icen][ib] - (fpercent/100.)*lptsmf[nj][icen][ib];      
    }else{
      smpt = recopt - (fmd[nj][icen]->Eval(recopt) + (fpercent/100.)*fmd[nj][icen]->Eval(recopt))*recopt;
      smpt = AfterBurnMean(nj,ic,smpt,recopt);
      if(recopt>120 && recopt<=240){
	ib=-1;
	ib = (int)(recopt - 120.)/40.;
	if(ib<0)tsmf+=0;
	else tsmf += afsmf[nj][icen][ib];
      }
    }
    smpt   = gRandom->Gaus(smpt,tsmf);
  }
  else if(strcmp(csys,"up")==0){
    tsmf  = fsm[nj][icen]->Eval(recopt) + (fpercent/100.)*fsm[nj][icen]->Eval(recopt);
    if(recopt<=90){
      ib = (int)(recopt - 30.)/20.;
      if(ib<0)return smpt;
      tsmf  = lptsmf[nj][icen][ib] + (fpercent/100.)*lptsmf[nj][icen][ib];      
    }else{
      smpt = recopt - (fmd[nj][icen]->Eval(recopt) - (fpercent/100.)*fmd[nj][icen]->Eval(recopt))*recopt;
      smpt = AfterBurnMean(nj,ic,smpt,recopt);
      if(recopt>120 && recopt<=240){
	ib=-1;
	ib = (int)(recopt - 120.)/40.;
	if(ib<0)tsmf+=0;
	else tsmf += afsmf[nj][icen][ib];
      }
    }
    smpt = gRandom->Gaus(smpt,tsmf);

  }
  else{
    //! Smearing currently doing from 30 GeV onwards 
    smpt = recopt - fmd[nj][icen]->Eval(recopt)*recopt;
    smpt = AfterBurnMean(nj,ic,smpt,recopt);

    if(recopt>90){
      tsmf = fsm[nj][icen]->Eval(recopt);
      if(recopt>120 && recopt<=240){
	ib=-1;
	ib = (int)(recopt - 120.)/40.;
	if(ib<0)tsmf+=0;
	else tsmf += afsmf[nj][icen][ib];
      }
    }else{
      ib = (int)(recopt - 30.)/20.;
      if(ib<0)return smpt;
      else tsmf = lptsmf[nj][icen][ib];
    }
    smpt = gRandom->Gaus(smpt,tsmf);
  }
  return smpt;
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
  int ib = GetPtBin(refpt);
  if(ib<0)return smpt; //! do not shift this
  else{
    smpt += smpt*amdiff[nj][ic][ib];
  }
  return smpt;
}
float GetPtBin(float smpt)
{
  for(int i=0;i<10;i++){
    if(smpt>=PT_BINS[i] && smpt<PT_BINS[i+1])return i;
  }
  return -1;
}
#endif
