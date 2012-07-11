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
float GetSmearedPtData(int /*nj*/,int /*ic*/,float /*recopt*/);
void LoadParameters();
int GetCBin(int /*bin*/);

const int NCEN=7;
const int KNJ =7;

//! Smearing function                                                                                                                                          
TF1 *fsm[KNJ][NCEN], *fmd[KNJ][NCEN];
//!                             0-5%    5-10%  10-30%  30-50%   50-70%  70-90%   pp 
double smearf[KNJ][NCEN][3]={ {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo                                                           
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2PF                                                               
                              {{19.4669,-0.0881832,0.000121327},
                               {17.8516,-0.0826758,0.000116252},
                               {16.0364,-0.0725557,0.000106632},
                               {12.0566,-0.059805,9.10766e-05},
                               {11.3832,-0.0680122,0.000114237},
                               {5.6704,-0.0423724,9.28418e-05},
                               {0,0,0}}, //! ak3PF                                                                                                               
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak4PF                                                               
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu2PF                                                             
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu3PF                                                             
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}} //! akPu4PF                                                              
};

//! for low pT <90 Gev/c 30-50,50-70 and 70-90                                                                                                                   
double lptsmf[KNJ][NCEN][3]={ {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo                                                           
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2PF                                                               
                              {{9.12286,8.54926,7.39564},//! 0-5%                                                                                                
                               {8.70894,7.61884,8.21884},//! 5-10%                                                                                               
                               {7.99014,6.75841,5.47648},//! 10-30%                                                                                              
                               {6.44246,5.28939,5.086377},//! 30-50%                                                                                             
                               {6.15183,4.59964,4.390092},//! 50-70%                                                                                             
                               {5.21342,4.52342,1.11234},//! 70-90%                                                                                              
                               {0,0,0}}, //! ak3PF                                                                                                               
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak4PF                                                               
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu2PF                                                             
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! akPu3PF                                                             
                              {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}} //! akPu4PF                                                              
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
  if(smpt>90)smpt = gRandom->Gaus(smpt,fsm[nj][icen]->Eval(refpt));
  else{
    int ib=-1;
    ib = (int)(smpt - 30.)/20.;
    if(ib<0)return recopt;
    else smpt = gRandom->Gaus(smpt,lptsmf[nj][icen][ib]);
  }
  return smpt;
}

float GetSmearedPtData(int nj,int ic,float recopt)
{
  //! Get Centrality bin
  //int icen = GetCBin(ic);
  int icen=ic;
  
  //! Smearing currently doing from 30 GeV onwards                                                                                                       
  float smpt = recopt - fmd[nj][icen]->Eval(recopt)*recopt;
  if(smpt>90)smpt = gRandom->Gaus(smpt,fsm[nj][icen]->Eval(recopt));
  else{
    int ib=-1;
    ib = (int)(smpt - 30.)/20.;
    if(ib<0)return recopt;
    else smpt = gRandom->Gaus(smpt,lptsmf[nj][icen][ib]);
  }
  return smpt;
}

int GetCBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 2.5% bins of cross section                                                                                                        
  //! in 0-39 bins                                                                                                                                               

  if(bin<2)ibin=0; //! 0-5%                                                                                                                                      
  else if(bin>=2 && bin<4)ibin=1;   //! 5-10%                                                                                                                    
  else if(bin>=4 && bin<12)ibin=2;  //! 10-30%                                                                                                                   
  else if(bin>=12&& bin<20)ibin=3;  //! 30-50%                                                                                                                   
  else if(bin>=20&& bin<28)ibin=4;  //! 50-70%                                                                                                                   
  else if(bin>=28&& bin<36)ibin=5;  //! 70-90%                                                                                                                   
  else ibin=5; //! 90-100% use the same 70-90% factors

  return ibin;
}

#endif
