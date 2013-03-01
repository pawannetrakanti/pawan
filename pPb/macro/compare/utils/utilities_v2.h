#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <utility>
#include <TROOT.h>
#include <TStyle.h>

#include "TFile.h"
#include "TCanvas.h"

#include "TH1F.h"
#include "TH1D.h"

#include "TH2F.h"
#include "TH2D.h"

#include "TF1.h"

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "TDirectory.h"
#include "TDirectoryFile.h"

#include "TChain.h"
#include "TGraph.h"
#include "TSpline.h"

#include "TKey.h"
#include "TMath.h"

#include "TLegend.h"
#include "TLine.h"

#include "TSpline.h"

using namespace std;


// In case any of these are called in this file =====================================
// In case any of these are called in this file =====================================    


TH1D* GetDummyHist1(Int_t bins,Float_t xmin, Float_t xmax, Float_t min, Float_t max,const char *ttl,const char *xttl,const char *yttl, Bool_t stat) {

   TH1D *dum;
   dum = new TH1D(ttl,"",bins,xmin,xmax);

   dum->SetMinimum(min);
   dum->SetMaximum(max);
   dum->SetStats(stat);

   dum->GetYaxis()->SetTitle(yttl);
   dum->GetYaxis()->CenterTitle();
   dum->GetXaxis()->SetTitle(xttl);
   dum->GetXaxis()->CenterTitle();

   return dum;

}


TH1D* GetDummyHist(Float_t xmin, Float_t xmax, Float_t min, Float_t max,const char *ttl,const char *xttl,const char *yttl, Bool_t stat) {

   TH1D *dum;
   dum = new TH1D(ttl,"",100,xmin,xmax);

   dum->SetMinimum(min);
   dum->SetMaximum(max);
   dum->SetStats(stat);

   dum->GetYaxis()->SetTitle(yttl);
   dum->GetYaxis()->CenterTitle();
   dum->GetXaxis()->SetTitle(xttl);
   dum->GetXaxis()->CenterTitle();

   return dum;

}

// convert hist to tgraph
TGraphErrors* TgraphItWithRange(TH1D* hist, double min, double max){

  int nbins = 0;

  for(int i=0;i<hist->GetNbinsX();i++){
    double xvalue = hist->GetBinCenter(i+1);
    if(xvalue>min && xvalue<max) nbins++;
  }

  TGraphErrors *tg;

  const int nlines = nbins;

  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];

  int counter = 0;
  for(int i = 0; i<hist->GetNbinsX(); i++ ){
    double xvalue = hist->GetBinCenter(i+1);
    if(xvalue>min && xvalue<max){
      pt[counter] = hist->GetBinCenter(i+1);
      xsec[counter] = hist->GetBinContent(i+1);
      xsecerr[counter] = hist->GetBinError(i+1);
      pterr[counter] = 0;
      counter++;
    }
  }

  tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  return tg;
}


TGraphErrors* TgraphIt(TH1D* hist){

  TGraphErrors *tg;
  int nbins = hist->GetNbinsX();

  const int nlines = nbins;

  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];

  for(int i = 0; i<nbins; i++ ){
    pt[i] = hist->GetBinCenter(i+1);
    xsec[i] = hist->GetBinContent(i+1);
    xsecerr[i] = hist->GetBinError(i+1);
    pterr[i] = 0;
  }

  tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  return tg;
}


TGraphErrors* TgraphIt(TH1F* hist){

  TGraphErrors *tg;
  int nbins = hist->GetNbinsX();

  const int nlines = nbins;

  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];

  for(int i = 0; i<nbins; i++ ){
    pt[i] = hist->GetBinCenter(i+1);
    xsec[i] = hist->GetBinContent(i+1);
    xsecerr[i] = hist->GetBinError(i+1);
    pterr[i] = 0;
  }

  tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  return tg;
}

TH1D* HistIt(TGraph* tg){

  int nbins = tg->GetN();
  const int nlines = nbins;

  double pt   [nlines], xsec[nlines], xsecerr[nlines];
  //double pterr[nlines], xsecerr[nlines];

  for(int i = 0; i<nbins; i++ ){
    tg->GetPoint(i,pt[i],xsec[i]);
    //pterr[i]   = tg->GetErrorX(i);
    xsecerr[i] = tg->GetErrorY(i);
  }

  double bins[nlines];
  bins[0] = 0.0;
  for(int i=0;i<nlines;i++){
    double dbin = pt[i] - bins[i];
    bins[i+1]   = pt[i] + dbin;
  }

  TH1D *hist = new TH1D("","",nbins,bins);
  
  for(int i = 0; i<nbins; i++ ){
    hist->SetBinContent(i+1,xsec[i]);
    hist->SetBinError  (i+1,xsecerr[i]);
  }
  
  return hist;
}

TH1D* ImprovedHistIt(TGraph* tg){

  int nbins = tg->GetN();
  const int nlines = nbins;

  double pt[nlines], xsec[nlines], xsecerr[nlines];
  //double pterr[nlines], xsecerr[nlines];

  for(int i = 0; i<nbins; i++ ){
    tg->GetPoint(i,pt[i],xsec[i]);
    //pterr[i]   = tg->GetErrorX(i);
    xsecerr[i] = tg->GetErrorY(i);
  }

  double bins[nlines];
  if(pt[0]!=0) bins[0] = pt[0] - 0.5*(pt[1]-pt[0]);
  else bins[0] = 0.0;

  int ncount=0;
  for(int i=0;i<nlines;i++){
    double dbin = pt[i] - bins[i];
    bins[i+1] = pt[i] + dbin;
    cout<<"pt = "<<pt[i]<<" dbin = "<<dbin<<" bin boundary = "<<bins[i]<<endl;
    ncount++;
  }

  bins[ncount+1] = bins[ncount] + (bins[ncount]-bins[ncount-1]); // final bin


  TH1D *hist = new TH1D("","",nbins,bins);

  for(int i = 0; i<nbins; i++ ){
    hist->SetBinContent(i+1,xsec[i]);
    hist->SetBinError  (i+1,xsecerr[i]);
  }

  return hist;
}


TH1D* DivideByNearBinRAA(TH1D* num, TH1D* den){

  int nNum = num->GetXaxis()->GetNbins();
  TH1D* ratio = (TH1D*) num->Clone("ratio");
  ratio->Reset();

  for(int i=0;i<nNum;i++){

    double xvalue = num->GetXaxis()->GetBinCenter(i+1);
    double denvalue = den->GetBinContent(den->GetXaxis()->FindBin(xvalue));
    if (denvalue==0) return 0; //                                                                                                                           
    double rat = num->GetBinContent(i+1)/denvalue;
    ratio->SetBinContent(i+1,rat);
    ratio->SetBinError(i+1,0.0001*rat); // need to be calculated correctly                                                                                  
  }
  return ratio;
}


TH1D* ratio_hist_to_tspline(TH1D* num, TSpline3* ts, double minx, double maxx){

  cout<<"[Ratio to fit with range used]"<<endl;

  TH1D *hRatio = (TH1D*) num->Clone("hRatio");

  int nbin = num->GetNbinsX();

  for(int i=0;i<nbin;i++){

    double cms_value = (double) ts->Eval(hRatio->GetBinCenter(i+1));

    double ratio = -999;
    double ratio_err = 0;

    if(hRatio->GetBinCenter(i+1)>minx && hRatio->GetBinCenter(i+1)<maxx){
      ratio = hRatio->GetBinContent(i+1)/cms_value;
      ratio_err =  (hRatio->GetBinError(i+1)/hRatio->GetBinContent(i+1))*ratio;
    }else{
      ratio = -999;
      ratio_err = 0.0;
    }

    hRatio->SetBinContent(i+1,ratio);
    hRatio->SetBinError(i+1,ratio_err);
  }

  return hRatio;

}

TGraphErrors* ratio_func_to_tspline(TGraphErrors* num, TSpline3* ts, double minx, double maxx){

  cout<<"[Ratio to fit used]"<<endl;

  TGraphErrors *tg;

  int nbin = num->GetN();

  const int nlines = nbin;
  double pt[nlines], pterr[nlines];
  double xsec[nlines], xsecerr[nlines];
  double ratio[nlines], ratioerr[nlines];

  for(int i=0;i<nbin;i++){

    num->GetPoint(i,pt[i],xsec[i]);
    xsecerr[i] = num->GetErrorY(i);
    double cms_value = (double) ts->Eval(pt[i]);

    if(pt[i]>minx && pt[i]<maxx){
      ratio[i] = xsec[i]/cms_value;
      ratioerr[i] = (xsecerr[i]/xsec[i])*ratio[i];
    }else{
      ratio[i] = -999;
      ratioerr[i] = 0.0;
    }

    pterr[i] = 0;
  }

  tg = new TGraphErrors(nlines,pt,ratio,pterr,ratioerr);
  return tg;
}


TGraphErrors* ratio_tg_to_hist(TGraphErrors* num, TH1D *hist, double minx, double maxx){

  cout<<"[Ratio to fit used]"<<endl;

  TGraphErrors *tg;

  int nbin = num->GetN();

  const int nlines = nbin;
  double pt[nlines], pterr[nlines];
  double xsec[nlines], xsecerr[nlines];
  double ratio[nlines], ratioerr[nlines];

  for(int i=0;i<nbin;i++){

    num->GetPoint(i,pt[i],xsec[i]);
    xsecerr[i] = num->GetErrorY(i);

    if(pt[i]>minx && pt[i]<maxx){
      double ptx= pt[i]+0.0001;
      int ptbin = hist->FindBin(ptx);
      double cms_value = (double) hist->GetBinContent(ptbin);
      ratio[i] = xsec[i]/cms_value;
      ratioerr[i] = (xsecerr[i]/xsec[i])*ratio[i];
      cout<<"pt = "<<pt[i]<<" pt bin = "<<hist->GetBinCenter(ptbin)<<endl;
      //cout<<"cms value = "<<cms_value<<" xsec[i] = "<<xsec[i]<<" and ratio = "<<ratio[i]<<endl;
    }else{
      ratio[i] = -999;
      ratioerr[i] = 0.0;
    }

    //cout<<"pt = "<<pt[i]<<" ratio = "<<ratio[i]<<endl;

    pterr[i] = 0;
  }

  tg = new TGraphErrors(nlines,pt,ratio,pterr,ratioerr);
  return tg;
}


TH1D* ratio_hist_to_func(TH1D* num, TF1* f3, double minx, double maxx){

  cout<<"[Ratio to fit used]"<<endl;

  TH1D *hRatio = (TH1D*) num->Clone("hRatio");

  int nbin = hRatio->GetNbinsX();

  for(int i=0;i<nbin;i++){
    double cms_value = (double) f3->Eval(hRatio->GetBinCenter(i+1));

    double ratio = -999;
    double ratio_err = 0;

    if(hRatio->GetBinCenter(i+1)>minx && hRatio->GetBinCenter(i+1)<maxx){
      ratio = hRatio->GetBinContent(i+1)/cms_value;
      ratio_err =  (hRatio->GetBinError(i+1)/hRatio->GetBinContent(i+1))*ratio;
    }else{
      ratio = -999;
      ratio_err = 0.0;
    }

    hRatio->SetBinContent(i+1,ratio);
    hRatio->SetBinError(i+1,ratio_err);
  }
  return hRatio;
}

TH1F* ratio_hist_to_func(TH1F* num, TF1* f3, double minx, double maxx){

  cout<<"[Ratio to fit used]"<<endl;

  TH1F *hRatio = (TH1F*) num->Clone("hRatio");

  int nbin = hRatio->GetNbinsX();

  for(int i=0;i<nbin;i++){
    double cms_value = (double) f3->Eval(hRatio->GetBinCenter(i+1));

    double ratio = -999;
    double ratio_err = 0;

    if(hRatio->GetBinCenter(i+1)>minx && hRatio->GetBinCenter(i+1)<maxx){
      ratio = hRatio->GetBinContent(i+1)/cms_value;
      ratio_err =  (hRatio->GetBinError(i+1)/hRatio->GetBinContent(i+1))*ratio;
    }else{
      ratio = -999;
      ratio_err = 0.0;
    }

    hRatio->SetBinContent(i+1,ratio);
    hRatio->SetBinError(i+1,ratio_err);
  }
  return hRatio;
}

TGraphErrors* ratio_tg_to_tg(TGraphErrors* num1, TGraphErrors *num2){
  cout<<"[Ratio to tg used]"<<endl;

  TGraphErrors *tg;

  int nbin = num1->GetN();

  const int nlines = nbin;
  double x1[nlines], x1err[nlines];
  double y1[nlines], y1err[nlines];
  double y2[nlines], y2err[nlines];

  double ratio[nlines], ratioerr[nlines];

  for(int i=0;i<nbin;i++){

    num1->GetPoint(i,x1[i],y1[i]);
    y1err[i] = num1->GetErrorY(i);

    num2->GetPoint(i,x1[i],y2[i]);
    y2err[i] = num2->GetErrorY(i);

    if(y2[i]!=0){
      ratio[i]    = y1[i]/y2[i];
      ratioerr[i] = ratio[i]*sqrt(pow(y1err[i]/y1[i],2) + pow(y2err[i]/y2[i],2));
    }else{
      ratio[i] = -999;
      ratioerr[i] = 0.0;
    }
  }

  tg = new TGraphErrors(nlines,x1,ratio,x1err,ratioerr);
  return tg;
}

void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge, const Float_t asyoffset) {
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];

  Float_t PadWidth =
     (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
                       (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
     (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
                         (1.0/(1.0-edge))+(Float_t)rows-2.0);

   Xlow[0] = leftOffset - asyoffset;
   Xup[0]  = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
     Xlow[i] = Xup[0] + (i-1)*PadWidth;
     Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
     Ylow[i] = Yup[rows-1] + ct*PadHeight;
     Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
     ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
     for(Int_t j=0;j<rows;j++) {
       canv->cd();
       padName = Form("p_%d_%d",i,j);
       pad[i][j] = new TPad(padName.Data(),padName.Data(),
                            Xlow[i],Ylow[j],Xup[i],Yup[j]);

       // this is hacked version to create aysmmetric pads around low 
       if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
       else pad[i][j]->SetLeftMargin(0);

       if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
       else pad[i][j]->SetRightMargin(0);

       if(j==0) pad[i][j]->SetTopMargin(edge);
       else pad[i][j]->SetTopMargin(0);

       //if(j==0) pad[i][j]->SetTopMargin(edge*3.5);
       //else pad[i][j]->SetTopMargin(0);


       if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
       else pad[i][j]->SetBottomMargin(0);

       pad[i][j]->Draw();
       pad[i][j]->cd();
       pad[i][j]->SetNumber(columns*j+i+1);
     }
   }
}


void makeMultiPanelCanvasAsy(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge, const Float_t asyoffset) {
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
   Float_t PadWidth =
     (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
                       (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
     (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
                         (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
     Xlow[i] = Xup[0] + (i-1)*PadWidth;
     Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;

   for(Int_t i=rows-2;i>0;i--) {
     Ylow[i] = Yup[rows-1] + ct*PadHeight;
     Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
     ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
     for(Int_t j=0;j<rows;j++) {
       canv->cd();
       padName = Form("p_%d_%d",i,j);

       // this is hacked version to create aysmmetric pads around low 
       if(j==0){
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
       }else{
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
       }

       if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
       else pad[i][j]->SetLeftMargin(0);

       if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
       else pad[i][j]->SetRightMargin(0);

       if(j==0) pad[i][j]->SetTopMargin(edge);
       else pad[i][j]->SetTopMargin(0);

       if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
       else pad[i][j]->SetBottomMargin(0);

       pad[i][j]->Draw();
       pad[i][j]->cd();
       pad[i][j]->SetNumber(columns*j+i+1);
     }
   }
}

void makeMultiPanelCanvasWithGap(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge, const Float_t asyoffset) {
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];

  Float_t PadWidth =
     (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
                       (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
     (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
                         (1.0/(1.0-edge))+(Float_t)rows-2.0);

   //PadHeight = 0.5*PadHeight;

   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
     Xlow[i] = Xup[0] + (i-1)*PadWidth;
     Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
     if(i==rows-2){
       Ylow[i] = Yup[rows-1] + ct*PadHeight;
       Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
     }else{
       Ylow[i] = Yup[rows-1] + ct*PadHeight;
       Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
       //Yup[i] = 0.2*Yup[i];
     }
     ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
     for(Int_t j=0;j<rows;j++) {
       canv->cd();
       padName = Form("p_%d_%d",i,j);
       //pad[i][j] = new TPad(padName.Data(),padName.Data(),
       //Xlow[i],Ylow[j],Xup[i],Yup[j]);
       // this is hacked version to create aysmmetric pads around low 
       if(j==0){
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
       }else{
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
       }


       if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
       else pad[i][j]->SetLeftMargin(0);

       if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
       else pad[i][j]->SetRightMargin(0);

       if(j==0) pad[i][j]->SetTopMargin(edge);
       //else pad[i][j]->SetTopMargin(0.01);
       else pad[i][j]->SetTopMargin(0.02);

       //if(j==0) pad[i][j]->SetTopMargin(edge*3.5);
       //else pad[i][j]->SetTopMargin(0.0);

       if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
       else pad[i][j]->SetBottomMargin(0.15);

       pad[i][j]->Draw();
       pad[i][j]->cd();
       pad[i][j]->SetNumber(columns*j+i+1);
     }
   }
}

#endif
