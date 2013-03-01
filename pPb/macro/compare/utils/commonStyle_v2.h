
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

#include "TLatex.h"

// color scheme for different color
Int_t color[13]={kViolet+2,kBlue,kAzure+6,kGreen-3,kOrange-5,kOrange-3,kOrange+4,kRed-3,kRed,kRed+2,kMagenta+1,kRed+1,kGreen+3};

// color scheme for fill area - fcolor with lcolor as a boundary
Int_t fcolor[5]={kRed-10,kBlue-10,kGreen-8,kOrange-4, kRed-6};
Int_t lcolor[5]={kRed+1,kBlue-3,kGreen+3,kOrange+3, kRed+2};

// marker scheme 
Int_t fmstyle[6] = {20,21,22,23,29,3};
Int_t emstyle[6] = {24,25,26,27,30,28};


void putAorB(const char *text){
  TLatex * texB;
  texB = new TLatex(0.2,0.88,text);
  texB->SetTextSize(0.045);
  texB->SetTextFont(42);
  texB->SetLineWidth(2);
  texB->SetNDC();
  texB->Draw();
}

void printCanvases(TCanvas * MyCanvas, const char * name, int log=0)
{
   MyCanvas->cd();
   MyCanvas->SetLogx(log);
   printf("canvas name: %s\n",name);

   // add some text labels  
   //double ndcX = 0.25;
   //double ndcY = 0.9;
   //TLatex * tex;

   //MyCanvas->Print(Form("./fig/%s.eps",name));
   MyCanvas->Print(Form("./fig/%s.gif",name));
   MyCanvas->Print(Form("./fig/%s.pdf",name));

}



void th1Style1( TGraph *g=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1)
{
   g->SetMarkerStyle(marker);
   if(mcol>95 || mcol<14){
      g->SetMarkerColor(mcol);
   }else{
      g->SetMarkerColor(color[mcol-14]);
   }
   g->SetMarkerSize(msize);
   if(lcol>95 || lcol<14){
      g->SetLineColor(lcol);
   }else{
      g->SetLineColor(color[lcol-14]);
   }
   g->SetLineWidth(lsize);
   g->SetLineStyle(lstyle);

   if(dStyle==1) g->Draw("PZsame");
   if(dStyle==2) g->Draw("Xlsame");
   if(dStyle==3) g->Draw("histsame");
   if(dStyle==4) g->Draw("plxsame");
}

void th1Style1( TH1 *h=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1)
{
   h->SetMarkerStyle(marker);
   if(mcol>95 || mcol<14){
      h->SetMarkerColor(mcol);
   }else{
      h->SetMarkerColor(color[mcol-14]);
   }
   h->SetMarkerSize(msize);
   if(lcol>95 || lcol<14){
      h->SetLineColor(lcol);
   }else{
      h->SetLineColor(color[lcol-14]);
   }
   h->SetLineWidth(lsize);
   h->SetLineStyle(lstyle);

   if(dStyle==1) h->Draw("PZsame");
   if(dStyle==2) h->Draw("Csame");
   if(dStyle==3) h->Draw("histsame");
   if(dStyle==4) h->Draw("PCsame");
}


void th1Style2( TH1 *h=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1, bool drawstat=false)
{
  h->SetMarkerStyle(marker);
  if(mcol>95 || mcol<14){
    h->SetMarkerColor(mcol);
  }else{
    h->SetMarkerColor(color[mcol-14]);
  }
  h->SetMarkerSize(msize);
  if(lcol>95 || lcol<14){
    h->SetLineColor(lcol);
  }else{
    h->SetLineColor(color[lcol-14]);
  }
  h->SetLineWidth(lsize);
  h->SetLineStyle(lstyle);

  if(drawstat){
    if(dStyle==1) h->Draw("PZsames");
    if(dStyle==2) h->Draw("Csames");
    if(dStyle==3) h->Draw("histsames");
    if(dStyle==4) h->Draw("PCsames");
  }else{
    if(dStyle==1) h->Draw("PZsame");
    if(dStyle==2) h->Draw("Csame");
    if(dStyle==3) h->Draw("histsame");
    if(dStyle==4) h->Draw("PCsame");
  }
}


void setPad1by2(TCanvas *can){

  can->cd();

  TPad *pp1=0, *pp2=0;

  pp1 = new TPad("p1","",0,0.34,1,1,0,0,0);
  pp1->SetBottomMargin(0.0);
  pp1->SetTopMargin(0.05*(1/0.72));

  pp1->Draw();
  pp1->cd();
  pp1->SetNumber(1);

  can->cd();

  pp2 = new TPad("pp2","",0,0.0,1,0.34,0,0,0);
  pp2->SetTopMargin(0);
  pp2->SetBottomMargin(0.14*(1/0.34));
  pp2->Draw();
  pp2->cd();
  pp2->SetNumber(2);

}


void setLowerPad1by2(TH1D *hlow){
  /*
  hlow->GetXaxis()->SetLabelSize(19);
  hlow->GetXaxis()->SetLabelFont(43);
  hlow->GetXaxis()->SetLabelOffset(0.06);
  hlow->GetXaxis()->SetTitleSize(23);
  hlow->GetXaxis()->SetTitleFont(43);
  hlow->GetXaxis()->SetTitleOffset(4.6);
  hlow->GetXaxis()->CenterTitle();

  hlow->GetYaxis()->SetLabelSize(16);
  hlow->GetYaxis()->SetLabelFont(43);
  hlow->GetYaxis()->SetTitleSize(18);
  hlow->GetYaxis()->SetTitleFont(43);
  hlow->GetYaxis()->SetTitleOffset(2.8);
  hlow->GetYaxis()->CenterTitle();

  hlow->GetYaxis()->SetNdivisions(906);
  */

  hlow->GetYaxis()->SetLabelFont(43);
  hlow->GetYaxis()->SetTitleFont(43);
  hlow->GetYaxis()->SetLabelSize(20);
  hlow->GetYaxis()->SetTitleSize(22);
  hlow->GetYaxis()->SetTitleOffset(2.1);

  hlow->GetXaxis()->SetLabelFont(43);
  hlow->GetXaxis()->SetTitleFont(43);
  hlow->GetXaxis()->SetLabelSize(20);
  hlow->GetXaxis()->SetTitleSize(22);
  hlow->GetXaxis()->SetTitleOffset(3.8);

  hlow->GetYaxis()->SetNdivisions(906);

}
void MakeHist(TH1 *histo,const char *xtitle,const char *ytitle)
{
  histo->SetStats(0);
  histo->SetTitle("");
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetXaxis()->CenterTitle(true);
  histo->GetXaxis()->SetTitleFont(42);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetXaxis()->SetTitleSize(0.07);
  histo->GetXaxis()->SetTitleOffset(1.15);
  histo->GetXaxis()->SetLabelSize(0.07);
  histo->GetXaxis()->SetLabelOffset(0.005);
  histo->GetXaxis()->SetNdivisions(507);

  histo->GetYaxis()->SetTitle(ytitle);
  histo->GetYaxis()->CenterTitle(true);
  histo->GetYaxis()->SetTitleSize(0.07);
  histo->GetYaxis()->SetTitleOffset(1.50);
  histo->GetYaxis()->SetLabelSize(0.07);
  histo->GetYaxis()->SetNdivisions(507);
}
void MakeHist2D(TH2 *histo,const char *xtitle,const char *ytitle)
{
  histo->SetStats(0);
  histo->SetTitle("");
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetXaxis()->CenterTitle(true);
  histo->GetXaxis()->SetTitleFont(42);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetXaxis()->SetTitleSize(0.06);
  histo->GetXaxis()->SetTitleOffset(1.15);
  histo->GetXaxis()->SetLabelSize(0.06);
  histo->GetXaxis()->SetLabelOffset(0.005);
  histo->GetXaxis()->SetNdivisions(507);

  histo->GetYaxis()->SetTitle(ytitle);
  histo->GetYaxis()->CenterTitle(true);
  histo->GetYaxis()->SetTitleSize(0.06);
  histo->GetYaxis()->SetTitleOffset(1.23);
  histo->GetYaxis()->SetLabelSize(0.06);
  histo->GetYaxis()->SetNdivisions(507);

  histo->GetZaxis()->SetTitleSize(0.06);
  histo->GetZaxis()->SetTitleOffset(1.50);
  histo->GetZaxis()->SetLabelSize(0.06);
  histo->GetZaxis()->SetNdivisions(507);

}

void LoadStyle()
{
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetErrorX(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetPadBorderSize(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetPalette(1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
}
void putCMSPrel(double x, double y, double size){
  TLatex *tex=0;
  tex = new TLatex(x,y,"CMS Preliminary");
  tex->SetTextSize(size);
  tex->SetLineWidth(2);
  tex->SetNDC();
  tex->Draw();
}
void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}

void drawText2(const char *text, float xp, float yp, int size,Color_t ccol){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(ccol);
  //tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
