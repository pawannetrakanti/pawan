#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>
#include "MultiCanvas.h"

using namespace std;

const double pi=acos(-1.);

int MatchEffAn(){

  float ketacut=2.0;
  const int ncen=7;
  const char *ksp  [ncen] = {"pbpb","pbpb","pbpb","pbpb","pbpb","pbpb","pp"};
  const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","pp"};
  //!                          0       1       2        3        4        5      6  

  bool iSave=true;
  bool iSigma=false;
  double min=80.;
  double max=425.;

  TFile *fin[ncen];
  for(int ic=0;ic<ncen;ic++){
    fin[ic] = new TFile(Form("rootfiles/MatchingEfficiency/Histo_%s_%s.root",ksp[ic],ccent[ic]));
  }

  const int knj=7;
  const char *calgo[knj]= {
    "icPu5Calo", 
    "ak2PF","ak3PF","ak4PF",
    "akPu2PF","akPu3PF","akPu4PF"
  };                                     

  const int minnj[4]={1,7, 13,19};
  const int maxnj[4]={7,13,19,24};
  
  const char *meta[2] = {"|#eta|<1.3","1.3<|#eta|<2.0"};

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
  gStyle->SetOptFit(1);
  gStyle->SetPadBorderMode(0);

  //! pt binning
  //double ptbins[] ={80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400};
  //const int bins  = sizeof(ptbins)/sizeof(Double_t) - 1;
  //const int nbins = bins;
  //std::cout<<"# of pt bins : " <<nbins<<std::endl;
  //std::cout<<std::endl;


  TH1F *hEffVsPt[knj][ncen], *hEffVsEta[knj][ncen], *hEffVsPhi[knj][ncen];

  int sty    [7] = {24,25,26,28,30,32,20};
  Color_t col[7] = {kRed,kBlue,kMagenta,kCyan,kGreen+1,kBlue+4,kRed+4};
  //return 0;

  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){
      int icen=ic;

      hEffVsPt[nj][ic] = (TH1F*)fin[ic]->Get(Form("hEffVsPt%s_%s",ksp[ic],calgo[nj]));
      hEffVsPt[nj][ic]->SetName(Form("hEffVsPt_%s_%s_%d",ksp[ic],calgo[nj],icen));
      hEffVsPt[nj][ic]->SetTitle("");
      hEffVsPt[nj][ic]->SetMaximum(1.0075);      
      hEffVsPt[nj][ic]->SetMinimum(0.964);      
      hEffVsPt[nj][ic]->GetXaxis()->SetRangeUser(min,max);
      hEffVsPt[nj][ic]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
      hEffVsPt[nj][ic]->GetXaxis()->CenterTitle(true);
      hEffVsPt[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hEffVsPt[nj][ic]->GetXaxis()->SetNoExponent();
      hEffVsPt[nj][ic]->GetXaxis()->SetTitleFont(42);
      hEffVsPt[nj][ic]->GetXaxis()->SetLabelFont(42);
      hEffVsPt[nj][ic]->GetXaxis()->SetTitleSize(0.06);
      hEffVsPt[nj][ic]->GetXaxis()->SetTitleOffset(1.15);
      hEffVsPt[nj][ic]->GetXaxis()->SetLabelSize(0.06);
      hEffVsPt[nj][ic]->GetXaxis()->SetLabelOffset(0.005);
      hEffVsPt[nj][ic]->GetXaxis()->SetNdivisions(505);
      hEffVsPt[nj][ic]->GetYaxis()->SetTitle("Jet Reconstrucntion Efficiency");
      hEffVsPt[nj][ic]->GetYaxis()->SetTitleSize(0.06);
      hEffVsPt[nj][ic]->GetYaxis()->SetTitleOffset(1.43);
      hEffVsPt[nj][ic]->GetYaxis()->SetLabelSize(0.06);
      hEffVsPt[nj][ic]->GetYaxis()->SetNdivisions(507);
      hEffVsPt[nj][ic]->GetYaxis()->SetDecimals(true);
      hEffVsPt[nj][ic]->GetYaxis()->SetTitleFont(42);
      hEffVsPt[nj][ic]->GetYaxis()->SetLabelFont(42);

      hEffVsEta[nj][ic] = (TH1F*)fin[ic]->Get(Form("hEffVsEta%s_%s",ksp[ic],calgo[nj]));
      hEffVsEta[nj][ic]->SetName(Form("hEffVsEta_%s_%s_%d",ksp[ic],calgo[nj],icen));
      hEffVsEta[nj][ic]->SetTitle("");
      hEffVsEta[nj][ic]->SetMaximum(1.0075);      
      hEffVsEta[nj][ic]->SetMinimum(0.964);      
      hEffVsEta[nj][ic]->GetXaxis()->SetRangeUser(-ketacut,ketacut);
      hEffVsEta[nj][ic]->GetXaxis()->SetTitle("GenJet #eta ");
      hEffVsEta[nj][ic]->GetXaxis()->CenterTitle(true);
      hEffVsEta[nj][ic]->GetXaxis()->SetTitleFont(42);
      hEffVsEta[nj][ic]->GetXaxis()->SetLabelFont(42);
      hEffVsEta[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hEffVsEta[nj][ic]->GetXaxis()->SetNoExponent();
      hEffVsEta[nj][ic]->GetXaxis()->SetTitleSize(0.06);
      hEffVsEta[nj][ic]->GetXaxis()->SetTitleOffset(0.98);
      hEffVsEta[nj][ic]->GetXaxis()->SetLabelSize(0.06);
      hEffVsEta[nj][ic]->GetXaxis()->SetNdivisions(507);
      hEffVsEta[nj][ic]->GetYaxis()->SetTitle("Jet Reconstruction Efficiency");    
      hEffVsEta[nj][ic]->GetYaxis()->SetTitleSize(0.06);
      hEffVsEta[nj][ic]->GetYaxis()->SetTitleOffset(1.41);
      hEffVsEta[nj][ic]->GetYaxis()->SetLabelSize(0.06);
      hEffVsEta[nj][ic]->GetYaxis()->SetNdivisions(507);
      hEffVsEta[nj][ic]->GetYaxis()->SetTitleFont(42);
      hEffVsEta[nj][ic]->GetYaxis()->SetLabelFont(42);
      hEffVsEta[nj][ic]->GetYaxis()->SetDecimals(true);

      hEffVsPhi[nj][ic] = (TH1F*)fin[ic]->Get(Form("hEffVsPhi%s_%s",ksp[ic],calgo[nj]));
      hEffVsPhi[nj][ic]->SetName(Form("hEffVsPhi_%s_%s_%d",ksp[ic],calgo[nj],icen));
      hEffVsPhi[nj][ic]->SetTitle("");
      hEffVsPhi[nj][ic]->SetMaximum(1.0075);      
      hEffVsPhi[nj][ic]->SetMinimum(0.964);      
      hEffVsPhi[nj][ic]->GetXaxis()->SetRangeUser(-pi,pi);
      hEffVsPhi[nj][ic]->GetXaxis()->SetTitle("GenJet #phi (radian) ");
      hEffVsPhi[nj][ic]->GetXaxis()->CenterTitle(true);
      hEffVsPhi[nj][ic]->GetXaxis()->SetTitleFont(42);
      hEffVsPhi[nj][ic]->GetXaxis()->SetLabelFont(42);
      hEffVsPhi[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hEffVsPhi[nj][ic]->GetXaxis()->SetNoExponent();
      hEffVsPhi[nj][ic]->GetXaxis()->SetTitleSize(0.06);
      hEffVsPhi[nj][ic]->GetXaxis()->SetTitleOffset(1.15);
      hEffVsPhi[nj][ic]->GetXaxis()->SetLabelSize(0.06);
      hEffVsPhi[nj][ic]->GetXaxis()->SetNdivisions(507);
      hEffVsPhi[nj][ic]->GetYaxis()->SetTitle("Jet Reconstruction Efficiency");    
      hEffVsPhi[nj][ic]->GetYaxis()->SetTitleSize(0.06);
      hEffVsPhi[nj][ic]->GetYaxis()->SetTitleOffset(1.37);
      hEffVsPhi[nj][ic]->GetYaxis()->SetLabelSize(0.06);
      hEffVsPhi[nj][ic]->GetYaxis()->SetNdivisions(507);
      hEffVsPhi[nj][ic]->GetYaxis()->SetTitleFont(42);
      hEffVsPhi[nj][ic]->GetYaxis()->SetLabelFont(42);
      hEffVsPhi[nj][ic]->GetYaxis()->SetDecimals(true);
    }
  }
  //return 0;

  const int pal[3]  = {4,5,6}; 
  //const int cbin[4] = {0,2,4,5,6};
  int ipad=0;
  int maxc=3;
  int maxr=1;

  TLine *line = new TLine(min,1.0,max+100,1.0);
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  TCanvas *c3 = new TCanvas("c3","pT Matching Efficiency",102,141,1399,497);
  //makeMultiPanelCanvas(c3,maxc,maxr,0.0,0.0,0.22,0.22,0.02);
  c3->Divide(maxc,maxr,0,0);
  TLegend *l3;
  l3 = new TLegend(0.2596222,0.2582326,0.4984747,0.7083629,NULL,"BRNDC");
  l3->SetHeader("PYTHIA (Z2) + HYDJET 1.8");
  l3->SetBorderSize(0);
  l3->SetTextFont(42);
  l3->SetTextSize(0.05);
  l3->SetLineColor(1);
  l3->SetLineStyle(1);
  l3->SetLineWidth(1);
  l3->SetFillColor(10);
  l3->SetFillStyle(1001);

  ipad=0;
  for(int nj=0;nj<3;nj++){
    //std::cout<<std::endl;

    c3->cd(++ipad);
    if(ipad==1)gPad->SetLeftMargin(0.17);
    else if(ipad==3)gPad->SetRightMargin(0.07);    
    gPad->SetBottomMargin(0.15);        

    for(int ic=0;ic<ncen;ic++){
      hEffVsPt[pal[nj]][ic]->SetMarkerStyle(sty[ic]);
      hEffVsPt[pal[nj]][ic]->SetMarkerColor(col[ic]);
      hEffVsPt[pal[nj]][ic]->SetLineColor(col[ic]);
      hEffVsPt[pal[nj]][ic]->SetMarkerSize(1.0);

      if(ic==0){
	hEffVsPt[pal[nj]][ic]->Draw("p");
      }else{
	hEffVsPt[pal[nj]][ic]->Draw("psame");
      }
      if(nj==0){
	l3->AddEntry(hEffVsPt[pal[nj]][ic],Form("%s",ccent[ic]),"p"); 
      }
    }//! ic 
    if(ipad==1){
      l3->Draw();
      TPaveText *pt0 = new TPaveText(0.6130501,0.883784,0.9029564,0.9487681,"brNDC");
      pt0->SetBorderSize(0);
      pt0->SetFillColor(10);
      pt0->SetTextFont(62);
      TText *text0 = pt0->AddText("CMS Preliminary");
      text0->SetTextSize(0.05);
      pt0->Draw();
      
      TPaveText *pt1 = new TPaveText(0.553961,0.474384,0.8401741,0.5826909,"brNDC");
      pt1->SetBorderSize(0);
      pt1->SetFillColor(10);
      pt1->SetTextFont(42);
      TText *text1 = pt1->AddText("Anti-k_{T} (R=0.2)");
      text1->SetTextSize(0.055);
      TText *text2 = pt1->AddText("Particle Flow Jets");
      text2->SetTextSize(0.055);
      pt1->Draw();
    }else if(ipad==2){      
      TPaveText *pt2 = new TPaveText(0.553961,0.474384,0.8401741,0.5826909,"brNDC");
      pt2->SetBorderSize(0);
      pt2->SetFillColor(10);
      pt2->SetTextFont(42);
      TText *text3 = pt2->AddText("Anti-k_{T} (R=0.3)");
      text3->SetTextSize(0.055);
      TText *text4 = pt2->AddText("Particle Flow Jets");
      text4->SetTextSize(0.055);
      pt2->Draw();
    }else if(ipad==3){
      TPaveText *pt3 = new TPaveText(0.553961,0.474384,0.8401741,0.5826909,"brNDC");
      pt3->SetBorderSize(0);
      pt3->SetFillColor(10);
      pt3->SetTextFont(42);
      TText *text5 = pt3->AddText("Anti-k_{T} (R=0.4)");
      text5->SetTextSize(0.055);
      TText *text6 = pt3->AddText("Particle Flow Jets");
      text6->SetTextSize(0.055);
      pt3->Draw();
    }
  }      
  if(iSave){
    c3->SaveAs("AN/Efficiency/JetRecoEff_Pt.png");
    c3->SaveAs("AN/Efficiency/JetRecoEff_Pt.pdf");
    c3->SaveAs("AN/Efficiency/JetRecoEff_Pt.C"  );
    c3->SaveAs("AN/Efficiency/JetRecoEff_Pt.eps");
  }


  TCanvas *c4 = new TCanvas("c4","eta Matching Efficiency",102,141,1399,497);
  //makeMultiPanelCanvas(c4,maxc,maxr,0.0,0.0,0.22,0.22,0.02);
  c4->Divide(maxc,maxr,0,0);

  TLegend *l4;
  l4 = new TLegend(0.2596222,0.2582326,0.4984747,0.7083629,NULL,"BRNDC");
  l4->SetHeader("PYTHIA (Z2) + HYDJET 1.8");
  l4->SetBorderSize(0);
  l4->SetTextFont(42);
  l4->SetTextSize(0.05);
  l4->SetLineColor(1);
  l4->SetLineStyle(1);
  l4->SetLineWidth(1);
  l4->SetFillColor(10);
  l4->SetFillStyle(1001);

  ipad=0;
  for(int nj=0;nj<3;nj++){
    //std::cout<<std::endl;

    c4->cd(++ipad);
    if(ipad==1)gPad->SetLeftMargin(0.17);
    else if(ipad==3)gPad->SetRightMargin(0.07);    
    gPad->SetBottomMargin(0.15);        

    for(int ic=0;ic<ncen;ic++){
      hEffVsEta[pal[nj]][ic]->SetMarkerStyle(sty[ic]);
      hEffVsEta[pal[nj]][ic]->SetMarkerColor(col[ic]);
      hEffVsEta[pal[nj]][ic]->SetLineColor(col[ic]);
      hEffVsEta[pal[nj]][ic]->SetMarkerSize(1.0);

      if(ic==0){
	hEffVsEta[pal[nj]][ic]->Draw("p");
      }else{
	hEffVsEta[pal[nj]][ic]->Draw("psame");
      }
      if(nj==0){
	l4->AddEntry(hEffVsEta[pal[nj]][ic],Form("%s",ccent[ic]),"p"); 
      }
    }//! ic 
    if(ipad==1){
      l4->Draw();
      TPaveText *pt0 = new TPaveText(0.6130501,0.883784,0.9029564,0.9487681,"brNDC");
      pt0->SetBorderSize(0);
      pt0->SetFillColor(10);
      pt0->SetTextFont(62);
      TText *text0 = pt0->AddText("CMS Preliminary");
      text0->SetTextSize(0.05);
      pt0->Draw();
      
      TPaveText *pt1 = new TPaveText(0.553961,0.474384,0.8401741,0.5826909,"brNDC");
      pt1->SetBorderSize(0);
      pt1->SetFillColor(10);
      pt1->SetTextFont(42);
      TText *text1 = pt1->AddText("Anti-k_{T} (R=0.2)");
      text1->SetTextSize(0.055);
      TText *text2 = pt1->AddText("Particle Flow Jets");
      text2->SetTextSize(0.055);
      pt1->Draw();
    }else if(ipad==2){      
      TPaveText *pt2 = new TPaveText(0.553961,0.474384,0.8401741,0.5826909,"brNDC");
      pt2->SetBorderSize(0);
      pt2->SetFillColor(10);
      pt2->SetTextFont(42);
      TText *text3 = pt2->AddText("Anti-k_{T} (R=0.3)");
      text3->SetTextSize(0.055);
      TText *text4 = pt2->AddText("Particle Flow Jets");
      text4->SetTextSize(0.055);
      pt2->Draw();
    }else if(ipad==3){
      TPaveText *pt3 = new TPaveText(0.553961,0.474384,0.8401741,0.5826909,"brNDC");
      pt3->SetBorderSize(0);
      pt3->SetFillColor(10);
      pt3->SetTextFont(42);
      TText *text5 = pt3->AddText("Anti-k_{T} (R=0.4)");
      text5->SetTextSize(0.055);
      TText *text6 = pt3->AddText("Particle Flow Jets");
      text6->SetTextSize(0.055);
      pt3->Draw();
    }
  }      
  if(iSave){
    c4->SaveAs("AN/Efficiency/JetRecoEff_Eta.png");
    c4->SaveAs("AN/Efficiency/JetRecoEff_Eta.pdf");
    c4->SaveAs("AN/Efficiency/JetRecoEff_Eta.C"  );
    c4->SaveAs("AN/Efficiency/JetRecoEff_Eta.eps");
  }


  TCanvas *c5 = new TCanvas("c5","phi Matching Efficiency",102,141,1399,497);
  //makeMultiPanelCanvas(c5,maxc,maxr,0.0,0.0,0.22,0.22,0.02);
  c5->Divide(maxc,maxr,0,0);

  TLegend *l5;
  l5 = new TLegend(0.2596222,0.2582326,0.4984747,0.7083629,NULL,"BRNDC");
  l5->SetHeader("PYTHIA (Z2) + HYDJET 1.8");
  l5->SetBorderSize(0);
  l5->SetTextFont(42);
  l5->SetTextSize(0.05);
  l5->SetLineColor(1);
  l5->SetLineStyle(1);
  l5->SetLineWidth(1);
  l5->SetFillColor(10);
  l5->SetFillStyle(1001);

  ipad=0;
  for(int nj=0;nj<3;nj++){
    //std::cout<<std::endl;

    c5->cd(++ipad);
    if(ipad==1)gPad->SetLeftMargin(0.17);
    else if(ipad==3)gPad->SetRightMargin(0.07);    
    gPad->SetBottomMargin(0.15);        

    for(int ic=0;ic<ncen;ic++){
      hEffVsPhi[pal[nj]][ic]->SetMarkerStyle(sty[ic]);
      hEffVsPhi[pal[nj]][ic]->SetMarkerColor(col[ic]);
      hEffVsPhi[pal[nj]][ic]->SetLineColor(col[ic]);
      hEffVsPhi[pal[nj]][ic]->SetMarkerSize(1.0);

      if(ic==0){
	hEffVsPhi[pal[nj]][ic]->Draw("p");
      }else{
	hEffVsPhi[pal[nj]][ic]->Draw("psame");
      }
      if(nj==0){
	l5->AddEntry(hEffVsPhi[pal[nj]][ic],Form("%s",ccent[ic]),"p"); 
      }
    }//! ic 
    if(ipad==1){
      l5->Draw();
      TPaveText *pt0 = new TPaveText(0.6130501,0.883784,0.9029564,0.9487681,"brNDC");
      pt0->SetBorderSize(0);
      pt0->SetFillColor(10);
      pt0->SetTextFont(62);
      TText *text0 = pt0->AddText("CMS Preliminary");
      text0->SetTextSize(0.05);
      pt0->Draw();
      
      TPaveText *pt1 = new TPaveText(0.553961,0.474384,0.8401741,0.5826909,"brNDC");
      pt1->SetBorderSize(0);
      pt1->SetFillColor(10);
      pt1->SetTextFont(42);
      TText *text1 = pt1->AddText("Anti-k_{T} (R=0.2)");
      text1->SetTextSize(0.055);
      TText *text2 = pt1->AddText("Particle Flow Jets");
      text2->SetTextSize(0.055);
      pt1->Draw();
    }else if(ipad==2){      
      TPaveText *pt2 = new TPaveText(0.553961,0.474384,0.8401741,0.5826909,"brNDC");
      pt2->SetBorderSize(0);
      pt2->SetFillColor(10);
      pt2->SetTextFont(42);
      TText *text3 = pt2->AddText("Anti-k_{T} (R=0.3)");
      text3->SetTextSize(0.055);
      TText *text4 = pt2->AddText("Particle Flow Jets");
      text4->SetTextSize(0.055);
      pt2->Draw();
    }else if(ipad==3){
      TPaveText *pt3 = new TPaveText(0.553961,0.474384,0.8401741,0.5826909,"brNDC");
      pt3->SetBorderSize(0);
      pt3->SetFillColor(10);
      pt3->SetTextFont(42);
      TText *text5 = pt3->AddText("Anti-k_{T} (R=0.4)");
      text5->SetTextSize(0.055);
      TText *text6 = pt3->AddText("Particle Flow Jets");
      text6->SetTextSize(0.055);
      pt3->Draw();
    }
  }      

  if(iSave){
    c5->SaveAs("AN/Efficiency/JetRecoEff_Phi.png");
    c5->SaveAs("AN/Efficiency/JetRecoEff_Phi.pdf");
    c5->SaveAs("AN/Efficiency/JetRecoEff_Phi.C"  );
    c5->SaveAs("AN/Efficiency/JetRecoEff_Phi.eps");
  }

  return 0;

}
