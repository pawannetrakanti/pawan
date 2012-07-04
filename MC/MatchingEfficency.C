#include <cstdlib>
#include <cmath>
#include <iostream>
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
#include <TColor.h>
#include "MultiCanvas.h"

using namespace std;

const double pi=acos(-1.);

/// set binomial errors for efficiency from two historgrams
void setBinomialErrors(TH1F* hEff,const TH1F* hEnum, const TH1F* hDenom);


int MatchingEfficency(int icen=0,const char *ksp="pp")
{
  const char *cdR="0.3";
  const double ketacut=2.0;
  double min=80.;
  double max=400.;
  bool iSave=false;

  //! Input files 
  TFile *fin     = 0;
  if(strcmp(ksp,"pp")==0)fin   = new TFile("input/Convolution/pp/Response_newHiForest_DJ_merged_pp_2012.root","r");  
  else fin   = new TFile("input/Convolution/pbpb/Response_newHiForest_DJ_merged_pbpb_2012.root","r");

  //if(strcmp(ksp,"pp")==0)fin = new TFile("input/NewCent/coarseptbins/withvtxcut/Response_newHiForest_DJ_merged_pp_2012.root","r");
  //else fin = new TFile("input/NewCent/coarseptbins/withvtxcut/Response_newHiForest_DJ_merged_pbpb_2012.root","r");


  std::cout<<"\t"<<std::endl
	   <<"Input file name   : "<<fin->GetName()<<"\t GenReco matching : "<<cdR<<std::endl
	   <<std::endl;
  //return 0;

  const int knj=7;
  const char *calgo[knj]= {
    "icPu5Calo", 
    "ak2PF","ak3PF","ak4PF",
    "akPu2PF","akPu3PF","akPu4PF"
  };                                     

  const char *ccent[7] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","pp"};
  


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
  gStyle->SetPadBorderMode(0);

  //! pt binning
  //double ptbins[] ={80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400};
  //const int bins  = sizeof(ptbins)/sizeof(Double_t) - 1;
  //const int nbins = bins;
  //std::cout<<"# of pt bins : " <<nbins<<std::endl;
  //std::cout<<std::endl;

  TH1F *hPtAll [knj], *hPtSel [knj];
  TH1F *hEtaAll[knj], *hEtaSel[knj];
  TH1F *hPhiAll[knj], *hPhiSel[knj];
  
  TH1F *hEffVsPt[knj], *hEffVsEta[knj], *hEffVsPhi[knj];

  for(int nj=0;nj<knj;nj++){
    std::cout<<"nj : "<<nj<<Form("\t %s",calgo[nj])<<std::endl;

    //int icen=0;
    
    //! pp /////////////////////////////
    hPtAll[nj] = (TH1F*)fin->Get(Form("hPtAll%d_%d",nj,icen));
    hPtAll[nj]->SetName(Form("hPtAll%s_%s",ksp,calgo[nj]));

    hEtaAll[nj] = (TH1F*)fin->Get(Form("hEtaAll%d_%d",nj,icen));
    hEtaAll[nj]->SetName(Form("hEtaAll%s_%s",ksp,calgo[nj]));

    hPhiAll[nj] = (TH1F*)fin->Get(Form("hPhiAll%d_%d",nj,icen));
    hPhiAll[nj]->SetName(Form("hPhiAll%s_%s",ksp,calgo[nj]));


    hPtSel[nj] = (TH1F*)fin->Get(Form("hPtSel%d_%d",nj,icen));
    hPtSel[nj]->SetName(Form("hPtSel%s_%s",ksp,calgo[nj]));

    hEtaSel[nj] = (TH1F*)fin->Get(Form("hEtaSel%d_%d",nj,icen));
    hEtaSel[nj]->SetName(Form("hEtaSel%s_%s",ksp,calgo[nj]));

    hPhiSel[nj] = (TH1F*)fin->Get(Form("hPhiSel%d_%d",nj,icen));
    hPhiSel[nj]->SetName(Form("hPhiSel%s_%s",ksp,calgo[nj]));

    //hEffVsPt[nj] = new TH1F(Form("hEffVsPtpp_%s",calgo[nj]),Form("hEffVsPtpp_%s",calgo[nj]),
    //			       hPtAll[nj]->GetNbinsX(),hPtAll[nj]->GetXaxis()->GetXmin(),hPtAll[nj]->GetXaxis()->GetXmax());

    hEffVsPt[nj] = new TH1F(Form("hEffVsPt%s_%s",ksp,calgo[nj]),Form("hEffVsPt%s_%s",ksp,calgo[nj]),
			    hPtSel[nj]->GetNbinsX(),hPtSel[nj]->GetXaxis()->GetXmin(),hPtSel[nj]->GetXaxis()->GetXmax());
    hEffVsPt[nj]->Sumw2();
    hEffVsPt[nj]->Divide(hPtSel[nj],hPtAll[nj],1.0,1.0,"B");
    //hEffVsPt[nj]->Divide(hPtSel[nj],hPtAll[nj],1.0,1.0);
    //setBinomialErrors(hEffVsPt[nj],hPtSel[nj],hPtAll[nj]);
    hEffVsPt[nj]->SetTitle("");
    hEffVsPt[nj]->GetXaxis()->SetRangeUser(min,max);
    hEffVsPt[nj]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    hEffVsPt[nj]->GetXaxis()->CenterTitle(true);
    hEffVsPt[nj]->GetXaxis()->SetTitleFont(42);
    hEffVsPt[nj]->GetXaxis()->SetLabelFont(42);
    hEffVsPt[nj]->GetXaxis()->SetMoreLogLabels();
    hEffVsPt[nj]->GetXaxis()->SetNoExponent();
    hEffVsPt[nj]->GetXaxis()->SetTitleSize(0.06);
    hEffVsPt[nj]->GetXaxis()->SetTitleOffset(1.15);
    hEffVsPt[nj]->GetXaxis()->SetLabelSize(0.06);
    hEffVsPt[nj]->GetXaxis()->SetNdivisions(507);
    hEffVsPt[nj]->GetYaxis()->SetTitle("Matching Efficiency");    
    hEffVsPt[nj]->GetYaxis()->SetTitleSize(0.06);
    hEffVsPt[nj]->GetYaxis()->SetTitleOffset(1.18);
    hEffVsPt[nj]->GetYaxis()->SetLabelSize(0.06);
    //hEffVsPt[nj]->GetYaxis()->SetNdivisions(507);
    hEffVsPt[nj]->GetYaxis()->SetDecimals(true);
    hEffVsPt[nj]->GetYaxis()->SetTitleFont(42);
    hEffVsPt[nj]->GetYaxis()->SetLabelFont(42);
    //hEffVsPt[nj]->SetMaximum(1.153);
    //hEffVsPt[nj]->SetMinimum(0.735);

    hEffVsPt[nj]->GetYaxis()->SetNdivisions(507);    
    hEffVsPt[nj]->SetMaximum(1.035);
    hEffVsPt[nj]->SetMinimum(0.835);


    hEffVsEta[nj] = new TH1F(Form("hEffVsEta%s_%s",ksp,calgo[nj]),Form("hEffVsEta%s_%s",ksp,calgo[nj]),
				hEtaAll[nj]->GetNbinsX(),hEtaAll[nj]->GetXaxis()->GetXmin(),hEtaAll[nj]->GetXaxis()->GetXmax());
    hEffVsEta[nj]->Sumw2();
    hEffVsEta[nj]->Divide(hEtaSel[nj],hEtaAll[nj],1.0,1.0,"B");
    //hEffVsEta[nj]->Divide(hEtaSel[nj],hEtaAll[nj],1.0,1.0);
    //setBinomialErrors(hEffVsEta[nj],hEtaSel[nj],hEtaAll[nj]);


    hEffVsEta[nj]->SetTitle("");
    hEffVsEta[nj]->GetXaxis()->SetRangeUser(-ketacut,ketacut);
    hEffVsEta[nj]->GetXaxis()->SetTitle("GenJet #eta ");
    hEffVsEta[nj]->GetXaxis()->CenterTitle(true);
    hEffVsEta[nj]->GetXaxis()->SetTitleFont(42);
    hEffVsEta[nj]->GetXaxis()->SetLabelFont(42);
    hEffVsEta[nj]->GetXaxis()->SetMoreLogLabels();
    hEffVsEta[nj]->GetXaxis()->SetNoExponent();
    hEffVsEta[nj]->GetXaxis()->SetTitleSize(0.06);
    hEffVsEta[nj]->GetXaxis()->SetTitleOffset(0.98);
    hEffVsEta[nj]->GetXaxis()->SetLabelSize(0.06);
    hEffVsEta[nj]->GetXaxis()->SetNdivisions(507);
    hEffVsEta[nj]->GetYaxis()->SetTitle("Matching Efficency");    
    hEffVsEta[nj]->GetYaxis()->SetTitleSize(0.06);
    hEffVsEta[nj]->GetYaxis()->SetTitleOffset(1.07);
    hEffVsEta[nj]->GetYaxis()->SetLabelSize(0.06);
    hEffVsEta[nj]->GetYaxis()->SetNdivisions(507);
    hEffVsEta[nj]->GetYaxis()->SetTitleFont(42);
    hEffVsEta[nj]->GetYaxis()->SetLabelFont(42);
    hEffVsEta[nj]->GetYaxis()->SetDecimals(true);

    hEffVsPhi[nj] = new TH1F(Form("hEffVsPhi%s_%s",ksp,calgo[nj]),Form("hEffVsPhi%s_%s",ksp,calgo[nj]),
				hPhiAll[nj]->GetNbinsX(),hPhiAll[nj]->GetXaxis()->GetXmin(),hPhiAll[nj]->GetXaxis()->GetXmax());
    hEffVsPhi[nj]->Sumw2();
    hEffVsPhi[nj]->Divide(hPhiSel[nj],hPhiAll[nj],1.0,1.0,"B");
    //hEffVsPhi[nj]->Divide(hPhiSel[nj],hPhiAll[nj],1.0,1.0);
    //setBinomialErrors(hEffVsPhi[nj],hPhiSel[nj],hPhiAll[nj]);


    hEffVsPhi[nj]->SetTitle("");
    hEffVsPhi[nj]->GetXaxis()->SetRangeUser(-pi,pi);
    hEffVsPhi[nj]->GetXaxis()->SetTitle("GenJet #phi (radian) ");
    hEffVsPhi[nj]->GetXaxis()->CenterTitle(true);
    hEffVsPhi[nj]->GetXaxis()->SetTitleFont(42);
    hEffVsPhi[nj]->GetXaxis()->SetLabelFont(42);
    hEffVsPhi[nj]->GetXaxis()->SetMoreLogLabels();
    hEffVsPhi[nj]->GetXaxis()->SetNoExponent();
    hEffVsPhi[nj]->GetXaxis()->SetTitleSize(0.06);
    hEffVsPhi[nj]->GetXaxis()->SetTitleOffset(1.15);
    hEffVsPhi[nj]->GetXaxis()->SetLabelSize(0.06);
    hEffVsPhi[nj]->GetXaxis()->SetNdivisions(507);
    hEffVsPhi[nj]->GetYaxis()->SetTitle("Matching Efficency");    
    hEffVsPhi[nj]->GetYaxis()->SetTitleSize(0.06);
    hEffVsPhi[nj]->GetYaxis()->SetTitleOffset(1.12);
    hEffVsPhi[nj]->GetYaxis()->SetLabelSize(0.06);
    hEffVsPhi[nj]->GetYaxis()->SetNdivisions(507);
    hEffVsPhi[nj]->GetYaxis()->SetTitleFont(42);
    hEffVsPhi[nj]->GetYaxis()->SetLabelFont(42);
    hEffVsPhi[nj]->GetYaxis()->SetDecimals(true);


    if(strcmp(ksp,"pp")==0 || icen>3){
	hEffVsPt[nj]->SetMaximum(1.015);
	hEffVsPt[nj]->SetMinimum(0.935);
	hEffVsEta[nj]->SetMaximum(1.015);
	hEffVsEta[nj]->SetMinimum(0.935);
	hEffVsPhi[nj]->SetMaximum(1.015);
	hEffVsPhi[nj]->SetMinimum(0.935);
      }else{
	
	hEffVsPt[nj]->SetMaximum(1.053);
	hEffVsPt[nj]->SetMinimum(0.535);
	hEffVsEta[nj]->SetMaximum(1.053);
	hEffVsEta[nj]->SetMinimum(0.535);
	hEffVsPhi[nj]->SetMaximum(1.053);
	hEffVsPhi[nj]->SetMinimum(0.535);

	if(icen>1){
	  hEffVsPt[nj]->SetMinimum(0.735);
	  hEffVsEta[nj]->SetMinimum(0.735);
	  hEffVsPhi[nj]->SetMinimum(0.735);
	}
      }

  }

  TFile *fout = new TFile(Form("rootfiles/MatchingEfficiency/Histo_%s_%s.root",ksp,ccent[icen]),"RECREATE");

  fout->cd();
  for(int nj=0;nj<knj;nj++){
    hEffVsPt[nj]->Write();
    hEffVsEta[nj]->Write();
    hEffVsPhi[nj]->Write();
  }
  fout->Close();
  return 0;

  int sty[6] = {24,25,26,28,30,32};
  //int col[6] = { 2, 4, 6, 8,48,44};
  Color_t col[7] = {kRed,kBlue,kMagenta,kCyan,kGreen+1,kBlue+4,kRed+4};

  int ipad=0;
  const int maxc=4;
  const int maxr=3;
  
  TLegend *l3[maxc];
  int ij=-1;
  int ic=0;

  //TCanvas *c3 = new TCanvas("c3","EffVsPt pp ",105,201,1208,452);
  //makeMultiPanelCanvas(c3,maxc,maxr,0.0,0.0,0.22,0.22,0.02);

  TCanvas *c3 = new TCanvas("c3", "Matching Eff pp ",114,69,1059,856);
  c3->Divide(maxc,maxr,0,0);

  TLine *line1 = new TLine(min,1.0,max,1.0);
  line1->SetLineWidth(2);
  line1->SetLineStyle(2);

  TLine *line2 = new TLine(-ketacut,1.0,ketacut,1.0);
  line2->SetLineWidth(2);
  line2->SetLineStyle(2);

  TLine *line3 = new TLine(-pi,1.0,pi,1.0);
  line3->SetLineWidth(2);
  line3->SetLineStyle(2);

  for(int nj=1;nj<knj;nj++){

    if((nj-1)%3==0){
      c3->cd(++ipad);
      if((ipad-1)%maxc==0){
	gPad->SetLeftMargin(0.14);
      }else if(ipad%maxc==0){
	gPad->SetRightMargin(0.02);
      }

      gPad->SetBottomMargin(0.15);
      gPad->SetLogx();
      
      ic=0;

      l3[++ij] = new TLegend(0.4025777,0.1958268,0.6523522,0.571499,NULL,"BRNDC");
      l3[ij]->SetHeader("");
      l3[ij]->SetBorderSize(0);
      l3[ij]->SetTextFont(42);
      l3[ij]->SetTextSize(0.06);
      l3[ij]->SetLineColor(1);
      l3[ij]->SetLineStyle(1);
      l3[ij]->SetLineWidth(1);
      l3[ij]->SetFillColor(10);
      l3[ij]->SetFillStyle(1001);
      l3[ij]->SetHeader("");			  

      hEffVsPt[nj]->SetMarkerStyle(sty[ic]);
      hEffVsPt[nj]->SetMarkerColor(col[ic]);
      hEffVsPt[nj]->SetLineColor(col[ic]);
      hEffVsPt[nj]->SetMarkerSize(1.0);
      hEffVsPt[nj]->Draw("p");
      l3[ij]->AddEntry(hEffVsPt[nj] ,calgo[nj],"p"); 
      l3[ij]->Draw();
      line1->Draw();

      c3->cd(ipad+maxc);
      if((ipad-1)% (ipad+maxc)==0){
	gPad->SetLeftMargin(0.14);
      }else if(ipad%4==0){
	gPad->SetRightMargin(0.02);
      }

      gPad->SetBottomMargin(0.15);
      hEffVsEta[nj]->SetMarkerStyle(sty[ic]);
      hEffVsEta[nj]->SetMarkerColor(col[ic]);
      hEffVsEta[nj]->SetLineColor(col[ic]);
      hEffVsEta[nj]->SetMarkerSize(1.0);
      hEffVsEta[nj]->Draw("p");
      line2->Draw();

      c3->cd(ipad+2*maxc);
      if((ipad-1)% (ipad+2*maxc)==0){
	gPad->SetLeftMargin(0.14);
      }else if(ipad%4==0){
	gPad->SetRightMargin(0.02);
      }
      gPad->SetBottomMargin(0.25);
      hEffVsPhi[nj]->SetMarkerStyle(sty[ic]);
      hEffVsPhi[nj]->SetMarkerColor(col[ic]);
      hEffVsPhi[nj]->SetLineColor(col[ic]);
      hEffVsPhi[nj]->SetMarkerSize(1.0);
      hEffVsPhi[nj]->Draw("p");
      line3->Draw();



      if(ipad+2*maxc>9){
	hEffVsPhi[nj]->GetXaxis()->SetTitleSize(0.07);
	hEffVsPhi[nj]->GetXaxis()->SetTitleOffset(1.03);
	hEffVsPhi[nj]->GetXaxis()->SetLabelSize(0.07);
	hEffVsPhi[nj]->GetXaxis()->SetLabelOffset(0.004);
      }

      ic++;

    }else{

      c3->cd(ipad);
      gPad->SetLogx();
      hEffVsPt[nj]->SetMarkerStyle(sty[ic]);
      hEffVsPt[nj]->SetMarkerColor(col[ic]);
      hEffVsPt[nj]->SetLineColor(col[ic]);
      hEffVsPt[nj]->SetMarkerSize(1.0);
      hEffVsPt[nj]->Draw("psame"); 
      l3[ij]->AddEntry(hEffVsPt[nj] ,calgo[nj],"p"); 
      line1->Draw();

      c3->cd(ipad+maxc);
      hEffVsEta[nj]->SetMarkerStyle(sty[ic]);
      hEffVsEta[nj]->SetMarkerColor(col[ic]);
      hEffVsEta[nj]->SetLineColor(col[ic]);
      hEffVsEta[nj]->SetMarkerSize(1.0);
      hEffVsEta[nj]->Draw("psame"); 
      line2->Draw();

      c3->cd(ipad+2*maxc);

      hEffVsPhi[nj]->SetMarkerStyle(sty[ic]);
      hEffVsPhi[nj]->SetMarkerColor(col[ic]);
      hEffVsPhi[nj]->SetLineColor(col[ic]);
      hEffVsPhi[nj]->SetMarkerSize(1.0);
      hEffVsPhi[nj]->Draw("psame"); 
      line3->Draw();

      ic++;
    }
  }

  /*  
  c3->cd(maxc);
  gPad->SetLogx();
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);

  hEffVsPt[0]->SetMarkerStyle(20);
  hEffVsPt[0]->SetMarkerColor(1);
  hEffVsPt[0]->SetLineColor(1);
  hEffVsPt[0]->SetMarkerSize(1.0);
  hEffVsPt[0]->Draw("p");

  //l3[++ij] = new TLegend(0.1792499,0.1958268,0.4290245,0.571499,NULL,"BRNDC");
  l3[++ij] = new TLegend(0.4025777,0.1958268,0.6523522,0.571499,NULL,"BRNDC");
  l3[ij]->SetHeader("");
  l3[ij]->SetBorderSize(0);
  l3[ij]->SetTextFont(42);
  l3[ij]->SetTextSize(0.06);
  l3[ij]->SetLineColor(1);
  l3[ij]->SetLineStyle(1);
  l3[ij]->SetLineWidth(1);
  l3[ij]->SetFillColor(10);
  l3[ij]->SetFillStyle(1001);
  l3[ij]->SetHeader("");		
  l3[ij]->AddEntry(hEffVsPt[0] ,calgo[0],"p"); 
  l3[ij]->Draw();
  line1->Draw();
  */
  c3->cd(1);
  TPaveText *pt3 = new TPaveText(0.3849466,0.571499,0.6729219,0.7106368,"brNDC");
  pt3->SetBorderSize(0);
  pt3->SetFillColor(10);
  pt3->SetTextFont(42);
  TText *text3 = pt3->AddText("CMS Simulation");
  text3->SetTextSize(0.06);
  TText *text4 = 0;
  if(strcmp(ksp,"pp")==0)text4=pt3->AddText("PYTHIA Z2");
  else text4=pt3->AddText("PYTHIA Z2 + HYDJET 1.8");
  text4->SetTextSize(0.06);
  pt3->Draw();

  c3->cd(2);
  TPaveText *pt4 = new TPaveText(0.402332,0.5436714,0.6885664,0.7210721,"brNDC");
  pt4->SetBorderSize(0);
  pt4->SetFillColor(10);
  pt4->SetTextFont(42);
  TText *text5 = pt4->AddText(Form("|#eta|<%0.0f, dr < %s",ketacut,cdR));
  text5->SetTextSize(0.07);
  TText *text6 = 0;
  if(strcmp(ksp,"pp")==0)text6=pt4->AddText("pp");
  else text6=pt4->AddText(Form("PbPb %s",ccent[icen]));
  text6->SetTextSize(0.07);
  pt4->Draw();

  /*
  c3->cd(2*maxc);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  hEffVsEta[0]->SetMarkerStyle(20);
  hEffVsEta[0]->SetMarkerColor(1);
  hEffVsEta[0]->SetLineColor(1);
  hEffVsEta[0]->SetMarkerSize(1.0);
  hEffVsEta[0]->Draw("p");
  line2->Draw();

  c3->cd(3*maxc);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.25);
  hEffVsPhi[0]->SetMarkerStyle(20);
  hEffVsPhi[0]->SetMarkerColor(1);
  hEffVsPhi[0]->SetLineColor(1);
  hEffVsPhi[0]->SetMarkerSize(1.0);
  hEffVsPhi[0]->Draw("p");
  line3->Draw();
  */

  if(iSave){
    if(strcmp(ksp,"pbpb")==0){
      c3->SaveAs(Form("AN/MatchingEfficiency_%s_%s.png",ksp,ccent[icen]));
      c3->SaveAs(Form("AN/MatchingEfficiency_%s_%s.pdf",ksp,ccent[icen]));
      c3->SaveAs(Form("AN/MatchingEfficiency_%s_%s.C",ksp,ccent[icen]));
      c3->SaveAs(Form("AN/MatchingEfficiency_%s_%s.eps",ksp,ccent[icen]));
    }else{
      c3->SaveAs(Form("AN/MatchingEfficiency_%s.png",ksp));
      c3->SaveAs(Form("AN/MatchingEfficiency_%s.pdf",ksp));
      c3->SaveAs(Form("AN/MatchingEfficiency_%s.C",ksp));
      c3->SaveAs(Form("AN/MatchingEfficiency_%s.eps",ksp));
    }
  }


  return 0;
}
void setBinomialErrors(TH1F* hEff,const TH1F* hEnum, const TH1F* hDenom)
{
  for (int i=1;i<=hEff->GetNbinsX();i++) {
    float nenum =hEnum ->GetBinContent(i);
    float ndenom=hDenom->GetBinContent(i);
    float eeff=(ndenom>0.0) ? sqrt(nenum/(ndenom*ndenom)*(1-nenum/ndenom)):0.0;
    hEff->SetBinError(i,eeff);
  }
}
