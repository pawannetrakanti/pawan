#include <string.h>
#include <iostream>
#include <TH3.h>
#include <TF1.h>

#include "utils/utilities_v2.h"
#include "utils/commonStyle_v2.h"
#include "utils/MultiCanvas.h"
using namespace std;
const double pi = TMath::Pi();
int plot()
{
  LoadStyle();

  bool iSave=true;
  const int nmult=6;
  const char *chf[nmult] = {"(0-100%)",
			    "0 < E_{T}^{HF}[|#eta| > 4]  #leq 20",
			    "20 < E_{T}^{HF}[|#eta| > 4] #leq 25",
			    "25 < E_{T}^{HF}[|#eta| > 4] #leq 30",
			    "30 < E_{T}^{HF}[|#eta| > 4] #leq 40",
			    "40 #leq E_{T}^{HF}[|#eta| > 4]"};

  const char *mhf[nmult] = {"0-100",
			    "0-20",
			    "20-25",
			    "25-30",
			    "30-40",
			    ">40"};
  

  TFile *fin[4];
  
  //! Kurt
  fin[0] = new TFile("input/Kurt/hist_output_etaLT3_akPu3PF_ForPawan.root","r");

  //! Pawan
  fin[1] = new TFile("input/Pawan/Pawan.root","r");

  //! Ruchi
  //fin[2] = new TFile("input/Ruchi/","r");

  //! Doga
  //fin[3] = new TFile("input/Doga/","r");
  const char *ctype[5] ={"eta1","eta2","dijeteta","pt2bypt1","delphi"};



  //! fit function for dphi (from gamma-jet paper)
  TF1 *func = new TF1("func","[0]*((TMath::Exp((x - TMath::Pi())/[1]))/((1.0-TMath::Exp(-TMath::Pi()/[1]))*[1]))",0,pi);
  func->SetParameters(0.01,0.2);
  func->SetNpx(2000);

  TH1D *heta1  [4][nmult];
  TH1D *heta2  [4][nmult];
  TH1D *hdjeta [4][nmult];
  TH1D *hpt2pt1[4][nmult];
  TH1D *hdphi  [4][nmult];

  TH1D *haveta1  [4];
  TH1D *haveta2  [4];
  TH1D *havdjeta [4];
  TH1D *havpt2pt1[4];
  TH1D *havdphi  [4];

  Color_t icol[4] = {kBlack,kRed,kBlue,kGreen+3};
  int     isty[4] = {24    ,25  ,28   ,30};



  for(int i=0; i<4;i ++){
    haveta1[i] = new TH1D(Form("haveta1_%d",i),"<eta_{1}>",nmult,-0.5,nmult-0.5);
    haveta1[i]->Sumw2();
    haveta1[i]->SetLineColor(icol[i]);
    haveta1[i]->SetMarkerColor(icol[i]);
    haveta1[i]->SetMarkerStyle(isty[i]);
    haveta1[i]->SetMarkerSize(1.2);
    haveta1[i]->GetXaxis()->SetLabelSize(0.03);

    haveta2[i] = new TH1D(Form("haveta2_%d",i),"<eta_{2}>",nmult,-0.5,nmult-0.5);
    haveta2[i]->Sumw2();
    haveta2[i]->SetLineColor(icol[i]);
    haveta2[i]->SetMarkerColor(icol[i]);
    haveta2[i]->SetMarkerStyle(isty[i]);
    haveta2[i]->SetMarkerSize(1.2);
    haveta2[i]->GetXaxis()->SetLabelSize(0.03);


    havdjeta[i] = new TH1D(Form("havdjeta_%d",i),"< (eta_{1} + eta_{2})/2 >",nmult,-0.5,nmult-0.5);
    havdjeta[i]->Sumw2();
    havdjeta[i]->SetLineColor(icol[i]);
    havdjeta[i]->SetMarkerColor(icol[i]);
    havdjeta[i]->SetMarkerStyle(isty[i]);
    havdjeta[i]->SetMarkerSize(1.2);
    //havdjeta[i]->GetXaxis()->SetLabelSize(0.03);

    havpt2pt1[i] = new TH1D(Form("havpt2pt1_%d",i),"< p_{T,2} / p_{T,1} >",nmult,-0.5,nmult-0.5);
    havpt2pt1[i]->Sumw2();
    havpt2pt1[i]->SetLineColor(icol[i]);
    havpt2pt1[i]->SetMarkerColor(icol[i]);
    havpt2pt1[i]->SetMarkerStyle(isty[i]);
    havpt2pt1[i]->SetMarkerSize(1.2);
    havpt2pt1[i]->GetXaxis()->SetLabelSize(0.03);

    havdphi[i] = new TH1D(Form("havdphi_%d",i),"#sigma(  #Delta#phi_{12} )",nmult,-0.5,nmult-0.5);
    havdphi[i]->Sumw2();
    havdphi[i]->SetLineColor(icol[i]);
    havdphi[i]->SetMarkerColor(icol[i]);
    havdphi[i]->SetMarkerStyle(isty[i]);
    havdphi[i]->SetMarkerSize(1.2);
    //havdphi[i]->GetXaxis()->SetLabelSize(0.03);

    for(int j=1;j<=nmult;j++){
      haveta1  [i]->GetXaxis()->SetBinLabel(j,mhf[j-1]);
      haveta2  [i]->GetXaxis()->SetBinLabel(j,mhf[j-1]);
      havdjeta [i]->GetXaxis()->SetBinLabel(j,mhf[j-1]);
      havpt2pt1[i]->GetXaxis()->SetBinLabel(j,mhf[j-1]);
      havdphi  [i]->GetXaxis()->SetBinLabel(j,mhf[j-1]);
    }

  }


  //! Read Kurt's file  
  for(int i=0; i<nmult;i++){ //! 0 inclusive and 5 central
    heta1  [0][i] = (TH1D*)fin[0]->Get(Form("eta_lead%d",i));
    //heta1  [0][i]->Scale(1./ heta1  [0][i]->Integral());
    heta2  [0][i] = (TH1D*)fin[0]->Get(Form("eta_sublead%d",i));
    //heta2  [0][i]->Scale(1./ heta2  [0][i]->Integral());
    hpt2pt1[0][i] = (TH1D*)fin[0]->Get(Form("mult%d_0",i));
    //hpt2pt1[0][i]->Scale(1./ hpt2pt1[0][i]->Integral());
    hdphi  [0][i] = (TH1D*)fin[0]->Get(Form("delphi%d_0",i));
    //hdphi  [0][i]->Scale(1./hdphi[0][i]->Integral());
    hdjeta [0][i] = (TH1D*)fin[0]->Get(Form("etaDistropPb_%d",i));
    hdjeta [0][i]->Scale(1./hdjeta[0][i]->Integral());
  }

  //! Read Pawan's file  
  int mbins[nmult-1] = {5,4,3,2,1};
  for(int i=1; i<nmult;i++){ //! 1 peripheral 5 central
    //cout<<i<<"\t"<<mbins[i-1]<<endl;
    heta1  [1][i] = (TH1D*)fin[1]->Get(Form("eta_data_akPu3PF_lead_%d",mbins[i-1]));
    heta2  [1][i] = (TH1D*)fin[1]->Get(Form("eta_data_akPu3PF_sublead_%d",mbins[i-1]));
    hpt2pt1[1][i] = (TH1D*)fin[1]->Get(Form("pt2pt1_data_akPu3PF_%d",mbins[i-1]));
    hdphi  [1][i] = (TH1D*)fin[1]->Get(Form("dphi_data_akPu3PF_%d",mbins[i-1]));
    hdjeta [1][i] = (TH1D*)fin[1]->Get(Form("djeta_data_akPu3PF_%d",mbins[i-1]));
  }
  heta1  [1][0] = (TH1D*)fin[1]->Get(Form("eta_data_akPu3PF_lead_%d"   ,0));
  heta2  [1][0] = (TH1D*)fin[1]->Get(Form("eta_data_akPu3PF_sublead_%d",0));
  hpt2pt1[1][0] = (TH1D*)fin[1]->Get(Form("pt2pt1_data_akPu3PF_%d"     ,0));
  hdphi  [1][0] = (TH1D*)fin[1]->Get(Form("dphi_data_akPu3PF_%d"       ,0));
  hdjeta [1][0] = (TH1D*)fin[1]->Get(Form("djeta_data_akPu3PF_%d"      ,0));

  /*
  //! Get inclusive 0
  heta1   [1][0] = (TH1D*)heta1  [1][1]->Clone("eta_data_akPu3PF_lead_0");
  heta2   [1][0] = (TH1D*)heta2  [1][1]->Clone("eta_data_akPu3PF_sublead_0");
  hpt2pt1 [1][0] = (TH1D*)hpt2pt1[1][1]->Clone("pt2pt1_data_akPu3PF_0");
  hdphi   [1][0] = (TH1D*)hdphi  [1][1]->Clone("dphi_data_akPu3PF_0");
  hdjeta  [1][0] = (TH1D*)hdjeta [1][1]->Clone("djeta_data_akPu3PF_0");

  for(int i=5;i<=2;i--){
    heta1   [1][0]->Add(heta1  [1][i]);
    heta2   [1][0]->Add(heta2  [1][i]);
    hpt2pt1 [1][0]->Add(hpt2pt1[1][i]);
    hdphi   [1][0]->Add(hdphi  [1][i]);
    hdjeta  [1][0]->Add(hdjeta [1][i]);
  }
  */
  //hdjeta[1][0]->Draw("p");
  //for(int j=1;j<nmult;j++)  heta2[1][j]->Draw("psame");
  //return 0;

  //! To read  Ruchi's files:
  int ist[5] = {21,22,1,0,3};
  for(int i=0;i<5;i++){
    for(int j=0;j<nmult;j++){
      fin[2] = new TFile(Form("input/Ruchi/rootfiles/%s/Histo_%d_%d.root",ctype[i],ist[i],(nmult-1)-j),"r");
      if     (i==0)heta1  [2][j] = (TH1D*)fin[2]->Get("h");
      else if(i==1)heta2  [2][j] = (TH1D*)fin[2]->Get("h");
      else if(i==2)hdjeta [2][j] = (TH1D*)fin[2]->Get("h");
      else if(i==3)hpt2pt1[2][j] = (TH1D*)fin[2]->Get("h");
      else if(i==4)hdphi  [2][j] = (TH1D*)fin[2]->Get("h");
    }
  }


  //! To read Doga's files:
  const char *cdoga[5] ={"leadetainhfbins.root","subetainhfbins.root","dijetetainhfbins.root","xinhfbins.root","dphiinhfbins.root"};
  const char *cname[5] ={"heta_data_mult"      ,"heta_data_mult"     ,"heta_data_mult"       ,"hx_data_mult"  ,"hdphi_data_mult"};
  const int min[nmult] = {0 ,0 ,20,25,30,40};  
  const int max[nmult] = {90,20,25,30,40,90};  
  for(int i=0;i<5;i++){
    fin[3] = new TFile(Form("input/Doga/%s",cdoga[i]),"r");
    for(int j=0;j<nmult;j++){
      if     (i==0)heta1  [3][j] = (TH1D*)fin[3]->Get(Form("%s%d_%d",cname[i],min[j],max[j]));
      else if(i==1)heta2  [3][j] = (TH1D*)fin[3]->Get(Form("%s%d_%d",cname[i],min[j],max[j]));
      else if(i==2)hdjeta [3][j] = (TH1D*)fin[3]->Get(Form("%s%d_%d",cname[i],min[j],max[j]));
      else if(i==3)hpt2pt1[3][j] = (TH1D*)fin[3]->Get(Form("%s%d_%d",cname[i],min[j],max[j]));
      else if(i==4)hdphi  [3][j] = (TH1D*)fin[3]->Get(Form("%s%d_%d",cname[i],min[j],max[j]));      
    }
  }


  for(int i=0;i<4;i++){
    cout<<i<<endl;
    //if(i==1)continue;
    for(int j=0;j<nmult;j++){
      //cout<<"i : "<<i<<"\t j : "<<j<<"\t eta 1 : "<<heta1  [i][j]->GetMean()<<endl;
      haveta1  [i]->SetBinContent(j+1,heta1  [i][j]->GetMean());
      haveta1  [i]->SetBinError  (j+1,heta1  [i][j]->GetMeanError());

      haveta2  [i]->SetBinContent(j+1,heta2  [i][j]->GetMean());
      haveta2  [i]->SetBinError  (j+1,heta2  [i][j]->GetMeanError());
      havdjeta [i]->SetBinContent(j+1,hdjeta [i][j]->GetMean());
      havdjeta [i]->SetBinError  (j+1,hdjeta [i][j]->GetMeanError());

      havpt2pt1[i]->SetBinContent(j+1,hpt2pt1[i][j]->GetMean());
      havpt2pt1[i]->SetBinError  (j+1,hpt2pt1[i][j]->GetMeanError());

      //! avdphi fit
      func->SetParameters(0.01,0.2);
      func->SetParLimits(0,0,1);
      func->SetParLimits(2,0.1,0.5);
      hdphi  [i][j]->Fit(func,"RQN","",2*pi/3.,pi);
      //cout<<"i : "<<i<<"\t j : "<<j<<" sigma : "<<func->GetParameter(1)<<"\t"<<func->GetParError(1)<<endl;
      havdphi[i]->SetBinContent(j+1,func->GetParameter(1));
      havdphi[i]->SetBinError  (j+1,func->GetParError(1));
    }
  }


  cout<<"eta1 : Kurt : "<<heta1[0][0]->GetNbinsX()<<"\t Pawan : "<<heta1[1][0]->GetNbinsX()<<"\t Ruchi  : "<<heta1[2][0]->GetNbinsX()<<"\t Doga : "<<heta1[3][0]->GetNbinsX()<<endl;
  cout<<"eta2 : Kurt : "<<heta2[0][0]->GetNbinsX()<<"\t Pawan : "<<heta2[1][0]->GetNbinsX()<<"\t Ruchi  : "<<heta2[2][0]->GetNbinsX()<<"\t Doga : "<<heta2[3][0]->GetNbinsX()<<endl;
  cout<<"ratio : Kurt : "<<hpt2pt1[0][0]->GetNbinsX()<<"\t Pawan : "<<hpt2pt1[1][0]->GetNbinsX()<<"\t Ruchi  : "<<hpt2pt1[2][0]->GetNbinsX()<<"\t Doga : "<<hpt2pt1[3][0]->GetNbinsX()<<endl;
  cout<<"dphi : Kurt : "<<hdphi[0][0]->GetNbinsX()<<"\t Pawan : "<<hdphi[1][0]->GetNbinsX()<<"\t Ruchi  : "<<hdphi[2][0]->GetNbinsX()<<"\t Doga : "<<hdphi[3][0]->GetNbinsX()<<endl;
  //return 0;

  //hpt2pt1[2][0]->Draw("p");
  //return 0;

  int ipad=0;
  const char *cpad[] = {"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)"};
  const float ketacut=3.0;

  TLegend *leg1 = new TLegend(0.7433258,0.2840909,0.9162209,0.5813342,"BRNDC");
  leg1->SetHeader("");
  leg1->SetFillStyle(0);
  leg1->SetFillColor(10);
  leg1->SetLineColor(10);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.04);
  leg1->AddEntry(haveta1 [0],"Kurt","p");
  leg1->AddEntry(haveta1 [1],"Pawan","p");
  leg1->AddEntry(haveta1 [2],"Ruchi","p");
  leg1->AddEntry(haveta1 [3],"Doga","p");

  TLegend *leg2 = new TLegend(0.2096901,0.2298806,0.6326831,0.5269267,"BRNDC");
  leg2->SetHeader("");
  leg2->SetFillStyle(0);
  leg2->SetFillColor(10);
  leg2->SetLineColor(10);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(heta1 [0][0],"Kurt","p");
  leg2->AddEntry(heta1 [1][0],"Pawan","p");
  leg2->AddEntry(heta1 [2][0],"Ruchi","p");
  leg2->AddEntry(heta1 [3][0],"Doga","p");

  TCanvas *c85a  = new TCanvas("c85a","<Lead Jet eta>",1165,424);
  makeMultiPanelCanvasWithGap(c85a,3,1,0.05,0.01,0.16,0.2,0.05,0.01); 
  //c85a->Divide(3,1,0,0);

  c85a->cd(++ipad);
  haveta1[0]->SetMaximum(0.02);
  haveta1[0]->SetMinimum(-0.53);
  MakeHist(haveta1[0],"E_{T}^{HF[|#eta|<4]}","< #eta >");
  haveta1[0]->GetYaxis()->SetTitleOffset(1.08);
  haveta1 [0]->Draw("p");
  haveta1 [1]->Draw("psame");
  haveta1 [2]->Draw("psame");
  haveta1 [3]->Draw("psame");
  leg1->Draw();
  drawText2("< #eta_{1} >",0.23,0.82,26,kBlack);
  drawText("pPb #sqrt{s_{NN}} = 5.02 TeV",0.33,0.88,22);
  drawText2("p_{T,1} > 120 GeV/c",0.30,0.75,19,kBlack);
  drawText2("p_{T,2} >  30 GeV/c",0.30,0.68,19,kBlack);
  drawText2("|#eta| < 3, #Delta#phi_{12} #geq #frac{2#pi}{3}",0.30,0.60,19,kBlack);

  //if(ipad==1){
  // drawText("CMS Preliminary",0.33,0.88,22);
  //drawText("ak3PF",0.20,0.80,21);
  // drawText2("< #eta_{1} >",0.53,0.80,26,kBlack);
  //}

  c85a->cd(++ipad);
  haveta2[0]->SetMaximum(0.02);
  haveta2[0]->SetMinimum(-0.53);
  MakeHist(haveta2[0],"E_{T}^{HF[|#eta|<4]}","< #eta >");
  haveta2[0]->GetYaxis()->SetTitleOffset(1.08);
  haveta2 [0]->Draw("p");
  haveta2 [1]->Draw("psame");
  haveta2 [2]->Draw("psame");
  haveta2 [3]->Draw("psame");
  drawText2("< #eta_{2} >",0.23,0.82,26,kBlack);

  c85a->cd(++ipad);
  havdjeta[0]->SetMaximum(0.02);
  havdjeta[0]->SetMinimum(-0.53);
  havdjeta[0]->GetYaxis()->SetTitleOffset(1.08);
  MakeHist(havdjeta[0],"E_{T}^{HF[|#eta|<4]}","< (#eta_{1} + #eta_{2})/2 >");
  havdjeta[0]->Draw("p");
  havdjeta[1]->Draw("psame");
  havdjeta[2]->Draw("psame");
  havdjeta[3]->Draw("psame");
  drawText2("< (#eta_{1} + <#eta_{2})/2 >",0.23,0.82,26,kBlack);

  if(iSave){
    c85a->SaveAs("giffiles/Aveta.gif");
    c85a->SaveAs("giffiles/Aveta.pdf");
    c85a->SaveAs("giffiles/Aveta.C");
    c85a->SaveAs("giffiles/Aveta.eps");
  }


  ipad=0;
  TCanvas *c85b  = new TCanvas("c85b","pt2/pT1 and dphi",968,466);
  c85b->Divide(2,1,0,0);

  c85b->cd(++ipad);
  gPad->SetTopMargin(0);
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.2);
  gPad->SetRightMargin(0.1);
  havpt2pt1[0]->SetMaximum(0.75);
  havpt2pt1[0]->SetMinimum(0.65);
  MakeHist(havpt2pt1[0],"E_{T}^{HF[|#eta|<4]}","< p_{T,2} / p_{T,1} >");
  havpt2pt1[0]->GetYaxis()->SetTitleOffset(1.03);
  havpt2pt1[0]->Draw("p");
  havpt2pt1[1]->Draw("psame");
  havpt2pt1[2]->Draw("psame");
  havpt2pt1[3]->Draw("psame");

  drawText("pPb #sqrt{s_{NN}} = 5.02 TeV",0.33,0.88,22);
  drawText2("p_{T,1} > 120 GeV/c",0.30,0.75,19,kBlack);
  drawText2("p_{T,2} >  30 GeV/c",0.30,0.68,19,kBlack);
  drawText2("|#eta| < 3, #Delta#phi_{12} #geq #frac{2#pi}{3}",0.30,0.60,19,kBlack);


  c85b->cd(++ipad);
  gPad->SetTopMargin(0);
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.2);
  gPad->SetRightMargin(0.01);
  havdphi[0]->SetMaximum(0.25);
  havdphi[0]->SetMinimum(0.15);

  MakeHist(havdphi[0],"E_{T}^{HF[|#eta|<4]}","#sigma ( #Delta#phi )");
  havdphi[0]->GetYaxis()->SetTitleOffset(1.03);
  havdphi[0]->Draw("p");
  havdphi[1]->Draw("psame");
  havdphi[2]->Draw("psame");
  havdphi[3]->Draw("psame");
  leg1->Draw();

  if(iSave){
    c85b->SaveAs("giffiles/Av_ptratio_dphi.gif");
    c85b->SaveAs("giffiles/Av_ptratio_dphi.pdf");
    c85b->SaveAs("giffiles/Av_ptratio_dphi.C");
    c85b->SaveAs("giffiles/Av_ptratio_dphi.eps");
  }





  ipad=0;
  TH1D *hDum;   
  float rymin = 1e-05;
  float rymax = 1.000;
  hDum = GetDummyHist(-ketacut,ketacut,rymin,rymax,"hDum","#eta_{jet}","Event Fraction",false);
  hDum->GetYaxis()->SetNdivisions(611); 
  hDum->GetYaxis()->SetLabelFont(43);
  hDum->GetYaxis()->SetTitleFont(43);
  hDum->GetYaxis()->SetLabelSize(20);
  hDum->GetYaxis()->SetTitleSize(22);
  hDum->GetYaxis()->SetTitleOffset(2.6);
  hDum->GetXaxis()->SetNdivisions(614); 
  hDum->GetXaxis()->SetLabelFont(43);
  hDum->GetXaxis()->SetTitleFont(43);
  hDum->GetXaxis()->SetLabelSize(20);
  hDum->GetXaxis()->SetTitleSize(22);
  hDum->GetXaxis()->SetTitleOffset(3.1);
  hDum->GetXaxis()->SetNoExponent();
  hDum->GetXaxis()->SetMoreLogLabels();
  hDum->SetLineColor(10);
  hDum->SetMarkerColor(10);
  hDum->SetMarkerSize(0);
  hDum->SetMarkerStyle(15);

  TH1D *hDuma = (TH1D*) hDum->Clone("hDuma");
  hDuma->GetXaxis()->SetTitleOffset(2.30);
  hDuma->GetXaxis()->SetTitle("#eta_{1}");
  hDuma->SetAxisRange(1e-05,130.386,"Y");
  ipad=0;
  TCanvas *c83a  = new TCanvas("c83a","Jet eta",1100,770);
  makeMultiPanelCanvasWithGap(c83a,3,2,0.01,0.01,0.16,0.2,0.04,0.04); 
  for(int im=0;im<nmult;im++){
    c83a->cd(++ipad);
    gPad->SetLogy();

    hDuma->Draw("");
    drawText2(chf[im],0.23,0.79,19,kBlack);
    if(ipad==1 || ipad==4)drawText(cpad[ipad-1],0.2,0.88,24);
    else drawText(cpad[ipad-1],0.08,0.88,24);

    if(ipad==2)leg2->Draw();

    if(ipad==3){
      drawText("pPb #sqrt{s_{NN}} = 5.02 TeV",0.33,0.88,22);
      drawText2("p_{T,1} > 120 GeV/c",0.30,0.45,19,kBlack);
      drawText2("p_{T,2} >  30 GeV/c",0.30,0.38,19,kBlack);
      drawText2("|#eta| < 3, #Delta#phi_{12} #geq #frac{2#pi}{3}",0.30,0.30,19,kBlack);
    }

    for(int i=0;i<4;i++){
      heta1[i][im]->SetLineColor(icol[i]);
      heta1[i][im]->SetMarkerColor(icol[i]);
      heta1[i][im]->SetMarkerStyle(isty[i]);
      heta1[i][im]->SetMarkerSize(1.2);
      heta1[i][im]->Draw("psame");
    }
  }

  if(iSave){
    c83a->SaveAs("giffiles/LeadingJetEta.gif");
    c83a->SaveAs("giffiles/LeadingJetEta.pdf");
    c83a->SaveAs("giffiles/LeadingJetEta.C");
    c83a->SaveAs("giffiles/LeadingJetEta.eps");
  }


  TH1D *hDumb = (TH1D*) hDum->Clone("hDumb");
  hDumb->GetXaxis()->SetTitleOffset(2.30);
  hDumb->GetXaxis()->SetTitle("#eta_{2}");
  hDumb->SetAxisRange(1e-05,130.386,"Y");
  ipad=0;
  TCanvas *c83b  = new TCanvas("c83b","Jet eta",1100,770);
  makeMultiPanelCanvasWithGap(c83b,3,2,0.01,0.01,0.16,0.2,0.04,0.04); 
  for(int im=0;im<nmult;im++){
    c83b->cd(++ipad);
    gPad->SetLogy();

    hDumb->Draw("");
    drawText2(chf[im],0.23,0.79,19,kBlack);
    if(ipad==1 || ipad==4)drawText(cpad[ipad-1],0.2,0.88,24);
    else drawText(cpad[ipad-1],0.08,0.88,24);

    if(ipad==2)leg2->Draw();
    if(ipad==3){
      drawText("pPb #sqrt{s_{NN}} = 5.02 TeV",0.33,0.88,22);
      drawText2("p_{T,1} > 120 GeV/c",0.30,0.45,19,kBlack);
      drawText2("p_{T,2} >  30 GeV/c",0.30,0.38,19,kBlack);
      drawText2("|#eta| < 3, #Delta#phi_{12} #geq #frac{2#pi}{3}",0.30,0.30,19,kBlack);
    }


    for(int i=0;i<4;i++){
      heta2[i][im]->SetLineColor(icol[i]);
      heta2[i][im]->SetMarkerColor(icol[i]);
      heta2[i][im]->SetMarkerStyle(isty[i]);
      heta2[i][im]->SetMarkerSize(1.2);
      heta2[i][im]->Draw("psame");
    }
  }
  if(iSave){
    c83b->SaveAs("giffiles/SubLeadingJetEta.gif");
    c83b->SaveAs("giffiles/SubLeadingJetEta.pdf");
    c83b->SaveAs("giffiles/SubLeadingJetEta.C");
    c83b->SaveAs("giffiles/SubLeadingJetEta.eps");
  }

  //return 0;

  TH1D *hDumc = (TH1D*) hDum->Clone("hDumc");
  hDumc->GetXaxis()->SetTitleOffset(2.30);
  hDumc->GetXaxis()->SetTitle("( #eta_{1} + #eta_{2} ) / 2");
  hDumc->SetAxisRange(0,0.6,"Y");
  ipad=0;

  TLegend *leg3 = new TLegend(0.6386408,0.4609164,0.9424808,0.7579625,"BRNDC");
  leg3->SetHeader("");
  leg3->SetFillStyle(0);
  leg3->SetFillColor(10);
  leg3->SetLineColor(10);
  leg3->SetTextFont(42);
  leg3->SetTextSize(0.04);
  leg3->AddEntry(hdjeta [0][0],"Kurt","p");
  leg3->AddEntry(hdjeta [1][0],"Pawan","p");
  leg3->AddEntry(hdjeta [2][0],"Ruchi","p");
  leg3->AddEntry(hdjeta [3][0],"Doga","p");


  TCanvas *c83c  = new TCanvas("c83c","DiJet eta",1100,770);
  makeMultiPanelCanvasWithGap(c83c,3,2,0.01,0.01,0.16,0.2,0.04,0.04); 
  for(int im=0;im<nmult;im++){
    c83c->cd(++ipad);
    //gPad->SetLogy();

    hDumc->Draw("");
    drawText2(chf[im],0.23,0.79,19,kBlack);
    if(ipad==1 || ipad==4)drawText(cpad[ipad-1],0.2,0.88,24);
    else drawText(cpad[ipad-1],0.08,0.88,24);

    if(ipad==2)leg3->Draw();
    if(ipad==3){
      drawText("pPb #sqrt{s_{NN}} = 5.02 TeV",0.33,0.88,22);
      drawText2("p_{T,1} > 120 GeV/c",0.30,0.70,19,kBlack);
      drawText2("p_{T,2} >  30 GeV/c",0.30,0.64,19,kBlack);
      drawText2("|#eta| < 3, #Delta#phi_{12} #geq #frac{2#pi}{3}",0.30,0.56,19,kBlack);
    }
    for(int i=0;i<4;i++){
      hdjeta[i][im]->SetLineColor(icol[i]);
      hdjeta[i][im]->SetMarkerColor(icol[i]);
      hdjeta[i][im]->SetMarkerStyle(isty[i]);
      hdjeta[i][im]->SetMarkerSize(1.2);
      hdjeta[i][im]->Draw("psame");
    }
  }
  if(iSave){
    c83c->SaveAs("giffiles/DiJetEta.gif");
    c83c->SaveAs("giffiles/DiJetEta.pdf");
    c83c->SaveAs("giffiles/DiJetEta.C");
    c83c->SaveAs("giffiles/DiJetEta.eps");
  }


  // pT2/pT1
  TH1D *hDum1 = GetDummyHist(0,1,0,0.5,"hDum1","p_{T,2} / p_{T,1}","Event Fraction",false);
  hDum1->GetYaxis()->SetNdivisions(611); 
  hDum1->GetYaxis()->SetLabelFont(43);
  hDum1->GetYaxis()->SetTitleFont(43);
  hDum1->GetYaxis()->SetLabelSize(20);
  hDum1->GetYaxis()->SetTitleSize(22);
  hDum1->GetYaxis()->SetTitleOffset(2.6);
  hDum1->GetXaxis()->SetNdivisions(614); 
  hDum1->GetXaxis()->SetLabelFont(43);
  hDum1->GetXaxis()->SetTitleFont(43);
  hDum1->GetXaxis()->SetLabelSize(20);
  hDum1->GetXaxis()->SetTitleSize(22);
  hDum1->GetXaxis()->SetTitleOffset(3.1);
  hDum1->GetXaxis()->SetNoExponent();
  hDum1->GetXaxis()->SetMoreLogLabels();
  hDum1->SetLineColor(10);
  hDum1->SetMarkerColor(10);
  hDum1->SetMarkerSize(0);
  hDum1->SetMarkerStyle(15);

  TH1D *hDumd = (TH1D*) hDum1->Clone("hDumd");
  hDumd->GetXaxis()->SetTitleOffset(2.30);
  hDumd->SetAxisRange(0,0.6,"Y");

  ipad=0;
  TCanvas *c83d  = new TCanvas("c83d","pT ratio",1100,770);
  makeMultiPanelCanvasWithGap(c83d,3,2,0.01,0.01,0.16,0.2,0.04,0.04); 
  for(int im=0;im<nmult;im++){
    c83d->cd(++ipad);
    //gPad->SetLogy();

    hDumd->Draw("");
    drawText2(chf[im],0.23,0.79,19,kBlack);
    if(ipad==1 || ipad==4)drawText(cpad[ipad-1],0.2,0.88,24);
    else drawText(cpad[ipad-1],0.08,0.88,24);

    if(ipad==2)leg3->Draw();
    if(ipad==3){
      drawText("pPb #sqrt{s_{NN}} = 5.02 TeV",0.33,0.88,22);
      drawText2("p_{T,1} > 120 GeV/c",0.30,0.70,19,kBlack);
      drawText2("p_{T,2} >  30 GeV/c",0.30,0.64,19,kBlack);
      drawText2("|#eta| < 3, #Delta#phi_{12} #geq #frac{2#pi}{3}",0.30,0.56,19,kBlack);
    }
    for(int i=0;i<4;i++){
      hpt2pt1[i][im]->SetLineColor(icol[i]);
      hpt2pt1[i][im]->SetMarkerColor(icol[i]);
      hpt2pt1[i][im]->SetMarkerStyle(isty[i]);
      hpt2pt1[i][im]->SetMarkerSize(1.2);
      hpt2pt1[i][im]->Draw("psame");
    }
  }
  if(iSave){
    c83d->SaveAs("giffiles/pT2bypT1.gif");
    c83d->SaveAs("giffiles/pT2bypT1.pdf");
    c83d->SaveAs("giffiles/pT2bypT1.C");
    c83d->SaveAs("giffiles/pT2bypT1.eps");
  }


  //! dphi
  TH1D *hDum2 = GetDummyHist(0,pi,1e-05,1.5,"hDum2","#Delta#phi_{1,2}","Event Fraction",false);
  hDum2->GetYaxis()->SetNdivisions(611); 
  hDum2->GetYaxis()->SetLabelFont(43);
  hDum2->GetYaxis()->SetTitleFont(43);
  hDum2->GetYaxis()->SetLabelSize(20);
  hDum2->GetYaxis()->SetTitleSize(22);
  hDum2->GetYaxis()->SetTitleOffset(2.6);
  hDum2->GetXaxis()->SetNdivisions(613); 
  hDum2->GetXaxis()->SetLabelFont(43);
  hDum2->GetXaxis()->SetTitleFont(43);
  hDum2->GetXaxis()->SetLabelSize(20);
  hDum2->GetXaxis()->SetTitleSize(22);
  hDum2->GetXaxis()->SetTitleOffset(3.1);
  hDum2->GetXaxis()->SetNoExponent();
  hDum2->GetXaxis()->SetMoreLogLabels();
  hDum2->SetLineColor(10);
  hDum2->SetMarkerColor(10);
  hDum2->SetMarkerSize(0);
  hDum2->SetMarkerStyle(15);



  TLegend *leg4 = new TLegend(0.05479128,0.4939215,0.2782031,0.7607129,"BRNDC");
  leg4->SetHeader("");
  leg4->SetFillStyle(0);
  leg4->SetFillColor(10);
  leg4->SetLineColor(10);
  leg4->SetTextFont(42);
  leg4->SetTextSize(0.04);
  leg4->AddEntry(hdphi [0][0],"Kurt","p");
  leg4->AddEntry(hdphi [1][0],"Pawan","p");
  leg4->AddEntry(hdphi [2][0],"Ruchi","p");
  leg4->AddEntry(hdphi [3][0],"Doga","p");


  TH1D *hDume = (TH1D*) hDum2->Clone("hDume");
  hDume->GetXaxis()->SetTitleOffset(2.30);

  ipad=0;
  TCanvas *c83e  = new TCanvas("c83e","dphi",1100,770);
  makeMultiPanelCanvasWithGap(c83e,3,2,0.01,0.01,0.16,0.2,0.04,0.04); 
  for(int im=0;im<nmult;im++){
    c83e->cd(++ipad);
    gPad->SetLogy();

    hDume->Draw("");
    drawText2(chf[im],0.23,0.79,19,kBlack);
    if(ipad==1 || ipad==4)drawText(cpad[ipad-1],0.2,0.88,24);
    else drawText(cpad[ipad-1],0.08,0.88,24);

    if(ipad==2)leg4->Draw();
    if(ipad==3){
      drawText("pPb #sqrt{s_{NN}} = 5.02 TeV",0.33,0.88,22);
      drawText2("p_{T,1} > 120 GeV/c",0.23,0.70,19,kBlack);
      drawText2("p_{T,2} >  30 GeV/c",0.23,0.64,19,kBlack);
    }
    for(int i=0;i<4;i++){
      hdphi[i][im]->SetLineColor(icol[i]);
      hdphi[i][im]->SetMarkerColor(icol[i]);
      hdphi[i][im]->SetMarkerStyle(isty[i]);
      hdphi[i][im]->SetMarkerSize(1.2);
      hdphi[i][im]->Draw("psame");
    }
  }
  if(iSave){
    c83e->SaveAs("giffiles/dphi.gif");
    c83e->SaveAs("giffiles/dphi.pdf");
    c83e->SaveAs("giffiles/dphi.C");
    c83e->SaveAs("giffiles/dphi.eps");
  }

  return 0;
}
