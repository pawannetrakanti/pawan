#include <TROOT.h>
#include <TMath.h>
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
double round_digits (double /*x*/, int /*y*/); 
void MakeSmooth(TH1 */*histo*/);

void MakeHistRMS(TH1 */*hRMS*/,float /*max*/,float /*min*/);
void MakeHistMean(TH1 */*Mean*/,float /*max*/,float /*min*/);

double xmin=80.;
double xmax=410.;
int kSmooth=80;

int ResponseAn(bool iFit=false){

  float ketacut=2.0;
  const int ncen=7;
  const char *ksp  [ncen] = {"pbpb","pbpb" ,"pbpb"  ,"pbpb"  ,"pbpb"  ,"pbpb"  ,"pp"};
  const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","pp"};

  double ptbins[] ={80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400};
  const int bins  = sizeof(ptbins)/sizeof(Double_t) - 1;
  const int nbins = bins;

  bool iSave=false;
  bool iSigma=false;
  //double max=265.;

  TFile *fin[ncen];
  for(int ic=0;ic<ncen;ic++){
    fin[ic] = new TFile(Form("rootfiles/JERJES/Histo_%s_%s.root",ksp[ic],ccent[ic]));
    cout<<"input file names : "<< fin[ic] ->GetName()<<endl;
  }

  const int knj=7;
  const char *calgo[knj]= {
    "icPu5", 
    "ak2PF","ak3PF","ak4PF",
    "akPu2PF","akPu3PF","akPu4PF"
  };                                     

  
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

  
  //const int nbins = bins;
  //std::cout<<"# of pt bins : " <<nbins<<std::endl;
  //std::cout<<std::endl;


  TH1F *hArM[knj][ncen], *hRMS[knj][ncen];
  TH1F *hArM_lead[knj][ncen], *hRMS_lead[knj][ncen];
  TH1F *hArM_slead[knj][ncen], *hRMS_slead[knj][ncen];
  TH1F *hArM_remain[knj][ncen], *hRMS_remain[knj][ncen];
  TH1F *hArM_genm[knj][ncen], *hRMS_genm[knj][ncen];

  TH1F *hArM_Raw[knj][ncen];  
  TH1F *hArM_eta[knj][ncen][2], *hRMS_eta[knj][ncen][2];

  int sty    [7] = {24,25,26,28,30,32,20};
  Color_t col[7] = {kRed,kBlue,kMagenta,kCyan,kGreen+1,kBlue+4,kRed+4};
  //return 0;

  for(int nj=0;nj<knj;nj++){
    cout<<"calgo : "<<calgo[nj]<<endl;
    for(int ic=0;ic<ncen;ic++){
      
      int icen=ic;

      //cout<<"\t  : "<<icen<<"\t ccent "<<ccent[icen]<<"\t"<<ksp[ic]<<"\t"<<Form("hArM_%s_%s_%d",ksp[ic],calgo[nj],icen)<<"\t"<<hArM[nj][ic]->GetName()<<endl;

      if(iFit){
	hArM[nj][ic] = (TH1F*)fin[ic]->Get(Form("hMean_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS[nj][ic] = (TH1F*)fin[ic]->Get(Form("hSigma_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hArM_Raw[nj][ic] = (TH1F*)fin[ic]->Get(Form("hMean_r_%s_%s_%d",ksp[ic],calgo[nj],icen));

	hArM_lead[nj][ic] = (TH1F*)fin[ic]->Get(Form("hMean_lead_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS_lead[nj][ic] = (TH1F*)fin[ic]->Get(Form("hSigma_lead_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hArM_slead[nj][ic] = (TH1F*)fin[ic]->Get(Form("hMean_slead_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS_slead[nj][ic] = (TH1F*)fin[ic]->Get(Form("hSigma_slead_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hArM_remain[nj][ic] = (TH1F*)fin[ic]->Get(Form("hMean_remain_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS_remain[nj][ic] = (TH1F*)fin[ic]->Get(Form("hSigma_remain_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hArM_genm[nj][ic] = (TH1F*)fin[ic]->Get(Form("hMean_genm_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS_genm[nj][ic] = (TH1F*)fin[ic]->Get(Form("hSigma_genm_%s_%s_%d",ksp[ic],calgo[nj],icen));

      }else{
	hArM[nj][ic] = (TH1F*)fin[ic]->Get(Form("hArM_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS[nj][ic] = (TH1F*)fin[ic]->Get(Form("hRMS_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hArM_Raw[nj][ic] = (TH1F*)fin[ic]->Get(Form("hArM_r_%s_%s_%d",ksp[ic],calgo[nj],icen));

	hArM_lead[nj][ic] = (TH1F*)fin[ic]->Get(Form("hArM_lead_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS_lead[nj][ic] = (TH1F*)fin[ic]->Get(Form("hRMS_lead_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hArM_slead[nj][ic] = (TH1F*)fin[ic]->Get(Form("hArM_slead_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS_slead[nj][ic] = (TH1F*)fin[ic]->Get(Form("hRMS_slead_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hArM_remain[nj][ic] = (TH1F*)fin[ic]->Get(Form("hArM_remain_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS_remain[nj][ic] = (TH1F*)fin[ic]->Get(Form("hRMS_remain_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hArM_genm[nj][ic] = (TH1F*)fin[ic]->Get(Form("hArM_genm_%s_%s_%d",ksp[ic],calgo[nj],icen));
	hRMS_genm[nj][ic] = (TH1F*)fin[ic]->Get(Form("hRMS_genm_%s_%s_%d",ksp[ic],calgo[nj],icen));
      }

      for(int ie=0;ie<2;ie++){
	hArM_eta[nj][ic][ie] = (TH1F*)fin[ic]->Get(Form("hArM_eta_%s_%s_%d_%d",ksp[ic],calgo[nj],icen,ie));
	hRMS_eta[nj][ic][ie] = (TH1F*)fin[ic]->Get(Form("hRMS_eta_%s_%s_%d_%d",ksp[ic],calgo[nj],icen,ie));
      }

    }
  }
  //return 0;


  /*
  //! Round off the digits
  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){
      MakeSmooth(hArM[nj][ic]);
      MakeSmooth(hRMS[nj][ic]);
      MakeSmooth(hArM_Raw[nj][ic]);
      for(int ie=0;ie<2;ie++){
	MakeSmooth(hArM_eta[nj][ic][ie]);
	MakeSmooth(hRMS_eta[nj][ic][ie]);
      }
    }
  }
  */


  TH1F *hsmf[knj][ncen-1];
  double smearf[knj][ncen][bins];
  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){
      for(int ix=0;ix<bins;ix++){
	smearf[nj][ic][ix]=0;
      }
    }
  }

  const int pal[3]    = {4,5,6}; 

  //! Smearing factors
  for(int nj=pal[0];nj<pal[2]+1;nj++){
    cout<<"For algo : "<<calgo[nj]<<endl;

    for(int ic=0;ic<ncen;ic++){

      if(ic>ncen-2)continue;
      //cout<<"Centrality : "<<ccent[ic]<<endl;
      hsmf[nj-3][ic] = new TH1F(Form("hsmf%d_%d",nj-3,ic),Form("Smearing factors %s %s",calgo[nj-3],ccent[ic]),bins,ptbins);
      hsmf[nj-3][ic]->Sumw2();
      hsmf[nj-3][ic]->SetMarkerStyle(20);
      hsmf[nj-3][ic]->SetMarkerColor(1);
      hsmf[nj-3][ic]->SetLineColor(1);
      hsmf[nj-3][ic]->SetMarkerSize(1.0);

      for(int ix=0;ix<hArM[nj][ic]->GetNbinsX();ix++){
	float sigpp   = hRMS[nj-3][ncen-1]->GetBinContent(ix+1);
	float esigpp  = hRMS[nj-3][ncen-1]->GetBinError(ix+1);
	float meanpp  = hArM[nj-3][ncen-1]->GetBinContent(ix+1);


	float sigpbpb  = hRMS[nj][ic]->GetBinContent(ix+1);
	float esigpbpb = hRMS[nj][ic]->GetBinError(ix+1);
	float meanpbpb = hArM[nj][ic]->GetBinContent(ix+1);

	float smf = (sqrt(pow(sigpbpb,2)/meanpbpb - pow(sigpp,2)/meanpp))*100.;
	//if(nj==5 && ic==5)cout<<"\t ix : "<<ix<<"\t pt : "<<ptbins[ix]<<"\t pbpb : "<<sigpbpb<<"\t pp : "<<sigpp<<"\t smf  "<<smf<<endl;

	if(sigpbpb > sigpp){
	  float err = (sqrt(pow(esigpbpb/sigpbpb,2) + pow(esigpp/sigpp,2)))*100.;
	  //cout<<"\t ix : "<<ix<<"\t pt : "<<ptbins[ix]<<"\t pbpb : "<<sigpbpb<<"\t pp : "<<sigpp<<"\t smf : "<<smf<<"\t "<<err<<endl;
	  hsmf  [nj-3][ic]->SetBinContent(ix,smf);
	  hsmf  [nj-3][ic]->SetBinError(ix,err);
	  smearf[nj-3][ic][ix] = smf;
	}else{
	  smearf[nj-3][ic][ix]=0;
	}
	if(ix==bins-1)cout<<smearf[nj-3][ic][ix];
	else cout<<smearf[nj-3][ic][ix]<<",";
      }
      cout<<endl;
    }
    cout<<endl;
  }


  /*
  //! Write out smearing factors in an array
  for(int nj=2;nj<3;nj++){
    for(int ic=0;ic<ncen-1;ic++){    
      for(int ix=0;ix<bins;ix++){
	if(ix==bins-1)cout<<smearf[nj][ic][ix];
	else cout<<smearf[nj][ic][ix]<<",";
      }
      cout<<endl;
      cout<<endl;
    }
    cout<<endl;
  }

  //return 0;
  */

  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){

      MakeHistRMS(hRMS[nj][ic],0.28,0.0001);
      MakeHistRMS(hRMS_lead[nj][ic],0.28,0.0001);
      MakeHistRMS(hRMS_slead[nj][ic],0.28,0.0001);
      MakeHistRMS(hRMS_remain[nj][ic],0.28,0.0001);
      MakeHistRMS(hRMS_genm[nj][ic],0.28,0.0001);

      MakeHistMean(hArM[nj][ic],1.068,0.94);
      MakeHistMean(hArM_lead[nj][ic],1.068,0.94);
      MakeHistMean(hArM_slead[nj][ic],1.068,0.94);
      MakeHistMean(hArM_remain[nj][ic],1.068,0.94);
      MakeHistMean(hArM_genm[nj][ic],1.068,0.94);
      
      MakeHistMean(hArM_Raw[nj][ic],1.058,0.758);
    }      
  }





  const char *algn[3] = {"Anti-k_{T}, PF, R = 0.2","Anti-k_{T}, PF, R = 0.3","Anti-k_{T}, PF, R = 0.4"};
  int ipad=0;
  int maxc=ncen-1;
  int maxr=2;

  //hsmf[5][0]->Draw("p");
  //return 0;

  ipad=0;
  TCanvas *c00 = new TCanvas("c00","Smf",102,141,1328,776);
  c00->Divide(3,2);
  for(int ic=ncen-2;ic>=0;ic--){
    c00->cd(++ipad);
    hsmf[2][ic]->Draw("p");
  }
  //return 0;



  //TLegend *l3[maxc];
  //int ij=-1;
  //int ic=0;


  TLine *line = new TLine(xmin,1.0,xmax+100,1.0);
  line->SetLineWidth(2);
  line->SetLineStyle(2);




  TCanvas *c3[3];
  TLegend *l3[3];
  for(int nj=0;nj<3;nj++){
    c3[nj] = new TCanvas(Form("c3_%d",nj),Form("%d JES JER",nj),97,118,1739,479);
    makeMultiPanelCanvas(c3[nj],maxc,maxr,0.0,0.0,0.22,0.22,0.02);
    //std::cout<<std::endl;
    ipad=0;

    l3[nj] = new TLegend(0.644011,0.6497811,0.8825236,0.9649781,NULL,"BRNDC");
    l3[nj]->SetHeader("");
    l3[nj]->SetBorderSize(0);
    l3[nj]->SetTextFont(42);
    l3[nj]->SetTextSize(0.09);
    l3[nj]->SetLineColor(1);
    l3[nj]->SetLineStyle(1);
    l3[nj]->SetLineWidth(1);
    l3[nj]->SetFillColor(10);
    l3[nj]->SetFillStyle(1001);
    l3[nj]->SetHeader("");			  

    for(int ic=ncen-2;ic>=0;ic--){


      hRMS[pal[nj]][ic]->SetMaximum(0.28);
      hRMS[pal[nj]][ic]->SetMinimum(0.0001);

      hRMS[pal[nj]][ic]->SetMarkerStyle(24);
      hRMS[pal[nj]][ic]->SetMarkerColor(1);
      hRMS[pal[nj]][ic]->SetLineColor(1);
      hRMS[pal[nj]][ic]->SetMarkerSize(1.3);

      //hArM[pal[nj]][ic]->SetMaximum(1.088);
      //hArM[pal[nj]][ic]->SetMinimum(0.288);

      hArM[pal[nj]][ic]->SetMarkerStyle(24);
      hArM[pal[nj]][ic]->SetMarkerColor(1);
      hArM[pal[nj]][ic]->SetLineColor(1);
      hArM[pal[nj]][ic]->SetMarkerSize(1.3);

      hArM_Raw[pal[nj]][ic]->SetMarkerStyle(25);
      hArM_Raw[pal[nj]][ic]->SetMarkerColor(1);
      hArM_Raw[pal[nj]][ic]->SetLineColor(1);
      hArM_Raw[pal[nj]][ic]->SetMarkerSize(1.3);
      
      c3[nj]->cd(++ipad);
      //gPad->SetLogx();
      hRMS[pal[nj]][ic]->Draw("p");

      hRMS[pal[nj]-3][ncen-1]->SetMarkerStyle(20);
      hRMS[pal[nj]-3][ncen-1]->SetMarkerColor(1);
      hRMS[pal[nj]-3][ncen-1]->SetLineColor(1);
      hRMS[pal[nj]-3][ncen-1]->SetMarkerSize(1.0);
      hRMS[pal[nj]-3][ncen-1]->Draw("psame");

      if(ipad==1){
	//! with pileup
	hRMS[pal[nj]][ncen-1]->SetMarkerStyle(30);
	hRMS[pal[nj]][ncen-1]->SetMarkerColor(2);
	hRMS[pal[nj]][ncen-1]->SetLineColor(2);
	hRMS[pal[nj]][ncen-1]->SetMarkerSize(1.2);
	hRMS[pal[nj]][ncen-1]->Draw("psame");
      }
      /*
      if(pal[nj]>3){
	//hRMS[pal[nj]-3][ncen-1]->SetMarkerStyle(20);
	//hRMS[pal[nj]-3][ncen-1]->SetMarkerColor(1);
	//hRMS[pal[nj]-3][ncen-1]->SetLineColor(1);
	//hRMS[pal[nj]-3][ncen-1]->SetMarkerSize(1.0);
	//hRMS[pal[nj]-3][ncen-1]->Draw("psame");
      }else{
	hRMS[pal[nj]][ncen-1]->SetMarkerStyle(20);
	hRMS[pal[nj]][ncen-1]->SetMarkerColor(1);
	hRMS[pal[nj]][ncen-1]->SetLineColor(1);
	hRMS[pal[nj]][ncen-1]->SetMarkerSize(1.0);
	hRMS[pal[nj]][ncen-1]->Draw("psame");
      }
      */

      if(ipad==6){
	l3[nj]->AddEntry(hRMS[pal[nj]-3][ncen-1],"pp","p"); 
	l3[nj]->AddEntry(hRMS[pal[nj]][ncen-1],"pp w pu","p"); 
	l3[nj]->AddEntry(hRMS[pal[nj]][0],"PbPb","p"); 
	l3[nj]->Draw();
      }
      if(ipad==1){
	TPaveText *pt3 = new TPaveText(0.4193406,0.7098186,0.7083456,0.8899312,"brNDC");
	pt3->SetBorderSize(0);
	pt3->SetFillColor(10);
	pt3->SetTextFont(42);
	TText *text3 = pt3->AddText("CMS Simulation");
	text3->SetTextSize(0.09);
	TText *text4 = pt3->AddText("PYTHIA Z2+HYDJET1.8");
	text4->SetTextSize(0.09);
	pt3->Draw();
      }
      
      TPaveText *pt = new TPaveText(0.2720047,0.03939961,0.5610097,0.2195122,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(10);
      pt->SetTextFont(42);
      TText *text = pt->AddText(Form("%s",ccent[ic]));
      text->SetTextSize(0.09);
      pt->Draw();

      if(ipad==2){
	TPaveText *pt1 = new TPaveText(0.1346399,0.6797999,0.9011061,0.8749218,"brNDC");
	pt1->SetBorderSize(0);
	pt1->SetFillColor(10);
	pt1->SetTextFont(42);
	//TText *text1 = pt1->AddText(Form("%s",calgo[pal[nj]]));
	TText *text1 = pt1->AddText(Form("%s",algn[nj]));
	text1->SetTextSize(0.09);
	TText *text2 = pt1->AddText("|#eta|<2.0");
	text2->SetTextSize(0.09);
	pt1->Draw();
      }


      c3[nj]->cd(ipad+(ncen-1));

      //cout<<"pad : "<<(ipad+(ncen-1))<<"\t "<<calgo[pal[nj]]<<"\t centrality : "<<ccent[ic]<<"\t name : "<<hArM[pal[nj]][ic]->GetName()<<endl;
      //gPad->SetLogx();
      hArM[pal[nj]][ic]->Draw("p");

      hArM[pal[nj]-3][ncen-1]->SetMarkerStyle(20);
      hArM[pal[nj]-3][ncen-1]->SetMarkerColor(1);
      hArM[pal[nj]-3][ncen-1]->SetLineColor(1);
      hArM[pal[nj]-3][ncen-1]->SetMarkerSize(1.0);
      hArM[pal[nj]-3][ncen-1]->Draw("psame");

      if(ipad==1){
	//! pp with pu
	hArM[pal[nj]][ncen-1]->SetMarkerStyle(30);
	hArM[pal[nj]][ncen-1]->SetMarkerColor(2);
	hArM[pal[nj]][ncen-1]->SetLineColor(2);
	hArM[pal[nj]][ncen-1]->SetMarkerSize(1.2);
	hArM[pal[nj]][ncen-1]->Draw("psame");
      }
      hArM_Raw[pal[nj]-3][ncen-1]->SetMarkerStyle(21);
      hArM_Raw[pal[nj]-3][ncen-1]->SetMarkerColor(1);
      hArM_Raw[pal[nj]-3][ncen-1]->SetLineColor(1);
      hArM_Raw[pal[nj]-3][ncen-1]->SetMarkerSize(1.0);
      //hArM_Raw[pal[nj]-3][ncen-1]->Draw("psame");
      //hArM_Raw[pal[nj]][ic]->Draw("psame");      

      if(ipad==6){
	TLatex *tex1 = new TLatex(288.8638,0.8372817,"Raw");
	tex1->SetTextFont(42);
	tex1->SetTextSize(0.08);
	tex1->SetLineWidth(2);
	tex1->Draw();
	TLatex *tex2 = new TLatex(286.0951,0.8069566,"Corrected");
	tex2->SetTextFont(42);
	tex2->SetTextSize(0.08);
	tex2->SetLineWidth(2);
	tex2->Draw();
	TMarker *marker1 = new TMarker(272.2515,0.848309,25);
	marker1->SetMarkerStyle(25);
	marker1->Draw();
	TMarker *marker2 = new TMarker(248.7175,0.848309,21);
	marker2->SetMarkerStyle(21);
	marker2->Draw();
	TMarker *marker3 = new TMarker(272.2515,0.8166055,24);
	marker3->SetMarkerStyle(24);
	marker3->Draw();
	TMarker *marker4 = new TMarker(248.7175,0.8166055,20);
	marker4->SetMarkerStyle(20);
	marker4->Draw();
      }

      /*
      if(pal[nj]>3){
	hArM[pal[nj]-3][ncen-1]->SetMarkerStyle(20);
	hArM[pal[nj]-3][ncen-1]->SetMarkerColor(1);
	hArM[pal[nj]-3][ncen-1]->SetLineColor(1);
	hArM[pal[nj]-3][ncen-1]->SetMarkerSize(1.0);
	hArM[pal[nj]-3][ncen-1]->Draw("psame");
      }else{
	hArM[pal[nj]][ncen-1]->SetMarkerStyle(20);
	hArM[pal[nj]][ncen-1]->SetMarkerColor(1);
	hArM[pal[nj]][ncen-1]->SetLineColor(1);
	hArM[pal[nj]][ncen-1]->SetMarkerSize(1.0);
	hArM[pal[nj]][ncen-1]->Draw("psame");
      }
      */
      line->Draw();
    }//! icen
    if(iSave){
      c3[nj]->SaveAs(Form("AN/JERJES/JERJES_%s.png",calgo[pal[nj]]));
      c3[nj]->SaveAs(Form("AN/JERJES/JERJES_%s.pdf",calgo[pal[nj]]));
      c3[nj]->SaveAs(Form("AN/JERJES/JERJES_%s.C",calgo[pal[nj]]));
      c3[nj]->SaveAs(Form("AN/JERJES/JERJES_%s.eps",calgo[pal[nj]]));
    }
  }//! algo
  return 0;


  //! Comparison
  TCanvas *c4[3];
  TLegend *l4[3];
  for(int nj=0;nj<3;nj++){
    c4[nj] = new TCanvas(Form("c4_%d",nj),Form("%d JES JER",nj),97,118,1739,479);
    makeMultiPanelCanvas(c4[nj],maxc+1,maxr,0.0,0.0,0.22,0.22,0.02);
    //std::cout<<std::endl;
    ipad=0;

    l4[nj] = new TLegend(0.3494706,0.5046904,0.8755845,0.9499687,NULL,"BRNDC");
    l4[nj]->SetHeader("");
    l4[nj]->SetBorderSize(0);
    l4[nj]->SetTextFont(42);
    l4[nj]->SetTextSize(0.09);
    l4[nj]->SetLineColor(1);
    l4[nj]->SetLineStyle(1);
    l4[nj]->SetLineWidth(1);
    l4[nj]->SetFillColor(10);
    l4[nj]->SetFillStyle(1001);
    l4[nj]->SetHeader("");			  

    for(int ic=ncen-1;ic>=0;ic--){
      
      c4[nj]->cd(++ipad);
      //gPad->SetLogx();

      if(ipad==1){ //! pp using akPF

	hRMS[pal[nj]-3][ic]->SetMarkerStyle(24);
	hRMS[pal[nj]-3][ic]->SetMarkerColor(1);
	hRMS[pal[nj]-3][ic]->SetLineColor(1);
	hRMS[pal[nj]-3][ic]->SetMarkerSize(1.0);

	
	hRMS_lead[pal[nj]-3][ic]->SetMarkerStyle(29);
	hRMS_lead[pal[nj]-3][ic]->SetMarkerColor(2);
	hRMS_lead[pal[nj]-3][ic]->SetLineColor(2);
	hRMS_lead[pal[nj]-3][ic]->SetMarkerSize(1.0);

	hRMS_slead[pal[nj]-3][ic]->SetMarkerStyle(30);
	hRMS_slead[pal[nj]-3][ic]->SetMarkerColor(2);
	hRMS_slead[pal[nj]-3][ic]->SetLineColor(2);
	hRMS_slead[pal[nj]-3][ic]->SetMarkerSize(1.0);

	hRMS_remain[pal[nj]-3][ic]->SetMarkerStyle(28);
	hRMS_remain[pal[nj]-3][ic]->SetMarkerColor(4);
	hRMS_remain[pal[nj]-3][ic]->SetLineColor(4);
	hRMS_remain[pal[nj]-3][ic]->SetMarkerSize(1.0);

	hRMS_genm[pal[nj]-3][ic]->SetMarkerStyle(27);
	hRMS_genm[pal[nj]-3][ic]->SetMarkerColor(8);
	hRMS_genm[pal[nj]-3][ic]->SetLineColor(8);
	hRMS_genm[pal[nj]-3][ic]->SetMarkerSize(1.0);
	
	hRMS[pal[nj]-3][ic]->Draw("p");	
	//hRMS_lead[pal[nj]-3][ic]->Draw("psame");
	//hRMS_slead[pal[nj]-3][ic]->Draw("psame");
	//hRMS_remain[pal[nj]-3][ic]->Draw("psame");
	hRMS_genm[pal[nj]-3][ic]->Draw("psame");

      }else{

	hRMS[pal[nj]][ic]->SetMarkerStyle(24);
	hRMS[pal[nj]][ic]->SetMarkerColor(1);
	hRMS[pal[nj]][ic]->SetLineColor(1);
	hRMS[pal[nj]][ic]->SetMarkerSize(1.0);

	
	hRMS_lead[pal[nj]][ic]->SetMarkerStyle(29);
	hRMS_lead[pal[nj]][ic]->SetMarkerColor(2);
	hRMS_lead[pal[nj]][ic]->SetLineColor(2);
	hRMS_lead[pal[nj]][ic]->SetMarkerSize(1.0);

	hRMS_slead[pal[nj]][ic]->SetMarkerStyle(30);
	hRMS_slead[pal[nj]][ic]->SetMarkerColor(2);
	hRMS_slead[pal[nj]][ic]->SetLineColor(2);
	hRMS_slead[pal[nj]][ic]->SetMarkerSize(1.0);

	hRMS_remain[pal[nj]][ic]->SetMarkerStyle(28);
	hRMS_remain[pal[nj]][ic]->SetMarkerColor(4);
	hRMS_remain[pal[nj]][ic]->SetLineColor(4);
	hRMS_remain[pal[nj]][ic]->SetMarkerSize(1.0);

	hRMS_genm[pal[nj]][ic]->SetMarkerStyle(27);
	hRMS_genm[pal[nj]][ic]->SetMarkerColor(8);
	hRMS_genm[pal[nj]][ic]->SetLineColor(8);
	hRMS_genm[pal[nj]][ic]->SetMarkerSize(1.0);

	hRMS       [pal[nj]][ic]->Draw("p");	
	//hRMS_lead  [pal[nj]][ic]->Draw("psame");
	//hRMS_slead [pal[nj]][ic]->Draw("psame");
	//hRMS_remain[pal[nj]][ic]->Draw("psame");
	hRMS_genm  [pal[nj]][ic]->Draw("psame");
	
      }

      if(ipad==ncen){
	l4[nj]->AddEntry(hRMS[pal[nj]][ic],"all","p"); 
	//l4[nj]->AddEntry(hRMS_lead[pal[nj]][ic],"leading","p"); 
	//l4[nj]->AddEntry(hRMS_slead[pal[nj]][ic],"sub-leading","p"); 
	//l4[nj]->AddEntry(hRMS_remain[pal[nj]][ic],"remaining","p"); 
	l4[nj]->AddEntry(hRMS_genm[pal[nj]][ic],"gen-matched","p"); 
	l4[nj]->Draw();
      }
      if(ipad==1){
	TPaveText *pt3 = new TPaveText(0.4193406,0.7098186,0.7083456,0.8899312,"brNDC");
	pt3->SetBorderSize(0);
	pt3->SetFillColor(10);
	pt3->SetTextFont(42);
	TText *text3 = pt3->AddText("CMS Simulation");
	text3->SetTextSize(0.09);
	//TText *text4 = pt3->AddText("PYTHIA Z2 + HYDJET1.8");
	//text4->SetTextSize(0.09);
	pt3->Draw();
      }
      
      TPaveText *pt = new TPaveText(0.2720047,0.03939961,0.5610097,0.2195122,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(10);
      pt->SetTextFont(42);
      TText *text = pt->AddText(Form("%s",ccent[ic]));
      text->SetTextSize(0.09);
      pt->Draw();

      if(ipad==2){
	TPaveText *pt1 = new TPaveText(0.1346399,0.6797999,0.9011061,0.8749218,"brNDC");
	pt1->SetBorderSize(0);
	pt1->SetFillColor(10);
	pt1->SetTextFont(42);
	//TText *text1 = pt1->AddText(Form("%s",calgo[pal[nj]]));
	TText *text1 = pt1->AddText(Form("%s",algn[nj]));
	text1->SetTextSize(0.09);
	TText *text2 = pt1->AddText("|#eta|<2.0");
	text2->SetTextSize(0.09);
	pt1->Draw();
      }


      c4[nj]->cd(ipad+ncen);


      if(ipad==1){ //! pp using akPF

	hArM[pal[nj]-3][ic]->SetMarkerStyle(24);
	hArM[pal[nj]-3][ic]->SetMarkerColor(1);
	hArM[pal[nj]-3][ic]->SetLineColor(1);
	hArM[pal[nj]-3][ic]->SetMarkerSize(1.0);

	
	hArM_lead[pal[nj]-3][ic]->SetMarkerStyle(29);
	hArM_lead[pal[nj]-3][ic]->SetMarkerColor(2);
	hArM_lead[pal[nj]-3][ic]->SetLineColor(2);
	hArM_lead[pal[nj]-3][ic]->SetMarkerSize(1.0);

	hArM_slead[pal[nj]-3][ic]->SetMarkerStyle(30);
	hArM_slead[pal[nj]-3][ic]->SetMarkerColor(2);
	hArM_slead[pal[nj]-3][ic]->SetLineColor(2);
	hArM_slead[pal[nj]-3][ic]->SetMarkerSize(1.0);

	hArM_remain[pal[nj]-3][ic]->SetMarkerStyle(28);
	hArM_remain[pal[nj]-3][ic]->SetMarkerColor(4);
	hArM_remain[pal[nj]-3][ic]->SetLineColor(4);
	hArM_remain[pal[nj]-3][ic]->SetMarkerSize(1.0);

	hArM_genm[pal[nj]-3][ic]->SetMarkerStyle(27);
	hArM_genm[pal[nj]-3][ic]->SetMarkerColor(8);
	hArM_genm[pal[nj]-3][ic]->SetLineColor(8);
	hArM_genm[pal[nj]-3][ic]->SetMarkerSize(1.0);

	hArM[pal[nj]-3][ic]->Draw("p");	
	//hArM_lead[pal[nj]-3][ic]->Draw("psame");
	//hArM_slead[pal[nj]-3][ic]->Draw("psame");
	//hArM_remain[pal[nj]-3][ic]->Draw("psame");
	hArM_genm[pal[nj]-3][ic]->Draw("psame");
	
      }else{

	hArM[pal[nj]][ic]->SetMarkerStyle(24);
	hArM[pal[nj]][ic]->SetMarkerColor(1);
	hArM[pal[nj]][ic]->SetLineColor(1);
	hArM[pal[nj]][ic]->SetMarkerSize(1.0);

	
	hArM_lead[pal[nj]][ic]->SetMarkerStyle(29);
	hArM_lead[pal[nj]][ic]->SetMarkerColor(2);
	hArM_lead[pal[nj]][ic]->SetLineColor(2);
	hArM_lead[pal[nj]][ic]->SetMarkerSize(1.0);

	hArM_slead[pal[nj]][ic]->SetMarkerStyle(30);
	hArM_slead[pal[nj]][ic]->SetMarkerColor(2);
	hArM_slead[pal[nj]][ic]->SetLineColor(2);
	hArM_slead[pal[nj]][ic]->SetMarkerSize(1.0);

	hArM_remain[pal[nj]][ic]->SetMarkerStyle(28);
	hArM_remain[pal[nj]][ic]->SetMarkerColor(4);
	hArM_remain[pal[nj]][ic]->SetLineColor(4);
	hArM_remain[pal[nj]][ic]->SetMarkerSize(1.0);

	hArM_genm[pal[nj]][ic]->SetMarkerStyle(27);
	hArM_genm[pal[nj]][ic]->SetMarkerColor(8);
	hArM_genm[pal[nj]][ic]->SetLineColor(8);
	hArM_genm[pal[nj]][ic]->SetMarkerSize(1.0);

	hArM       [pal[nj]][ic]->Draw("p");	
	//hArM_lead  [pal[nj]][ic]->Draw("psame");
	//hArM_slead [pal[nj]][ic]->Draw("psame");
	//hArM_remain[pal[nj]][ic]->Draw("psame");
	hArM_genm  [pal[nj]][ic]->Draw("psame");

      }
      line->Draw();
    }//! icen
    if(iSave){
      c4[nj]->SaveAs(Form("AN/JERJES/JERJES_Summary_%s.png",calgo[pal[nj]]));
      c4[nj]->SaveAs(Form("AN/JERJES/JERJES_Summary_%s.pdf",calgo[pal[nj]]));
      c4[nj]->SaveAs(Form("AN/JERJES/JERJES_Summary_%s.C",calgo[pal[nj]]));
      c4[nj]->SaveAs(Form("AN/JERJES/JERJES_Summary_%s.eps",calgo[pal[nj]]));
    }
  }//! algo

  return 0;
}
double round_digits (double x, int y) 
{ 
  return (floor(x*TMath::Power(10,y)+0.5)/TMath::Power(10,y)); 
}
void MakeSmooth(TH1 *histo){
  
  histo->Smooth(kSmooth);
  for(int ix=1;ix<=histo->GetNbinsX();ix++){
    double val = histo->GetBinContent(ix);
    if(val<0 || val==0)continue;
    double err = histo->GetBinError(ix);
    histo->SetBinContent(ix,round_digits(val,5));
    histo->SetBinError(ix,round_digits(err,5));
  }
}
void MakeHistRMS(TH1 *h1,float ymax,float ymin)
{
  h1->SetTitle("");
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetTitle("#sigma (RecoJet p_{T} / GenJet p_{T})");    
  h1->GetYaxis()->SetTitleSize(0.08);
  h1->GetYaxis()->SetTitleOffset(1.26);
  h1->GetYaxis()->SetLabelSize(0.08);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetDecimals(true);
  
}
void MakeHistMean(TH1 *h1,float ymax,float ymin)
{
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->SetTitle("");
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetLabelOffset(0.005);
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(1.50);
  h1->GetYaxis()->SetLabelSize(0.07);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelFont(42);
  
}
