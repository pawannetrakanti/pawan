#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TMath.h>
#include  "MultiCanvas.h"
#include "../Utils/utilities_v2.h";
#include "../Utils/commonStyle_v2.h"
void drawText(const char *text, float xp, float yp, int size);
void putCMSPrel(double x, double y, double size);
void drawText(const char *text, float xp, float yp, int size);
void drawText2(const char *text, float xp, float yp, int size);
TH1F *functionHist(TF1 */*f*/, TH1F* /*h*/,char */*fHistname*/);
using namespace std;
int spectra_wsyserr(const int nx=1)
{
  bool iSave=false;
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

  const int knj=3;
  const int kin=2;
  const int ncen=7;
  const int nvtx=4;

  TFile *fin[kin];
  //const char *kfile[kin] = {"HLT_HIJet80_v1_HiForest2_pbpb_pt80_2012.root","HLT_Jet60_v1_ppSmearing_HiForest2_pp_pt80_2012.root"}; //! BinByBin based
  const char *kfile[kin] = {"HLT_HIJet80_v1_HiForest2_pbpb_pt80_2012.root","HLT_Jet60_v1_ppSmearing_Fit_HiForest2_pp_pt80_2012.root"}; //! Fit based
  //const char *cinf[kin]  = {"PbPb","pp"};
  Color_t icol [ncen] = {kRed,kBlue,kMagenta,kCyan,kGreen+1,kBlue+4,kRed+4};
  const int sty [ncen] = {24,25,26,28,27,30,20};
  const int ssty[ncen] = {20,21,22,34,33,29,20};

  const char *calgo[knj] = {"Anti-k_{T}, PF, R = 0.2","Anti-k_{T}, PF, R = 0.3","Anti-k_{T}, PF, R = 0.4"};

  //! Luminosity 150 inv. mub  nn x-section : 7650  mb pbpb
  //! Luminosity 231 inv. nb   pp x-section : 64+-4 mb pp
  //!                      pbpb pp
  float lumi[kin] = {129.,212.};
  float xsec[kin] = {7.65,64.};
  float csca[ncen]= {2,2,8,8,8,8,1};  


  const char *ccen  [ncen]={"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","pp" };
  const float npart [ncen]={356  ,261    ,187     ,108     ,42       ,8.8     ,  2  };
  const float ncoll[ncen]   ={1660 ,1310  ,745 ,251 ,62.8 ,10.8    , 1   };

  double effc       [ncen]={0.998,0.998  ,0.995   ,0.991   ,0.987    ,0.977   ,0.966};

  //! be caureful 0-5:pbpb central to peripheral
  //! 6:pp

  TH1F *hjetrpt[knj][ncen];
  TH1F *hjetjpt[knj][ncen];
  TH1F *hjetspt[knj][ncen];
  TH1F *hjetspt_Syslo[knj][ncen], *hjetspt_Syshi[knj][ncen];
  TH1F *hjetraa[knj][ncen];
  TH1F *hjetraa_sm[knj][ncen], *hjetraa_sm_Sys[knj][ncen], *hjetraa_sm_Syslo[knj][ncen], *hjetraa_sm_Syshi[knj][ncen];

  TF1 *fNoise = new TF1("f","1+0.3*0.16*abs(1-([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x))");
  fNoise->SetParameters(0.9521,0.001105,-9.397e-6,3.32e-8,-5.618e-11);

  TH1F *hRatio_sm[knj][ncen];
  TH1F *hjetjptvz_pp[knj][nvtx],  *hjetjptvz_pbpb[knj][nvtx];

  for (int it=0;it<kin; it++){
    fin[it] = new TFile(Form("../Output/%s",kfile[it]),"r");
    cout<<"\t name : "<<fin[it]->GetName()<<endl;
  }

  /*
  TH1F *hvz_pp = (TH1F*)fin[1]->Get("hvz");
  hvz_pp->SetName("hVz_pp");
  hvz_pp->SetLineColor(1);
  hvz_pp->SetMarkerColor(1);
  hvz_pp->SetMarkerStyle(25);
  hvz_pp->Scale(1./hvz_pp->Integral());

  TH1F *hvz_pbpb = (TH1F*)fin[0]->Get("hvz");
  hvz_pbpb->SetName("hVz_pbpb");
  hvz_pbpb->SetLineColor(2);
  hvz_pbpb->SetMarkerColor(2);
  hvz_pbpb->SetMarkerStyle(24);
  hvz_pbpb->Scale(1./hvz_pbpb->Integral());

  hvz_pp->Draw("p");
  hvz_pbpb->Draw("psame");
  return 0;
  */

  int mt=-1;
  for(int nj=0;nj<knj;nj++){

    //! vertex dependence
    for(int iv=0;iv<nvtx;iv++){
      hjetjptvz_pp[nj][iv] = (TH1F*)fin[1]->Get(Form("hjetjptvz%d_%d",nj,iv));
      hjetjptvz_pp[nj][iv]->SetName(Form("hjetjptvz_pp_%d_%d",nj,iv));

      hjetjptvz_pp[nj][iv]->SetMarkerStyle(sty[iv]);
      hjetjptvz_pp[nj][iv]->SetMarkerColor(icol[iv]);
      hjetjptvz_pp[nj][iv]->SetLineColor(icol[iv]);

      hjetjptvz_pbpb[nj][iv] = (TH1F*)fin[0]->Get(Form("hjetjptvz%d_%d",nj,iv));
      hjetjptvz_pbpb[nj][iv]->SetName(Form("hjetjptvz_pbpb_%d_%d",nj,iv));

      hjetjptvz_pbpb[nj][iv]->SetMarkerStyle(sty[iv]);
      hjetjptvz_pbpb[nj][iv]->SetMarkerColor(icol[iv]);
      hjetjptvz_pbpb[nj][iv]->SetLineColor(icol[iv]);

      for(int ix=1;ix<=hjetjptvz_pp[nj][iv]->GetNbinsX();ix++){
	float width = hjetjptvz_pp[nj][iv]->GetBinWidth(ix);
	float yield = hjetjptvz_pp[nj][iv]->GetBinContent(ix)/width/lumi[1]/xsec[1]/1000000.;
	float err   = sqrt(hjetjptvz_pp[nj][iv]->GetBinContent(ix))/width/lumi[1]/xsec[1]/1000000.;
	hjetjptvz_pp[nj][iv]->SetBinContent(ix,yield);
	hjetjptvz_pp[nj][iv]->SetBinError(ix,err);
      }

      for(int ix=1;ix<=hjetjptvz_pbpb[nj][iv]->GetNbinsX();ix++){
	float width = hjetjptvz_pbpb[nj][iv]->GetBinWidth(ix);
	float yield = hjetjptvz_pbpb[nj][iv]->GetBinContent(ix)/width/lumi[0]/xsec[0]/1000000.;
	float err   = sqrt(hjetjptvz_pbpb[nj][iv]->GetBinContent(ix))/width/lumi[0]/xsec[0]/1000000.;
	hjetjptvz_pbpb[nj][iv]->SetBinContent(ix,yield);
	hjetjptvz_pbpb[nj][iv]->SetBinError(ix,err);
      }
    }

    mt=-1;
    for(int ic=0;ic<ncen;ic++){
      if(ic<ncen-1){
	//! lead lead
	mt=0;
	hjetrpt[nj][ic] = (TH1F*)fin[mt]->Get(Form("hjetrpt%d_%d",nj,ic));
	hjetjpt[nj][ic] = (TH1F*)fin[mt]->Get(Form("hjetjpt%d_%d",nj,ic));    

      }else{
	
	//! pp
	mt=1;
	hjetrpt[nj][ic] = (TH1F*)fin[mt]->Get(Form("hjetrpt%d_%d",nj,ic));
	hjetjpt[nj][ic] = (TH1F*)fin[mt]->Get(Form("hjetjpt%d_%d",nj,ic));    

	for(int ik=0;ik<ncen;ik++){

	  //! pp smeared pt 
	  hjetspt[nj][ik] = (TH1F*)fin[mt]->Get(Form("hjetspt%d_%d",nj,ik));    	
	  hjetspt[nj][ik]->SetMarkerStyle(sty[ik]);
	  hjetspt[nj][ik]->SetMarkerColor(icol[ik]);
	  hjetspt[nj][ik]->SetLineColor(icol[ik]);
	  hjetspt[nj][ik]->SetMarkerSize(1.3);

	  //! pp Sys lo smeared pt 
	  hjetspt_Syslo[nj][ik] = (TH1F*)fin[mt]->Get(Form("hjetspt_Syslo%d_%d",nj,ik));    	
	  hjetspt_Syslo[nj][ik]->SetMarkerStyle(sty[ik]);
	  hjetspt_Syslo[nj][ik]->SetMarkerColor(icol[ik]);
	  hjetspt_Syslo[nj][ik]->SetLineColor(icol[ik]);
	  hjetspt_Syslo[nj][ik]->SetMarkerSize(1.3);

	  //! pp Sys hi smeared pt 
	  hjetspt_Syshi[nj][ik] = (TH1F*)fin[mt]->Get(Form("hjetspt_Syshi%d_%d",nj,ik));    	
	  hjetspt_Syshi[nj][ik]->SetMarkerStyle(sty[ik]);
	  hjetspt_Syshi[nj][ik]->SetMarkerColor(icol[ik]);
	  hjetspt_Syshi[nj][ik]->SetLineColor(icol[ik]);
	  hjetspt_Syshi[nj][ik]->SetMarkerSize(1.3);


	  //cout<<"smeared : "<<ik<<"\t"<<hjetspt[nj][ik]->Integral()<<endl;
	  for(int ix=1; ix<=hjetspt[nj][ik]->GetNbinsX(); ix++) {
	    //if(nj==1 && ik==0)cout<< yield <<"\t "<<endl;

	    float width = hjetspt[nj][ik]->GetBinWidth(ix);
	    float yield = hjetspt[nj][ik]->GetBinContent(ix)/width/lumi[mt]/xsec[mt]/effc[ncen-1]/1000000.;
	    float err   = sqrt(hjetspt[nj][ik]->GetBinContent(ix))/width/lumi[mt]/xsec[mt]/effc[ncen-1]/1000000.;
	    hjetspt[nj][ik]->SetBinContent(ix,yield);
	    hjetspt[nj][ik]->SetBinError(ix,err);

	    width = hjetspt_Syslo[nj][ik]->GetBinWidth(ix);
	    yield = hjetspt_Syslo[nj][ik]->GetBinContent(ix)/width/lumi[mt]/xsec[mt]/effc[ik]/1000000.;
	    err   = sqrt(hjetspt_Syslo[nj][ik]->GetBinContent(ix))/width/lumi[mt]/xsec[mt]/effc[ncen-1]/1000000.;
	    hjetspt_Syslo[nj][ik]->SetBinContent(ix,yield);
	    hjetspt_Syslo[nj][ik]->SetBinError(ix,err);

	    width = hjetspt_Syshi[nj][ik]->GetBinWidth(ix);
	    yield = hjetspt_Syshi[nj][ik]->GetBinContent(ix)/width/lumi[mt]/xsec[mt]/effc[ik]/1000000.;
	    err   = sqrt(hjetspt_Syshi[nj][ik]->GetBinContent(ix))/width/lumi[mt]/xsec[mt]/effc[ncen-1]/1000000.;
	    hjetspt_Syshi[nj][ik]->SetBinContent(ix,yield);
	    hjetspt_Syshi[nj][ik]->SetBinError(ix,err);
	  }
	  hRatio_sm[nj][ik] = (TH1F*)hjetspt[nj][ik]->Clone(Form("hRatio_sm_%d_%d",nj,ik));
	}//! ik
      }

      hjetrpt[nj][ic]->SetMarkerStyle(sty[ic]);
      hjetrpt[nj][ic]->SetMarkerColor(icol[ic]);
      hjetrpt[nj][ic]->SetLineColor(icol[ic]);

      hjetjpt[nj][ic]->SetMarkerStyle(sty[ic]);
      hjetjpt[nj][ic]->SetMarkerColor(icol[ic]);
      hjetjpt[nj][ic]->SetLineColor(icol[ic]);
	

      for(int ix=1; ix<=hjetjpt[nj][ic]->GetNbinsX(); ix++) {
	float width = hjetjpt[nj][ic]->GetBinWidth(ix);
	float yield = hjetjpt[nj][ic]->GetBinContent(ix)/width/lumi[mt]/xsec[mt]/effc[ic]/1000000.;
	float err   = sqrt(hjetjpt[nj][ic]->GetBinContent(ix))/width/lumi[mt]/xsec[mt]/effc[ic]/1000000.;
	hjetjpt[nj][ic]->SetBinContent(ix,yield);
	hjetjpt[nj][ic]->SetBinError(ix,err);

	width = hjetrpt[nj][ic]->GetBinWidth(ix);
	yield = hjetrpt[nj][ic]->GetBinContent(ix)/width/lumi[mt]/xsec[mt]/effc[ic]/1000000.;
	err   = sqrt(hjetrpt[nj][ic]->GetBinContent(ix))/width/lumi[mt]/xsec[mt]/effc[ic]/1000000.;
	hjetrpt[nj][ic]->SetBinContent(ix,yield);
	hjetrpt[nj][ic]->SetBinError(ix,err);
      }    
      
      hjetraa[nj][ic] = (TH1F*)hjetjpt[nj][ic]->Clone(Form("hjetraa_%d_%d",nj,ic));
      hjetraa[nj][ic]->SetName(Form("hjetraa_%d_%d",nj,ic));

      hjetraa_sm[nj][ic] = (TH1F*)hjetjpt[nj][ic]->Clone(Form("hjetraa_sm_%d_%d",nj,ic));
      hjetraa_sm[nj][ic]->SetName(Form("hjetraa_sm_%d_%d",nj,ic));
      hjetraa_sm[nj][ic]->SetMarkerStyle(ssty[ic]);
      hjetraa_sm[nj][ic]->SetMarkerColor(icol[ic]);
      hjetraa_sm[nj][ic]->SetLineColor(icol[ic]);


      hjetraa_sm_Syslo[nj][ic] = (TH1F*)hjetjpt[nj][ic]->Clone(Form("hjetraa_sm_Syslo_%d_%d",nj,ic));
      hjetraa_sm_Syslo[nj][ic]->SetName(Form("hjetraa_sm_Syslo_%d_%d",nj,ic));

      hjetraa_sm_Syshi[nj][ic] = (TH1F*)hjetjpt[nj][ic]->Clone(Form("hjetraa_sm_Syshi_%d_%d",nj,ic));
      hjetraa_sm_Syshi[nj][ic]->SetName(Form("hjetraa_sm_Syshi_%d_%d",nj,ic));

    }//! ic 

    //! Calculate Raa
    for(int ic=0;ic<ncen-1;ic++){
      float scaf = 1./ncoll[ic]/0.025/csca[ic];
      hjetraa[nj][ic]->Scale(scaf);
      hjetraa[nj][ic]->Divide(hjetjpt[nj][ncen-1]);
    }

    double jecsys[ncen]={0.10,0.08,0.06,0.05,0.04,0.02,0.00};
    double ppsysm[ncen]={0.04,0.03,0.03,0.04,0.04,0.01,0.00};
    //! smeared pp Raa
    for(int ic=0;ic<ncen-1;ic++){
      float scaf = 1./ncoll[ic]/0.025/csca[ic];
      hjetraa_sm[nj][ic]->Scale(scaf);
      hjetraa_sm[nj][ic]->Divide(hjetspt[nj][ic]);

      hjetraa_sm_Syshi[nj][ic]->Scale(scaf);
      hjetraa_sm_Syshi[nj][ic]->Divide(hjetspt_Syshi[nj][ic]);

      hjetraa_sm_Syslo[nj][ic]->Scale(scaf);
      hjetraa_sm_Syslo[nj][ic]->Divide(hjetspt_Syslo[nj][ic]);
      
      hjetraa_sm_Sys[nj][ic] = (TH1F*)hjetraa_sm[nj][ic]->Clone(Form("hjetraa_sm_Sys_%d_%d",nj,ic));

      for(int ix=1;ix<=hjetraa_sm[nj][ic]->GetNbinsX();ix++){
	double lo     = fabs(hjetraa_sm_Syslo[nj][ic]->GetBinContent(ix)-hjetraa_sm[nj][ic]->GetBinContent(ix))/hjetraa_sm[nj][ic]->GetBinContent(ix);
	double hi     = fabs(hjetraa_sm_Syshi[nj][ic]->GetBinContent(ix)-hjetraa_sm[nj][ic]->GetBinContent(ix))/hjetraa_sm[nj][ic]->GetBinContent(ix);

	double sysnoise = fNoise->Integral(hjetraa_sm[nj][ic]->GetBinLowEdge(ix),hjetraa_sm[nj][ic]->GetBinLowEdge(ix+1))/hjetraa_sm[nj][ic]->GetBinWidth(ix);
	double syserr   = 1+sqrt(pow(lo,2) + pow(hi,2) + pow(0.01,2) + pow(0.01,2) + pow(jecsys[ic],2) + (pow(sysnoise,2)-1) + pow(ppsysm[ic],2));
	//double syserr = 1+sqrt(pow(lo,2) + pow(hi,2) + pow(ppsysm[ic],2)); //!
	//if(nj==1)cout<<ic<<"\t lo : "<<lo<<"\t hi : "<<hi<<"\t syserr : "<<syserr<<endl;
	hjetraa_sm_Sys[nj][ic]->SetBinContent(ix,syserr);
	hjetraa_sm_Sys[nj][ic]->SetBinError(ix,0);
      }
    }

    for(int ic=0;ic<ncen-1;ic++){
      //! Ratio
      hRatio_sm[nj][ic]->Divide(hjetjpt[nj][ncen-1]);
    }
  }//! nj 



  //! Calculate Raa
  const char *cname[knj]={"akPu2PF","akPu3PF","akPu4PF"};
  const char *mcen [ncen]={"0-5","5-10","10-30","30-50","50-70","70-90","pp" };
  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen-1;ic++){
      ofstream outf(Form("txtfiles/%s/Raa_BinByBinSmear_%s_%s.txt",cname[nj],cname[nj],mcen[ic]));
      for(int ix=1;ix<=hjetraa_sm[nj][ic]->GetNbinsX();ix++){
	double pt = hjetraa_sm[nj][ic]->GetBinCenter(ix);
	if(pt<100 || pt>400)continue;
	double raa  = hjetraa_sm[nj][ic]->GetBinContent(ix);      
	double eraa = hjetraa_sm[nj][ic]->GetBinError(ix);      
	double syserr = (hjetraa_sm_Sys[nj][ic]->GetBinContent(ix)-1)*raa;
	//if(nj==1)cout<<pt<<"\t"<<raa<<"\t"<<eraa<<"\t"<<syserr<<"\t"<<0<<"\t"<<ccen[ic]<<endl;
	outf<<pt<<"\t"<<raa<<"\t"<<eraa<<"\t"<<syserr<<"\t"<<0<<"\t"<<ccen[ic]<<endl;
	//if(nj==1 && ic==0)cout<<pt<<"\t"<<raa<<"\t"<<eraa<<"\t"<<syserr<<"\t"<<0<<"\t"<<ccen[ic]<<endl;
      }
      outf.close();
    }
  }
  return 0;


  int ipad=0;
  double xmin=94.;
  double xmax=410.;
  double ymin=5.04e-14;
  double ymax=5.04e+04;

  int ranfac[ncen] ={1e+06,1e+05,1e+04,1e+03,1e+02,1e+01,1e+00};

  TLine *line = new TLine(xmin,1,xmax,1);
  line->SetLineStyle(2);
  line->SetLineWidth(1);


  ipad=0;
  TCanvas *c1[knj];
  
  float ptmax = 300.0;
  float rymin = 0.0;
  float rymax = 2.0;
  
  TH1D *hDum, *hDum_1, *hDum2, *hDum2;    
  hDum = GetDummyHist(100.0,ptmax,rymin,rymax,"Jet p_{T} (GeV/c)","R_{AA}",false);
  hDum->GetYaxis()->SetNdivisions(610); 
  hDum->GetYaxis()->SetLabelFont(43);
  hDum->GetYaxis()->SetTitleFont(43);
  hDum->GetYaxis()->SetLabelSize(20);
  hDum->GetYaxis()->SetTitleSize(22);
  hDum->GetYaxis()->SetTitleOffset(2.6);
  hDum->GetXaxis()->SetLabelFont(43);
  hDum->GetXaxis()->SetTitleFont(43);
  hDum->GetXaxis()->SetLabelSize(20);
  hDum->GetXaxis()->SetTitleSize(22);
  hDum->GetXaxis()->SetTitleOffset(3.1);
  hDum->GetXaxis()->SetNoExponent();
  hDum->GetXaxis()->SetMoreLogLabels();

  // to show in up and down
  hDum->GetXaxis()->SetTitleOffset(2.4);
  hDum_1 = (TH1D*) hDum->Clone("hDum_1");
  hDum_1->SetAxisRange(0,2.0,"Y");
  ipad=0;
  for(int nj=1;nj<knj-1;nj++){
    c1[nj] = new TCanvas(Form("c1_%d",nj),Form("%s Jet Raa ",calgo[nj]),1100,770);
    makeMultiPanelCanvasWithGap(c1[nj],3,2,0.01,0.01,0.16,0.2,0.04,0.04); 
    ipad=0;	 

    TLegend *l07 = new TLegend(0.03393951,0.8019693,0.4599113,0.9559932,"BRNDC");
    l07->SetTextSize(0.05);
    l07->SetHeader("");
    l07->SetFillColor(10);
    l07->SetLineColor(10);
    l07->SetTextFont(42);
    l07->AddEntry(hjetraa[1][0],"PbPb/pp","p");  
    l07->AddEntry(hjetraa_sm[1][0],"PbPb/pp smeared","p");  

    for(int ic=ncen-2;ic>=0;ic--){
      TLine *line = new TLine(100.1,1,ptmax,1);
      line->SetLineStyle(2), line->SetLineWidth(2);
      
      c1[nj]->cd(++ipad);
      hDum->Draw("");

      drawText(ccen[ic],0.648173,0.8459761,22);
      if(ipad==1){
	drawText2(calgo[nj],0.2,0.85,17);  
      }
      if(ipad==3)l07->Draw();

      hjetraa[nj][ic]->SetMaximum(2.0);
      hjetraa[nj][ic]->SetMinimum(0.0);
      hjetraa_sm[nj][ic]->SetMaximum(2.0);
      hjetraa_sm[nj][ic]->SetMinimum(0.0);

      hjetraa_sm_Syshi[nj][ic]->SetMaximum(2.0);
      hjetraa_sm_Syshi[nj][ic]->SetMinimum(0.0);

      hjetraa_sm_Syslo[nj][ic]->SetMaximum(2.0);
      hjetraa_sm_Syslo[nj][ic]->SetMinimum(0.0);


      for(int ix=1;ix<=hjetraa_sm[nj][ic]->GetNbinsX();ix++){
	double val = hjetraa_sm[nj][ic]->GetBinContent(ix);
	double err1= hjetraa_sm_Sys[nj][ic]->GetBinContent(ix)-1;

	if(ic==0)cout<<"pt : "<<hjetraa_sm[nj][ic]->GetBinCenter(ix)<<"\t"<<val<<"\t"<<val*(1-err1)<<"\t"<<val*(1+err1)<<endl;
	TBox *b = new TBox(hjetraa_sm[nj][ic]->GetBinLowEdge(ix),val*(1-err1),hjetraa_sm[nj][ic]->GetBinLowEdge(ix+1),val*(1+err1));
	b->SetFillColor(kGray);
	b->SetFillStyle(1001);
	b->SetLineColor(kGray);
	b->Draw();
      }	
      th1Style1(hjetraa[nj][ic],1,20,1.2,1,1.5,1,1);    
      th1Style1(hjetraa_sm[nj][ic],1,24,1.2,1,1.5,1,1);    
      line->Draw();
      hDum->Draw("same");
    }
  }
  //return 0;
  /*
  TFile *fout = new TFile("Data_ppSmearingfit.root","RECREATE");
  fout->cd();
  for(int ic=0;ic<ncen-1;ic++){

    hjetjpt  [1][ic]->SetName(Form("hjetpt_akPu3PF_pbpb_%d",ic));
    hjetjpt  [1][ic]->SetTitle(Form("hjetpt_akPu3PF_pbpb_%s",ccen[ic]));
    hjetjpt  [1][ic]->Write();

    hjetspt  [1][ic]->SetName(Form("hjetspt_ak3PF_pp_smeared_%d",ic));
    hjetspt  [1][ic]->SetTitle(Form("hjetspt_ak3PF_pp_smeared_%s",ccen[ic]));
    hjetspt  [1][ic]->Write();


    hRatio_sm[1][ic]->SetName(Form("hRatio_pp_smearedbyunsmeared_ak3PF_%d",ic));
    hRatio_sm[1][ic]->SetTitle(Form("hRatio_pp_smearedbyunsmeared_ak3PF_%s",ccen[ic]));
    hRatio_sm[1][ic]->Write();

    hjetraa[1][ic]->SetName(Form("hjetraa_akPu3PF_%d",ic));
    hjetraa[1][ic]->SetTitle(Form("hjetraa_nosm_akPu3PF_%s",ccen[ic]));
    hjetraa[1][ic]->Write();

    hjetraa_sm[1][ic]->SetName(Form("hjetraa_sm_akPu3PF_%d",ic));
    hjetraa_sm[1][ic]->SetTitle(Form("hjetraa_sm_akPu3PF_%s",ccen[ic]));
    hjetraa_sm[1][ic]->Write();

    hjetraa_sm_Sys[nj][ic]->SetName(Form("hjetraa_sm_Sys_akPu3PF_%d",ic));
    hjetraa_sm_Sys[nj][ic]->SetTitle(Form("hjetraa_sm_Sys_akPu3PF_%d",ic));
  }
  hjetjpt[1][ncen-1]->SetName("hjetpt_ak3PF_pp");
  hjetjpt[1][ncen-1]->SetTitle("hjetpt_ak3PF_pp");
  hjetjpt[1][ncen-1]->Write();
  fout->Close();
  return 0;
  */



  /*
  const char *cvtx[nvtx] = {"-15 < vz (cm) < -7.5","-7.5 < vz (cm) < 0","0 < vz(cm) < 7.5","7.5 < vz(cm) < 15"};
  ipad=0;
  TCanvas *c4 = new TCanvas("c4","Vtx Jet Spectra",72,173,1620,770);
  makeMultiPanelCanvas(c4,3,1,0.0,0.0,0.22,0.22,0.02);
  ipad=0;
  for(int nj=0;nj<knj;nj++){
    c4->cd(++ipad);
    gPad->SetLogy();

    TLegend *l01 = 0;
    if(ipad==1){
      l01=new TLegend(0.2471927,0.7988418,0.4006227,0.9565841,"BRNDC");
      l01->SetTextSize(0.035);
    }else{
      l01=new TLegend(0.05246723,0.7989595,0.2067625,0.9573739,"BRNDC");
      l01->SetTextSize(0.045);
    }
    l01->SetHeader("pp");
    l01->SetFillColor(10);
    l01->SetLineColor(10);
    l01->SetTextFont(42);
    
    for(int iv=0;iv<nvtx;iv++){

      hjetjptvz_pp[nj][iv]->SetTitle("");
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetRangeUser(xmin,xmax);
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitle("p_{T}^{Jet} (GeV/c)");
      hjetjptvz_pp[nj][iv]->GetXaxis()->CenterTitle(true);
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetLabelFont(42);
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleFont(42);
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetMoreLogLabels();
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetNoExponent();
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleSize(0.05);
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleOffset(0.89);
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetLabelSize(0.05);
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetNdivisions(509);
      hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleOffset(1.13);
      hjetjptvz_pp[nj][iv]->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{dN}{dp_{T}}");    
      hjetjptvz_pp[nj][iv]->GetYaxis()->CenterTitle(true);
      hjetjptvz_pp[nj][iv]->GetYaxis()->SetLabelFont(42);
      hjetjptvz_pp[nj][iv]->GetYaxis()->SetTitleFont(42);
      hjetjptvz_pp[nj][iv]->GetYaxis()->SetTitleSize(0.05);
      hjetjptvz_pp[nj][iv]->GetYaxis()->SetTitleOffset(1.70);
      hjetjptvz_pp[nj][iv]->GetYaxis()->SetLabelSize(0.05);
      hjetjptvz_pp[nj][iv]->GetYaxis()->SetNdivisions(507);
      

      //hjetjptvz_pp[nj][iv]->Scale(ranfac[iv]);
      hjetjptvz_pp[nj][iv]->SetMinimum(2.345e-13);
      hjetjptvz_pp[nj][iv]->SetMaximum(2.345e-06);

      if(nj==0 && iv==0){
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetLabelOffset(0.001);
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleOffset(1.13);
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetLabelSize(0.04);
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleSize(0.04);
	hjetjptvz_pp[nj][iv]->GetYaxis()->SetLabelSize(0.04);
	hjetjptvz_pp[nj][iv]->GetYaxis()->SetTitleSize(0.04);
	hjetjptvz_pp[nj][iv]->GetYaxis()->SetTitleOffset(1.86);
      }else if (nj==1 && iv==0){
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetLabelOffset(-0.007);
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetLabelSize(0.05);
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleSize(0.05);
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleOffset(0.94);
      }else if(nj==2 &&iv==0){
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetLabelOffset(-0.007);
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetLabelSize(0.05);
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleSize(0.05);
	hjetjptvz_pp[nj][iv]->GetXaxis()->SetTitleOffset(0.92);
      }


      if(iv==0)hjetjptvz_pp[nj][iv]->Draw("p");
      else  hjetjptvz_pp[nj][iv]->Draw("psame");

      l01->AddEntry(hjetjptvz_pp[nj][iv],Form("%s",cvtx[iv]),"p");
    }

    l01->Draw();

    TPaveText *pt01   = 0;
    if(ipad==0){
      pt01 = new TPaveText(0.6207541,0.8368412,0.924797,0.9349893,"brNDC");
    }else{
      pt01 = new TPaveText(0.5803194,0.807569,0.883496,0.9384331,"brNDC");
    }
    pt01->SetBorderSize(0);
    pt01->SetFillColor(10);
    pt01->SetTextFont(42);
    TText *tx01 = pt01->AddText(Form("#splitline{CMS Preliminary}{%s}",calgo[nj]));
    if(ipad==1){
      tx01->SetTextSize(0.035);
    }else{
      tx01->SetTextSize(0.045);
    }
    pt01->Draw();
    
    if(ipad==1){
      TPaveText *pt02   = new TPaveText(0.6714279,0.7111428,0.9585795,0.8024033,"brNDC");
      pt02->SetBorderSize(0);
      pt02->SetFillColor(10);
      pt02->SetTextFont(42);
      TText *tx03 = pt02->AddText(Form("#int L^{pp} dt   = %0.0f nb^{-1}",lumi[1]));
      tx03->SetTextSize(0.028);
      pt02->Draw();
    }
  }//! nj


  ipad=0;
  TCanvas *c5 = new TCanvas("c5","Vtx pbpb Jet Spectra",72,173,1620,770);
  makeMultiPanelCanvas(c5,3,1,0.0,0.0,0.22,0.22,0.02);
  ipad=0;
  for(int nj=0;nj<knj;nj++){
    c5->cd(++ipad);
    gPad->SetLogy();

    TLegend *l01 = 0;
    if(ipad==1){
      l01=new TLegend(0.2471927,0.7988418,0.4006227,0.9565841,"BRNDC");
      l01->SetTextSize(0.035);
    }else{
      l01=new TLegend(0.05246723,0.7989595,0.2067625,0.9573739,"BRNDC");
      l01->SetTextSize(0.045);
    }
    l01->SetHeader("PbPb 0-90%");
    l01->SetFillColor(10);
    l01->SetLineColor(10);
    l01->SetTextFont(42);
    
    for(int iv=0;iv<nvtx;iv++){

      hjetjptvz_pbpb[nj][iv]->SetTitle("");
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetRangeUser(xmin,xmax);
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitle("p_{T}^{Jet} (GeV/c)");
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->CenterTitle(true);
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetLabelFont(42);
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleFont(42);
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetMoreLogLabels();
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetNoExponent();
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleSize(0.05);
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleOffset(0.89);
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetLabelSize(0.05);
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetNdivisions(509);
      hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleOffset(1.13);
      hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{dN}{dp_{T}}");    
      hjetjptvz_pbpb[nj][iv]->GetYaxis()->CenterTitle(true);
      hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetLabelFont(42);
      hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetTitleFont(42);
      hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetTitleSize(0.05);
      hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetTitleOffset(1.70);
      hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetLabelSize(0.05);
      hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetNdivisions(507);
      

      //hjetjptvz_pbpb[nj][iv]->Scale(ranfac[iv]);
      hjetjptvz_pbpb[nj][iv]->SetMinimum(2.345e-13);
      hjetjptvz_pbpb[nj][iv]->SetMaximum(2.345e-03);

      if(nj==0 && iv==0){
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetLabelOffset(0.001);
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleOffset(1.13);
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetLabelSize(0.04);
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleSize(0.04);
	hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetLabelSize(0.04);
	hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetTitleSize(0.04);
	hjetjptvz_pbpb[nj][iv]->GetYaxis()->SetTitleOffset(1.86);
      }else if (nj==1 && iv==0){
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetLabelOffset(-0.007);
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetLabelSize(0.05);
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleSize(0.05);
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleOffset(0.94);
      }else if(nj==2 &&iv==0){
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetLabelOffset(-0.007);
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetLabelSize(0.05);
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleSize(0.05);
	hjetjptvz_pbpb[nj][iv]->GetXaxis()->SetTitleOffset(0.92);
      }


      if(iv==0)hjetjptvz_pbpb[nj][iv]->Draw("p");
      else  hjetjptvz_pbpb[nj][iv]->Draw("psame");

      l01->AddEntry(hjetjptvz_pbpb[nj][iv],Form("%s",cvtx[iv]),"p");
    }

    l01->Draw();

    TPaveText *pt01   = 0;
    if(ipad==0){
      pt01 = new TPaveText(0.6207541,0.8368412,0.924797,0.9349893,"brNDC");
    }else{
      pt01 = new TPaveText(0.5803194,0.807569,0.883496,0.9384331,"brNDC");
    }
    pt01->SetBorderSize(0);
    pt01->SetFillColor(10);
    pt01->SetTextFont(42);
    TText *tx01 = pt01->AddText(Form("#splitline{CMS Preliminary}{%s}",calgo[nj]));
    if(ipad==1){
      tx01->SetTextSize(0.035);
    }else{
      tx01->SetTextSize(0.045);
    }
    pt01->Draw();
    
    if(ipad==1){
      TPaveText *pt02   = new TPaveText(0.6714279,0.7111428,0.9585795,0.8024033,"brNDC");
      pt02->SetBorderSize(0);
      pt02->SetFillColor(10);
      pt02->SetTextFont(42);
      TText *tx04 = pt02->AddText(Form("#int L^{PbPb} dt = %0.0f #mub^{-1}",lumi[0]));
      tx04->SetTextSize(0.028);
      pt02->Draw();
    }
  }//! nj
  */


  

  ipad=0;
  TCanvas *c0 = new TCanvas("c0","Jet Spectra",72,173,1620,770);
  makeMultiPanelCanvas(c0,3,1,0.0,0.0,0.22,0.22,0.02);
  for(int nj=0;nj<knj;nj++){
    c0->cd(++ipad);
    gPad->SetLogy();

    TLegend *l01 = 0;
    TLegend *l00 = 0;
    if(ipad==1){
      l01=new TLegend(0.2471927,0.7988418,0.4006227,0.9565841,"BRNDC");
      l01->SetTextSize(0.035);
      l00 = new TLegend(0.3990562,0.8316755,0.5869716,0.9539301,"BRNDC");
      l00->SetTextSize(0.035);
    }else{
      l01=new TLegend(0.05246723,0.7989595,0.2067625,0.9573739,"BRNDC");
      l01->SetTextSize(0.045);
      l00 = new TLegend(0.2419526,0.8333974,0.431438,0.955652,"BRNDC");
      l00->SetTextSize(0.045);
    }
    l01->SetHeader("");
    l01->SetFillColor(10);
    l01->SetLineColor(10);
    l01->SetTextFont(42);

    l00->SetHeader("");
    l00->SetFillColor(10);
    l00->SetLineColor(10);
    l00->SetTextFont(42);

    //if(ipad==1)gPad->SetLeftMargin(0.21);
    //if(ipad==3)gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.13);
    
    for(int ic=0;ic<ncen;ic++){
      hjetjpt[nj][ic]->SetTitle("");
      hjetjpt[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hjetjpt[nj][ic]->GetXaxis()->SetTitle("p_{T}^{Jet} (GeV/c)");
      hjetjpt[nj][ic]->GetXaxis()->CenterTitle(true);
      hjetjpt[nj][ic]->GetXaxis()->SetLabelFont(42);
      hjetjpt[nj][ic]->GetXaxis()->SetTitleFont(42);
      hjetjpt[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hjetjpt[nj][ic]->GetXaxis()->SetNoExponent();
      hjetjpt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
      hjetjpt[nj][ic]->GetXaxis()->SetTitleOffset(0.89);
      hjetjpt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
      hjetjpt[nj][ic]->GetXaxis()->SetNdivisions(509);
      hjetjpt[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
      hjetjpt[nj][ic]->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{dN}{dp_{T}}");    
      hjetjpt[nj][ic]->GetYaxis()->CenterTitle(true);
      hjetjpt[nj][ic]->GetYaxis()->SetLabelFont(42);
      hjetjpt[nj][ic]->GetYaxis()->SetTitleFont(42);
      hjetjpt[nj][ic]->GetYaxis()->SetTitleSize(0.05);
      hjetjpt[nj][ic]->GetYaxis()->SetTitleOffset(1.70);
      hjetjpt[nj][ic]->GetYaxis()->SetLabelSize(0.05);
      hjetjpt[nj][ic]->GetYaxis()->SetNdivisions(507);


      hjetjpt[nj][ic]->Scale(ranfac[ic]);
      hjetjpt[nj][ic]->SetMinimum(ymin);
      hjetjpt[nj][ic]->SetMaximum(ymax);

      if(nj==0){
	hjetjpt[nj][ic]->GetXaxis()->SetLabelOffset(0.001);
	hjetjpt[nj][ic]->GetXaxis()->SetLabelSize(0.04);
	hjetjpt[nj][ic]->GetXaxis()->SetTitleSize(0.04);
	hjetjpt[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
	hjetjpt[nj][ic]->GetYaxis()->SetLabelSize(0.04);
	hjetjpt[nj][ic]->GetYaxis()->SetTitleSize(0.04);
	hjetjpt[nj][ic]->GetYaxis()->SetTitleOffset(1.86);
      }else if (nj==1){
	hjetjpt[nj][ic]->GetXaxis()->SetLabelOffset(-0.007);
	hjetjpt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
	hjetjpt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
	hjetjpt[nj][ic]->GetXaxis()->SetTitleOffset(0.94);
      }else if(nj==2){
	hjetjpt[nj][ic]->GetXaxis()->SetLabelOffset(-0.007);
	hjetjpt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
	hjetjpt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
	hjetjpt[nj][ic]->GetXaxis()->SetTitleOffset(0.92);
      }


      if(ic==0)hjetjpt[nj][ic]->Draw("p");
      else  hjetjpt[nj][ic]->Draw("psame");

      if(ic<4)l01->AddEntry(hjetjpt[nj][ic],Form("%s",ccen[ic]),"p");
      else l00->AddEntry(hjetjpt[nj][ic],Form("%s",ccen[ic]),"p");
    }

    l01->Draw();
    l00->Draw();

    TPaveText *pt01   = 0;
    if(ipad==0){
      pt01 = new TPaveText(0.6207541,0.8368412,0.924797,0.9349893,"brNDC");
    }else{
      pt01 = new TPaveText(0.5803194,0.807569,0.883496,0.9384331,"brNDC");
    }
    pt01->SetBorderSize(0);
    pt01->SetFillColor(10);
    pt01->SetTextFont(42);
    TText *tx01 = pt01->AddText(Form("#splitline{CMS Preliminary}{%s}",calgo[nj]));
    if(ipad==1){
      tx01->SetTextSize(0.035);
    }else{
      tx01->SetTextSize(0.045);
    }
    pt01->Draw();
    
    if(ipad==1){
      TPaveText *pt02   = new TPaveText(0.6714279,0.7111428,0.9585795,0.8024033,"brNDC");
      pt02->SetBorderSize(0);
      pt02->SetFillColor(10);
      pt02->SetTextFont(42);
      TText *tx03 = pt02->AddText(Form("#int L^{pp} dt   = %0.0f nb^{-1}",lumi[1]));
      TText *tx04 = pt02->AddText(Form("#int L^{PbPb} dt = %0.0f #mub^{-1}",lumi[0]));
      tx03->SetTextSize(0.028);
      tx04->SetTextSize(0.028);
      pt02->Draw();
    }
  }//! nj


  ipad=0;
  TCanvas *c6 = new TCanvas("c6","Raw Jet Spectra",72,173,1620,770);
  makeMultiPanelCanvas(c6,3,1,0.0,0.0,0.22,0.22,0.02);
  for(int nj=0;nj<knj;nj++){
    c6->cd(++ipad);
    gPad->SetLogy();

    TLegend *l01 = 0;
    TLegend *l00 = 0;
    if(ipad==1){
      l01=new TLegend(0.2471927,0.7988418,0.4006227,0.9565841,"BRNDC");
      l01->SetTextSize(0.035);
      l00 = new TLegend(0.3990562,0.8316755,0.5869716,0.9539301,"BRNDC");
      l00->SetTextSize(0.035);
    }else{
      l01=new TLegend(0.05246723,0.7989595,0.2067625,0.9573739,"BRNDC");
      l01->SetTextSize(0.045);
      l00 = new TLegend(0.2419526,0.8333974,0.431438,0.955652,"BRNDC");
      l00->SetTextSize(0.045);
    }
    l01->SetHeader("");
    l01->SetFillColor(10);
    l01->SetLineColor(10);
    l01->SetTextFont(42);

    l00->SetHeader("");
    l00->SetFillColor(10);
    l00->SetLineColor(10);
    l00->SetTextFont(42);

    //if(ipad==1)gPad->SetLeftMargin(0.21);
    //if(ipad==3)gPad->SetRightMargin(0.01);
    //gPad->SetBottomMargin(0.13);
    
    for(int ic=0;ic<ncen;ic++){
      hjetrpt[nj][ic]->SetTitle("");
      hjetrpt[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hjetrpt[nj][ic]->GetXaxis()->SetTitle("p_{T}^{Jet} (GeV/c)");
      hjetrpt[nj][ic]->GetXaxis()->CenterTitle(true);
      hjetrpt[nj][ic]->GetXaxis()->SetLabelFont(42);
      hjetrpt[nj][ic]->GetXaxis()->SetTitleFont(42);
      hjetrpt[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hjetrpt[nj][ic]->GetXaxis()->SetNoExponent();
      hjetrpt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
      hjetrpt[nj][ic]->GetXaxis()->SetTitleOffset(0.89);
      hjetrpt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
      hjetrpt[nj][ic]->GetXaxis()->SetNdivisions(509);
      hjetrpt[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
      hjetrpt[nj][ic]->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{dN}{dp_{T}}");    
      hjetrpt[nj][ic]->GetYaxis()->CenterTitle(true);
      hjetrpt[nj][ic]->GetYaxis()->SetLabelFont(42);
      hjetrpt[nj][ic]->GetYaxis()->SetTitleFont(42);
      hjetrpt[nj][ic]->GetYaxis()->SetTitleSize(0.05);
      hjetrpt[nj][ic]->GetYaxis()->SetTitleOffset(1.70);
      hjetrpt[nj][ic]->GetYaxis()->SetLabelSize(0.05);
      hjetrpt[nj][ic]->GetYaxis()->SetNdivisions(507);


      hjetrpt[nj][ic]->Scale(ranfac[ic]);
      hjetrpt[nj][ic]->SetMinimum(ymin);
      hjetrpt[nj][ic]->SetMaximum(ymax);

      if(nj==0){
	hjetrpt[nj][ic]->GetXaxis()->SetLabelOffset(0.001);
	hjetrpt[nj][ic]->GetXaxis()->SetLabelSize(0.04);
	hjetrpt[nj][ic]->GetXaxis()->SetTitleSize(0.04);
	hjetrpt[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
	hjetrpt[nj][ic]->GetYaxis()->SetLabelSize(0.04);
	hjetrpt[nj][ic]->GetYaxis()->SetTitleSize(0.04);
	hjetrpt[nj][ic]->GetYaxis()->SetTitleOffset(1.86);
      }else if (nj==1){
	hjetrpt[nj][ic]->GetXaxis()->SetLabelOffset(-0.007);
	hjetrpt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
	hjetrpt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
	hjetrpt[nj][ic]->GetXaxis()->SetTitleOffset(0.94);
      }else if(nj==2){
	hjetrpt[nj][ic]->GetXaxis()->SetLabelOffset(-0.007);
	hjetrpt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
	hjetrpt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
	hjetrpt[nj][ic]->GetXaxis()->SetTitleOffset(0.92);
      }


      if(ic==0)hjetrpt[nj][ic]->Draw("p");
      else  hjetrpt[nj][ic]->Draw("psame");

      if(ic<4)l01->AddEntry(hjetrpt[nj][ic],Form("%s",ccen[ic]),"p");
      else l00->AddEntry(hjetrpt[nj][ic],Form("%s",ccen[ic]),"p");
    }

    l01->Draw();
    l00->Draw();

    TPaveText *pt01   = 0;
    if(ipad==0){
      pt01 = new TPaveText(0.6207541,0.8368412,0.924797,0.9349893,"brNDC");
    }else{
      pt01 = new TPaveText(0.5803194,0.807569,0.883496,0.9384331,"brNDC");
    }
    pt01->SetBorderSize(0);
    pt01->SetFillColor(10);
    pt01->SetTextFont(42);
    TText *tx01 = pt01->AddText(Form("#splitline{CMS Preliminary}{%s}",calgo[nj]));
    if(ipad==1){
      tx01->SetTextSize(0.035);
    }else{
      tx01->SetTextSize(0.045);
    }
    pt01->Draw();
    
    if(ipad==1){
      TPaveText *pt02   = new TPaveText(0.6714279,0.7111428,0.9585795,0.8024033,"brNDC");
      pt02->SetBorderSize(0);
      pt02->SetFillColor(10);
      pt02->SetTextFont(42);
      TText *tx03 = pt02->AddText(Form("#int L^{pp} dt   = %0.0f nb^{-1}",lumi[1]));
      TText *tx04 = pt02->AddText(Form("#int L^{PbPb} dt = %0.0f #mub^{-1}",lumi[0]));
      tx03->SetTextSize(0.028);
      tx04->SetTextSize(0.028);
      pt02->Draw();
    }
  }//! nj
  

  ipad=0;
  TCanvas *c2[knj];
  for(int nj=0;nj<knj;nj++){
    ipad=0;
    c2[nj] = new TCanvas(Form("c2_%d",nj),Form("%s Smeared Jet Spectra",calgo[nj]),50,164,1742,541);
    makeMultiPanelCanvas(c2[nj],ncen-1,2,0.0,0.0,0.22,0.22,0.02);
    TLegend *l07 = new TLegend(0.6002655,0.6005989,0.8107852,0.8023582,"BRNDC");
    l07->SetTextSize(0.07);
    l07->SetHeader("");
    l07->SetFillColor(10);
    l07->SetLineColor(10);
    l07->SetTextFont(42);
    l07->AddEntry(hjetspt[nj][0],"pp Smeared","p");  
    l07->AddEntry(hjetrpt[nj][ncen-1],"pp","p");  

    for(int ic=ncen-2;ic>=0;ic--){
      hjetspt[nj][ic]->SetMaximum(5.435e-06);
      hjetspt[nj][ic]->SetMinimum(5.435e-12);
      //hjetspt[nj][ic]->SetTitle(ccen[ic]);
      hjetspt[nj][ic]->SetTitle("");
      hjetspt[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hjetspt[nj][ic]->GetXaxis()->SetTitle("p_{T}^{Jet} (GeV/c)");
      hjetspt[nj][ic]->GetXaxis()->CenterTitle(true);
      hjetspt[nj][ic]->GetXaxis()->SetLabelFont(42);
      hjetspt[nj][ic]->GetXaxis()->SetTitleFont(42);
      hjetspt[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hjetspt[nj][ic]->GetXaxis()->SetNoExponent();
      hjetspt[nj][ic]->GetXaxis()->SetTitleSize(0.07);
      hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(0.89);
      hjetspt[nj][ic]->GetXaxis()->SetLabelSize(0.07);
      hjetspt[nj][ic]->GetXaxis()->SetNdivisions(509);
      hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
      hjetspt[nj][ic]->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{dN}{dp_{T}}");    
      hjetspt[nj][ic]->GetYaxis()->CenterTitle(true);
      hjetspt[nj][ic]->GetYaxis()->SetLabelFont(42);
      hjetspt[nj][ic]->GetYaxis()->SetTitleFont(42);
      hjetspt[nj][ic]->GetYaxis()->SetTitleSize(0.08);
      hjetspt[nj][ic]->GetYaxis()->SetTitleOffset(1.22);
      hjetspt[nj][ic]->GetYaxis()->SetLabelSize(0.08);
      hjetspt[nj][ic]->GetYaxis()->SetNdivisions(507);

      if(nj==0){
	hjetspt[nj][ic]->GetXaxis()->SetLabelOffset(0.010);
	hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
	hjetspt[nj][ic]->GetYaxis()->SetTitleSize(0.05);
	hjetspt[nj][ic]->GetYaxis()->SetTitleOffset(1.86);
      }else if (nj==1){
	hjetspt[nj][ic]->GetXaxis()->SetLabelOffset(0.010);
	hjetspt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
	hjetspt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
	hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(0.94);
      }else if(nj==2){
	hjetspt[nj][ic]->GetXaxis()->SetLabelOffset(0.010);
	hjetspt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
	hjetspt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
	hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(0.92);
      }
    
      c2[nj]->cd(++ipad);
      gPad->SetLogy();
      hjetspt[nj][ic]->SetMarkerStyle(30);
      hjetspt[nj][ic]->SetMarkerColor(2);
      hjetspt[nj][ic]->SetLineColor(2);
      hjetspt[nj][ic]->SetMarkerSize(1.2);

      hjetspt[nj][ncen-1]->SetMarkerStyle(20);
      hjetspt[nj][ncen-1]->SetMarkerColor(1);
      hjetspt[nj][ncen-1]->SetLineColor(1);
      hjetspt[nj][ncen-1]->SetMarkerSize(1.0);

      hjetspt[nj][ic]->Draw("p");
      hjetjpt[nj][ncen-1]->Draw("psame");
      
      if(ipad==1){
	TPaveText *pt03 = new TPaveText(0.4186168,0.837257,0.7240939,0.9340231,"brNDC");
	pt03->SetBorderSize(0);
	pt03->SetFillColor(10);
	pt03->SetTextFont(42);
	TText *tx03 = pt03->AddText(Form("%s",calgo[nj]));
	tx03->SetTextSize(0.09);
	pt03->Draw();
	l07->Draw();
      }

      if(ipad==2){
	TPaveText *pt02   = new TPaveText(0.299003,0.7053031,0.5854784,0.8680462,"brNDC");
	pt02->SetBorderSize(0);
	pt02->SetFillColor(10);
	pt02->SetTextFont(42);
	TText *tx03 = pt02->AddText(Form("#int L^{pp} dt   = %0.0f nb^{-1}",lumi[1]));
	tx03->SetTextSize(0.08);
	pt02->Draw();
      }

      TPaveText *pt04 = new TPaveText(0.2517358,0.06752635,0.5572128,0.1642925,"brNDC");
      pt04->SetBorderSize(0);
      pt04->SetFillColor(10);
      pt04->SetTextFont(42);
      TText *tx04 = pt04->AddText(ccen[ic]);
      tx04->SetTextSize(0.09);
      pt04->Draw();

      c2[nj]->cd(ipad+(ncen-1));
      hRatio_sm[nj][ic]->SetMaximum(1.48);
      hRatio_sm[nj][ic]->SetMinimum(0.58);
      //hRatio_sm[nj][ic]->SetTitle(ccen[ic]);
      hRatio_sm[nj][ic]->SetTitle("");
      hRatio_sm[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hRatio_sm[nj][ic]->GetXaxis()->SetTitle("p_{T}^{Jet} (GeV/c)");
      hRatio_sm[nj][ic]->GetXaxis()->CenterTitle(true);
      hRatio_sm[nj][ic]->GetXaxis()->SetLabelFont(42);
      hRatio_sm[nj][ic]->GetXaxis()->SetTitleFont(42);
      hRatio_sm[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hRatio_sm[nj][ic]->GetXaxis()->SetNoExponent();
      hRatio_sm[nj][ic]->GetXaxis()->SetTitleSize(0.07);
      hRatio_sm[nj][ic]->GetXaxis()->SetTitleOffset(0.89);
      hRatio_sm[nj][ic]->GetXaxis()->SetLabelSize(0.07);
      hRatio_sm[nj][ic]->GetXaxis()->SetNdivisions(509);
      hRatio_sm[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
      hRatio_sm[nj][ic]->GetYaxis()->SetTitle("Ratio");    
      hRatio_sm[nj][ic]->GetYaxis()->CenterTitle(true);
      hRatio_sm[nj][ic]->GetYaxis()->SetLabelFont(42);
      hRatio_sm[nj][ic]->GetYaxis()->SetTitleFont(42);
      hRatio_sm[nj][ic]->GetYaxis()->SetTitleSize(0.07);
      hRatio_sm[nj][ic]->GetYaxis()->SetTitleOffset(1.48);
      hRatio_sm[nj][ic]->GetYaxis()->SetLabelSize(0.07);
      hRatio_sm[nj][ic]->GetYaxis()->SetNdivisions(507);

      hRatio_sm[nj][ic]->SetMarkerStyle(30);
      hRatio_sm[nj][ic]->SetMarkerColor(2);
      hRatio_sm[nj][ic]->SetLineColor(2);
      hRatio_sm[nj][ic]->SetMarkerSize(1.2);
      hRatio_sm[nj][ic]->Draw("p");
      line->Draw();
    }//! ic
  }//! nj loop

  if(iSave){
    c2[1]->SaveAs("../mc/Smearing/NewForest/pbpb/AN/Smearing/Data_akPu3PF_Fit_ppSmeared.png");
    c2[1]->SaveAs("../mc/Smearing/NewForest/pbpb/AN/Smearing/Data_akPu3PF_Fit_ppSmeared.eps");
    c2[1]->SaveAs("../mc/Smearing/NewForest/pbpb/AN/Smearing/Data_akPu3PF_Fit_ppSmeared.C");
    c2[1]->SaveAs("../mc/Smearing/NewForest/pbpb/AN/Smearing/Data_akPu3PF_Fit_ppSmeared.pdf");
  }

  /*
  ipad=0;
  TCanvas *c2 = new TCanvas("c2","Smeared Jet Spectra",24,97,1744,866);
  makeMultiPanelCanvas(c2,ncen-1,3,0.0,0.0,0.22,0.22,0.02);
  TLegend *l07 = new TLegend(0.6554558,0.5520004,0.8673488,0.756689,"BRNDC");
  l07->SetTextSize(0.07);
  l07->SetHeader("");
  l07->SetFillColor(10);
  l07->SetLineColor(10);
  l07->SetTextFont(42);
  l07->AddEntry(hjetspt[0][0],"pp Smeared","p");  
  l07->AddEntry(hjetrpt[0][ncen-1],"pp","p");  
  

  for(int ic=0;ic<ncen-1;ic++){
    for(int nj=0;nj<knj;nj++){
      hjetspt[nj][ic]->SetMaximum(1e-06);
      hjetspt[nj][ic]->SetMinimum(1e-12);
      hjetspt[nj][ic]->SetTitle(ccen[ic]);
      hjetspt[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hjetspt[nj][ic]->GetXaxis()->SetTitle("p_{T}^{Jet} (GeV/c)");
      hjetspt[nj][ic]->GetXaxis()->CenterTitle(true);
      hjetspt[nj][ic]->GetXaxis()->SetLabelFont(42);
      hjetspt[nj][ic]->GetXaxis()->SetTitleFont(42);
      hjetspt[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hjetspt[nj][ic]->GetXaxis()->SetNoExponent();
      hjetspt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
      hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(0.89);
      hjetspt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
      hjetspt[nj][ic]->GetXaxis()->SetNdivisions(509);
      hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
      hjetspt[nj][ic]->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{dN}{dp_{T}}");    
      hjetspt[nj][ic]->GetYaxis()->CenterTitle(true);
      hjetspt[nj][ic]->GetYaxis()->SetLabelFont(42);
      hjetspt[nj][ic]->GetYaxis()->SetTitleFont(42);
      hjetspt[nj][ic]->GetYaxis()->SetTitleSize(0.05);
      hjetspt[nj][ic]->GetYaxis()->SetTitleOffset(1.70);
      hjetspt[nj][ic]->GetYaxis()->SetLabelSize(0.05);
      hjetspt[nj][ic]->GetYaxis()->SetNdivisions(507);

      if(nj==0){
	hjetspt[nj][ic]->GetXaxis()->SetLabelOffset(0.001);
	hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
	hjetspt[nj][ic]->GetYaxis()->SetTitleSize(0.05);
	hjetspt[nj][ic]->GetYaxis()->SetTitleOffset(1.86);
      }else if (nj==1){
	hjetspt[nj][ic]->GetXaxis()->SetLabelOffset(-0.007);
	hjetspt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
	hjetspt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
	hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(0.94);
      }else if(nj==2){
	hjetspt[nj][ic]->GetXaxis()->SetLabelOffset(-0.007);
	hjetspt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
	hjetspt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
	hjetspt[nj][ic]->GetXaxis()->SetTitleOffset(0.92);
      }
    }

    c2->cd(++ipad);
    gPad->SetLogy();
    hjetspt[0][ic]->Draw("p");
    hjetjpt[0][ncen-1]->Draw("psame");

    if(ipad==1){
      TPaveText *pt03 = new TPaveText(0.2483553,0.8654867,0.5526316,0.9646018,"brNDC");
      pt03->SetBorderSize(0);
      pt03->SetFillColor(10);
      pt03->SetTextFont(42);
      TText *tx03 = pt03->AddText(Form("%s",calgo[0]));
      tx03->SetTextSize(0.06);
      pt03->Draw();
      l07->Draw();
    }
    c2->cd(ipad+(ncen-1));
    gPad->SetLogy();


    hjetspt[1][ic]->Draw("p");
    hjetjpt[1][ncen-1]->Draw("psame");

    if(ipad==1){
      TPaveText *pt03 = new TPaveText(0.2483553,0.8654867,0.5526316,0.9646018,"brNDC");
      pt03->SetBorderSize(0);
      pt03->SetFillColor(10);
      pt03->SetTextFont(42);
      TText *tx03 = pt03->AddText(Form("%s",calgo[1]));
      tx03->SetTextSize(0.06);
      pt03->Draw();
    }
    c2->cd(ipad+2*(ncen-1));
    gPad->SetLogy();

    hjetspt[1][ic]->Draw("p");
    hjetjpt[1][ncen-1]->Draw("psame");

    if(ipad==1){
      TPaveText *pt03 = new TPaveText(0.2483553,0.8654867,0.5526316,0.9646018,"brNDC");
      pt03->SetBorderSize(0);
      pt03->SetFillColor(10);
      pt03->SetTextFont(42);
      TText *tx03 = pt03->AddText(Form("%s",calgo[2]));
      tx03->SetTextSize(0.06);
      pt03->Draw();
    }
  }

  ipad=0;
  TCanvas *c3 = new TCanvas("c3","Ratio Smeared Jet Spectra",24,97,1744,866);
  makeMultiPanelCanvas(c3,ncen-1,3,0.0,0.0,0.22,0.22,0.02);
  for(int ic=0;ic<ncen-1;ic++){
    for(int nj=0;nj<knj;nj++){
      hRatio_sm[nj][ic]->SetMaximum(2.45);
      hRatio_sm[nj][ic]->SetMinimum(-0.05);
      hRatio_sm[nj][ic]->SetTitle(ccen[ic]);
      hRatio_sm[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hRatio_sm[nj][ic]->GetXaxis()->SetTitle("p_{T}^{Jet} (GeV/c)");
      hRatio_sm[nj][ic]->GetXaxis()->CenterTitle(true);
      hRatio_sm[nj][ic]->GetXaxis()->SetLabelFont(42);
      hRatio_sm[nj][ic]->GetXaxis()->SetTitleFont(42);
      hRatio_sm[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hRatio_sm[nj][ic]->GetXaxis()->SetNoExponent();
      hRatio_sm[nj][ic]->GetXaxis()->SetTitleSize(0.05);
      hRatio_sm[nj][ic]->GetXaxis()->SetTitleOffset(0.89);
      hRatio_sm[nj][ic]->GetXaxis()->SetLabelSize(0.05);
      hRatio_sm[nj][ic]->GetXaxis()->SetNdivisions(509);
      hRatio_sm[nj][ic]->GetXaxis()->SetTitleOffset(1.13);
      hRatio_sm[nj][ic]->GetYaxis()->SetTitle("ppSmeared / pp");    
      hRatio_sm[nj][ic]->GetYaxis()->CenterTitle(true);
      hRatio_sm[nj][ic]->GetYaxis()->SetLabelFont(42);
      hRatio_sm[nj][ic]->GetYaxis()->SetTitleFont(42);
      hRatio_sm[nj][ic]->GetYaxis()->SetTitleSize(0.05);
      hRatio_sm[nj][ic]->GetYaxis()->SetTitleOffset(1.70);
      hRatio_sm[nj][ic]->GetYaxis()->SetLabelSize(0.05);
      hRatio_sm[nj][ic]->GetYaxis()->SetNdivisions(507);
    }

    c3->cd(++ipad);
    hRatio_sm[0][ic]->Draw("p");
    line->Draw();
    if(ipad==1){
      TPaveText *pt03 = new TPaveText(0.2483553,0.8654867,0.5526316,0.9646018,"brNDC");
      pt03->SetBorderSize(0);
      pt03->SetFillColor(10);
      pt03->SetTextFont(42);
      TText *tx03 = pt03->AddText(Form("%s",calgo[0]));
      tx03->SetTextSize(0.06);
      pt03->Draw();
    }

    c3->cd(ipad+(ncen-1));
    hRatio_sm[1][ic]->Draw("p");
    line->Draw();
    if(ipad==1){
      TPaveText *pt03 = new TPaveText(0.2483553,0.8654867,0.5526316,0.9646018,"brNDC");
      pt03->SetBorderSize(0);
      pt03->SetFillColor(10);
      pt03->SetTextFont(42);
      TText *tx03 = pt03->AddText(Form("%s",calgo[1]));
      tx03->SetTextSize(0.06);
      pt03->Draw();
    }

    c3->cd(ipad+2*(ncen-1));
    hRatio_sm[2][ic]->Draw("p");
    line->Draw();
    if(ipad==1){
      TPaveText *pt03 = new TPaveText(0.2483553,0.8654867,0.5526316,0.9646018,"brNDC");
      pt03->SetBorderSize(0);
      pt03->SetFillColor(10);
      pt03->SetTextFont(42);
      TText *tx03 = pt03->AddText(Form("%s",calgo[2]));
      tx03->SetTextSize(0.06);
      pt03->Draw();
    }

  }
  */
  return 0;
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

void drawText2(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
