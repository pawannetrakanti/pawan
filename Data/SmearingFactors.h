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
TF1 *fReWe; //! only for 0-10% and 50-100% and ak3PF pp jets


//! For pp we are using jet algorithm w/o pileup subtraction.
//! For PbPb we are suing jet algorithm w pileup subtraction



//!                             0-5%    5-10%  10-30%  30-50%   50-70%  70-90%   pp 
double resol[KNJ][NCEN][3]={ {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo
			     {
			       {0.0353701 ,1.56192 ,4.39179}, //! 0-5%
			       {0.0302473 ,1.5641  ,3.16197}, //! 5-10%
			       {0.0495822 ,1.35557 ,5.26617}, //! 10-30%
			       {0.04471   ,1.36655 ,3.54741}, //! 30-50%
			       {0.0539697 ,1.27212 ,4.05015}, //! 50-70%
			       {0.0609397 ,1.16975 ,4.71212}, //! 70-90%
			       {0.054108  ,1.32498 ,2.58301} //! pp
			     }, //! ak2PF
			     {
			       {0.0318181 ,1.49645 ,5.87603},//! 0-5% ak3PF
			       {0.0429251 ,1.38785 ,6.1478}, //! 5-10%  
			       {0.0466294 ,1.32239 ,5.0623}, //! 10-30%  
			       {0.0450478 ,1.26727 ,3.54234}, //! 30-50%  
			       {0.0525918 ,1.15779 ,4.05004}, //! 50-70%  
			       {0.0587943 ,1.04695 ,4.48503}, //! 70-90%	  
			       {0.053795  ,1.19538 ,0.00106437}  //!  pp 
			     },//! ak3PF
			     {
			       {0.0279665 ,1.4998  ,7.61268}, //! 0-5%
			       {0.0439427 ,1.35271 ,7.72843}, //! 5-10%
			       {0.0378725 ,1.35547 ,6.05271}, //! 10-30%
			       {0.0459776 ,1.22252 ,5.105},   //! 30-50%
			       {0.0504345 ,1.1264  ,4.83181}, //! 50-70% 
			       {0.055947  ,1.01899 ,5.01181}, //! 70-90%
			       {0.0534341 ,1.13202 ,0.00135424} //! pp
			     }, //! ak4PF
			     {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2Calo
			     {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak3Calo
			     {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}  //! ak4Calo
};


double scale[KNJ][NCEN][3]={ {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! icpu5calo
                             {
			       {1.00317  ,-3.76187  ,20.4157},  //! 0-5%
			       {1.00225  ,-3.10271  ,-9.67917}, //! 5-10%
			       {1.00165  ,-3.1616   ,23.7143},  //! 10-30%
			       {0.997761 ,-1.19316  ,-42.3774}, //! 30-50%
			       {1.00065  ,-1.41227  ,-6.27402}, //! 50-70%
			       {0.998528 ,-0.245675 ,-21.9381}, //! 70-90%
			       {1,0,0}  //! pp
			     }, //! ak2PF
                             {
			       {1.00624 ,-4.00334 ,143.305}, //! 0-5% 
			       {1.00363 ,-3.49866 ,101.918}, //! 5-10%  
			       {1.00196 ,-3.24098 ,71.7601}, //! 10-30%  
			       {1.00297 ,-3.29939 ,80.2865}, //! 30-50%  
			       {1.00281 ,-2.77505 ,58.3106}, //! 50-70%  
			       {0.99972 ,-1.06665 ,8.24232}, //! 70-90%   
			       {1.00000 , 0.00000 ,0.0000} //! pp
			     }, //! ak3PF	  
                             {
			       {1.00904 ,-1.67615 ,262.56},  //! 0-5%
			       {1.00588 ,-1.20147 ,179.757}, //! 5-10%
			       {1.00333 ,-1.98672 ,144.324}, //! 10-30%
			       {1.00193 ,-2.28906 ,92.6395}, //! 30-50%
			       {1.0025  ,-2.93797 ,91.361},  //! 50-70%
			       {0.99993 ,-1.53522 ,14.3988}, //! 70-90%
			       {1,0,0}  //! pp
			     }, //! ak4PF
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak2Calo
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}, //! ak3Calo
                             {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}  //! ak4Calo
};

//! after burner for mean diff                       
double amdiff[KNJ][NCEN][33]={ {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
				{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}}, //! icpu5calo
                               {
				 {0.0239412,0.01243,0.0143796,0.000505626,-0.00131786,0.0141358,-0.00442976,0.00213772,-0.00114882,0.000402868,0.0018369,0.00950521,0.00381398,-0.00156349,-0.000285923,0.00643075,0.00118148,0.000638068,0.0106224,-0.0024339,-0.00566161,0.0171672,0.002267,6.9201e-05,-0.00178653,-0.00182015,0.00878799,0.00130862,0.00111222,-0.00181752,0.00138253,0.00340289,0.00552464}, //! 0-5%
				 {0.0180368,0.00527436,0.00911248,-0.00301123,0.0152886,0.0123391,-0.0171759,-0.00536734,0.00135446,-0.00401354,0.0036822,0.00432181,0.00690466,0.00711948,0.00195038,-0.000469983,0.0110755,-0.0067665,0.00263351,-0.00372314,-0.00377512,0.0130012,0.00827807,-0.00378406,-0.00551111,0.00506836,0.0125707,-0.00443441,0.00259203,0.00498235,0.00615728,-0.00267005,0.0045433}, //! 5-10%
				 {-0.00298578,-0.00197583,0.0117723,0.00480539,0.00675088,0.00555003,-0.00230259,-0.00404507,-0.00529146,-0.00280112,-0.00271034,0.00290549,0.00630212,0.0109406,-0.000270367,0.0075264,0.0113877,0.00128502,0.00512886,-0.00858867,-0.0066945,0.0159115,-0.000799716,0.00842524,-0.0100992,0.00487304,0.00603622,-0.00134212,0.00139964,-0.00596696,-0.00298631,0.00672066,0.00374007}, //! 10-30%
				 {0.0159886,0.0117351,0.0111744,0.0079124,0.0018732,0.00147218,0.00596976,-0.000104249,-0.00558257,-0.00185359,0.00210446,0.00421846,0.00360948,0.00897342,0.000481665,0.0020231,0.0128304,0.000904322,0.00705063,-0.00379092,-0.00811958,0.010891,4.31538e-05,0.00839055,-0.00302792,0.00858408,0.0102143,0.000877321,-0.00516325,-0.00359613,-0.00322837,0.00372231,0.00214446}, //! 30-50%
				 {0.00907046,0.00839227,0.00966918,0.00977111,0.00452095,-0.0026018,0.0012759,0.00421834,-0.00492626,-0.000146508,0.00156122,0.00139749,0.00457841,0.0119306,-0.000428915,0.00253266,0.0116203,-0.00200152,0.00180191,0.00216228,-0.00886106,0.0135353,0.00320607,0.00971735,-0.00658941,0.00722986,0.00130761,0.00301373,-0.00343943,-0.00172973,0.00235015,0.00183892,0.00305146}, //! 50-70%
				 {0.0066371,0.00577068,0.0124451,0.00754064,0.00320292,0.00212824,0.00348872,-0.00321072,-0.0115995,-0.00157541,-0.00335181,0.011375,0.00554776,0.0119193,0.00823766,0.00602472,0.00506467,-0.00307351,0.00216413,0.00133753,-0.00281668,0.0130697,0.00464362,0.00522619,0.00135881,-0.000138402,0.00343561,0.00129431,0.00398809,-0.00539625,7.24792e-05,-6.98566e-05,-0.0019787}, //! 70-90%
				 {0,0,0,0,0,0,0,0,0,0}  //! pp
			       }, //! ak2PF
                               {
				 {-0.0957503,-0.0366432,-0.00239193,-0.00111985,0.00406867,0.0129364,-0.000904024,0.00290388,0.000446916,0.014991,0.00441539,0.0102446,7.05719e-05,-0.0119217,-0.00174958,-0.000980139,0.00811815,0.00433797,0.00240582,-0.00321567,-0.00202787,0.00674027,0.000503898,3.97563e-05,0.00317174,0.000580072,0.00688094,0.00317574,-0.000688612,-0.00319052,0.0101621,0.00435334,0.00570643},//! 0-5%
				 {-0.105525,-0.0316963,-0.00745112,-0.00189686,0.00352544,0.00309801,-0.00600392,-0.0121201,0.00319189,0.0059008,0.000226438,0.00754851,0.00883555,-0.00923926,0.0188677,-0.00176173,0.00513905,0.00692612,0.00295156,-0.00486153,-0.00638604,-0.0018059,0.00472921,-0.00532085,-0.000431478,0.00263381,0.00789589,0.000686705,0.00072521,0.00738513,0.00688565,-0.00357831,0.00438887},//! 5-10%
				 {-0.0516163,-0.00960845,-0.00638938,0.00395107,0.00338864,0.00748146,0.00172335,0.000775218,-0.0170699,-0.00944746,-0.00334263,-0.000624597,0.00544232,-0.000833631,0.0170863,0.00166899,0.0073728,0.0096156,0.00508732,-0.0119,-0.00246751,0.0054419,0.00258797,-0.00117266,-0.000197947,-0.00679743,0.00264871,0.00251925,-0.000521064,-0.00134331,0.0041877,-0.000174522,0.004668},//! 10-30%
				 {-0.0296276,-0.0102345,-0.00182575,0.00346065,-0.000382781,0.00336576,0.00171405,0.000446558,-0.00463557,-0.00347996,-0.000578165,0.000472009,0.00233442,-0.00186092,0.01033,0.00327784,0.00852287,0.00558299,0.00221169,-0.00364447,-0.00356013,0.000103354,0.00445223,0.000342906,-0.00261116,0.00263554,0.0114818,0.00180966,-0.00305378,-0.0103016,0.00147146,0.00506431,0.00155097},//! 30-50%
				 {-0.031548,-0.00918972,-0.00248939,-0.000692368,0.000444174,0.00123191,0.0023036,0.000940561,-0.00794965,-0.00277925,-0.00262636,0.000455439,0.00198525,-0.00150663,0.0165374,0.000724912,0.00977331,0.00610918,-0.00144237,-0.00246561,-0.00128436,0.00214124,0.00414962,0.0038625,1.8537e-05,-0.000538766,0.00695884,0.00556296,-0.00397176,-0.00839287,0.00705308,0.000589132,0.00280684},//! 50-70%
				 {-0.018268,-0.00804049,0.000407219,0.00383049,0.00274009,-0.00042206,-0.000943899,-0.00426054,-0.00859624,-0.00543714,-0.00393707,0.000441194,0.00764209,-0.00200003,0.0171912,0.00492692,0.00745159,0.00339669,0.00219131,-0.0038476,0.00133491,0.00189805,0.00218612,0.00452137,0.00087893,-0.00327235,0.00359207,0.00202942,0.00302291,-0.00834006,0.00341874,-0.00160706,-0.00174284},//! 70-90%
				 {0,0,0,0,0,0,0,0,0,0}}, //!pp //! ak3PF
                               {
				 {-0.0127124,-0.00598019,0.0171051,0.00119317,-0.00346321,0.00694638,-0.00738776,0.00591868,-0.00169849,0.00939095,0.003739,0.0082643,-0.00253433,-0.0124138,0.000649095,-0.000511825,0.00873768,-1.03712e-05,0.00254047,-0.00194591,-0.00522834,0.00782293,0.000170052,5.25117e-05,0.0036006,-0.00140399,0.00635427,0.00329226,-0.000954032,-0.00322616,0.0107534,0.00533843,0.00623101},//! 0-5%
				 {-0.0301312,-0.00641578,0.0107445,0.00492507,0.00387597,0.00327742,-0.0112711,-0.0110382,0.00213593,0.00252652,-0.00120747,0.00749069,0.0046578,-0.0117459,0.0176239,-0.00368863,0.0012641,0.00351328,0.00133175,-0.00312757,-0.00560999,-0.00240052,0.00250256,-0.00371259,0.00106132,-0.00011462,0.00878227,0.000581324,-0.000164211,0.00783348,0.00587869,-0.00223535,0.00570685},//! 5-10%
				 {-0.0199876,-0.00111932,0.0117522,0.00957417,0.0051834,0.0051254,-0.0019294,0.00262815,-0.0177857,-0.00988227,-0.00379318,-0.00124007,0.00637007,0.000297964,0.0170801,-0.000648499,0.00726557,0.00887084,0.00441515,-0.0105748,-0.0074079,0.00457609,0.00218737,-0.00267065,0.000511229,-0.00580484,0.00424474,0.00351477,-0.00273573,-0.0018695,0.00455701,0.000540853,0.0050236},//! 10-30%
				 {-0.037591,-0.0087679,0.0113178,0.00869817,0.000976503,0.00380665,-0.000551999,8.15988e-05,-0.00668401,-0.00397885,-0.00093776,0.00108021,0.00158101,-0.00303012,0.0101748,0.00378764,0.0084042,0.00755316,0.00213617,-0.00399554,-0.0038929,0.000941753,0.00513417,0.001037,-0.00230628,0.00285375,0.0116087,0.00216866,-0.00278056,-0.0100698,0.00179732,0.00547415,0.00184155},//! 30-50%
				 {-0.0346615,-0.00257123,0.0121688,0.00340623,0.00191236,0.000513136,-1.82986e-05,0.00200391,-0.00782925,-0.00323248,-0.00229836,0.0013119,0.00223309,-0.00133055,0.0167521,0.000724912,0.00977373,0.00609052,-0.00153118,-0.0025562,-0.00151426,0.00198376,0.00398856,0.00363183,-0.000365019,-0.000891984,0.00612265,0.00515777,-0.00431561,-0.00892586,0.00650108,0.000134349,0.00234717}, //! 50-70%
				 {-0.0302156,-0.00475556,0.0109398,0.00750124,0.00389731,-1.66893e-05,-0.00231707,-0.00292313,-0.00947666,-0.00548315,-0.00363231,0.000988305,0.00792807,-0.00175643,0.0174814,0.00530344,0.0075919,0.00369561,0.00294036,-0.00367105,0.00138885,0.00202137,0.00243568,0.00452137,0.00087893,-0.00327235,0.00359207,0.00202942,0.00302291,-0.00834006,0.00341874,-0.00160706,-0.00183988},//! 70-90%
				 {0,0,0,0,0,0,0,0,0,0}  //! pp
			       }, //! ak4PF

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
      for(int im=0;im<3;im++){
	fresol[nj][i]->SetParameter(im,resol[nj][i][im]);
	fscale[nj][i]->SetParameter(im,scale[nj][i][im]);
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
  float rpbpb = fresol[nj][icen]->Eval(refpt);   //! PbPb
  
  //! Calculate the smearing factor
  float smf = 0;
  if(rpp>rpbpb){
    smf=0;
  }else{
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*refpt;
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
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*refpt;
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
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*refpt;
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

  if(recopt<60)return recopt;


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
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*recopt;
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

  if(recopt<60)return recopt;

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
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*recopt;
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

  if(recopt<60)return recopt;

  float smpt = recopt;
  //! Get resolutions
  float rpp   = fresol[nj][NCEN-1]->Eval(recopt); //! pp
  float rpbpb = fresol[nj][icen]->Eval(recopt);     //! PbPb
  
  //! Calculate the smearing factor
  float smf = 0;
  if(rpp>rpbpb){
    smf=0;
  }else{
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*recopt;
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
  if(recopt<60)return recopt;

  int icen = GetCBin(hibin);
  //! Get the jet energy scale
  float mpbpb = fscale[nj][icen]->Eval(recopt);

  //! Now correct scale first
  float corr_recopt = recopt/mpbpb;
  return corr_recopt;
}
float GetReWeight(int nj,int ic,float smpt)
{
  //! Only to be used for pp data 

  float rewe=1;
  if(nj!=2 || smpt<60)return rewe;

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
  if(nj!=2 || smpt<60)return rewe;

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
    smf = sqrt(pow(rpbpb,2) - pow(rpp,2))*recopt;
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
