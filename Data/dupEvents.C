#include <iostream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <iomanip>   
#include <vector>
#include <algorithm>

#include <TFile.h>
#include <TNtuple.h>

using namespace std;

bool cmp_by_run(const std::pair<float, float>& firE,const std::pair<float, float>& secE)
{
  return firE.first < secE.first;
}

int dupEvents()
{
  
  TFile *f1 = new TFile("input/ntuple/ntuple_2013_ppJet40.root","r");
  TNtuple *n1 = (TNtuple*)f1->Get("ntjet");
  TFile *f2 = new TFile("input/ntuple/ntuple_2013_ppJet80.root","r");
  TNtuple *n2 = (TNtuple*)f2->Get("ntjet");

  float r1, e1;
  std::vector < std::pair < float, float > >evts1;
  n1->SetBranchAddress("run",&r1);
  n1->SetBranchAddress("evt",&e1);

  float r2, e2;
  std::vector < std::pair < float, float > >evts2;
  n2->SetBranchAddress("run",&r2);
  n2->SetBranchAddress("evt",&e2);



  std::cout<<"# of entries in jet40 : " <<n1->GetEntries()<<std::endl;
  for(int i=0; i<n1->GetEntries(); i++){
    n1->GetEntry(i);
    evts1.push_back(std::make_pair(r1,e1));    
  }
  std::cout<<"# of entries in jet80 : " <<n2->GetEntries()<<std::endl;
  for(int i=0; i<n2->GetEntries(); i++){
    n2->GetEntry(i);
    evts2.push_back(std::make_pair(r2,e2));    
  }
  f1->Close();
  f2->Close();

  std::vector < std::pair < float, float > > dupEvts;
  std::vector < std::pair < float, float > >::iterator it; 

  // for(it = evts1.begin(); it != evts1.end(); ++it){
  //   std::cout<< " r1 : " << (*it).first << " e1 : "<< (*it).second << std::endl;
  // }

  std::sort(evts1.begin(), evts1.end(), cmp_by_run);
  std::sort(evts2.begin(), evts2.end(), cmp_by_run);

  // for(it = evts2.begin(); it != evts2.end(); ++it){
  //   std::cout<< " r2 : " << (*it).first << " e2 : "<< (*it).second << std::endl;
  // }


  std::set_intersection(evts1.begin(), evts1.end(), evts2.begin(), evts2.end(), std::back_inserter(dupEvts));
  for(it = dupEvts.begin(); it != dupEvts.end(); ++it){
    std::cout<< " run : " << (*it).first << " event : " << std::setprecision(11) << (*it).second << std::endl;
  }

  return 0;

}
