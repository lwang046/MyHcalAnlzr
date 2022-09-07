#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH3.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <assert.h> 

using namespace std;

int main(int argc, char *argv[])
{
  if(argc!=3){
    cerr << "Invalid arguments provided, correct format is: ./main finame foname\n";
    exit(0);
  }

  string finame = argv[1];
  string foname = argv[2];

  vector<float> lumi = {9.6497, 9.6496, 9.6495, 9.6494, 9.6493, 9.6492, 9.6491, 9.649, 9.165, 8.585, 7.950, 6.7885, 6.788, 6.434, 6.333, 5.950, 5.442, 5.245, 4.837, 4.408, 3.643, 3.312, 2.817, 2.786, 2.077, 1.820, 1.487, 1.104, 0.484, 0.178, 0.110, 0.034, 0.001};
  vector<string> runid = {"358277", "358222", "358179", "358160", "358101", "358087", "357996", "357968", "357845", "357787", "357743", "357646", "357622", "357564", "357501", "357456", "357415", "357337", "357287", "357142", "357008", "356958", "356926" , "356836", "356646", "356590", "356538", "356457", "356115", "355882", "355710", "355538", "355079"};

  //vector<float> lumi = {8.24, 8.25, 8.26, 8.27, 8.28, 8.29, 8.30, 8.31};
  //vector<string> runid = {"358030", "358092", "358152", "358155", "358164", "358209", "358253", "358283"};

  reverse(lumi.begin(), lumi.end());
  reverse(runid.begin(), runid.end());

  assert(lumi.size()==runid.size());

  int nruns = lumi.size();

  TFile *f = new TFile((finame+".root").c_str(), "read");
  TNtuple* qiedigi = (TNtuple*)f->Get("MyHcalAnlzr/qiedigi");
  int ntot = qiedigi->GetEntries();
  std::cout << "Reading in input file, total " << ntot << " Entries." << std::endl;
  float RunNum, LumiNum, EvtNum, ieta, iphi, depth, sumADC, type, shunt;
  qiedigi->SetBranchAddress("RunNum", &RunNum);
  qiedigi->SetBranchAddress("LumiNum", &LumiNum);
  qiedigi->SetBranchAddress("EvtNum", &EvtNum);
  qiedigi->SetBranchAddress("ieta", &ieta);
  qiedigi->SetBranchAddress("iphi", &iphi);
  qiedigi->SetBranchAddress("depth", &depth);
  qiedigi->SetBranchAddress("sumADC", &sumADC);
  qiedigi->SetBranchAddress("type", &type);
  qiedigi->SetBranchAddress("shunt", &shunt);

  TFile *ofile = new TFile(("hist_"+finame+"_"+foname+"_digi.root").c_str(), "recreate");

  std::cout << "Creating histograms..." << std::endl;

  TH1F ****** histarray;
  histarray = new TH1F*****[nruns];
  for(int n=0; n<nruns; n++){
    histarray[n] = new TH1F****[2];
    for(int t=0; t<2; t++){
      histarray[n][t] = new TH1F***[58];
      for(int i=0; i<58; i++){
        histarray[n][t][i] = new TH1F**[72];
        for(int j=0; j<72; j++){
          histarray[n][t][i][j] = new TH1F*[7];
          for(int k=0; k<7; k++){
            histarray[n][t][i][j][k] = new TH1F(("hist_run"+runid.at(n)+"_sipmT"+to_string(t)+"_ieta"+to_string(i<=29?i-29:i-28)+"_iphi"+to_string(j+1)+"_depth"+to_string(k+1)).c_str(), "Pedestal per Channel; ADC; Entries", 1600, 0, 40);
          }
        }
      }
    }
  }

  TH1F **** pedMean;
  TH1F **** pedRMS;
  pedMean = new TH1F***[2];
  pedRMS = new TH1F***[2];
  for(int nh=0; nh<2; nh++){
    string det = nh==0?"HB":"HE";
    pedMean[nh]=new TH1F**[2];
    pedRMS[nh]=new TH1F**[2];
    for(int nt=0; nt<2; nt++){
      pedMean[nh][nt]=new TH1F*[nruns];
      pedRMS[nh][nt]=new TH1F*[nruns];
      for(int n=0; n<nruns; n++){
        pedMean[nh][nt][n] = new TH1F((det+"_sipmT"+to_string(nt)+"pedMean_run"+runid.at(n)).c_str(), "Pedestal Mean; ADC; Entries", 1600, 0, 40);
        pedRMS[nh][nt][n] = new TH1F((det+"_sipmT"+to_string(nt)+"pedRMS_run"+runid.at(n)).c_str(), "Pedestal RMS; ADC; Entries", 1600, 0, 40);
      }
    }
  }

  std::cout << "Looping over input events..." << std::endl;

  for(int i=0; i<ntot; i++){
    if(i%50000000==0) std::cout << i << "-th entry." << std::endl;
    qiedigi->GetEntry(i);

    if(find(runid.begin(), runid.end(), to_string((int)RunNum))==runid.end()) continue;

    if(shunt!=6.0) continue;

    int runidx = distance(runid.begin(), find(runid.begin(), runid.end(), to_string((int)RunNum)) );
    int ietaidx = ieta<0?ieta+29:ieta+28;
    int iphiidx = iphi-1;
    int depthidx = depth-1;
    int sipmidx = -1;
    if(type==3||type==5) sipmidx=0;
    else if(type==4||type==6) sipmidx=1;

    //std::cout << runidx << ", " << sipmidx << ", " << ietaidx << ", " << iphiidx << ", " << depthidx << std::endl;

    histarray[runidx][sipmidx][ietaidx][iphiidx][depthidx]->Fill(sumADC/8.0);
 

  }


  std::cout << "Postprocessing histograms..." << std::endl;

  TH3F* h3 = new TH3F("h3", "h3", 59, -29, 30, 72, 0, 72, 7, 1, 8);

  for(int n=0; n<nruns; n++){
    for(int t=0; t<2; t++){
      for(int i=0; i<58; i++){
        for(int j=0; j<72; j++){
          for(int k=0; k<7; k++){
            if(histarray[n][t][i][j][k]->GetEntries()>0){
              int nh=-1;
              if((i>=14 && i<=43) || (i==13 && k<3) || (i==44 && k<3)) nh=0;
              else if((i<=12 || i>=45 ) || (i==13 && k>=3) || (i==44 && k>=3)) nh=1;
              pedMean[nh][t][n]->Fill(histarray[n][t][i][j][k]->GetMean());
              pedRMS[nh][t][n]->Fill(histarray[n][t][i][j][k]->GetRMS());
              if(t==1 && nh==0 && histarray[n][t][i][j][k]->GetMean() < 5.6) h3->Fill(i<29?i-29:i-28, j, k+1); 
            }
          }
        }
      }
    }
  }


  vector<float> MeanofMeanPedVal[2][2];
  vector<float> RMSofMeanPedVal[2][2];
  vector<float> MeanofRMSPedVal[2][2];
  for(int nh=0; nh<2; nh++){
    for(int nt=0; nt<2; nt++){
      for(int n=0; n<nruns; n++){
        MeanofMeanPedVal[nh][nt].push_back(pedMean[nh][nt][n]->GetMean());
        RMSofMeanPedVal[nh][nt].push_back(pedMean[nh][nt][n]->GetRMS());
        MeanofRMSPedVal[nh][nt].push_back(pedRMS[nh][nt][n]->GetMean());
      }
    }
  }


  TGraph *** MeanofMeanvsLumi;
  TGraph *** RMSofMeanvsLumi;
  TGraph *** MeanofRMSvsLumi;
  MeanofMeanvsLumi = new TGraph**[2];
  RMSofMeanvsLumi = new TGraph**[2];
  MeanofRMSvsLumi = new TGraph**[2];
  for(int nh=0; nh<2; nh++){
    string det = nh==0?"HB":"HE";
    MeanofMeanvsLumi[nh]=new TGraph*[2];
    RMSofMeanvsLumi[nh]=new TGraph*[2];
    MeanofRMSvsLumi[nh]=new TGraph*[2];
    for(int nt=0; nt<2; nt++){
      MeanofMeanvsLumi[nh][nt]=new TGraph(lumi.size(), &(lumi[0]), &(MeanofMeanPedVal[nh][nt][0]));
      RMSofMeanvsLumi[nh][nt]=new TGraph(lumi.size(), &(lumi[0]), &(RMSofMeanPedVal[nh][nt][0]));
      MeanofRMSvsLumi[nh][nt]=new TGraph(lumi.size(), &(lumi[0]), &(MeanofRMSPedVal[nh][nt][0]));
      MeanofMeanvsLumi[nh][nt]->SetTitle(det+"Mean_of_pedMean_sipmT"+nt);
      RMSofMeanvsLumi[nh][nt]->SetTitle(det+"RMS_of_pedMean_sipmT"+nt);
      MeanofRMSvsLumi[nh][nt]->SetTitle(det+"Mean_of_pedRMS_sipmT"+nt);
      MeanofMeanvsLumi[nh][nt]->SetName(det+"Mean_of_pedMean_sipmT"+nt);
      RMSofMeanvsLumi[nh][nt]->SetName(det+"RMS_of_pedMean_sipmT"+nt);
      MeanofRMSvsLumi[nh][nt]->SetName(det+"Mean_of_pedRMS_sipmT"+nt);
    }
  }



  std::cout << "Saving results..." << std::endl;

  ofile->cd();
  h3->Write();

  for(int nh=0; nh<2; nh++){
    for(int nt=0; nt<2; nt++){
      MeanofMeanvsLumi[nh][nt]->Write();
      RMSofMeanvsLumi[nh][nt]->Write();
      MeanofRMSvsLumi[nh][nt]->Write();
    }
  }

  for(int nh=0; nh<2; nh++){
    for(int nt=0; nt<2; nt++){
      for(int n=0; n<nruns; n++){
        pedMean[nh][nt][n]->Write();
        pedRMS[nh][nt][n]->Write();
      }
    }
  }

  ofile->Close();


  std::cout << "End Job." << std::endl;
  return 0;

}


