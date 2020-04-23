////Analyzer_etau.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom Analyzer_etau analyze
//
//To run, assuming this is compiled to an executable named 'analyze':
//$ ./analyze /hdfs/store/user/jmadhusu/LatestNtuples/ /afs/hep.wisc.edu/user/ms/CMSSW_9_4_4/src/2017_analysis/etau/output.root -1 10000
//./analyze /hdfs/store/user/jmadhusu/MonoHiggs_MC2017/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/crab_ZZZ/180603_185329/0000/ /afs/hep.wisc.edu/user/ms/CMSSW_9_4_4/src/2017_analysis/analyzer/output.root -1 10000
//Runs over every event in the folder LatestNtuples, reporting progress every 10000 events
//and storing the resulting histograms in the file output.root.
//
//To plot, for example, single photon trigger efficiency as a function of photon pt:
//$ root -l
//root[0] TFile *f = new TFile("output.root");
//root[1] TGraphAsymmErrors *efficiency = new TGraphAsymmErrors((TH1F*)f->Get("Photon_Et_300_2"),(TH1F*)f->Get("Photon_Et_300_1"));
//root[2] efficiency->Draw("AP")
//root[3] efficiency->SetTitle("Single photon trigger efficiency")
//root[4] efficiency->GetXaxis()->SetTitle("Photon p_{T}")
//root[5] efficiency->Draw("AP")
//
#define PrintEvents_cxx
#include "Analyzer_etau.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <set>
#include "TMath.h" //M_PI is in TMath
#include "TRandom3.h"
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>

#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <stdio.h>
using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{

  std::string input = *(argv + 1);
  std::string output = *(argv + 2);


	Long64_t maxEvents = atof(argv[3]);
	if (maxEvents < -1LL)
	{
		std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
		return 1;
	}
	int reportEvery = atof(argv[4]);
	if (reportEvery < 1)
	{
		std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
		return 1;
	}
	std::string sample = argv[5];

	TFile *f_Double = new TFile(input.c_str());
	TH1F* nbevt = (TH1F*) f_Double->Get("isData" );
	Double_t ngen =0;
        //if (sample=="Data_2017"){  ngen = nbevt->GetBinContent(2);}
	//else {  ngen = nbevt->GetBinContent(1);}
	
	ngen = nbevt->GetBinContent(2);
	Long_t Luminosity = 44980.0;

	ofstream myfile;
	myfile.open ("Event_number_initial.txt", ios::app);
	myfile <<" Sample name " << ":"<<std::setw(10)<< " Number of events " <<std::setw(10)<<  " luminosity " <<std::setw(10)<<  "\n";
	myfile <<" "<< sample  << ":"<<std::setw(10)<< ngen <<std::setw(10)<< Luminosity <<std::setw(10)<<"\n";
	myfile.close();
       
}


  
 
