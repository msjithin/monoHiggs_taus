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
#define Analyzer_etau_cxx
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

using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{
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
	std::string SampleName = argv[5];
	Analyzer_etau t(argv[1],argv[2]);
	t.Loop(maxEvents,reportEvery, SampleName);

        Long_t nEvents = h_Events_level->Integral();
        Long_t Luminosity = 44980.0;
	TH1F * h1 = new TH1F(“h1”,“h1 title” , 20, 0, 4);
	h1 = (TH1F*)f.Get("ggNtuplizer/hEvents");
	nEvents = h1->GetEntries();
        ofstream myfile;
        myfile.open ("Event_number_initial.txt", ios::app);
        myfile <<" Sample name " << ":"<<std::setw(10)<< " Number of events " <<std::setw(10)<<  " luminosity " <<std::setw(10)<<  "\n";
        myfile <<" "<< sample  << ":"<<std::setw(10)<< nEvents <<std::setw(10)<< Luminosity <<std::setw(10)<<"\n";
        myfile.close();





	return 0;
}

void Analyzer_etau::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
{



}

void Analyzer_etau::BookHistos(const char* file2)
{

}

//Fill the sequential histos at a particular spot in the sequence
void Analyzer_etau::fillHistos(int histoNumber, double event_weight,int eleIndex, int tauIndex)
{
	//*********** fill electrons  ***********
// 		h_electron_En[histoNumber]->Fill((eleEn->at(eleIndex)),event_weight);
// 		h_electron_Pt[histoNumber]->Fill((elePt->at(eleIndex)),event_weight);		
// 		h_electron_eta[histoNumber]->Fill(eleEta->at(eleIndex),event_weight);
// 		h_electron_SCEta[histoNumber]->Fill(eleSCEta->at(eleIndex),event_weight);
// 		h_electron_phi[histoNumber]->Fill(elePhi->at(eleIndex),event_weight);
// 		h_electron_SCPhi[histoNumber]->Fill(eleSCPhi->at(eleIndex),event_weight);
// 	
// 	//*********** fill taus  ***********
// 
// 		h_tau_En[histoNumber]->Fill((tauEnergy->at(tauIndex)),event_weight);
// 		h_tau_Pt[histoNumber]->Fill((tauPt->at(tauIndex)),event_weight);
// 		h_tau_eta[histoNumber]->Fill(tauEta->at(tauIndex),event_weight);
// 		h_tau_phi[histoNumber]->Fill(tauPhi->at(tauIndex),event_weight);
// 		
// 	//*********** fill met  ***********		
// 		h_pfMET[histoNumber]->Fill(pfMET,event_weight);
// 		double dPhi_etau = DeltaPhi(elePhi->at(eleIndex),tauPhi->at(tauIndex));
// 		h_dPhi[histoNumber]->Fill(dPhi_etau,event_weight);
// 
// 
// 	//*********** fill rest  ***********
// 
// 		float mT_eMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
// 		h_Mt[histoNumber]->Fill(mT_eMet,event_weight);
// 		TLorentzVector myTau; 
// 		myTau.SetPtEtaPhiE(tauPt->at(tauIndex),tauEta->at(tauIndex),tauPhi->at(tauIndex), tauEnergy->at(tauIndex));		
// 		TLorentzVector myEle; 
// 		myEle.SetPtEtaPhiE(elePt->at(eleIndex),eleEta->at(eleIndex),elePhi->at(eleIndex), eleEn->at(eleIndex));
// 		double visMass_etau = VisMass_F(myTau, myEle);
// 		h_VisibleMass[histoNumber]->Fill(visMass_etau,event_weight);
// 		double HiggsPt = pTvecsum_F(elePt->at(eleIndex),tauPt->at(tauIndex),elePhi->at(eleIndex),tauPhi->at(tauIndex) );
// 		h_HiggsPt[histoNumber]->Fill(HiggsPt,event_weight);
// 		h_nVtx[histoNumber]->Fill(nVtx);
// 		h_nEvents[histoNumber]->Fill(isData);
// 		//		if(jetsize>0){
// 		//h_leadingJetPt[histoNumber]->Fill(jetPt->at(),event_weight);
// 		//h_leadingJetPt_300[histoNumber]->Fill(jetPt->at(),event_weight);
// 		//h_leadingJetEta[histoNumber]->Fill(jetEta->at(),event_weight);
// 		//}
// 	    
}




//---------------------------------------------------                                                                                                                                
// get a electron candiate based on pt eta and isolation                                                                                                                               
//----------------------------------------------------                                                                                                                               

std::vector<int> Analyzer_etau::getEleCand(double elePtCut, double eleEtaCut){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over electrons                                                                                                                                                             
  for(int iEle=0;iEle<nEle;iEle++)
    {
      bool kinematic = false;
      if( (*elePt)[iEle] > elePtCut  && fabs((*eleEta)[iEle])< eleEtaCut ) kinematic = true;
      bool electronId =false;
      if( eleIDbit->at(iEle)>>3&1==1) electronId =true;
      bool relative_iso = false;
     
      float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));

      if( relEleIso < 0.10 ) relative_iso = true;
      if(electronId && kinematic && relative_iso){
	tmpCand.push_back(iEle);
      }                                                                                                                                                              
    }                                                                                                                                                                

  return tmpCand;

}

std::vector<int> Analyzer_etau::getTauCand(double tauPtCut, double tauEtaCut){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons                                                                                                                                                             
  for(int iTau=0;iTau<nTau;iTau++)
    {
      bool kinematic = false;
      bool tauId = false;
      bool decayModeCut = true;
      bool tauIsolation = false;

      if( tauPt->at(iTau) > tauPtCut  && fabs( tauEta->at(iTau)< tauEtaCut ) )kinematic = true;
      if( tauByMVA6TightElectronRejection->at(iTau) == 1 && tauByLooseMuonRejection3->at(iTau) == 1) tauId = true;
      if( (tauByTightIsolationMVArun2v1DBoldDMwLT->at(iTau)==1)) tauIsolation = true;
      if( tauDecayMode->at(iTau)==1 || tauDecayMode->at(iTau)==3 ) decayModeCut = true;
  
      if(tauId && kinematic && tauIsolation && decayModeCut){
	tmpCand.push_back(iTau);
      }                                                                                                                                                              
    }                                                                                                                                                                
  
  return tmpCand;

}

bool Analyzer_etau::passSingleEleTrg(Long_t nEntries){

  int tmpCand =0;
  bool eleTrgPassed= false;
  for(int ientry=0; ientry<nEntries ;ientry++)
    {
      if (eleFiredSingleTrgs->at(ientry)>>5&1 ==1) tmpCand++ ;
    }                                                                                                                                                       

  if(tmpCand>0)eleTrgPassed= true;
   return eleTrgPassed;
}




//Veto failed if a electron is found that passes Loose Muon ID, Loose Muon Isolation, and muPtcut, and does not overlap the candidate photon within dR of 0.5                                                                                                                                                                    

double Analyzer_etau::dR(int ele_index, int tau_index)
{
  double deltaeta = abs(eleEta->at(ele_index) - tauEta->at(tau_index));
  double electron_Phi = elePhi->at(ele_index);
  double tau_Phi = tauPhi->at(tau_index);

  double deltaphi = DeltaPhi(electron_Phi, tau_Phi);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}



double Analyzer_etau::DeltaPhi(double phi1, double phi2)
//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
{
  double pi = TMath::Pi();
  double dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}

float Analyzer_etau::TMass_F(float LepPt, float LepPhi , float met, float metPhi) {
    return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}

float Analyzer_etau::TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met) {
  float totalTMass = (a + b+ met).M();
  return totalTMass;
}


float Analyzer_etau::VisMass_F(TLorentzVector a, TLorentzVector b){
  float visibleMass = (a + b).M();
  return visibleMass;
}

float Analyzer_etau::pTvecsum_F(float pt1, float pt2, float phi1, float phi2) {
  float pt_vecSum = sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2), 2) + pow(pt1*sin(phi1) + pt2*sin(phi2), 2));
  return pt_vecSum;
}

bool Analyzer_etau::passBjetVeto()
{
  int tmpCand = 0;
  bool veto = false;
  for(int iJets=0; iJets<nJet ; iJets++){
    if(jetCSV2BJetTags->at(iJets) > 0.8484) tmpCand++;
  }
  if(tmpCand>0) veto = true;
}
