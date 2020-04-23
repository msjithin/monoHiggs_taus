#define Analyzer_tautau_cxx
#include "Analyzer_tautau.h"
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


int main(int argc, const char* argv[]){
   string SampleName = argv[5];
   Long64_t maxEvents = atof(argv[3]);
   int reportEvery = atof(argv[4]);
   
   Analyzer_tautau t(argv[1], argv[2]);
   t.Loop(maxEvents, reportEvery, SampleName);
   return 0;
}


void Analyzer_tautau::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
{
   if (fChain == 0) return;

   // Initialize the vectors for the two tau candidates
   vector<int> tau1Cand;
   vector<int> tau2Cand;
   tau1Cand.clear();
   tau2Cand.clear();

   // Load up the file
   TString sample = TString(SampleName);

   // Make some histograms right quick
   TH1F* h_Events_level = new TH1F("Events_level","Events at each selection level", 15,0,30);
   TH1F* h_Cutflow= new TH1F("Cutflow","Events at each level of selection",6,0,12);

   // Making some variables that are used for the cutflow
   double numberOfEvents, nMETFiltersPassed, nSingleTrgPassed, nGoodTau1Passed, nGoodTau2Passed, nGoodPairPassed, nDeltaRPassed, nPassedThirdLeptonVeto, nPassedBJetVeto, nPassedMT, nPassedHiggsPtcut, nPassVisibleMasscut, nPassedMETcut;
   numberOfEvents= nMETFiltersPassed= nSingleTrgPassed= nGoodTau1Passed= nGoodTau2Passed= nGoodPairPassed= nDeltaRPassed= nPassedThirdLeptonVeto= nPassedBJetVeto= nPassedMT= nPassedHiggsPtcut= nPassVisibleMasscut= nPassedMETcut=0;


   // Here please include the histograms that are needed for weighting things
   auto tauSFs = TauTriggerSFs2017("tauTriggerEfficiencies2017.root","ditau","2017", "tight","MVAv2");

   TFile* prefiring_file = TFile::Open("L1PrefiringMaps_new.root");
   TH2F* prefiring_photon = (TH2F*)prefiring_file->Get("L1prefiring_photonptvseta_2017BtoF");
   TH2F* prefiring_jet = (TH2F*)prefiring_file->Get("L1prefiring_jetptvseta_2017BtoF");


   Long64_t nentries = fChain->GetEntries();
   Long64_t nentriesToCheck = nentries;
   if (maxEvents != -1LL && nentries > maxEvents) nentriesToCheck = maxEvents;
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++) {

      event_.clear();
      event_info.clear();

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // if (jentry%4!=0) continue;
 
      // start the event grooming and selection process here

      //stuff to avoid negative weights
      double inspected_event_weight=1.0;
      if (isData==0) fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight/fabs(genWeight) : inspected_event_weight = 0.0;      
      
      // this is where the weighting stuff gets going, declare everything one that will be modified later
      double event_weight=1.0;

      // make the tau1 and tau2 selection that will be used later
      int skip_me = -1;
      tau1Cand = getTau1Cand(55, 2.1);
      if (tau1Cand.size()>0) skip_me=tau1Cand[0];
      tau2Cand = getTau2Cand(45, 2.1, skip_me);


      // make an entry for the number of total events
      numberOfEvents+=event_weight;

      //start selection process
      
      //met filter
      if (metFilters==0){
         //fix the weight
         if (isData==0) fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;
         nMETFiltersPassed+=event_weight;

         // trigger
         if( (HLTEleMuX>>29&1 == 1) or (HLTEleMuX>>30&1 == 1) or (HLTEleMuX>>31&1 == 1) ){
            nSingleTrgPassed+=event_weight;
            
            // first tau
            if (tau1Cand.size()>0){
               if (isData==0) event_weight = event_weight*0.89; //Single tau ID
               nGoodTau1Passed+=event_weight;

               //second tau
               if (tau2Cand.size()>0){
                  if (isData==0) event_weight = event_weight*0.89; // single tau ID
                  nGoodTau2Passed+=event_weight;


                  // do the weighting and the scale factor stuff here
                  if (isData==0) event_weight = event_weight*tauSFs.getTriggerScaleFactor(tauPt->at(tau1Cand[0]), tauEta->at(tau1Cand[0]), tauPhi->at(tau1Cand[0]), tauDecayMode->at(tau1Cand[0]) );
                  if (isData==0) event_weight = event_weight*tauSFs.getTriggerScaleFactor(tauPt->at(tau2Cand[0]), tauEta->at(tau2Cand[0]), tauPhi->at(tau1Cand[0]), tauDecayMode->at(tau2Cand[0]));

                  if (isData==0) event_weight = event_weight * prefiring_weight(prefiring_photon, prefiring_jet);

                  if (isData==1 && event_weight != 1.0) cout<<"$$$$$$ WARNING DATA BEING WEIGHTED!"<<endl;

                  if (event_weight<0.0) cout << "WARNING WARNING TOU FUCKED SOMETHING UP WARNING WANRING WARNING"<<endl;

                  // opposite sign selection
                  if (tauCharge->at(tau1Cand[0])*tauCharge->at(tau2Cand[0]) < 0){
                     nGoodPairPassed+=event_weight;
                     fillHistos(0,event_weight,tau1Cand[0],tau2Cand[0]);

                     // dR selection
                     if ((dR(tau1Cand[0], tau2Cand[0])>0.3)){
                        nDeltaRPassed+=event_weight;
                        fillHistos(1,event_weight,tau1Cand[0],tau2Cand[0]);

                        // third lepton veto
                        if (thirdLeptonVeto()==false){
                           nPassedThirdLeptonVeto+=event_weight;
                           fillHistos(2,event_weight,tau1Cand[0],tau2Cand[0]);

                           // bjet veto
                           if (passBjetVeto() == true){
                              nPassedBJetVeto+=event_weight;
                              fillHistos(3,event_weight,tau1Cand[0],tau2Cand[0]);

                              //tranverse mass cut
                              if ((TMass_F(tauPt->at(tau1Cand[0]),tauPhi->at(tau1Cand[0]),pfMET,pfMETPhi))>50){
                                 nPassedMT+=event_weight;
                                 fillHistos(4,event_weight,tau1Cand[0],tau2Cand[0]);

                                 //higgs pt cut
                                 if ( (pTvecsum_F(tauPt->at(tau1Cand[0]), tauPt->at(tau2Cand[0]), tauPhi->at(tau1Cand[0]), tauPhi->at(tau2Cand[0])))>65){
                                    nPassedHiggsPtcut+=event_weight;
                                    fillHistos(5,event_weight,tau1Cand[0],tau2Cand[0]);

                                    // visible mass cut
                                    TLorentzVector myTau1;
                                    TLorentzVector myTau2;
                                    TLorentzVector myMET;
                                    myTau1.SetPtEtaPhiE(tauPt->at(tau1Cand[0]), tauEta->at(tau1Cand[0]), tauPhi->at(tau1Cand[0]), tauEnergy->at(tau1Cand[0]));
                                    myTau2.SetPtEtaPhiE(tauPt->at(tau2Cand[0]), tauEta->at(tau2Cand[0]), tauPhi->at(tau2Cand[0]), tauEnergy->at(tau2Cand[0]));
                                    myMET.SetPtEtaPhiE(pfMET, 0, pfMETPhi, pfMET);
                                    if ( VisMass_F(myTau1,myTau2) < 125){
                                       nPassVisibleMasscut+=event_weight;
                                       fillHistos(6,event_weight,tau1Cand[0],tau2Cand[0]);

                                       // MET cut
                                       if ((pfMET > 105)){
                                          nPassedMETcut+=event_weight;
                                          fillHistos(7,event_weight,tau1Cand[0],tau2Cand[0]);

                                          if ((TotTMass_F(myTau1, myTau2, myMET) > 260)){
                                             fillHistos(8,event_weight,tau1Cand[0],tau2Cand[0]);
              
   
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }


      h_Events_level->SetBinContent(1, numberOfEvents);
      h_Events_level->SetBinContent(2, nMETFiltersPassed);
      h_Events_level->SetBinContent(3, nSingleTrgPassed);
      h_Events_level->SetBinContent(4, nGoodTau1Passed);
      h_Events_level->SetBinContent(5, nGoodPairPassed);
      h_Events_level->SetBinContent(6, nDeltaRPassed);
      h_Events_level->SetBinContent(7, nPassedThirdLeptonVeto);
      h_Events_level->SetBinContent(8, nPassedBJetVeto);
      h_Events_level->SetBinContent(9, nPassedMT);
      h_Events_level->SetBinContent(10, nPassedHiggsPtcut);
      h_Events_level->SetBinContent(11, nPassVisibleMasscut);
      h_Events_level->SetBinContent(12, nPassedMETcut);

      h_Cutflow->SetBinContent(1, nGoodTau1Passed);
      h_Cutflow->SetBinContent(2, nGoodTau2Passed);
      h_Cutflow->SetBinContent(3, nGoodPairPassed);
      h_Cutflow->SetBinContent(4, nDeltaRPassed);
      h_Cutflow->SetBinContent(5, nPassedThirdLeptonVeto);
      h_Cutflow->SetBinContent(6, nPassedBJetVeto);

      tree->Fill();

      if (jentry%reportEvery==0){
         std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
      }
   }


   cout<<"initial value: "<<h_Events_level->GetBinContent(1)<<endl;
   cout<<"met filter: "<<h_Events_level->GetBinContent(2)<<endl;
   cout<<"single trigger: "<<h_Events_level->GetBinContent(3)<<endl;
   cout<<"first tau: "<<h_Events_level->GetBinContent(4)<<endl;
   cout<<"second tau: "<<h_Events_level->GetBinContent(5)<<endl;
   cout<<"good pair: "<<h_Events_level->GetBinContent(6)<<endl;
   cout<<"delta R: "<<h_Events_level->GetBinContent(7)<<endl;
   cout<<"Third lepton veto: "<<h_Events_level->GetBinContent(8)<<endl;
   cout<<"Bjet veto: "<<h_Events_level->GetBinContent(9)<<endl;
   cout<<"passed MT: "<<h_Events_level->GetBinContent(10)<<endl;
   cout<<"passed higgs pt: "<<h_Events_level->GetBinContent(11)<<endl;
   cout<<"passed visible mass: "<<h_Events_level->GetBinContent(12)<<endl;
   cout<<"passed MET: "<<h_Events_level->GetBinContent(13)<<endl;



}



void Analyzer_tautau::BookHistos(const char* file2){

   fileName = new TFile(file2, "RECREATE");
   tree = new TTree("ADD","ADD");
   tree->Branch("event_","std::vector<unsigned int>",&event_);
   tree->Branch("event_info","std::vector<double>",&event_info);
   fileName->cd();


   Float_t PtBins[21]={0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,110, 120,130, 140,150, 160,170, 180,190, 200.};
   Float_t MetBins[14]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180, 200, 300, 400, 500};
   Float_t TrMassBins[17]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180, 200, 225, 250,275, 300,400,500 };


   //Now make a loop that gets the histograms ready for action
   for (int i=0;i<10;i++){
      
      char ptbins[100];
      sprintf(ptbins, "_%d", i);
      std::string histname(ptbins);
      h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(), "nVtx",40,0,40);h_nVtx[i]->Sumw2();
      h_nEvents[i] = new TH1F(("nEvents"+histname).c_str(), "nEvents",3,0,3);h_nEvents[i]->Sumw2();
      h_genWeight[i] = new TH1F(("genWeight"+histname).c_str(), "genWeight",10, -10.0, 10.0);h_genWeight[i]->Sumw2();
      h_genHT[i] = new TH1F(("genHT"+histname).c_str(), "genHT",20,PtBins);h_genHT[i]->Sumw2();
      
      // first lepton
      h_tau1_En[i] = new TH1F(("Tau1_En"+histname).c_str(), "Tau1_En",20,PtBins);h_tau1_En[i]->Sumw2();
      h_tau1_Pt[i] = new TH1F(("Tau1_Pt"+histname).c_str(), "Tau1_Pt",20,PtBins);h_tau1_Pt[i]->Sumw2();
      h_tau1_eta[i] = new TH1F(("Tau1_eta"+histname).c_str(), "Tau1_eta",20,-3.0, 3.0);h_tau1_eta[i]->Sumw2();
      h_tau1_phi[i] = new TH1F(("Tau1_phi"+histname).c_str(), "Tau1_phi", 21,-3.14,3.14);h_tau1_phi[i]->Sumw2();


      // second lepton
      h_tau2_En[i] = new TH1F(("Tau2_En"+histname).c_str(), "Tau2_En",20,PtBins);h_tau2_En[i]->Sumw2();
      h_tau2_Pt[i] = new TH1F(("Tau2_Pt"+histname).c_str(), "Tau2_Pt",20,PtBins);h_tau2_Pt[i]->Sumw2();
      h_tau2_eta[i] = new TH1F(("Tau2_eta"+histname).c_str(), "Tau2_eta",20,-3.0, 3.0);h_tau2_eta[i]->Sumw2();
      h_tau2_phi[i] = new TH1F(("Tau2_phi"+histname).c_str(), "Tau2_phi", 21,-3.14,3.14);h_tau2_phi[i]->Sumw2();



      // MET stuff
      h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "pfMET",13,MetBins);h_pfMET[i]->Sumw2();
      h_pfMET_phi[i] = new TH1F(("pfMET_phi"+histname).c_str(), "pfMET_phi",21,-3.14,3.14);h_pfMET_phi[i]->Sumw2();
      h_dPhi[i] = new TH1F(("h_dPhi"+histname).c_str(),"h_dPhi",20,0,3.15);h_dPhi[i]->Sumw2();
      h_dR[i] = new TH1F(("h_dR"+histname).c_str(),"h_dR",20,0,3.15);h_dR[i]->Sumw2();
      h_nJet[i] = new TH1F(("nJet"+histname).c_str(), "nJet",20,0,20);h_nJet[i]->Sumw2();
      h_leadingJetPt[i] = new TH1F(("leadingJetPt"+histname).c_str(),"leadingJetPt",30,20,1000);h_leadingJetPt[i]->Sumw2();
      h_leadingJetEta[i] = new TH1F(("h_leadingJetEta"+histname).c_str(),"h_leadingJetEta",40,-1.4442,1.4442);h_leadingJetEta[i]->Sumw2();



      // Other misc variables
      h_Mt[i]= new TH1F(("Mt"+histname).c_str(),"MT",15,TrMassBins);h_Mt[i]->Sumw2();
      h_totTrMass[i]= new TH1F(("TotalTrMass"+histname).c_str(),"TotalTrMass",15,TrMassBins);h_totTrMass[i]->Sumw2();
      h_VisibleMass[i]= new TH1F(("VisibleMass"+histname).c_str(),"VisibleMass",40,0,200);h_VisibleMass[i]->Sumw2();        
      h_HiggsPt[i]= new TH1F(("HiggsPt"+histname).c_str(),"HiggsPt",20,PtBins);h_HiggsPt[i]->Sumw2();
      h_tauIso[i]= new TH1F(("Tau_iso"+histname).c_str(),"Tau_iso", 12, 0.0, 1.2);h_tauIso[i]->Sumw2();
      h_tauMass[i]= new TH1F(("Tau_mass"+histname).c_str(),"Tau_mass", 18, 0.0, 1.8);h_tauMass[i]->Sumw2();
      h_tauDecayMode[i]= new TH1F(("Tau_Decay_Mode"+histname).c_str(),"Tau_Decay_Mode", 11, 0.0, 11.0);h_tauDecayMode[i]->Sumw2();
  

   }
}





void Analyzer_tautau::fillHistos(int histoNumber, double event_weight, int tau1Index, int tau2Index){

   // fill the first lepton
   h_tau1_En[histoNumber]->Fill(tauEnergy->at(tau1Index), event_weight);
   h_tau1_Pt[histoNumber]->Fill(tauPt->at(tau1Index), event_weight);
   h_tau1_eta[histoNumber]->Fill(tauEta->at(tau1Index), event_weight);
   h_tau1_phi[histoNumber]->Fill(tauPhi->at(tau1Index), event_weight);

   // fill the second lepton
   h_tau2_En[histoNumber]->Fill(tauEnergy->at(tau2Index), event_weight);
   h_tau2_Pt[histoNumber]->Fill(tauPt->at(tau2Index), event_weight);
   h_tau2_eta[histoNumber]->Fill(tauEta->at(tau2Index), event_weight);
   h_tau2_phi[histoNumber]->Fill(tauPhi->at(tau2Index), event_weight);

   // fill the MET thi
   h_pfMET[histoNumber]->Fill(pfMET,event_weight);
   h_pfMET_phi[histoNumber]->Fill(pfMETPhi,event_weight);
   double dPhi_tautau = DeltaPhi(tauPhi->at(tau1Index),tauPhi->at(tau2Index));
   h_dPhi[histoNumber]->Fill(dPhi_tautau,event_weight);
   double dR_tautau = dR(tau1Index,tau2Index);
   h_dR[histoNumber]->Fill(dR_tautau,event_weight);


   // misc other things
   float mT_eMet = TMass_F((tauPt->at(tau1Index)),(tauPhi->at(tau1Index)),pfMET,pfMETPhi  );
   h_Mt[histoNumber]->Fill(mT_eMet,event_weight);
   TLorentzVector myTau1;
   myTau1.SetPtEtaPhiE(tauPt->at(tau1Index),tauEta->at(tau1Index),tauPhi->at(tau1Index), tauEnergy->at(tau1Index));
   TLorentzVector myTau2;
   myTau2.SetPtEtaPhiE(tauPt->at(tau2Index),tauEta->at(tau2Index),tauPhi->at(tau2Index), tauEnergy->at(tau2Index));
   double visMass_tautau = VisMass_F(myTau1, myTau2);
   h_VisibleMass[histoNumber]->Fill(visMass_tautau,event_weight);
   TLorentzVector myMet;
   myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
   float tot_tr_mass = TotTMass_F(myTau1, myTau2, myMet );
   h_totTrMass[histoNumber]->Fill(tot_tr_mass,event_weight);
   double HiggsPt = pTvecsum_F(tauPt->at(tau1Index),tauPt->at(tau2Index),tauPhi->at(tau1Index),tauPhi->at(tau2Index) );
   h_HiggsPt[histoNumber]->Fill(HiggsPt,event_weight);
   h_nVtx[histoNumber]->Fill(nVtx,event_weight);
   h_nEvents[histoNumber]->Fill(isData,event_weight);
   h_nJet[histoNumber]->Fill(nJet ,event_weight);

   h_genHT[histoNumber]->Fill(genHT,event_weight);
   h_genWeight[histoNumber]->Fill(genWeight,event_weight);

}



std::vector<int> Analyzer_tautau::getTau1Cand(double tauPtCut, double tauEtaCut){

   std::vector<int> tmpCand;
   tmpCand.clear();
   
   for (int iTau=0; iTau<nTau; iTau++){
      bool kinematic = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;

      if( tauPt->at(iTau) > tauPtCut  && fabs( tauEta->at(iTau))< tauEtaCut && taudz->at(iTau)<0.2 )kinematic = true;
      if( tauByLooseMuonRejection3->at(iTau) == 1 && tauByMVA6VLooseElectronRejection->at(iTau) == 1) tauId = true;
      if( (taubyMediumIsolationMVArun2017v2DBoldDMwLT2017->at(iTau)==1)) tauIsolation = true;
      if( tauDecayMode->at(iTau)==0 || tauDecayMode->at(iTau)==1 || tauDecayMode->at(iTau)==10 ) decayModeCut = true;

      if(tauId==true && kinematic==true && tauIsolation==true && decayModeCut==true) tmpCand.push_back(iTau);
   }
   return tmpCand;
}


std::vector<int> Analyzer_tautau::getTau2Cand(double tauPtCut, double tauEtaCut, int skip_me){

   std::vector<int> tmpCand;
   tmpCand.clear();

   for (int iTau=0; iTau<nTau; iTau++){
      if (iTau==skip_me) continue;

      bool kinematic = false;
      bool tauId = false;
      bool decayModeCut = false;
      bool tauIsolation = false;

      if( tauPt->at(iTau) > tauPtCut  && fabs( tauEta->at(iTau))< tauEtaCut )kinematic = true;
      if( tauByLooseMuonRejection3->at(iTau) == 1 && tauByMVA6VLooseElectronRejection->at(iTau) == 1) tauId = true;
      if( (taubyTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau)==1)) tauIsolation = true;
      if( tauDecayMode->at(iTau)==0 || tauDecayMode->at(iTau)==1 || tauDecayMode->at(iTau)==10 ) decayModeCut = true;

      if(tauId==true && kinematic==true && tauIsolation==true && decayModeCut==true) tmpCand.push_back(iTau);
   }
   return tmpCand;
}



bool Analyzer_tautau::thirdLeptonVeto(){

   int tmpCandEle = 0;
   int tmpCandMu = 0;
   bool thirdLepVeto=false;

   //Loop over electrons
   for(int iEle=0; iEle < nEle;iEle++){

      bool kinematicEle, electronId, relativeIsoEle;
      kinematicEle = electronId = relativeIsoEle = false;

      if( (*elePt)[iEle] > 10.0  && fabs((*eleEta)[iEle])< 2.5 && (*eleD0)[iEle] < 0.045 && (*eleDz)[iEle] < 0.2 ) kinematicEle = true;
      if( eleIDbit->at(iEle)>>0&1==1) electronId =true;

      float relEleIso = ( elePFChIso->at(iEle) + max( elePFNeuIso->at(iEle) + elePFPhoIso->at(iEle) - 0.5 *elePFPUIso->at(iEle) , 0.0 )) / (elePt->at(iEle));
      if( relEleIso < 0.3 ) relativeIsoEle = true;

      if(electronId==true && kinematicEle==true && relativeIsoEle==true) tmpCandEle++;
   }

   // l(op of the muons
   for(int iMu=0; iMu<nMu; iMu++){

      bool kinematicMu, muonId, relativeIsoMu;
      kinematicMu = muonId = relativeIsoMu = false;

      if( (*muPt)[iMu] > 10.0  && fabs((*muEta)[iMu])< 2.4 && (*muD0)[iMu] < 0.045 && (*muDz)[iMu] < 0.2 ) kinematicMu = true;
      if( muIDbit->at(iMu)>>0&1==1) muonId =true;

      float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
      if( relMuIso < 0.3 ) relativeIsoMu = true;

      if(muonId==true && kinematicMu==true && relativeIsoMu==true) tmpCandMu++;
   }

   if (tmpCandEle>0 || tmpCandMu>0) thirdLepVeto=true;
   return thirdLepVeto;
}


bool Analyzer_tautau::passBjetVeto(){
   int tmpCand=0;
   bool veto = true;
   bool foundBjet = false;
   for (int iJets=0; iJets<nJet; iJets++){
      if(jetCSV2BJetTags->at(iJets) > 0.8838) tmpCand++;
   }
   if(tmpCand>0){ veto = false; foundBjet = true; }
   return veto;
}


bool Analyzer_tautau::check_jets(){
   int tmpCand=0;
   bool found_event=false;
   for (int iJets=0; iJets<nJet; iJets++){
      if ( (jetPt->at(iJets) > 100) && (fabs(jetEta->at(iJets) < 3.0 )) &&  (fabs(jetEta->at(iJets) >2.25 ))){
         found_event=true;
      }
   }
   return found_event;
}


double Analyzer_tautau::prefiring_weight(TH2F *prefiring_photon, TH2F *prefiring_jet){
    double weight = 1.0;
    double map_val_ele = 0.0;
    int picked_ele = -1;
    double map_val_jet = 0.0;
    int x_bin;
    int y_bin;
    for (int iEle=0; iEle<nEle; iEle++){
        if (fabs(eleEta->at(iEle)) > 2.0){
            x_bin = prefiring_photon->GetXaxis()->FindBin( eleEta->at(iEle) );
            y_bin = prefiring_photon->GetYaxis()->FindBin( elePt->at(iEle) );
            map_val_ele = prefiring_photon->GetBinContent(x_bin, y_bin);
            picked_ele = iEle;
        }
    }
    double dEta;
    double dPhi;
    double dr;
    for (int iJet=0; iJet<nJet; iJet++){
        if (fabs(jetEta->at(iJet)) > 2.0 ){
            if (picked_ele > -1){
                dEta = fabs(eleEta->at(picked_ele)) - fabs(jetEta->at(iJet));
                dPhi = DeltaPhi(elePhi->at(picked_ele), jetPhi->at(iJet));
                dr = sqrt(dEta*dEta + dPhi*dPhi);
                if (dr<0.3) continue;
            }
            x_bin = prefiring_jet->GetXaxis()->FindBin( jetEta->at(iJet) );
            y_bin = prefiring_jet->GetYaxis()->FindBin( jetPt->at(iJet) );
            map_val_jet = prefiring_jet->GetBinContent( x_bin, y_bin);
        }
    }
    weight = (1-map_val_ele)*(1-map_val_jet);
    return weight; 
}


//All the little functions that make the various variables go under here

float Analyzer_tautau::TMass_F(float LepPt, float LepPhi, float met, float metPhi){
   return sqrt(pow(LepPt + met, 2) - pow(LepPt* cos(LepPhi) + met * cos(metPhi), 2) - pow(LepPt * sin(LepPhi) + met * sin(metPhi), 2));
}


float Analyzer_tautau::TotTMass_F(TLorentzVector a, TLorentzVector b, TLorentzVector met){
   float totalTMass = ( a + b + met).M();
   return totalTMass;
}

float Analyzer_tautau::VisMass_F(TLorentzVector a, TLorentzVector b){
  float visibleMass = (a + b).M();
  return visibleMass;
}

float Analyzer_tautau::pTvecsum_F(float pt1, float pt2, float phi1, float phi2){
   float pt_vecSum = sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2), 2) + pow(pt1*sin(phi1) + pt2*sin(phi2), 2));
  return pt_vecSum;
}


double Analyzer_tautau::dR(int tau1_index, int tau2_index){
   double deltaeta = abs(tauEta->at(tau1_index) - tauEta->at(tau2_index));
   double deltaphi = DeltaPhi(tauPhi->at(tau1_index), tauPhi->at(tau2_index));
   double deltaR = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
   return deltaR;
}


double Analyzer_tautau::DeltaPhi(double phi1, double phi2){
   double pi = TMath::Pi();
   double dphi = phi1-phi2;
   if (dphi>pi) dphi = 2.0*pi - dphi;
   if (dphi<= -1*pi) dphi = 2.0*pi+dphi;
   return fabs(dphi);
}



















