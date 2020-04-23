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

  std::string SampleName = argv[5];

  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
    {
      std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
      return 1;
    }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
    {
      std::cout<<"Please enter a valid value for reportEvery (parameter 4) "<<std::endl;
      return 1;
    }
  //std::string SampleName = argv[5];
  Analyzer_etau t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery, SampleName);
  return 0;
}

void Analyzer_etau::Loop(Long64_t maxEvents, int reportEvery, string SampleName)
{

   if (fChain == 0) return;
   int nTotal;
   nTotal = 0;
   int nInspected;
   nInspected = 0; 
   double nInspected_genWeighted; 
   nInspected_genWeighted = 0.0; 
   
   bool debug=false;  double netWeight = 1.0;
   int nBackupTriggerEvents, nBTMediumEvents, nBTMediumHLTsinglePhoEvents, nEffPhoptden, nEffPhoptnum, nEffMETden, nEffMETnum;
   nBackupTriggerEvents = nBTMediumEvents = nBTMediumHLTsinglePhoEvents = nEffPhoptden = nEffPhoptnum = nEffMETden = nEffMETnum = 0;
   double nHLTPassed,n_eventWeight, nSingleTrgPassed, nGoodElectronPassed, nElectronPtPassed, nGoodTauPassed, nTauPtPassed,numberOfEvents,nMETPassed, nDPhiPassed, nqcdden,nqcdnum, nMETFiltersPassed,nLeptonVetoPassed,nPassedBjetVeto,nNoisyCrystals,nDPhiJetMETPassed, nGoodETauPassed,nDeltaRPassed,nPassedThirdLepVeto, nPassedHiggsPtcut, nPassedVisibleMasscut, nPassedMETcut, nFinal_afterSelections,nGoodElectronPassed_qcd ,nGoodTauPassed_qcd, nGoodETauPassed_qcd,nDeltaRPassed_qcd,nPassedThirdLepVeto_qcd, nPassedHiggsPtcut_qcd, nPassedVisibleMasscut_qcd, nPassedMETcut_qcd, nFinal_afterSelections_qcd, nPassedBjetVeto_qcd,nMtPassed,nPasstottrmass, nPassJetsSelection, nL1PrefiringPassed;
   nHLTPassed = n_eventWeight = nSingleTrgPassed = nGoodElectronPassed = nElectronPtPassed = nGoodTauPassed = nTauPtPassed= numberOfEvents = nMETPassed = nDPhiPassed = nqcdden= nqcdnum=nMETFiltersPassed= nLeptonVetoPassed=nPassedBjetVeto=nNoisyCrystals=nDPhiJetMETPassed= nGoodETauPassed = nDeltaRPassed= nPassedThirdLepVeto=nPassedHiggsPtcut=nMtPassed=nPassedVisibleMasscut=nPassedMETcut=nFinal_afterSelections=nGoodElectronPassed_qcd=nGoodTauPassed_qcd=nGoodETauPassed_qcd = nDeltaRPassed_qcd= nPassedThirdLepVeto_qcd=nPassedHiggsPtcut_qcd=nPassedVisibleMasscut_qcd=nPassedMETcut_qcd=nFinal_afterSelections_qcd=nPassedBjetVeto_qcd=nPasstottrmass= nPassJetsSelection=nL1PrefiringPassed=0;

   std::vector<int> eleCand;
   eleCand.clear();
   
   std::vector<int> tauCand;
   tauCand.clear();


   TString sample = TString(SampleName);

   Float_t met_bins[16]={0.0, 15, 30, 45, 60, 90, 105, 120, 135, 150,175., 190., 250., 400., 600.0,800.0};
   TH1F* h_Events_level= new TH1F("Events_level","Events at each level of selection",15,0,30);
   TH1F* h_Events_level_qcd= new TH1F("Events_level_qcd","Events at each level of selection",15,0,30);

   TH1F* h_Cutflow= new TH1F("Cutflow","Events at each level of selection",6,0,12);
   TH1F* h_Cutflow_qcd= new TH1F("Cutflow_qcd","Events at each level of selection",6,0,12);

   TH1F* h_metfilter= new TH1F("metfilter","metfilter values",15,0,30);

   TH1F* h_MET_0= new TH1F("MET_0","MET before any selection",15,met_bins );  h_MET_0->Sumw2();
   TH1F* h_MET_1= new TH1F("MET_1","MET after met filters",15,met_bins ); h_MET_1->Sumw2();
   TH1F* h_MET_2= new TH1F("MET_2","MET after single ele trg",15,met_bins ); h_MET_2->Sumw2();
   TH1F* h_MET_3= new TH1F("MET_3","MET after ele selection",15,met_bins ); h_MET_3->Sumw2();
   TH1F* h_MET_4= new TH1F("MET_4","MET after tau selection",15,met_bins ); h_MET_4->Sumw2();

   //TH1F* h_genWeight = new TH1F((sample+"genWeight").c_str(), "genWeight",3000, -30000.0, 30000.0);h_genWeight[i]->Sumw2();
   //TH1F* h_genHT = new TH1F((sample+"genHT").c_str(), "genHT",15,PtBins);h_genHT[i]->Sumw2();
   TH1F* h_isData= new TH1F("isData","isData",4,0.0,2.0 );  h_isData->Sumw2();
   TH1F* h_insEvents= new TH1F("insEvents","insEvents",2,0.0,2.0 );  h_insEvents->Sumw2();

   int iphi = 41;
   int ieta = 5;

   //TFile *f_eleReconstrucSF_lowpt=new TFile("egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root");
   TFile *f_eleReconstrucSF_highpt=new TFile("egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root");
   TFile *f_eleIDeffSF=new TFile("egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root");
   
   //   TH2F *h_eleRecoSF_lowpt=(TH2F*) f_eleReconstrucSF_lowpt->Get("EGamma_SF2D");
   TH2F *h_eleRecoSF_highpt=(TH2F*) f_eleReconstrucSF_highpt->Get("EGamma_SF2D");
   TH2F *h_eleIDSF=(TH2F*) f_eleIDeffSF->Get("EGamma_SF2D");

   TFile* prefiring_file = TFile::Open("L1PrefiringMaps_new.root");
   TH2F* prefiring_photon = (TH2F*)prefiring_file->Get("L1prefiring_photonptvseta_2017BtoF");
   TH2F* prefiring_jet = (TH2F*)prefiring_file->Get("L1prefiring_jetptvseta_2017BtoF");

   //   TFile *f_kfactors = new TFile("kfactors.root");
   //TH1D *ewkCorrection = (TH1D*)f_kfactors->Get("EWKcorr/W");
   //TH1D *NNLOCorrection = (TH1D*)f_kfactors->Get("WJets_LO/inv_pt");

   //std::cout<<"Pont 1" << endl;
   Long64_t nentries = fChain->GetEntries();
   std::cout<<"Coming in: "<<std::endl;
   std::cout<<"nentries:"<<nentries<<std::endl;
   //Look at up to maxEvents events, or all if maxEvents == -1.
   Long64_t nentriesToCheck = nentries;
	if (maxEvents != -1LL && nentries > maxEvents)
		nentriesToCheck = maxEvents;
	nTotal = nentriesToCheck;
	Long64_t nbytes = 0, nb = 0;

	std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
	TStopwatch sw;
	sw.Start();
	for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
	{
	  
	  event_.clear();
	  event_info.clear();
	  
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  //if (jentry < 90000) continue;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;
	  double inspected_event_weight = 1.0; 
	  fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight/fabs(genWeight) : inspected_event_weight = 0.0;
	  nInspected_genWeighted += inspected_event_weight;  
	  nInspected += 1; 
          h_insEvents->SetBinContent(1, nInspected_genWeighted);
	  
	  //=1.0 for real data
	  double event_weight=1.0;
	  double eleRecoSF_corr=1.0;
          double eleEffSF_corr=1.0;
	  
	  double eletrgsf = 1.0;
	  double sf_tauID = 1.0;
	  double sf_Zvtx = 0.991;
	  
	  double EWK_corrected_weight=1.0;
	  double NNLO_weight = 1.0;
	  double kfactor = 1.0;

	  /*	  int bosonPID;
	  double bosonPt;
	  bool Wfound = false;
	  //check which mc particle is W boson                                                                                               
	  for(int i=0; i<nMC;i++){
	    if((*mcPID)[i] == 24){
	      Wfound=true;
	      bosonPID = (*mcPID)[i];
	      bosonPt = (*mcPt)[i];
	    }
	  }
	  */
	  eleCand = getEleCand(40,2.1); // Electron pt, eta cuts
	  tauCand = getTauCand(20,2.3); // Tau pt, eta cuts

	  if (debug==true && jentry%reportEvery == 0){ std::cout<<"Sample name ="<< sample <<endl;}
	  
	  
	  //h_metfilter->Fill(metFilters, event_weight);
	  //std::cout<<"event_weight (PA)=  "<< event_weight<<std::endl;
	  numberOfEvents+=event_weight;
         
	  h_MET_0->Fill(pfMET, event_weight);
	  h_isData->Fill(isData, event_weight);
	  if(metFilters==0) 
	    {
	      fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0;

	      nMETFiltersPassed+=event_weight;
	      h_metfilter->Fill(metFilters, event_weight);
	      //std::cout<<"event_weight (PB)=  "<< event_weight<<std::endl;

	      h_MET_1->Fill(pfMET, event_weight);
	      if (debug==true && jentry%reportEvery == 0){ std::cout<<"This works! Met filters" <<endl;}
	      
	      
	      if(HLTEleMuX>>3&1 == 1) // Single electron trigger
		{
		  nSingleTrgPassed+=event_weight;
		  if (debug==true && jentry%1000 == 0) std::cout<<"This works! Single electron trigger" <<endl;
		  h_MET_2->Fill(pfMET, event_weight);
		  
		  event_weight = event_weight*prefiring_weight(prefiring_photon, prefiring_jet);
		  nL1PrefiringPassed+=event_weight; 

		  //std::cout<<"event_weight (PC) =  "<< event_weight<<std::endl;


		  /*	  if (sample=="WJetsToLNu" || sample=="WJetsToLNu_inc" ||sample=="WJetsToLNu_HT-100To200" || sample=="WJetsToLNu_HT-1200To2500" || sample==" WJetsToLNu_HT-200To400"  || sample=="WJetsToLNu_HT-2500ToInf" || sample==" WJetsToLNu_HT-400To600" || sample=="WJetsToLNu_HT-600To800" || sample=="WJetsToLNu_HT-800To1200" || sample=="WJetsToLNu_ext"  || sample=="WJetsToLNu_inc_ext" || sample=="W1JetsToLNu"  || sample=="W2JetsToLNu"  || sample=="W3JetsToLNu"  || sample=="W4JetsToLNu")
		    {
		  EWK_corrected_weight = 1.0*(ewkCorrection->GetBinContent(ewkCorrection->GetXaxis()->FindBin(bosonPt)));
		  NNLO_weight = 1.0*(NNLOCorrection->GetBinContent(NNLOCorrection->GetXaxis()->FindBin(bosonPt)));
		  if(EWK_corrected_weight!=0 && NNLO_weight!=0)
		    kfactor = (EWK_corrected_weight/NNLO_weight);
		  else
		    kfactor=1.21;
		  event_weight*=kfactor;
		  sf_tauID =1.0;
		    }
		  */
		  if(eleCand.size() >0) 
		    {
		      nGoodElectronPassed+=event_weight;
		      if (debug==true &&jentry%1000 == 0) std::cout<<"This works! eleCand" <<endl;
		      
		      //std::cout<<"event_weight (PD) =  "<< event_weight<<std::endl;
		      h_MET_3->Fill(pfMET, event_weight);
		      // if (jentry%1000 == 0) std::cout<<"This works! P 2" <<endl;
		      if( tauCand.size() >0 )
			{
			  nGoodTauPassed+=event_weight;
			  h_MET_4->Fill(pfMET, event_weight);
			  

			  if (debug==true ) std::cout<<"event_weight (PE) =  "<< event_weight<<std::endl;
			  
			  eleRecoSF_corr=h_eleRecoSF_highpt->GetBinContent(h_eleRecoSF_highpt->GetXaxis()->FindBin(eleSCEta->at(eleCand[0])),h_eleRecoSF_highpt->GetYaxis()->FindBin(elePt->at(eleCand[0])));
			  if (debug==true ) std::cout<<"eleRecoSF_corr =  "<< eleRecoSF_corr<<std::endl;
			  eleEffSF_corr=h_eleIDSF->GetBinContent(h_eleIDSF->GetXaxis()->FindBin(eleSCEta->at(eleCand[0])),h_eleIDSF->GetYaxis()->FindBin(elePt->at(eleCand[0])));
			  if (debug==true ) std::cout<<"eleEffSF_corr =  "<< eleEffSF_corr<<std::endl;
			  
			  if (debug==true ) std::cout<<"This works line 269 "<<std::endl;
			  eletrgsf = EletriggerSF(elePt->at(eleCand[0]), eleEta->at(eleCand[0]));
			  if (debug==true ) std::cout<<"eletrgsf =  "<< eletrgsf << "Line 271"<<std::endl;

			  event_weight=event_weight*eleRecoSF_corr*eleEffSF_corr*sf_Zvtx*sf_tauID*eletrgsf;
			  //if (debug==true ) std::cout<<"event_weight (AFTER SF) =  "<< event_weight<<std::endl;

			  if (debug==true ) std::cout<<"event_weight =  "<< event_weight<<" event number = "<<jentry <<std::endl;
			  netWeight = event_weight;
			  if( (eleCharge->at(eleCand[0]))*(tauCharge->at(tauCand[0])) < 0 ) // opposite charge condition
			    {
			      if (debug==true ) std::cout<<"This works! Charge criteria" <<endl;
			      nGoodETauPassed+=event_weight;
			      fillHistos(0,event_weight,eleCand[0], tauCand[0]);
			      //std::cout<<"event_weight (AFTER SF 2) =  "<< event_weight<<std::endl;

			      /*
			      if ( (eleEta->at(eleCand[0]) < 0.0 && eleEta->at(eleCand[0]) > -0.3) || (eleEta->at(eleCand[0]) < -1.5 && eleEta->at(eleCand[0]) > -1.8))
				{
				  std::cout<<"event_weight (in eta range) =  "<< event_weight<<std::endl;
				
				  std::cout<<"inspected_event_weight (in eta range) =  "<<inspected_event_weight<<std::endl;
				}
			      */
			      double deltaR = dR(eleCand[0],tauCand[0]);
			      if(deltaR > 0.3 )
				{
				  if (debug==true && jentry%5 == 0) std::cout<<"This works! dR cut" <<endl;
				  nDeltaRPassed+=event_weight;
				  fillHistos(1,event_weight,eleCand[0], tauCand[0]);
				  //std::cout<<"event_weight (AFTER SF 3) =  "<< event_weight<<std::endl;

				  if( thirdLeptonVeto()==true )
				  {
				    nPassedThirdLepVeto+=event_weight;
				    if (debug==true && jentry%5 == 0) std::cout<<"This works! third lepon veto" <<endl;
				    fillHistos(2,event_weight,eleCand[0], tauCand[0]);
				    if( passBjetVeto() == true)
				      {
					nPassedBjetVeto+=event_weight;
					if (debug==true && jentry%5 == 0) std::cout<<"This works! bJetVeto" <<endl;
					fillHistos(3,event_weight,eleCand[0], tauCand[0]);
					if (debug==true && jentry%5 == 0) std::cout<<"This works! bJetVeto **** after fill" <<endl;
					
					
					float eleMet_trMass = TMass_F((elePt->at(eleCand[0])),(elePhi->at(eleCand[0])),pfMET,pfMETPhi  );
					if( (eleMet_trMass>50) )
					  {
					    fillHistos(4,event_weight,eleCand[0], tauCand[0]);
					    nMtPassed+=event_weight;
					    double HiggsPt = pTvecsum_F(elePt->at(eleCand[0]),tauPt->at(tauCand[0]),elePhi->at(eleCand[0]),tauPhi->at(tauCand[0]));
					    if( (HiggsPt>65) )
					      {
						nPassedHiggsPtcut+=event_weight;
						fillHistos(5,event_weight,eleCand[0], tauCand[0]);
						
						TLorentzVector myTau;
						myTau.SetPtEtaPhiE(tauPt->at(tauCand[0]),tauEta->at(tauCand[0]),tauPhi->at(tauCand[0]), tauEnergy->at(tauCand[0]));
						TLorentzVector myEle;
						myEle.SetPtEtaPhiE(elePt->at(eleCand[0]),eleEta->at(eleCand[0]),elePhi->at(eleCand[0]), eleEn->at(eleCand[0]));
						double visMass_etau = VisMass_F(myTau, myEle);
						
						if( (visMass_etau < 125) ) 
						  {
						    nPassedVisibleMasscut+=event_weight;
						    fillHistos(6,event_weight,eleCand[0], tauCand[0]);
						    
						    
						    if( (pfMET > 105) )
						      {
							nPassedMETcut+=event_weight;
							fillHistos(7,event_weight,eleCand[0], tauCand[0]);
							
							if (debug==true && jentry%5 == 0) std::cout<<"******************* This works! All selections passed ***************" <<endl;
							nFinal_afterSelections+=event_weight;
							
							TLorentzVector myTau_1; 
							myTau_1.SetPtEtaPhiE(tauPt->at(tauCand[0]),tauEta->at(tauCand[0]),tauPhi->at(tauCand[0]), tauEnergy->at(tauCand[0]));		
							TLorentzVector myEle_1; 
							myEle_1.SetPtEtaPhiE(elePt->at(eleCand[0]),eleEta->at(eleCand[0]),elePhi->at(eleCand[0]), eleEn->at(eleCand[0]));
							TLorentzVector myMet_1;
							myMet_1.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
							float tot_tr_mass = TotTMass_F(myEle_1, myTau_1, myMet_1 );
							
							if (tot_tr_mass > 260)
							  {
							    nPasstottrmass+=event_weight;
							    fillHistos(8,event_weight,eleCand[0], tauCand[0]);
							    if (check_jets()==false)
							      {
								nPassJetsSelection+=event_weight;
								fillHistos(9,event_weight,eleCand[0], tauCand[0]);
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
	    }
	      
	
	  //h_Events_level= new TH1F("Events_level","Events at each level of selection",20,0,20);
	  h_Events_level->SetBinContent(1, nInspected_genWeighted);
	  h_Events_level->SetBinContent(2, nMETFiltersPassed);
          h_Events_level->SetBinContent(3, nSingleTrgPassed);          
	  h_Events_level->SetBinContent(4, nGoodElectronPassed);
	  h_Events_level->SetBinContent(5, nGoodTauPassed);	  
          h_Events_level->SetBinContent(6, nGoodETauPassed);
          h_Events_level->SetBinContent(7, nDeltaRPassed);
          h_Events_level->SetBinContent(8, nPassedThirdLepVeto); 
	  h_Events_level->SetBinContent(9, nPassedBjetVeto);
	  h_Events_level->SetBinContent(10, nMtPassed);
	  h_Events_level->SetBinContent(11, nPassedHiggsPtcut); 
	  h_Events_level->SetBinContent(12, nPassedVisibleMasscut); 
	  h_Events_level->SetBinContent(13, nPassedMETcut); 
	  
	  //          h_Events_level->SetBinContent(13, nFinal_afterSelections);
          h_Cutflow->SetBinContent(1, nGoodElectronPassed);
          h_Cutflow->SetBinContent(2, nGoodTauPassed);
          h_Cutflow->SetBinContent(3, nGoodETauPassed);
          h_Cutflow->SetBinContent(4, nDeltaRPassed);
          h_Cutflow->SetBinContent(5, nPassedThirdLepVeto);
          h_Cutflow->SetBinContent(6, nPassedBjetVeto);
         
        
	 
	  tree->Fill();
	  
	  if (jentry%reportEvery == 0)
	    {
	      std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
	    }
	}
	//h_Events_level->Write();
	std::cout.setf( std::ios::fixed, std:: ios::floatfield );
	if((nentriesToCheck-1)%reportEvery != 0)
		std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
	sw.Stop();
	std::cout<<"All events checked."<<std::endl;

	std::cout<<"*******************************************"<<std::endl;
	std::cout<<"******************Jithin's original*************************"<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"Initial entries "<<numberOfEvents<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"Inspected genWeightd "<<nInspected_genWeighted<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"METFiltersPassed "<<nMETFiltersPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"SingleTrgPassed "<<nSingleTrgPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"nL1PrefiringPassed "<<nL1PrefiringPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"GoodElectronPassed "<<nGoodElectronPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"GoodTauPassed "<<nGoodTauPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"opp charge "<<nGoodETauPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"DeltaRPassed "<<nDeltaRPassed<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"PassedThirdLepVeto "<<nPassedThirdLepVeto<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"PassedBjetVeto "<<nPassedBjetVeto<<std::endl;
	std::cout<<"*******************************************"<<std::endl;

	std::cout<<"*******************************************"<<std::endl;
	std::cout<<std::setw(20) <<std::right <<"Number of events inspected: " << nInspected <<std::endl;
	std::cout<<std::setw(20) <<std::right << "Number of events inspected (minus negative gen. weights): " << nInspected_genWeighted << std::endl; 
	std::cout<<std::setw(20) <<std::right << " (net weight with SF): " << netWeight << std::endl;
	std::cout<<"*******************************************"<<std::endl;
	
	f_eleReconstrucSF_highpt->Close();
	f_eleIDeffSF->Close();
	//f_kfactors->Close();

}

void Analyzer_etau::BookHistos(const char* file2)
{
	fileName = new TFile(file2, "RECREATE");
	tree = new TTree("ADD","ADD");
	tree->Branch("event_","std::vector<unsigned int>",&event_);
	tree->Branch("event_info","std::vector<double>",&event_info);
	fileName->cd();
	
	Float_t PtBins[21]={0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,110, 120,130, 140,150, 160,170, 180,190, 200.};
	Float_t MetBins[15]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180., 200, 300., 400., 600.0,800.0};
	Float_t TrMassBins[24]={0.0, 20, 40, 60, 80, 100, 120, 140, 160,180., 200, 220, 240,260,280,300.,320,340,360,380, 400., 600.0,800.0, 1000.0};


	//Set up the histos to be filled with method fillHistos
	for(int i=0; i<20; i++)
	{
		char ptbins[100];
		sprintf(ptbins, "_%d", i);
		std::string histname(ptbins);
		h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(), "nVtx",40,0,40);h_nVtx[i]->Sumw2();
		h_nEvents[i] = new TH1F(("nEvents"+histname).c_str(), "nEvents",3,0,3);h_nEvents[i]->Sumw2();
                h_genWeight[i] = new TH1F(("genWeight"+histname).c_str(), "genWeight",10, -10.0, 10.0);h_genWeight[i]->Sumw2();
                h_genHT[i] = new TH1F(("genHT"+histname).c_str(), "genHT",20,PtBins);h_genHT[i]->Sumw2();
		//**************  electrons  **************
		h_electron_En[i] = new TH1F(("Electron_En"+histname).c_str(), "Electron_En",20,PtBins);h_electron_En[i]->Sumw2();
		h_electron_Pt[i] = new TH1F(("Electron_Pt"+histname).c_str(), "Electron_Pt",20,PtBins);h_electron_Pt[i]->Sumw2();
		h_electron_eta[i] = new TH1F(("Electron_eta"+histname).c_str(), "Electron_eta",20,-3.0, 3.0);h_electron_eta[i]->Sumw2();
		h_electron_SCEta[i] = new TH1F(("Electron_SCeta"+histname).c_str(), "Electron_SCeta",20, -3.0, 3.0);h_electron_SCEta[i]->Sumw2();
		h_electron_phi[i] = new TH1F(("Electron_phi"+histname).c_str(), "Electron_phi", 21,-3.14,3.14);h_electron_phi[i]->Sumw2();
		h_electron_SCPhi[i] = new TH1F(("electron_SCphi"+histname).c_str(), "electron_SCphi", 21,-3.14,3.14);h_electron_SCPhi[i]->Sumw2();
		h_electron_IDbit[i] = new TH1F(("electron_ID_bit"+histname).c_str(), "electron_ID_bit",8,0,8);h_electron_IDbit[i]->Sumw2();
		//**************  tau  **************
		h_tau_En[i] = new TH1F(("Tau_En"+histname).c_str(), "Tau_En",20,PtBins);h_tau_En[i]->Sumw2();
		h_tau_Pt[i] = new TH1F(("Tau_Pt"+histname).c_str(), "Tau_Pt",20,PtBins);h_tau_Pt[i]->Sumw2();
		h_tau_eta[i] = new TH1F(("Tau_eta"+histname).c_str(), "Tau_eta",20,-3.0, 3.0);h_tau_eta[i]->Sumw2();

		h_tau_phi[i] = new TH1F(("Tau_phi"+histname).c_str(), "Tau_phi", 21,-3.14,3.14);h_tau_phi[i]->Sumw2();
		//**************  Met  **************
		h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "pfMET",14,MetBins);h_pfMET[i]->Sumw2();
		h_pfMET_phi[i] = new TH1F(("pfMET_phi"+histname).c_str(), "pfMET_phi",21,-3.14,3.14);h_pfMET_phi[i]->Sumw2();

		h_dPhi[i] = new TH1F(("h_dPhi"+histname).c_str(),"h_dPhi",20,0,3.15);h_dPhi[i]->Sumw2();
                h_dR[i] = new TH1F(("h_dR"+histname).c_str(),"h_dR",20,0,3.15);h_dR[i]->Sumw2();

		h_nJet[i] = new TH1F(("nJet"+histname).c_str(), "nJet",20,0,20);h_nJet[i]->Sumw2();
		h_leadingJetPt[i] = new TH1F(("leadingJetPt"+histname).c_str(),"leadingJetPt",30,20,1000);h_leadingJetPt[i]->Sumw2();

		h_leadingJetEta[i] = new TH1F(("h_leadingJetEta"+histname).c_str(),"h_leadingJetEta",40,-1.4442,1.4442);h_leadingJetEta[i]->Sumw2();

		//**************  rest  **************

		h_Mt[i]= new TH1F(("Mt"+histname).c_str(),"MT",23,TrMassBins);h_Mt[i]->Sumw2();
		h_totTrMass[i]= new TH1F(("TotalTrMass"+histname).c_str(),"TotalTrMass",23,TrMassBins);h_totTrMass[i]->Sumw2();

		h_VisibleMass[i]= new TH1F(("VisibleMass"+histname).c_str(),"VisibleMass",40,0,200);h_VisibleMass[i]->Sumw2();
		h_HiggsPt[i]= new TH1F(("HiggsPt"+histname).c_str(),"HiggsPt",20,PtBins);h_HiggsPt[i]->Sumw2();
		h_electronIso[i]= new TH1F(("Electron iso"+histname).c_str(),"Electron iso",40,0.0,1.0);h_electronIso[i]->Sumw2();
		h_tauIso[i]= new TH1F(("Tau iso"+histname).c_str(),"Tau iso",10,-2.0,2.0);h_tauIso[i]->Sumw2();


	}
}

//Fill the sequential histos at a particular spot in the sequence
void Analyzer_etau::fillHistos(int histoNumber, double event_weight,int eleIndex, int tauIndex)
{
	//*********** fill electrons  ***********
		h_electron_En[histoNumber]->Fill((eleEn->at(eleIndex)),event_weight);
		h_electron_Pt[histoNumber]->Fill((elePt->at(eleIndex)),event_weight);		
		h_electron_eta[histoNumber]->Fill(eleEta->at(eleIndex),event_weight);
		h_electron_SCEta[histoNumber]->Fill(eleSCEta->at(eleIndex),event_weight);
		h_electron_phi[histoNumber]->Fill(elePhi->at(eleIndex),event_weight);
		h_electron_SCPhi[histoNumber]->Fill(eleSCPhi->at(eleIndex),event_weight);
		h_electron_IDbit[histoNumber]->Fill( eleIDbit->at(eleIndex)>>3&1,event_weight);
	//*********** fill taus  ***********

		h_tau_En[histoNumber]->Fill((tauEnergy->at(tauIndex)),event_weight);
		h_tau_Pt[histoNumber]->Fill((tauPt->at(tauIndex)),event_weight);
		h_tau_eta[histoNumber]->Fill(tauEta->at(tauIndex),event_weight);
		h_tau_phi[histoNumber]->Fill(tauPhi->at(tauIndex),event_weight);
		
	//*********** fill met  ***********		
		h_pfMET[histoNumber]->Fill(pfMET,event_weight);
		h_pfMET_phi[histoNumber]->Fill(pfMETPhi,event_weight);

		double dPhi_etau = DeltaPhi(elePhi->at(eleIndex),tauPhi->at(tauIndex));
		h_dPhi[histoNumber]->Fill(dPhi_etau,event_weight);
		double dR_etau = dR(eleIndex,tauIndex);
		h_dR[histoNumber]->Fill(dR_etau,event_weight);


	//*********** fill rest  ***********

		float mT_eMet = TMass_F((elePt->at(eleIndex)),(elePhi->at(eleIndex)),pfMET,pfMETPhi  );
		h_Mt[histoNumber]->Fill(mT_eMet,event_weight);
		TLorentzVector myTau; 
		myTau.SetPtEtaPhiE(tauPt->at(tauIndex),tauEta->at(tauIndex),tauPhi->at(tauIndex), tauEnergy->at(tauIndex));		
		TLorentzVector myEle; 
		myEle.SetPtEtaPhiE(elePt->at(eleIndex),eleEta->at(eleIndex),elePhi->at(eleIndex), eleEn->at(eleIndex));
		double visMass_etau = VisMass_F(myTau, myEle);
		h_VisibleMass[histoNumber]->Fill(visMass_etau,event_weight);
		TLorentzVector myMet;
                myMet.SetPtEtaPhiE(pfMET ,0,pfMETPhi,pfMET);
		float tot_tr_mass = TotTMass_F(myEle, myTau, myMet );
		h_totTrMass[histoNumber]->Fill(tot_tr_mass,event_weight);
		double HiggsPt = pTvecsum_F(elePt->at(eleIndex),tauPt->at(tauIndex),elePhi->at(eleIndex),tauPhi->at(tauIndex) );
		h_HiggsPt[histoNumber]->Fill(HiggsPt,event_weight);
		h_nVtx[histoNumber]->Fill(nVtx,event_weight);
		h_nEvents[histoNumber]->Fill(isData,event_weight);
		
		h_nJet[histoNumber]->Fill(nJet ,event_weight);
		//		if(jetsize>0){
		//h_leadingJetPt[histoNumber]->Fill(jetPt->at(),event_weight);
		//h_leadingJetPt_300[histoNumber]->Fill(jetPt->at(),event_weight);
		//h_leadingJetEta[histoNumber]->Fill(jetEta->at(),event_weight);
		//}
		h_genHT[histoNumber]->Fill(genHT,event_weight);
		h_genWeight[histoNumber]->Fill(genWeight,event_weight);
		float rel_eleIso = ( elePFChIso->at(eleIndex) + max( elePFNeuIso->at(eleIndex) + elePFPhoIso->at(eleIndex) - 0.5 *elePFPUIso->at(eleIndex) , 0.0 )) / (elePt->at(eleIndex));
		h_electronIso[histoNumber]->Fill(rel_eleIso,event_weight);
		h_tauIso[histoNumber]->Fill(tauByTightIsolationMVArun2v1DBoldDMwLT->at(tauIndex),event_weight);
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
      bool decayModeCut = false;
      bool tauIsolation = false;

      if( tauPt->at(iTau) > tauPtCut  && fabs( tauEta->at(iTau))< tauEtaCut && taudz->at(iTau)<0.2 )kinematic = true;
      if( tauByMVA6TightElectronRejection->at(iTau) == 1 && tauByLooseMuonRejection3->at(iTau) == 1) tauId = true;
      if( (taubyTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau)==1)) tauIsolation = true;
      if( tauDecayMode->at(iTau)==0 || tauDecayMode->at(iTau)==1 || tauDecayMode->at(iTau)==10 ) decayModeCut = true;
      
      if(tauId==true && kinematic==true && tauIsolation==true && decayModeCut==true)
	//if(tauId==true && kinematic==true && decayModeCut==true)
	{
	  tmpCand.push_back(iTau);
	}                                                                                                                                                              
    }                                                                                                                                                                
  
  return tmpCand;

}


bool Analyzer_etau::thirdLeptonVeto(){

  //  std::vector<int> tmpCand;
  //tmpCand.clear();
  int tmpCand = 0;
  bool thirdLep = false;
  bool thirdLepVeto=true;
   //Loop over muons                                                                                                                       

   for(int iMu=0; iMu < nMu;iMu++)
     {
       bool kinematic = false;
       if( (*muPt)[iMu] > 10.0  && fabs((*muEta)[iMu])< 2.4 && (*muD0)[iMu] < 0.045 && (*muDz)[iMu] < 0.2 ) kinematic = true;
       bool muonId =false;
       if( muIDbit->at(iMu)>>0&1==1) muonId =true;
       bool relative_iso = false;

       float relMuIso = ( muPFChIso->at(iMu) + max( muPFNeuIso->at(iMu) + muPFPhoIso->at(iMu) - 0.5 *muPFPUIso->at(iMu) , 0.0 )) / (muPt->at(iMu));
       if( relMuIso < 0.3 ) relative_iso = true;
       if(muonId==true && kinematic==true && relative_iso==true){
	 //tmpCand.push_back(iMu);
	 tmpCand++;
       }                                                                                                                                                    
     }                                                                                                                      
   if(tmpCand > 0){ thirdLep = true; thirdLepVeto=false;}
   //   if( thirdLep = true ) thirdLepVeto=false;
   //if( tmpCand = 0 ) thirdLepVeto=true;

   return thirdLepVeto;

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
  double dphi = phi1-phi2;
  if(dphi>pi) dphi = 2.0*pi - dphi;
  if(dphi<= -1*pi) dphi =  2.0*pi +dphi;
  return fabs(dphi);
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
  bool veto = true;
  bool foundBjet = false;
  for(int iJets=0; iJets<nJet ; iJets++){
    if(jetCSV2BJetTags->at(iJets) > 0.8838) tmpCand++;
    // CSV B jet tag for selecting bJets is medium WP (jetCSV2BJetTags > 0.8838.)
  }
  if(tmpCand>0){ veto = false; foundBjet = true; } 
  return veto;
}
//from https://lathomas.web.cern.ch/lathomas/HLTDQM/Full2017/SingleElectron/2017BtoF/elewptight35_ptvseta_ele_den_ptvseta.pdf
float Analyzer_etau::EletriggerSF(float pt, float eta){
  float sf = 1.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 0.8)
    {
      if(pt > 40.0  && pt < 50) sf = 0.79;
      if(pt > 50.0  && pt < 60) sf = 0.82;
      if(pt > 60.0  && pt < 100) sf = 0.85;
      if(pt > 100.0  && pt < 150) sf = 0.87;
      if(pt > 150.0  && pt < 200) sf = 0.88;
      if(pt > 200) sf = 0.89;
    }
  if(fabs(eta) >= 0.8   && fabs(eta) < 1.442 ) 
    {
      if(pt > 40.0  && pt < 50) sf = 0.77;
      if(pt > 50.0  && pt < 60) sf = 0.81;
      if(pt > 60.0  && pt < 100) sf = 0.85;
      if(pt > 100.0  && pt < 150) sf = 0.87;
      if(pt > 150.0  && pt < 300) sf = 0.89;
      if(pt > 300.0) sf = 0.87;
    }
    
  if(fabs(eta) >= 1.442   && fabs(eta) < 1.56 ) 
    {
      if(pt > 40.0  && pt < 50) sf = 0.73;
      if(pt > 50.0  && pt < 60) sf = 0.75;
      if(pt > 60.0  && pt < 100) sf = 0.76;
      if(pt > 100.0  && pt < 150) sf = 0.72;
      if(pt > 150.0  && pt < 300) sf = 0.78;
      if(pt > 300.0) sf = 0.67;


    }
  if(fabs(eta) >= 1.56   && fabs(eta) < 2.0 ) 
    {
      if(pt > 40.0  && pt < 50) sf = 0.80;
      if(pt > 50.0  && pt < 60) sf = 0.84;
      if(pt > 60.0  && pt < 100) sf = 0.87;
      if(pt > 100.0  && pt < 150) sf = 0.88;
      if(pt > 150.0  && pt < 300) sf = 0.89;
      if(pt > 300.0) sf = 0.87;


    }

  if(fabs(eta) >= 2.0   && fabs(eta) < 2.5 ) 

    {

      if(pt > 40.0  && pt < 50) sf = 0.73;
      if(pt > 50.0  && pt < 60) sf = 0.78;
      if(pt > 60.0  && pt < 100) sf = 0.83;
      if(pt > 100.0  && pt < 150) sf = 0.86;
      if(pt > 150.0  && pt < 300) sf = 0.89;
      if(pt > 300.0) sf = 0.86;

    }


  return sf;
}
bool Analyzer_etau::check_jets()
{
  int tmpCand=0;
  bool found_event = false;
  for (int iJets=0; iJets < nJet; iJets++)
    {
      if ( (jetPt->at(iJets)>100) && (fabs(jetEta->at(iJets)) < 3.0) && (fabs(jetEta->at(iJets) >2.25)))
        {
          tmpCand++;
        }
    }
  if(tmpCand>0)  found_event=true;
  return found_event;
}
double Analyzer_etau::prefiring_weight(TH2F* prefiring_photon, TH2F* prefiring_jet)
{
  std::vector<int> tmpCand_ele;
  tmpCand_ele.clear();
  std::vector<int> tmpCand_jet;
  tmpCand_jet.clear();

  double weight = 1.0;
  double map_val_ele = 0.0;
  int picked_ele = -1;
  double map_val_jet = 0.0;
  int x_bin;
  int y_bin;
  for (int iEle=0; iEle<nEle; iEle++)
    {
      if(fabs(eleEta->at(iEle)) > 2.0 && eleIDbit->at(iEle)>>2&1==1)
	{
	  x_bin = prefiring_photon->GetXaxis()->FindBin( eleEta->at(iEle));
	  y_bin = prefiring_photon->GetYaxis()->FindBin( elePt->at(iEle));
	  map_val_ele = prefiring_photon->GetBinContent(x_bin,y_bin);
	  picked_ele = iEle;
	  tmpCand_ele.push_back(iEle);
	}
    }
  double dEta;
  double dPhi;
  double dr;
  int x_bin_j;
  int y_bin_j;
  //for (int ele_i : tmpCand_ele) 
    {
      //picked_ele = ele_i;
      for (int iJet = 0; iJet<nJet; iJet++)
	{
	  if (fabs(jetEta->at(iJet))> 2.0)
	    {
	      if (picked_ele != -1)
		{
		  dEta = fabs(eleEta->at(picked_ele) - jetEta->at(iJet));
		  dPhi = DeltaPhi(elePhi->at(picked_ele), jetPhi->at(iJet));
		  dr = sqrt (dEta*dEta + dPhi*dPhi);
		  if (dr < 0.3) continue; 
		}
	      x_bin_j = prefiring_jet->GetXaxis()->FindBin( jetEta->at(iJet));
	      y_bin_j = prefiring_jet->GetYaxis()->FindBin( jetPt->at(iJet));
	      map_val_jet = prefiring_jet->GetBinContent(x_bin_j, y_bin_j);
	      tmpCand_jet.push_back(iJet);
	    }
	}
    }
  //std::cout<<ele_i<<" electrons  and "<< jet_i<< " jets "<< endl;
  //for (int x : tmpCand_ele ) std::cout<<"tmpCand_ele elements  "<< x<< endl;
  weight = (1-map_val_ele)*(1-map_val_jet);
  //if (weight != 1)  std::cout<<"weigh = " << weight << endl;

  return weight;
  
}
