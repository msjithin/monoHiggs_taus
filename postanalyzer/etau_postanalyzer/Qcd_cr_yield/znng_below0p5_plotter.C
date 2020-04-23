#include <fstream>
#include <vector>
#include <iomanip>
#include "TFile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TColor.h"

void plot_DM(string histname_string, string sample_name, Float_t xsec, int nevents_total, bool isNLO)
{
  TString histname = TString(histname_string+"_3"); // Photon+MET: 0, Dphi(photon,MET): 1, lepton veto: 2, Dphi(jets,MET): 3
  TString histname_JESUp = TString(histname_string+"_7"); // Photon+MET: 4, Dphi(photon,MET): 5, lepton veto: 6, Dphi(jets,MET): 7
  TString histname_JESDown = TString(histname_string+"_11"); // Photon+MET: 8, Dphi(photon,MET): 9, lepton veto: 10, Dphi(jets,MET): 11
  TString histname_PESUp = TString(histname_string+"_15"); // Photon+MET: 12, Dphi(photon,MET): 13, lepton veto: 14, Dphi(jets,MET): 15
  TString histname_PESDown = TString(histname_string+"_19"); // Photon+MET: 16, Dphi(photon,MET): 17, lepton veto: 18, Dphi(jets,MET): 19
  
  Float_t int_lumi = 36814.0;
  Float_t scale_factor = 0.984; // Flat pixel seed veto SF
  Float_t frac_below0p5 = 1.0;
  
  double photon_scale_factor_unc = 0.009;

  TFile *f_DM = new TFile(TString("ZnnG_JESPES_"+sample_name+".root"));
  TH1F* histo_DM = (TH1F*)((TH1F*)f_DM->Get(histname))->Clone(TString("histo_"+sample_name));
  TH1F* histo_DM_JESUp = (TH1F*)((TH1F*)f_DM->Get(histname_JESUp))->Clone(TString("histo_"+sample_name+"_JESUp"));
  TH1F* histo_DM_JESDown = (TH1F*)((TH1F*)f_DM->Get(histname_JESDown))->Clone(TString("histo_"+sample_name+"_JESDown"));
  TH1F* histo_DM_PESUp = (TH1F*)((TH1F*)f_DM->Get(histname_PESUp))->Clone(TString("histo_"+sample_name+"_PESUp"));
  TH1F* histo_DM_PESDown = (TH1F*)((TH1F*)f_DM->Get(histname_PESDown))->Clone(TString("histo_"+sample_name+"_PESDown"));
  const int nBins = histo_DM->GetXaxis()->GetNbins();
  histo_DM->Scale(int_lumi*scale_factor*frac_below0p5*xsec/nevents_total);
  histo_DM_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*xsec/nevents_total);
  histo_DM_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*xsec/nevents_total);
  histo_DM_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*xsec/nevents_total);
  histo_DM_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*xsec/nevents_total);
  Float_t int_DM = histo_DM->Integral()+histo_DM->GetBinContent(0)+histo_DM->GetBinContent(nBins+1);
  Float_t int_DM_JESUp = histo_DM_JESUp->Integral()+histo_DM_JESUp->GetBinContent(0)+histo_DM_JESUp->GetBinContent(nBins+1);
  Float_t int_DM_JESDown = histo_DM_JESDown->Integral()+histo_DM_JESDown->GetBinContent(0)+histo_DM_JESDown->GetBinContent(nBins+1);
  Float_t int_DM_PESUp = histo_DM_PESUp->Integral()+histo_DM_PESUp->GetBinContent(0)+histo_DM_PESUp->GetBinContent(nBins+1);
  Float_t int_DM_PESDown = histo_DM_PESDown->Integral()+histo_DM_PESDown->GetBinContent(0)+histo_DM_PESDown->GetBinContent(nBins+1);
  double jeserr_DM = (fabs(int_DM_JESUp-int_DM)+fabs(int_DM_JESDown-int_DM))/2.0;
  double peserr_DM = (fabs(int_DM_PESUp-int_DM)+fabs(int_DM_PESDown-int_DM))/2.0;
  Float_t err_DM = 0.0;
  if(int_DM > 0.0)
    err_DM = sqrt(int_DM*int_DM*((xsec*int_lumi*scale_factor*frac_below0p5-int_DM)/(nevents_total*int_DM))+(jeserr_DM*jeserr_DM)+(peserr_DM*peserr_DM));
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_DM->GetBinContent(i);
    double jesup = histo_DM_JESUp->GetBinContent(i);
    double jesdown = histo_DM_JESDown->GetBinContent(i);
    double pesup = histo_DM_PESUp->GetBinContent(i);
    double pesdown = histo_DM_PESDown->GetBinContent(i);
    // cout<<"DM: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((xsec*int_lumi*scale_factor*frac_below0p5-int_bin)/(nevents_total*int_bin));
    histo_DM->SetBinError(i,err_bin);
  }

  if (histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3")
  {
    TFile* f_DM_histos;
    if(histname == "Photon_Et_range_3"){
      if(isNLO)
        f_DM_histos = new TFile("DM_NLO_histos_below0p5_Pt.root", "UPDATE");
      else
        f_DM_histos = new TFile("DM_LO_histos_below0p5_Pt.root", "UPDATE");
    }
    else if(histname == "pfMET_3"){
      if(isNLO)
        f_DM_histos = new TFile("DM_NLO_histos_below0p5_MET.root", "UPDATE");
      else
        f_DM_histos = new TFile("DM_LO_histos_below0p5_MET.root", "UPDATE");
    }
    else if(histname == "Mt_3"){
      if(isNLO)
        f_DM_histos = new TFile("DM_NLO_histos_below0p5_Mt.root", "UPDATE");
      else
        f_DM_histos = new TFile("DM_LO_histos_below0p5_Mt.root", "UPDATE");
    }
    f_DM_histos->cd();
    histo_DM->Write();
    histo_DM_JESUp->Write();
    histo_DM_JESDown->Write();
    histo_DM_PESUp->Write();
    histo_DM_PESDown->Write();
    f_DM_histos->Close();
  }
}

void plot(string histname_string, Double_t leg_xoffset, Double_t leg_yoffset, TString xaxis_title, TString plotname)
{
  TString histname = TString(histname_string+"_3"); // Photon+MET: 0, Dphi(photon,MET): 1, lepton veto: 2, Dphi(jets,MET): 3
  TString histname_JESUp = TString(histname_string+"_7"); // Photon+MET: 4, Dphi(photon,MET): 5, lepton veto: 6, Dphi(jets,MET): 7
  TString histname_JESDown = TString(histname_string+"_11"); // Photon+MET: 8, Dphi(photon,MET): 9, lepton veto: 10, Dphi(jets,MET): 11
  TString histname_PESUp = TString(histname_string+"_15"); // Photon+MET: 12, Dphi(photon,MET): 13, lepton veto: 14, Dphi(jets,MET): 15
  TString histname_PESDown = TString(histname_string+"_19"); // Photon+MET: 16, Dphi(photon,MET): 17, lepton veto: 18, Dphi(jets,MET): 19
  TString histname_data = TString(histname_string+"_17"); // Photon: 12, Photon+MET: 13, Dphi(photon,MET): 14, lepton veto: 15, noisy crystal: 16, Dphi(jets,MET): 17
  TString histname_wenu_errUp = TString(histname_string+"_19"); // Dphi(jets,MET): 19
  TString histname_wenu_errDown = TString(histname_string+"_20"); // Dphi(jets,MET): 20
  TString histname_wenu_unweighted = TString(histname_string+"_21"); // Dphi(jets,MET): 21
  TString histname_qcd_sidebandUp = TString(histname_string+"_19");
  TString histname_qcd_sidebandDown = TString(histname_string+"_20");
  TString histname_qcd_METUp = TString(histname_string+"_21");
  TString histname_qcd_METDown = TString(histname_string+"_22");
  TString histname_qcd_binningUp = TString(histname_string+"_23");
  TString histname_qcd_binningDown = TString(histname_string+"_24");
  TString histname_qcd_sieieLeft = TString(histname_string+"_25");
  TString histname_qcd_sieieRight = TString(histname_string+"_26");
  TString histname_qcd_templateUp = TString(histname_string+"_27");
  TString histname_qcd_templateDown = TString(histname_string+"_28");
  TString histname_qcd_unweighted = TString(histname_string+"_29");
  TString histname_uncorrected = TString(histname_string+"_22");
  
  std::vector<TH1F*> histo_vector;
  histo_vector.clear();
  
  Float_t int_lumi = 36814.0;
  Float_t scale_factor = 0.984; // Flat pixel seed veto SF
  Float_t frac_below0p5 = 1.0/3.14159;
  
  double photon_scale_factor_unc = 0.009; // Flat pixel seed veto uncertainty
  
  double total_background = 0.0;
  double background_unc_sumsquares = 0.0;
  
  TCanvas *c = new TCanvas("c", "canvas",700,640);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  float t_m = 0.08; //top margin
  float b_m = 0.4; //botton margin
  float l_m = 0.09; //left margin
  float r_m = 0.05; //right margin
  c->SetTopMargin(t_m);
  c->SetBottomMargin(b_m);
  c->SetLeftMargin(l_m);
  c->SetRightMargin(r_m);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->cd();
  
  //DEBUG
  // cout<<"Setup complete"<<endl;

  TFile *f_data = new TFile("ZnnG_data_below0p5_all.root");
  // Switch this name after unblinding
  // TH1F* histo_data = (TH1F*)((TH1F*)f_data->Get(histname_data))->Clone("data_obs");
  TH1F* histo_data = (TH1F*)((TH1F*)f_data->Get(histname_data))->Clone("histo_data");
  const int nBins = histo_data->GetXaxis()->GetNbins();
  Float_t int_data = histo_data->Integral()+histo_data->GetBinContent(0)+histo_data->GetBinContent(nBins+1);
  histo_data->SetLineWidth(3);
  histo_data->SetLineColor(kWhite);
  histo_data->SetMarkerStyle(kFullSquare);
  histo_data->SetMarkerColor(kWhite);
  histo_vector.push_back(histo_data);
  
  //DEBUG
  // cout<<"Data got"<<endl;
  
  // Now that nBins has been specified,
  // initialize binned systematic shift containers to the appropriate length
  std::vector<double> jesup_shift;
  jesup_shift.clear();
  std::vector<double> jesdown_shift;
  jesdown_shift.clear();
  std::vector<double> pesup_shift;
  pesup_shift.clear();
  std::vector<double> pesdown_shift;
  pesdown_shift.clear();
  std::vector<double> renup_shift;
  renup_shift.clear();
  std::vector<double> rendown_shift;
  rendown_shift.clear();
  std::vector<double> facup_shift;
  facup_shift.clear();
  std::vector<double> facdown_shift;
  facdown_shift.clear();
  std::vector<double> pdfup_shift;
  pdfup_shift.clear();
  std::vector<double> pdfdown_shift;
  pdfdown_shift.clear();
  std::vector<double> xsec_shift_ZNuNuG;
  xsec_shift_ZNuNuG.clear();
  std::vector<double> xsec_shift_WG;
  xsec_shift_WG.clear();
  std::vector<double> xsec_shift_ZLLG;
  xsec_shift_ZLLG.clear();
  std::vector<double> syst_shiftUp_jetfake;
  syst_shiftUp_jetfake.clear();
  std::vector<double> syst_shiftDown_jetfake;
  syst_shiftDown_jetfake.clear();
  std::vector<double> syst_shiftUp_elefake;
  syst_shiftUp_elefake.clear();
  std::vector<double> syst_shiftDown_elefake;
  syst_shiftDown_elefake.clear();
  for(int i = 1; i <= nBins; i++){
    jesup_shift.push_back(0);
    jesdown_shift.push_back(0);
    pesup_shift.push_back(0);
    pesdown_shift.push_back(0);
    renup_shift.push_back(0);
    rendown_shift.push_back(0);
    facup_shift.push_back(0);
    facdown_shift.push_back(0);
    pdfup_shift.push_back(0);
    pdfdown_shift.push_back(0);
    xsec_shift_ZNuNuG.push_back(0);
    xsec_shift_WG.push_back(0);
    xsec_shift_ZLLG.push_back(0);
    syst_shiftUp_jetfake.push_back(0);
    syst_shiftDown_jetfake.push_back(0);
    syst_shiftUp_elefake.push_back(0);
    syst_shiftDown_elefake.push_back(0);
  }
  
  TFile *phoSFfile = new TFile("scaleFactor_ptalt.root");
  TH1F* phoSF = (TH1F*)phoSFfile->Get("sf_truth");
  TH1F* phoSF_up = (TH1F*)phoSF->Clone("phoSF_up");
  TH1F* phoSF_down = (TH1F*)phoSF->Clone("phoSF_down");
  for(int i = 1; i <= nBins; i++){
    phoSF_up->SetBinContent(i, phoSF->GetBinContent(i) + phoSF->GetBinError(i));
    phoSF_down->SetBinContent(i, phoSF->GetBinContent(i) - phoSF->GetBinError(i));
  }
  
  TFile *f_jetfake = new TFile("ZnnG_qcd_below0p5_all.root");
  TH1F* histo_jetfake = (TH1F*)((TH1F*)f_jetfake->Get(histname_data))->Clone("histo_jetfake");
  TH1F* histo_jetfake_sidebandUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sidebandUp))->Clone("histo_jetfake_sidebandUp");
  TH1F* histo_jetfake_sidebandDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sidebandDown))->Clone("histo_jetfake_sidebandDown");
  TH1F* histo_jetfake_METUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_METUp))->Clone("histo_jetfake_METUp");
  TH1F* histo_jetfake_METDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_METDown))->Clone("histo_jetfake_METDown");
  TH1F* histo_jetfake_binningUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_binningUp))->Clone("histo_jetfake_binningUp");
  TH1F* histo_jetfake_binningDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_binningDown))->Clone("histo_jetfake_binningDown");
  TH1F* histo_jetfake_sieieLeft = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sieieLeft))->Clone("histo_jetfake_sieieLeft");
  TH1F* histo_jetfake_sieieRight = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sieieRight))->Clone("histo_jetfake_sieieRight");
  TH1F* histo_jetfake_templateUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_templateUp))->Clone("histo_jetfake_templateUp");
  TH1F* histo_jetfake_templateDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_templateDown))->Clone("histo_jetfake_templateDown");
  TH1F* histo_jetfake_unweighted = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_unweighted))->Clone("histo_jetfake_unweighted");
  // Scale to appropriate lumi until we unblind
  histo_jetfake->Scale(int_lumi/12950.0);
  histo_jetfake_sidebandUp->Scale(int_lumi/12950.0);
  histo_jetfake_sidebandDown->Scale(int_lumi/12950.0);
  histo_jetfake_METUp->Scale(int_lumi/12950.0);
  histo_jetfake_METDown->Scale(int_lumi/12950.0);
  histo_jetfake_binningUp->Scale(int_lumi/12950.0);
  histo_jetfake_binningDown->Scale(int_lumi/12950.0);
  histo_jetfake_sieieLeft->Scale(int_lumi/12950.0);
  histo_jetfake_sieieRight->Scale(int_lumi/12950.0);
  histo_jetfake_templateUp->Scale(int_lumi/12950.0);
  histo_jetfake_templateDown->Scale(int_lumi/12950.0);
  histo_jetfake_unweighted->Scale(int_lumi/12950.0);
  TH1F* histo_jetfake_errUp = (TH1F*)histo_jetfake->Clone("histo_jetfake_errUp");
  TH1F* histo_jetfake_errDown = (TH1F*)histo_jetfake->Clone("histo_jetfake_errDown");
  Float_t int_jetfake = histo_jetfake->Integral()+histo_jetfake->GetBinContent(0)+histo_jetfake->GetBinContent(nBins+1);
  Float_t max_int_jetfake = 0.0;
  Float_t min_int_jetfake = 0.0;
  Float_t stat_jetfake = 0.0;
  for(int i = 1; i <= nBins; i++){
    double int_bin_jetfake = histo_jetfake->GetBinContent(i);
    double int_bin_jetfake_sidebandUp = histo_jetfake_sidebandUp->GetBinContent(i);
    double int_bin_jetfake_sidebandDown = histo_jetfake_sidebandDown->GetBinContent(i);
    double int_bin_jetfake_METUp = histo_jetfake_METUp->GetBinContent(i);
    double int_bin_jetfake_METDown = histo_jetfake_METDown->GetBinContent(i);
    double int_bin_jetfake_binningUp = histo_jetfake_binningUp->GetBinContent(i);
    double int_bin_jetfake_binningDown = histo_jetfake_binningDown->GetBinContent(i);
    double int_bin_jetfake_sieieLeft = histo_jetfake_sieieLeft->GetBinContent(i);
    double int_bin_jetfake_sieieRight = histo_jetfake_sieieRight->GetBinContent(i);
    double int_bin_jetfake_templateUp = histo_jetfake_templateUp->GetBinContent(i);
    double int_bin_jetfake_templateDown = histo_jetfake_templateDown->GetBinContent(i);
    double int_bin_jetfake_unweighted = histo_jetfake_unweighted->GetBinContent(i);
    double ints_bin[] = {int_bin_jetfake, int_bin_jetfake_sidebandUp, int_bin_jetfake_sidebandDown, int_bin_jetfake_METUp, int_bin_jetfake_METDown, int_bin_jetfake_binningUp, int_bin_jetfake_binningDown, int_bin_jetfake_sieieLeft, int_bin_jetfake_sieieRight, int_bin_jetfake_templateUp, int_bin_jetfake_templateDown};
    double max_int_bin = *max_element(ints_bin, ints_bin+11);
    double min_int_bin = *min_element(ints_bin, ints_bin+11);
    histo_jetfake_errUp->SetBinContent(i, max_int_bin);
    histo_jetfake_errDown->SetBinContent(i, min_int_bin);
    max_int_jetfake += max_int_bin;
    min_int_jetfake += min_int_bin;
    syst_shiftUp_jetfake[i-1] = max_int_jetfake-int_bin_jetfake;
    syst_shiftDown_jetfake[i-1] = min_int_jetfake-int_bin_jetfake;
    double stat_bin_jetfake = 0.0;
    if (int_bin_jetfake_unweighted > 0)
      stat_bin_jetfake = int_bin_jetfake/sqrt(int_bin_jetfake_unweighted);
    histo_jetfake->SetBinError(i, stat_bin_jetfake);
    stat_jetfake += stat_bin_jetfake*stat_bin_jetfake;
  }
  Float_t syst_jetfake = TMath::Max(max_int_jetfake-int_jetfake, int_jetfake-min_int_jetfake);
  stat_jetfake = sqrt(stat_jetfake);
  Float_t err_jetfake = sqrt(syst_jetfake*syst_jetfake + stat_jetfake*stat_jetfake);
  total_background += int_jetfake;
  background_unc_sumsquares += err_jetfake*err_jetfake;
  histo_jetfake->SetFillColor(TColor::GetColor("#FFFFCC"));
  histo_vector.push_back(histo_jetfake);
  
  //DEBUG
  // cout<<"Jetfakes got"<<endl;
  
  TFile *f_spikes = new TFile("SpikesTemplates.root");
  TH1F* histo_spikes = (TH1F*)f_spikes->Get("sscutPt");
  if(histname == "pfMET_3")
    histo_spikes = (TH1F*)f_spikes->Get("sscutMet");
  else if(histname == "Mt_3")
    histo_spikes = (TH1F*)f_spikes->Get("sscutMT");
  histo_spikes = (TH1F*)histo_spikes->Clone("histo_spikes"); // Change the name
  Float_t int_spikes = histo_spikes->Integral()+histo_spikes->GetBinContent(0)+histo_spikes->GetBinContent(nBins+1);
  histo_spikes->Scale(23.90*frac_below0p5/int_spikes);
  int_spikes = histo_spikes->Integral()+histo_spikes->GetBinContent(0)+histo_spikes->GetBinContent(nBins+1);
  Float_t err_spikes = 1.00*int_spikes;
  for(int i = 1; i <= nBins; i++){
    double int_bin_spikes = histo_spikes->GetBinContent(i);
    histo_spikes->SetBinError(i,1.00*int_bin_spikes);
  }
  total_background += int_spikes;
  background_unc_sumsquares += err_spikes*err_spikes;
  histo_spikes->SetFillColor(kOrange+10);
  histo_vector.push_back(histo_spikes);
  
  TFile *f_elefake = new TFile("ZnnG_wenu_below0p5_all.root");
  TH1F* histo_elefake = (TH1F*)((TH1F*)f_elefake->Get(histname_data))->Clone("histo_elefake");
  TH1F* histo_elefake_errUp = (TH1F*)((TH1F*)f_elefake->Get(histname_wenu_errUp))->Clone("histo_elefake_errUp");
  TH1F* histo_elefake_errDown = (TH1F*)((TH1F*)f_elefake->Get(histname_wenu_errDown))->Clone("histo_elefake_errDown");
  TH1F* histo_elefake_unweighted = (TH1F*)((TH1F*)f_elefake->Get(histname_wenu_unweighted))->Clone("histo_elefake_unweighted");
  // Scale to appropriate luminosity until we unblind
  histo_elefake->Scale(int_lumi/12950.0);
  histo_elefake_errUp->Scale(int_lumi/12950.0);
  histo_elefake_errDown->Scale(int_lumi/12950.0);
  histo_elefake_unweighted->Scale(int_lumi/12950.0);
  Float_t int_elefake = histo_elefake->Integral()+histo_elefake->GetBinContent(0)+histo_elefake->GetBinContent(nBins+1);
  Float_t int_elefake_errUp = histo_elefake_errUp->Integral()+histo_elefake_errUp->GetBinContent(0)+histo_elefake_errUp->GetBinContent(nBins+1);
  Float_t int_elefake_errDown = histo_elefake_errDown->Integral()+histo_elefake_errDown->GetBinContent(0)+histo_elefake_errDown->GetBinContent(nBins+1);
  Float_t max_int_elefake = 0.0;
  Float_t min_int_elefake = 0.0;
  Float_t stat_elefake = 0.0;
  for(int i = 1; i <= nBins; i++){
    double int_bin_elefake = histo_elefake->GetBinContent(i);
    double int_bin_elefake_errUp = histo_elefake_errUp->GetBinContent(i);
    double int_bin_elefake_errDown = histo_elefake_errDown->GetBinContent(i);
    double int_bin_unweighted = histo_elefake_unweighted->GetBinContent(i);
    syst_shiftUp_elefake[i-1] = int_bin_elefake_errUp-int_bin_elefake;
    syst_shiftDown_elefake[i-1] = int_bin_elefake_errDown-int_bin_elefake;
    max_int_elefake += int_bin_elefake_errUp;
    min_int_elefake += int_bin_elefake_errDown;
    double stat_bin_elefake = 0.0;
    if (int_bin_unweighted > 0)
      stat_bin_elefake = int_bin_elefake/sqrt(int_bin_unweighted);
    histo_elefake->SetBinError(i,stat_bin_elefake);
    stat_elefake += stat_bin_elefake*stat_bin_elefake;
  }
  Float_t syst_elefake = TMath::Max(max_int_elefake-int_elefake, int_elefake-min_int_elefake);
  stat_elefake = sqrt(stat_elefake);
  Float_t err_elefake = sqrt(syst_elefake*syst_elefake + stat_elefake*stat_elefake);
  total_background += int_elefake;
  background_unc_sumsquares += err_elefake*err_elefake;
  // histo_elefake->SetFillColor(kBlue-7);
  histo_elefake->SetFillColor(kBlue-8);
  histo_vector.push_back(histo_elefake);
  
  //DEBUG
  // cout<<"elefakes got"<<endl;
  
  TFile *f_bhalo = new TFile("ZnnG_bhalo_below0p5_all.root");
  TH1F* histo_bhalo = (TH1F*)((TH1F*)f_bhalo->Get(histname_data))->Clone("histo_bhalo");
  // Scale to appropriate luminosity until we unblind
  histo_bhalo->Scale(int_lumi/12950.0);
  Float_t int_bhalo = histo_bhalo->Integral()+histo_bhalo->GetBinContent(0)+histo_bhalo->GetBinContent(nBins+1);
  histo_bhalo->Scale(7.3*frac_below0p5/int_bhalo);
  int_bhalo = histo_bhalo->Integral()+histo_bhalo->GetBinContent(0)+histo_bhalo->GetBinContent(nBins+1);
  Float_t err_bhalo = (4.7/7.3)*int_bhalo;
  for(int i = 1; i <= nBins; i++){
    double int_bin_bhalo = histo_bhalo->GetBinContent(i);
    histo_bhalo->SetBinError(i,(4.7/7.3)*int_bin_bhalo);
  }
  total_background += int_bhalo;
  background_unc_sumsquares += err_bhalo*err_bhalo;
  histo_bhalo->SetFillColor(12);
  // histo_bhalo->SetFillColor(kSpring-6);
  // histo_bhalo->SetFillColor(kViolet-1);
  // histo_bhalo->SetFillColor(kAzure-9);
  histo_vector.push_back(histo_bhalo);
  
  //DEBUG
  // cout<<"beam halo got"<<endl;
  
  TFile* f_ZNuNuG = new TFile("ZnnG_JESPES_ZNuNuGJets.root");
  TH1F* histo_ZNuNuG = (TH1F*)((TH1F*)f_ZNuNuG->Get(histname))->Clone("histo_ZNuNuG");
  TH1F* histo_ZNuNuG_JESUp = (TH1F*)((TH1F*)f_ZNuNuG->Get(histname_JESUp))->Clone("histo_ZNuNuG_JESUp");
  TH1F* histo_ZNuNuG_JESDown = (TH1F*)((TH1F*)f_ZNuNuG->Get(histname_JESDown))->Clone("histo_ZNuNuG_JESDown");
  TH1F* histo_ZNuNuG_PESUp = (TH1F*)((TH1F*)f_ZNuNuG->Get(histname_PESUp))->Clone("histo_ZNuNuG_PESUp");
  TH1F* histo_ZNuNuG_PESDown = (TH1F*)((TH1F*)f_ZNuNuG->Get(histname_PESDown))->Clone("histo_ZNuNuG_PESDown");
  TH1F* histo_ZNuNuG_uncorrected = (TH1F*)((TH1F*)f_ZNuNuG->Get(histname_uncorrected))->Clone("histo_ZNuNuG_uncorrected");
  TFile* f_ZNuNuG_ext = new TFile("ZnnG_JESPES_ZNuNuGJets_ext.root");
  TH1F* histo_ZNuNuG_ext = (TH1F*)((TH1F*)f_ZNuNuG_ext->Get(histname))->Clone("histo_ZNuNuG_ext");
  TH1F* histo_ZNuNuG_JESUp_ext = (TH1F*)((TH1F*)f_ZNuNuG_ext->Get(histname_JESUp))->Clone("histo_ZNuNuG_JESUp_ext");
  TH1F* histo_ZNuNuG_JESDown_ext = (TH1F*)((TH1F*)f_ZNuNuG_ext->Get(histname_JESDown))->Clone("histo_ZNuNuG_JESDown_ext");
  TH1F* histo_ZNuNuG_PESUp_ext = (TH1F*)((TH1F*)f_ZNuNuG_ext->Get(histname_PESUp))->Clone("histo_ZNuNuG_PESUp_ext");
  TH1F* histo_ZNuNuG_PESDown_ext = (TH1F*)((TH1F*)f_ZNuNuG_ext->Get(histname_PESDown))->Clone("histo_ZNuNuG_PESDown_ext");
  TH1F* histo_ZNuNuG_uncorrected_ext = (TH1F*)((TH1F*)f_ZNuNuG_ext->Get(histname_uncorrected))->Clone("histo_ZNuNuG_uncorrected_ext");
  histo_ZNuNuG->SetStats(0);
  histo_ZNuNuG->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0); // 329894(standard) + 1710471(ext) = 2040365
  histo_ZNuNuG_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG_uncorrected->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0); // 329894(standard) + 1710471(ext) = 2040365
  histo_ZNuNuG_JESUp_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG_JESDown_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG_PESUp_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG_PESDown_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG_uncorrected_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.18213/2040365.0);
  histo_ZNuNuG->Add(histo_ZNuNuG_ext);
  histo_ZNuNuG_JESUp->Add(histo_ZNuNuG_JESUp_ext);
  histo_ZNuNuG_JESDown->Add(histo_ZNuNuG_JESDown_ext);
  histo_ZNuNuG_PESUp->Add(histo_ZNuNuG_PESUp_ext);
  histo_ZNuNuG_PESDown->Add(histo_ZNuNuG_PESDown_ext);
  histo_ZNuNuG_uncorrected->Add(histo_ZNuNuG_uncorrected_ext);
  Float_t int_ZNuNuG = histo_ZNuNuG->Integral()+histo_ZNuNuG->GetBinContent(0)+histo_ZNuNuG->GetBinContent(nBins+1);
  Float_t int_ZNuNuG_JESUp = histo_ZNuNuG_JESUp->Integral()+histo_ZNuNuG_JESUp->GetBinContent(0)+histo_ZNuNuG_JESUp->GetBinContent(nBins+1);
  Float_t int_ZNuNuG_JESDown = histo_ZNuNuG_JESDown->Integral()+histo_ZNuNuG_JESDown->GetBinContent(0)+histo_ZNuNuG_JESDown->GetBinContent(nBins+1);
  Float_t int_ZNuNuG_PESUp = histo_ZNuNuG_PESUp->Integral()+histo_ZNuNuG_PESUp->GetBinContent(0)+histo_ZNuNuG_PESUp->GetBinContent(nBins+1);
  Float_t int_ZNuNuG_PESDown = histo_ZNuNuG_PESDown->Integral()+histo_ZNuNuG_PESDown->GetBinContent(0)+histo_ZNuNuG_PESDown->GetBinContent(nBins+1);
  Float_t int_ZNuNuG_uncorrected = histo_ZNuNuG_uncorrected->Integral()+histo_ZNuNuG_uncorrected->GetBinContent(0)+histo_ZNuNuG_uncorrected->GetBinContent(nBins+1);
  double jeserr_ZNuNuG = (fabs(int_ZNuNuG_JESUp-int_ZNuNuG)+fabs(int_ZNuNuG_JESDown-int_ZNuNuG))/2.0;
  double peserr_ZNuNuG = (fabs(int_ZNuNuG_PESUp-int_ZNuNuG)+fabs(int_ZNuNuG_PESDown-int_ZNuNuG))/2.0;
  double xsecerr_ZNuNuG = fabs(int_ZNuNuG_uncorrected-int_ZNuNuG);
  TH1F* histo_ZNuNuG_RenUp;
  TH1F* histo_ZNuNuG_RenDown;
  TH1F* histo_ZNuNuG_FacUp;
  TH1F* histo_ZNuNuG_FacDown;
  TH1F* histo_ZNuNuG_PDFUp;
  TH1F* histo_ZNuNuG_PDFDown;
  TH1F* histo_ZNuNuG_RenUp_ext;
  TH1F* histo_ZNuNuG_RenDown_ext;
  TH1F* histo_ZNuNuG_FacUp_ext;
  TH1F* histo_ZNuNuG_FacDown_ext;
  TH1F* histo_ZNuNuG_PDFUp_ext;
  TH1F* histo_ZNuNuG_PDFDown_ext;
  double renerr_ZNuNuG = 0.0;
  double facerr_ZNuNuG = 0.0;
  double pdferr_ZNuNuG = 0.0;
  if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3"){
    string xaxis_variable = "XAXIS_VARIABLE_NOT_SET";
    if (histname == "Photon_Et_range_3") xaxis_variable = "Pt";
    else if (histname == "pfMET_3") xaxis_variable = "MET";
    else if (histname == "Mt_3") xaxis_variable = "Mt";
    TString filename = TString("histos_"+xaxis_variable+"_ZnnG_pdfscale_ZNuNuGJets.root");
    TString filename_ext = TString("histos_"+xaxis_variable+"_ZnnG_pdfscale_ZNuNuGJets_ext.root");
    TFile* f_pdfscale_ZNuNuG = new TFile(filename);
    histo_ZNuNuG_RenUp = (TH1F*)f_pdfscale_ZNuNuG->Get("h_ZnnG_pdfscale_ZNuNuGJets_renUp");
    histo_ZNuNuG_RenDown = (TH1F*)f_pdfscale_ZNuNuG->Get("h_ZnnG_pdfscale_ZNuNuGJets_renDown");
    histo_ZNuNuG_FacUp = (TH1F*)f_pdfscale_ZNuNuG->Get("h_ZnnG_pdfscale_ZNuNuGJets_facUp");
    histo_ZNuNuG_FacDown = (TH1F*)f_pdfscale_ZNuNuG->Get("h_ZnnG_pdfscale_ZNuNuGJets_facDown");
    histo_ZNuNuG_PDFUp = (TH1F*)f_pdfscale_ZNuNuG->Get("h_ZnnG_pdfscale_ZNuNuGJets_pdfUp");
    histo_ZNuNuG_PDFDown = (TH1F*)f_pdfscale_ZNuNuG->Get("h_ZnnG_pdfscale_ZNuNuGJets_pdfDown");
    TFile* f_pdfscale_ZNuNuG_ext = new TFile(filename_ext);
    histo_ZNuNuG_RenUp_ext = (TH1F*)f_pdfscale_ZNuNuG_ext->Get("h_ZnnG_pdfscale_ZNuNuGJets_ext_renUp");
    histo_ZNuNuG_RenDown_ext = (TH1F*)f_pdfscale_ZNuNuG_ext->Get("h_ZnnG_pdfscale_ZNuNuGJets_ext_renDown");
    histo_ZNuNuG_FacUp_ext = (TH1F*)f_pdfscale_ZNuNuG_ext->Get("h_ZnnG_pdfscale_ZNuNuGJets_ext_facUp");
    histo_ZNuNuG_FacDown_ext = (TH1F*)f_pdfscale_ZNuNuG_ext->Get("h_ZnnG_pdfscale_ZNuNuGJets_ext_facDown");
    histo_ZNuNuG_PDFUp_ext = (TH1F*)f_pdfscale_ZNuNuG_ext->Get("h_ZnnG_pdfscale_ZNuNuGJets_ext_pdfUp");
    histo_ZNuNuG_PDFDown_ext = (TH1F*)f_pdfscale_ZNuNuG_ext->Get("h_ZnnG_pdfscale_ZNuNuGJets_ext_pdfDown");
    histo_ZNuNuG_RenUp->Add(histo_ZNuNuG_RenUp_ext);
    histo_ZNuNuG_RenDown->Add(histo_ZNuNuG_RenDown_ext);
    histo_ZNuNuG_FacUp->Add(histo_ZNuNuG_FacUp_ext);
    histo_ZNuNuG_FacDown->Add(histo_ZNuNuG_FacDown_ext);
    histo_ZNuNuG_PDFUp->Add(histo_ZNuNuG_PDFUp_ext);
    histo_ZNuNuG_PDFDown->Add(histo_ZNuNuG_PDFDown_ext);
    histo_ZNuNuG_RenUp->Scale(frac_below0p5);
    histo_ZNuNuG_RenDown->Scale(frac_below0p5);
    histo_ZNuNuG_FacUp->Scale(frac_below0p5);
    histo_ZNuNuG_FacDown->Scale(frac_below0p5);
    histo_ZNuNuG_PDFUp->Scale(frac_below0p5);
    histo_ZNuNuG_PDFDown->Scale(frac_below0p5);
    Float_t int_ZNuNuG_RenUp = histo_ZNuNuG_RenUp->Integral()+histo_ZNuNuG_RenUp->GetBinContent(0)+histo_ZNuNuG_RenUp->GetBinContent(nBins+1);
    Float_t int_ZNuNuG_RenDown = histo_ZNuNuG_RenDown->Integral()+histo_ZNuNuG_RenDown->GetBinContent(0)+histo_ZNuNuG_RenDown->GetBinContent(nBins+1);
    Float_t int_ZNuNuG_FacUp = histo_ZNuNuG_FacUp->Integral()+histo_ZNuNuG_FacUp->GetBinContent(0)+histo_ZNuNuG_FacUp->GetBinContent(nBins+1);
    Float_t int_ZNuNuG_FacDown = histo_ZNuNuG_FacDown->Integral()+histo_ZNuNuG_FacDown->GetBinContent(0)+histo_ZNuNuG_FacDown->GetBinContent(nBins+1);
    Float_t int_ZNuNuG_PDFUp = histo_ZNuNuG_PDFUp->Integral()+histo_ZNuNuG_PDFUp->GetBinContent(0)+histo_ZNuNuG_PDFUp->GetBinContent(nBins+1);
    Float_t int_ZNuNuG_PDFDown = histo_ZNuNuG_PDFDown->Integral()+histo_ZNuNuG_PDFDown->GetBinContent(0)+histo_ZNuNuG_PDFDown->GetBinContent(nBins+1);
    renerr_ZNuNuG = (fabs(int_ZNuNuG_RenUp-int_ZNuNuG)+fabs(int_ZNuNuG_RenDown-int_ZNuNuG))/2.0;
    facerr_ZNuNuG = (fabs(int_ZNuNuG_FacUp-int_ZNuNuG)+fabs(int_ZNuNuG_FacDown-int_ZNuNuG))/2.0;
    pdferr_ZNuNuG = (fabs(int_ZNuNuG_PDFUp-int_ZNuNuG)+fabs(int_ZNuNuG_PDFDown-int_ZNuNuG))/2.0;
  }
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_ZNuNuG->GetBinContent(i);
    double jesup = histo_ZNuNuG_JESUp->GetBinContent(i);
    double jesdown = histo_ZNuNuG_JESDown->GetBinContent(i);
    double pesup = histo_ZNuNuG_PESUp->GetBinContent(i);
    double pesdown = histo_ZNuNuG_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    double uncorrected = histo_ZNuNuG_uncorrected->GetBinContent(i);
    xsec_shift_ZNuNuG[i-1] = fabs(uncorrected-int_bin);
    if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3"){
      double renup = histo_ZNuNuG_RenUp->GetBinContent(i);
      double rendown = histo_ZNuNuG_RenDown->GetBinContent(i);
      double facup = histo_ZNuNuG_FacUp->GetBinContent(i);
      double facdown = histo_ZNuNuG_FacDown->GetBinContent(i);
      double pdfup = histo_ZNuNuG_PDFUp->GetBinContent(i);
      double pdfdown = histo_ZNuNuG_PDFDown->GetBinContent(i);
      renup_shift[i-1] += renup-int_bin;
      rendown_shift[i-1] += rendown-int_bin;
      facup_shift[i-1] += facup-int_bin;
      facdown_shift[i-1] += facdown-int_bin;
      pdfup_shift[i-1] += pdfup-int_bin;
      pdfdown_shift[i-1] += pdfdown-int_bin;
    }
    // cout<<"ZNuNuG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((0.18213*int_lumi*scale_factor*frac_below0p5-int_bin)/(2040365.0*int_bin));
    histo_ZNuNuG->SetBinError(i,err_bin);
  }
  Float_t err_ZNuNuG = 0.0;
  if(int_ZNuNuG > 0.0){
    err_ZNuNuG = sqrt(int_ZNuNuG*int_ZNuNuG*((photon_scale_factor_unc*photon_scale_factor_unc)+(0.18213*int_lumi*scale_factor*frac_below0p5-int_ZNuNuG)/(2040365.0*int_ZNuNuG))+(xsecerr_ZNuNuG*xsecerr_ZNuNuG)+(jeserr_ZNuNuG*jeserr_ZNuNuG)+(peserr_ZNuNuG*peserr_ZNuNuG)+(renerr_ZNuNuG*renerr_ZNuNuG)+(facerr_ZNuNuG*facerr_ZNuNuG)+(pdferr_ZNuNuG*pdferr_ZNuNuG));
    if(histname == "Photon_Et_range_3"){
      cout<<"ZnnG:"<<int_ZNuNuG<<endl;
      cout<<"stat err = "<<sqrt(((0.18213*int_lumi*scale_factor*frac_below0p5-int_ZNuNuG)/(2040365.0*int_ZNuNuG)))<<endl;
      cout<<"jes err = "<<jeserr_ZNuNuG/int_ZNuNuG<<endl;
      cout<<"pes err = "<<peserr_ZNuNuG/int_ZNuNuG<<endl;
      cout<<"ren err = "<<renerr_ZNuNuG/int_ZNuNuG<<endl;
      cout<<"fac err = "<<facerr_ZNuNuG/int_ZNuNuG<<endl;
      cout<<"pdf err = "<<pdferr_ZNuNuG/int_ZNuNuG<<endl;
      cout<<"pho scale factor err = "<<photon_scale_factor_unc<<endl;
      cout<<"NNLO xsec err = "<<xsecerr_ZNuNuG/int_ZNuNuG<<endl;
      cout<<"Total systematic = "<<sqrt(pow(jeserr_ZNuNuG,2)+pow(peserr_ZNuNuG,2)+pow(renerr_ZNuNuG,2)+pow(facerr_ZNuNuG,2)+pow(pdferr_ZNuNuG,2)+pow(photon_scale_factor_unc*int_ZNuNuG,2)+pow(xsecerr_ZNuNuG,2))/int_ZNuNuG<<endl;
    }
  }
  total_background += int_ZNuNuG;
  background_unc_sumsquares += err_ZNuNuG*err_ZNuNuG;
  // histo_ZNuNuG->SetFillColor(kOrange);
  histo_ZNuNuG->SetFillColor(kOrange-4);
  histo_vector.push_back(histo_ZNuNuG);
  
  //DEBUG
  // cout<<"ZNuNuG got"<<endl;
  
  TFile* f_WG = new TFile("ZnnG_JESPES_WGJets.root");
  TH1F* histo_WG = (TH1F*)((TH1F*)f_WG->Get(histname))->Clone("histo_WG");
  TH1F* histo_WG_JESUp = (TH1F*)((TH1F*)f_WG->Get(histname_JESUp))->Clone("histo_WG_JESUp");
  TH1F* histo_WG_JESDown = (TH1F*)((TH1F*)f_WG->Get(histname_JESDown))->Clone("histo_WG_JESDown");
  TH1F* histo_WG_PESUp = (TH1F*)((TH1F*)f_WG->Get(histname_PESUp))->Clone("histo_WG_PESUp");
  TH1F* histo_WG_PESDown = (TH1F*)((TH1F*)f_WG->Get(histname_PESDown))->Clone("histo_WG_PESDown");
  TH1F* histo_WG_uncorrected = (TH1F*)((TH1F*)f_WG->Get(histname_uncorrected))->Clone("histo_WG_uncorrected");
  TFile* f_WG_ext = new TFile("ZnnG_JESPES_WGJets_ext.root");
  TH1F* histo_WG_ext = (TH1F*)((TH1F*)f_WG_ext->Get(histname))->Clone("histo_WG_ext");
  TH1F* histo_WG_JESUp_ext = (TH1F*)((TH1F*)f_WG_ext->Get(histname_JESUp))->Clone("histo_WG_JESUp_ext");
  TH1F* histo_WG_JESDown_ext = (TH1F*)((TH1F*)f_WG_ext->Get(histname_JESDown))->Clone("histo_WG_JESDown_ext");
  TH1F* histo_WG_PESUp_ext = (TH1F*)((TH1F*)f_WG_ext->Get(histname_PESUp))->Clone("histo_WG_PESUp_ext");
  TH1F* histo_WG_PESDown_ext = (TH1F*)((TH1F*)f_WG_ext->Get(histname_PESDown))->Clone("histo_WG_PESDown_ext");
  TH1F* histo_WG_uncorrected_ext = (TH1F*)((TH1F*)f_WG_ext->Get(histname_uncorrected))->Clone("histo_WG_uncorrected_ext");
  histo_WG->SetStats(0);
  histo_WG->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0); // 502176(standard) + 1852305(ext) = 2354481
  histo_WG_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG_uncorrected->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0); // 502176(standard) + 1852305(ext) = 2354481
  histo_WG_JESUp_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG_JESDown_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG_PESUp_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG_PESDown_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG_uncorrected_ext->Scale(int_lumi*scale_factor*frac_below0p5*0.65521/2354481.0);
  histo_WG->Add(histo_WG_ext);
  histo_WG_JESUp->Add(histo_WG_JESUp_ext);
  histo_WG_JESDown->Add(histo_WG_JESDown_ext);
  histo_WG_PESUp->Add(histo_WG_PESUp_ext);
  histo_WG_PESDown->Add(histo_WG_PESDown_ext);
  histo_WG_uncorrected->Add(histo_WG_uncorrected_ext);
  Float_t int_WG = histo_WG->Integral()+histo_WG->GetBinContent(0)+histo_WG->GetBinContent(nBins+1);
  Float_t int_WG_JESUp = histo_WG_JESUp->Integral()+histo_WG_JESUp->GetBinContent(0)+histo_WG_JESUp->GetBinContent(nBins+1);
  Float_t int_WG_JESDown = histo_WG_JESDown->Integral()+histo_WG_JESDown->GetBinContent(0)+histo_WG_JESDown->GetBinContent(nBins+1);
  Float_t int_WG_PESUp = histo_WG_PESUp->Integral()+histo_WG_PESUp->GetBinContent(0)+histo_WG_PESUp->GetBinContent(nBins+1);
  Float_t int_WG_PESDown = histo_WG_PESDown->Integral()+histo_WG_PESDown->GetBinContent(0)+histo_WG_PESDown->GetBinContent(nBins+1);
  Float_t int_WG_uncorrected = histo_WG_uncorrected->Integral()+histo_WG_uncorrected->GetBinContent(0)+histo_WG_uncorrected->GetBinContent(nBins+1);
  double jeserr_WG = (fabs(int_WG_JESUp-int_WG)+fabs(int_WG_JESDown-int_WG))/2.0;
  double peserr_WG = (fabs(int_WG_PESUp-int_WG)+fabs(int_WG_PESDown-int_WG))/2.0;
  double xsecerr_WG = fabs(int_WG_uncorrected-int_WG);
  TH1F* histo_WG_RenUp;
  TH1F* histo_WG_RenDown;
  TH1F* histo_WG_FacUp;
  TH1F* histo_WG_FacDown;
  TH1F* histo_WG_PDFUp;
  TH1F* histo_WG_PDFDown;
  TH1F* histo_WG_RenUp_ext;
  TH1F* histo_WG_RenDown_ext;
  TH1F* histo_WG_FacUp_ext;
  TH1F* histo_WG_FacDown_ext;
  TH1F* histo_WG_PDFUp_ext;
  TH1F* histo_WG_PDFDown_ext;
  double renerr_WG = 0.0;
  double facerr_WG = 0.0;
  double pdferr_WG = 0.0;
  if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3"){
    string xaxis_variable = "XAXIS_VARIABLE_NOT_SET";
    if (histname == "Photon_Et_range_3") xaxis_variable = "Pt";
    else if (histname == "pfMET_3") xaxis_variable = "MET";
    else if (histname == "Mt_3") xaxis_variable = "Mt";
    TString filename = TString("histos_"+xaxis_variable+"_ZnnG_pdfscale_WGJets.root");
    TString filename_ext = TString("histos_"+xaxis_variable+"_ZnnG_pdfscale_WGJets_ext.root");
    TFile* f_pdfscale_WG = new TFile(filename);
    histo_WG_RenUp = (TH1F*)f_pdfscale_WG->Get("h_ZnnG_pdfscale_WGJets_renUp");
    histo_WG_RenDown = (TH1F*)f_pdfscale_WG->Get("h_ZnnG_pdfscale_WGJets_renDown");
    histo_WG_FacUp = (TH1F*)f_pdfscale_WG->Get("h_ZnnG_pdfscale_WGJets_facUp");
    histo_WG_FacDown = (TH1F*)f_pdfscale_WG->Get("h_ZnnG_pdfscale_WGJets_facDown");
    histo_WG_PDFUp = (TH1F*)f_pdfscale_WG->Get("h_ZnnG_pdfscale_WGJets_pdfUp");
    histo_WG_PDFDown = (TH1F*)f_pdfscale_WG->Get("h_ZnnG_pdfscale_WGJets_pdfDown");
    TFile* f_pdfscale_WG_ext = new TFile(filename_ext);
    histo_WG_RenUp_ext = (TH1F*)f_pdfscale_WG_ext->Get("h_ZnnG_pdfscale_WGJets_ext_renUp");
    histo_WG_RenDown_ext = (TH1F*)f_pdfscale_WG_ext->Get("h_ZnnG_pdfscale_WGJets_ext_renDown");
    histo_WG_FacUp_ext = (TH1F*)f_pdfscale_WG_ext->Get("h_ZnnG_pdfscale_WGJets_ext_facUp");
    histo_WG_FacDown_ext = (TH1F*)f_pdfscale_WG_ext->Get("h_ZnnG_pdfscale_WGJets_ext_facDown");
    histo_WG_PDFUp_ext = (TH1F*)f_pdfscale_WG_ext->Get("h_ZnnG_pdfscale_WGJets_ext_pdfUp");
    histo_WG_PDFDown_ext = (TH1F*)f_pdfscale_WG_ext->Get("h_ZnnG_pdfscale_WGJets_ext_pdfDown");
    histo_WG_RenUp->Add(histo_WG_RenUp_ext);
    histo_WG_RenDown->Add(histo_WG_RenDown_ext);
    histo_WG_FacUp->Add(histo_WG_FacUp_ext);
    histo_WG_FacDown->Add(histo_WG_FacDown_ext);
    histo_WG_PDFUp->Add(histo_WG_PDFUp_ext);
    histo_WG_PDFDown->Add(histo_WG_PDFDown_ext);
    histo_WG_RenUp->Scale(frac_below0p5);
    histo_WG_RenDown->Scale(frac_below0p5);
    histo_WG_FacUp->Scale(frac_below0p5);
    histo_WG_FacDown->Scale(frac_below0p5);
    histo_WG_PDFUp->Scale(frac_below0p5);
    histo_WG_PDFDown->Scale(frac_below0p5);
    Float_t int_WG_RenUp = histo_WG_RenUp->Integral()+histo_WG_RenUp->GetBinContent(0)+histo_WG_RenUp->GetBinContent(nBins+1);
    Float_t int_WG_RenDown = histo_WG_RenDown->Integral()+histo_WG_RenDown->GetBinContent(0)+histo_WG_RenDown->GetBinContent(nBins+1);
    Float_t int_WG_FacUp = histo_WG_FacUp->Integral()+histo_WG_FacUp->GetBinContent(0)+histo_WG_FacUp->GetBinContent(nBins+1);
    Float_t int_WG_FacDown = histo_WG_FacDown->Integral()+histo_WG_FacDown->GetBinContent(0)+histo_WG_FacDown->GetBinContent(nBins+1);
    Float_t int_WG_PDFUp = histo_WG_PDFUp->Integral()+histo_WG_PDFUp->GetBinContent(0)+histo_WG_PDFUp->GetBinContent(nBins+1);
    Float_t int_WG_PDFDown = histo_WG_PDFDown->Integral()+histo_WG_PDFDown->GetBinContent(0)+histo_WG_PDFDown->GetBinContent(nBins+1);
    renerr_WG = (fabs(int_WG_RenUp-int_WG)+fabs(int_WG_RenDown-int_WG))/2.0;
    facerr_WG = (fabs(int_WG_FacUp-int_WG)+fabs(int_WG_FacDown-int_WG))/2.0;
    pdferr_WG = (fabs(int_WG_PDFUp-int_WG)+fabs(int_WG_PDFDown-int_WG))/2.0;
  }
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_WG->GetBinContent(i);
    double jesup = histo_WG_JESUp->GetBinContent(i);
    double jesdown = histo_WG_JESDown->GetBinContent(i);
    double pesup = histo_WG_PESUp->GetBinContent(i);
    double pesdown = histo_WG_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    double uncorrected = histo_WG_uncorrected->GetBinContent(i);
    xsec_shift_WG[i-1] = fabs(uncorrected-int_bin);
    if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3"){
      double renup = histo_WG_RenUp->GetBinContent(i);
      double rendown = histo_WG_RenDown->GetBinContent(i);
      double facup = histo_WG_FacUp->GetBinContent(i);
      double facdown = histo_WG_FacDown->GetBinContent(i);
      double pdfup = histo_WG_PDFUp->GetBinContent(i);
      double pdfdown = histo_WG_PDFDown->GetBinContent(i);
      renup_shift[i-1] += renup-int_bin;
      rendown_shift[i-1] += rendown-int_bin;
      facup_shift[i-1] += facup-int_bin;
      facdown_shift[i-1] += facdown-int_bin;
      pdfup_shift[i-1] += pdfup-int_bin;
      pdfdown_shift[i-1] += pdfdown-int_bin;
    }
    // cout<<"WG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((0.65521*int_lumi*scale_factor*frac_below0p5-int_bin)/(2354481.0*int_bin));
    histo_WG->SetBinError(i,err_bin);
  }
  Float_t err_WG = 0.0;
  if(int_WG > 0.0){
    err_WG = sqrt(int_WG*int_WG*((photon_scale_factor_unc*photon_scale_factor_unc)+(0.65521*int_lumi*scale_factor*frac_below0p5-int_WG)/(2354481.0*int_WG))+(xsecerr_WG*xsecerr_WG)+(jeserr_WG*jeserr_WG)+(peserr_WG*peserr_WG)+(renerr_WG*renerr_WG)+(facerr_WG*facerr_WG)+(pdferr_WG*pdferr_WG));
    if(histname == "Photon_Et_range_3"){
      cout<<"WG:"<<int_WG<<endl;
      cout<<"stat err = "<<sqrt(((0.65521*int_lumi*scale_factor*frac_below0p5-int_WG)/(2354481.0*int_WG)))<<endl;
      cout<<"jes err = "<<jeserr_WG/int_WG<<endl;
      cout<<"pes err = "<<peserr_WG/int_WG<<endl;
      cout<<"ren err = "<<renerr_WG/int_WG<<endl;
      cout<<"fac err = "<<facerr_WG/int_WG<<endl;
      cout<<"pdf err = "<<pdferr_WG/int_WG<<endl;
      cout<<"pho scale factor err = "<<photon_scale_factor_unc<<endl;
      cout<<"NNLO xsec err = "<<xsecerr_WG/int_WG<<endl;
      cout<<"Total systematic = "<<sqrt(pow(jeserr_WG,2)+pow(peserr_WG,2)+pow(renerr_WG,2)+pow(facerr_WG,2)+pow(pdferr_WG,2)+pow(photon_scale_factor_unc*int_WG,2)+pow(xsecerr_WG,2))/int_WG<<endl;
    }
  }
  total_background += int_WG;
  background_unc_sumsquares += err_WG*err_WG;
  // histo_WG->SetFillColor(kAzure+10);
  histo_WG->SetFillColor(kRed-6);
  histo_vector.push_back(histo_WG);
  
  //DEBUG
  // cout<<"WG got"<<endl;
  
  TFile *f_GJets_40to100 = new TFile(TString("ZnnG_JESPES_GJets_HT-40To100.root"));
  TFile *f_GJets_100to200 = new TFile(TString("ZnnG_JESPES_GJets_HT-100To200.root"));
  TFile *f_GJets_200to400 = new TFile(TString("ZnnG_JESPES_GJets_HT-200To400.root"));
  TFile *f_GJets_400to600 = new TFile(TString("ZnnG_JESPES_GJets_HT-400To600.root"));
  TFile *f_GJets_600toInf = new TFile(TString("ZnnG_JESPES_GJets_HT-600ToInf.root"));
  TH1F* histo_GJets_40to100 = (TH1F*)f_GJets_40to100->Get(histname);
  TH1F* histo_GJets_40to100_JESUp = (TH1F*)f_GJets_40to100->Get(histname_JESUp);
  TH1F* histo_GJets_40to100_JESDown = (TH1F*)f_GJets_40to100->Get(histname_JESDown);
  TH1F* histo_GJets_40to100_PESUp = (TH1F*)f_GJets_40to100->Get(histname_PESUp);
  TH1F* histo_GJets_40to100_PESDown = (TH1F*)f_GJets_40to100->Get(histname_PESDown);
  TH1F* histo_GJets_100to200 = (TH1F*)f_GJets_100to200->Get(histname);
  TH1F* histo_GJets_100to200_JESUp = (TH1F*)f_GJets_100to200->Get(histname_JESUp);
  TH1F* histo_GJets_100to200_JESDown = (TH1F*)f_GJets_100to200->Get(histname_JESDown);
  TH1F* histo_GJets_100to200_PESUp = (TH1F*)f_GJets_100to200->Get(histname_PESUp);
  TH1F* histo_GJets_100to200_PESDown = (TH1F*)f_GJets_100to200->Get(histname_PESDown);
  TH1F* histo_GJets_200to400 = (TH1F*)f_GJets_200to400->Get(histname);
  TH1F* histo_GJets_200to400_JESUp = (TH1F*)f_GJets_200to400->Get(histname_JESUp);
  TH1F* histo_GJets_200to400_JESDown = (TH1F*)f_GJets_200to400->Get(histname_JESDown);
  TH1F* histo_GJets_200to400_PESUp = (TH1F*)f_GJets_200to400->Get(histname_PESUp);
  TH1F* histo_GJets_200to400_PESDown = (TH1F*)f_GJets_200to400->Get(histname_PESDown);
  TH1F* histo_GJets_400to600 = (TH1F*)f_GJets_400to600->Get(histname);
  TH1F* histo_GJets_400to600_JESUp = (TH1F*)f_GJets_400to600->Get(histname_JESUp);
  TH1F* histo_GJets_400to600_JESDown = (TH1F*)f_GJets_400to600->Get(histname_JESDown);
  TH1F* histo_GJets_400to600_PESUp = (TH1F*)f_GJets_400to600->Get(histname_PESUp);
  TH1F* histo_GJets_400to600_PESDown = (TH1F*)f_GJets_400to600->Get(histname_PESDown);
  TH1F* histo_GJets_600toInf = (TH1F*)f_GJets_600toInf->Get(histname);
  TH1F* histo_GJets_600toInf_JESUp = (TH1F*)f_GJets_600toInf->Get(histname_JESUp);
  TH1F* histo_GJets_600toInf_JESDown = (TH1F*)f_GJets_600toInf->Get(histname_JESDown);
  TH1F* histo_GJets_600toInf_PESUp = (TH1F*)f_GJets_600toInf->Get(histname_PESUp);
  TH1F* histo_GJets_600toInf_PESDown = (TH1F*)f_GJets_600toInf->Get(histname_PESDown);  
  histo_GJets_40to100->Scale(int_lumi*scale_factor*frac_below0p5*20730.0/4269126.0);
  histo_GJets_40to100_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*20730.0/4269126.0);
  histo_GJets_40to100_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*20730.0/4269126.0);
  histo_GJets_40to100_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*20730.0/4269126.0);
  histo_GJets_40to100_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*20730.0/4269126.0);
  histo_GJets_100to200->Scale(int_lumi*scale_factor*frac_below0p5*9226.0/5131808.0);
  histo_GJets_100to200_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*9226.0/5131808.0);
  histo_GJets_100to200_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*9226.0/5131808.0);
  histo_GJets_100to200_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*9226.0/5131808.0);
  histo_GJets_100to200_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*9226.0/5131808.0);
  histo_GJets_200to400->Scale(int_lumi*scale_factor*frac_below0p5*2300.0/10036339.0);
  histo_GJets_200to400_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*2300.0/10036339.0);
  histo_GJets_200to400_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*2300.0/10036339.0);
  histo_GJets_200to400_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*2300.0/10036339.0);
  histo_GJets_200to400_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*2300.0/10036339.0);
  histo_GJets_400to600->Scale(int_lumi*scale_factor*frac_below0p5*277.4/2435892.0);
  histo_GJets_400to600_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*277.4/2435892.0);
  histo_GJets_400to600_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*277.4/2435892.0);
  histo_GJets_400to600_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*277.4/2435892.0);
  histo_GJets_400to600_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*277.4/2435892.0);
  histo_GJets_600toInf->Scale(int_lumi*scale_factor*frac_below0p5*93.38/2117687.0);
  histo_GJets_600toInf_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*93.38/2117687.0);
  histo_GJets_600toInf_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*93.38/2117687.0);
  histo_GJets_600toInf_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*93.38/2117687.0);
  histo_GJets_600toInf_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*93.38/2117687.0);
  Float_t int_GJets_40to100 = histo_GJets_40to100->Integral()+histo_GJets_40to100->GetBinContent(0)+histo_GJets_40to100->GetBinContent(nBins+1);
  Float_t int_GJets_40to100_JESUp = histo_GJets_40to100_JESUp->Integral()+histo_GJets_40to100_JESUp->GetBinContent(0)+histo_GJets_40to100_JESUp->GetBinContent(nBins+1);
  Float_t int_GJets_40to100_JESDown = histo_GJets_40to100_JESDown->Integral()+histo_GJets_40to100_JESDown->GetBinContent(0)+histo_GJets_40to100_JESDown->GetBinContent(nBins+1);
  Float_t int_GJets_40to100_PESUp = histo_GJets_40to100_PESUp->Integral()+histo_GJets_40to100_PESUp->GetBinContent(0)+histo_GJets_40to100_PESUp->GetBinContent(nBins+1);
  Float_t int_GJets_40to100_PESDown = histo_GJets_40to100_PESDown->Integral()+histo_GJets_40to100_PESDown->GetBinContent(0)+histo_GJets_40to100_PESDown->GetBinContent(nBins+1);
  double jeserr_GJets_40to100 = (fabs(int_GJets_40to100_JESUp-int_GJets_40to100)+fabs(int_GJets_40to100_JESDown-int_GJets_40to100))/2.0;
  double peserr_GJets_40to100 = (fabs(int_GJets_40to100_PESUp-int_GJets_40to100)+fabs(int_GJets_40to100_PESDown-int_GJets_40to100))/2.0;
  Float_t int_GJets_100to200 = histo_GJets_100to200->Integral()+histo_GJets_100to200->GetBinContent(0)+histo_GJets_100to200->GetBinContent(nBins+1);
  Float_t int_GJets_100to200_JESUp = histo_GJets_100to200_JESUp->Integral()+histo_GJets_100to200_JESUp->GetBinContent(0)+histo_GJets_100to200_JESUp->GetBinContent(nBins+1);
  Float_t int_GJets_100to200_JESDown = histo_GJets_100to200_JESDown->Integral()+histo_GJets_100to200_JESDown->GetBinContent(0)+histo_GJets_100to200_JESDown->GetBinContent(nBins+1);
  Float_t int_GJets_100to200_PESUp = histo_GJets_100to200_PESUp->Integral()+histo_GJets_100to200_PESUp->GetBinContent(0)+histo_GJets_100to200_PESUp->GetBinContent(nBins+1);
  Float_t int_GJets_100to200_PESDown = histo_GJets_100to200_PESDown->Integral()+histo_GJets_100to200_PESDown->GetBinContent(0)+histo_GJets_100to200_PESDown->GetBinContent(nBins+1);
  double jeserr_GJets_100to200 = (fabs(int_GJets_100to200_JESUp-int_GJets_100to200)+fabs(int_GJets_100to200_JESDown-int_GJets_100to200))/2.0;
  double peserr_GJets_100to200 = (fabs(int_GJets_100to200_PESUp-int_GJets_100to200)+fabs(int_GJets_100to200_PESDown-int_GJets_100to200))/2.0;
  Float_t int_GJets_200to400 = histo_GJets_200to400->Integral()+histo_GJets_200to400->GetBinContent(0)+histo_GJets_200to400->GetBinContent(nBins+1);
  Float_t int_GJets_200to400_JESUp = histo_GJets_200to400_JESUp->Integral()+histo_GJets_200to400_JESUp->GetBinContent(0)+histo_GJets_200to400_JESUp->GetBinContent(nBins+1);
  Float_t int_GJets_200to400_JESDown = histo_GJets_200to400_JESDown->Integral()+histo_GJets_200to400_JESDown->GetBinContent(0)+histo_GJets_200to400_JESDown->GetBinContent(nBins+1);
  Float_t int_GJets_200to400_PESUp = histo_GJets_200to400_PESUp->Integral()+histo_GJets_200to400_PESUp->GetBinContent(0)+histo_GJets_200to400_PESUp->GetBinContent(nBins+1);
  Float_t int_GJets_200to400_PESDown = histo_GJets_200to400_PESDown->Integral()+histo_GJets_200to400_PESDown->GetBinContent(0)+histo_GJets_200to400_PESDown->GetBinContent(nBins+1);
  double jeserr_GJets_200to400 = (fabs(int_GJets_200to400_JESUp-int_GJets_200to400)+fabs(int_GJets_200to400_JESDown-int_GJets_200to400))/2.0;
  double peserr_GJets_200to400 = (fabs(int_GJets_200to400_PESUp-int_GJets_200to400)+fabs(int_GJets_200to400_PESDown-int_GJets_200to400))/2.0;
  Float_t int_GJets_400to600 = histo_GJets_400to600->Integral()+histo_GJets_400to600->GetBinContent(0)+histo_GJets_400to600->GetBinContent(nBins+1);
  Float_t int_GJets_400to600_JESUp = histo_GJets_400to600_JESUp->Integral()+histo_GJets_400to600_JESUp->GetBinContent(0)+histo_GJets_400to600_JESUp->GetBinContent(nBins+1);
  Float_t int_GJets_400to600_JESDown = histo_GJets_400to600_JESDown->Integral()+histo_GJets_400to600_JESDown->GetBinContent(0)+histo_GJets_400to600_JESDown->GetBinContent(nBins+1);
  Float_t int_GJets_400to600_PESUp = histo_GJets_400to600_PESUp->Integral()+histo_GJets_400to600_PESUp->GetBinContent(0)+histo_GJets_400to600_PESUp->GetBinContent(nBins+1);
  Float_t int_GJets_400to600_PESDown = histo_GJets_400to600_PESDown->Integral()+histo_GJets_400to600_PESDown->GetBinContent(0)+histo_GJets_400to600_PESDown->GetBinContent(nBins+1);
  double jeserr_GJets_400to600 = (fabs(int_GJets_400to600_JESUp-int_GJets_400to600)+fabs(int_GJets_400to600_JESDown-int_GJets_400to600))/2.0;
  double peserr_GJets_400to600 = (fabs(int_GJets_400to600_PESUp-int_GJets_400to600)+fabs(int_GJets_400to600_PESDown-int_GJets_400to600))/2.0;
  Float_t int_GJets_600toInf = histo_GJets_600toInf->Integral()+histo_GJets_600toInf->GetBinContent(0)+histo_GJets_600toInf->GetBinContent(nBins+1);
  Float_t int_GJets_600toInf_JESUp = histo_GJets_600toInf_JESUp->Integral()+histo_GJets_600toInf_JESUp->GetBinContent(0)+histo_GJets_600toInf_JESUp->GetBinContent(nBins+1);
  Float_t int_GJets_600toInf_JESDown = histo_GJets_600toInf_JESDown->Integral()+histo_GJets_600toInf_JESDown->GetBinContent(0)+histo_GJets_600toInf_JESDown->GetBinContent(nBins+1);
  Float_t int_GJets_600toInf_PESUp = histo_GJets_600toInf_PESUp->Integral()+histo_GJets_600toInf_PESUp->GetBinContent(0)+histo_GJets_600toInf_PESUp->GetBinContent(nBins+1);
  Float_t int_GJets_600toInf_PESDown = histo_GJets_600toInf_PESDown->Integral()+histo_GJets_600toInf_PESDown->GetBinContent(0)+histo_GJets_600toInf_PESDown->GetBinContent(nBins+1);
  double jeserr_GJets_600toInf = (fabs(int_GJets_600toInf_JESUp-int_GJets_600toInf)+fabs(int_GJets_600toInf_JESDown-int_GJets_600toInf))/2.0;
  double peserr_GJets_600toInf = (fabs(int_GJets_600toInf_PESUp-int_GJets_600toInf)+fabs(int_GJets_600toInf_PESDown-int_GJets_600toInf))/2.0;
  Float_t err_GJets_40to100 = 0.0;
  Float_t err_GJets_100to200 = 0.0;
  Float_t err_GJets_200to400 = 0.0;
  Float_t err_GJets_400to600 = 0.0;
  Float_t err_GJets_600toInf = 0.0;
  if(int_GJets_40to100>0.0)
    err_GJets_40to100 = sqrt(int_GJets_40to100*int_GJets_40to100*((photon_scale_factor_unc*photon_scale_factor_unc)+(20730.0*int_lumi*scale_factor*frac_below0p5-int_GJets_40to100)/(4269126.0*int_GJets_40to100))+(jeserr_GJets_40to100*jeserr_GJets_40to100));
  if(int_GJets_100to200>0.0)
    err_GJets_100to200 = sqrt(int_GJets_100to200*int_GJets_100to200*((photon_scale_factor_unc*photon_scale_factor_unc)+(9226.0*int_lumi*scale_factor*frac_below0p5-int_GJets_100to200)/(5131808.0*int_GJets_100to200))+(jeserr_GJets_100to200*jeserr_GJets_100to200));
  if(int_GJets_200to400>0.0)
    err_GJets_200to400 = sqrt(int_GJets_200to400*int_GJets_200to400*((photon_scale_factor_unc*photon_scale_factor_unc)+(2300.0*int_lumi*scale_factor*frac_below0p5-int_GJets_200to400)/(10036339.0*int_GJets_200to400))+(jeserr_GJets_200to400*jeserr_GJets_200to400));
  if(int_GJets_400to600>0.0)
    err_GJets_400to600 = sqrt(int_GJets_400to600*int_GJets_400to600*((photon_scale_factor_unc*photon_scale_factor_unc)+(277.4*int_lumi*scale_factor*frac_below0p5-int_GJets_400to600)/(2435892.0*int_GJets_400to600))+(jeserr_GJets_400to600*jeserr_GJets_400to600));
  if(int_GJets_600toInf>0.0)
    err_GJets_600toInf = sqrt(int_GJets_600toInf*int_GJets_600toInf*((photon_scale_factor_unc*photon_scale_factor_unc)+(93.38*int_lumi*scale_factor*frac_below0p5-int_GJets_600toInf)/(2117687.0*int_GJets_600toInf))+(jeserr_GJets_600toInf*jeserr_GJets_600toInf));
  Float_t int_sum_GJets = int_GJets_40to100+int_GJets_100to200+int_GJets_200to400+int_GJets_400to600+int_GJets_600toInf;
  Float_t err_sum_GJets = sqrt(err_GJets_40to100*err_GJets_40to100+err_GJets_100to200*err_GJets_100to200+err_GJets_200to400*err_GJets_200to400+err_GJets_400to600*err_GJets_400to600+err_GJets_600toInf*err_GJets_600toInf);
  total_background += int_sum_GJets;
  background_unc_sumsquares += err_sum_GJets*err_sum_GJets;
  for(int i = 1; i <= nBins; i++){
    //Skipping jes and pes for now
    double int_bin_40to100 = histo_GJets_40to100->GetBinContent(i);
    double int_bin_100to200 = histo_GJets_100to200->GetBinContent(i);
    double int_bin_200to400 = histo_GJets_200to400->GetBinContent(i);
    double int_bin_400to600 = histo_GJets_400to600->GetBinContent(i);
    double int_bin_600toInf = histo_GJets_600toInf->GetBinContent(i);
    double err_bin_40to100 = 0.0;
    double err_bin_100to200 = 0.0;
    double err_bin_200to400 = 0.0;
    double err_bin_400to600 = 0.0;
    double err_bin_600toInf = 0.0;
    if(int_bin_40to100>0.0)
      err_bin_40to100 = int_bin_40to100*sqrt((20730.0*int_lumi*scale_factor*frac_below0p5-int_bin_40to100)/(4269126.0*int_bin_40to100));
    if(int_bin_100to200>0.0)
      err_bin_100to200 = int_bin_100to200*sqrt((9226.0*int_lumi*scale_factor*frac_below0p5-int_bin_100to200)/(5131808.0*int_bin_100to200));
    if(int_bin_200to400>0.0)
      err_bin_200to400 = int_bin_200to400*sqrt((2300.0*int_lumi*scale_factor*frac_below0p5-int_bin_200to400)/(10036339.0*int_bin_200to400));
    if(int_bin_400to600>0.0)
      err_bin_400to600 = int_bin_400to600*sqrt((277.4*int_lumi*scale_factor*frac_below0p5-int_bin_400to600)/(2435892.0*int_bin_400to600));
    if(int_bin_600toInf>0.0)
      err_bin_600toInf = int_bin_600toInf*sqrt((93.38*int_lumi*scale_factor*frac_below0p5-int_bin_600toInf)/(2117687.0*int_bin_600toInf));
    histo_GJets_40to100->SetBinError(i,err_bin_40to100);
    histo_GJets_100to200->SetBinError(i,err_bin_100to200);
    histo_GJets_200to400->SetBinError(i,err_bin_200to400);
    histo_GJets_400to600->SetBinError(i,err_bin_400to600);
    histo_GJets_600toInf->SetBinError(i,err_bin_600toInf);
  }
  TH1F* histo_GJets_40toInf = (TH1F*)histo_GJets_40to100->Clone("histo_GJets");
  TH1F* histo_GJets_40toInf_JESUp = (TH1F*)histo_GJets_40to100->Clone("histo_GJets_JESUp");
  TH1F* histo_GJets_40toInf_JESDown = (TH1F*)histo_GJets_40to100->Clone("histo_GJets_JESDown");
  TH1F* histo_GJets_40toInf_PESUp = (TH1F*)histo_GJets_40to100->Clone("histo_GJets_PESUp");
  TH1F* histo_GJets_40toInf_PESDown = (TH1F*)histo_GJets_40to100->Clone("histo_GJets_PESDown");
  histo_GJets_40toInf->Add(histo_GJets_100to200);
  histo_GJets_40toInf->Add(histo_GJets_200to400);
  histo_GJets_40toInf->Add(histo_GJets_400to600);
  histo_GJets_40toInf->Add(histo_GJets_600toInf);
  histo_GJets_40toInf_JESUp->Add(histo_GJets_100to200_JESUp);
  histo_GJets_40toInf_JESUp->Add(histo_GJets_200to400_JESUp);
  histo_GJets_40toInf_JESUp->Add(histo_GJets_400to600_JESUp);
  histo_GJets_40toInf_JESUp->Add(histo_GJets_600toInf_JESUp);
  histo_GJets_40toInf_JESDown->Add(histo_GJets_100to200_JESDown);
  histo_GJets_40toInf_JESDown->Add(histo_GJets_200to400_JESDown);
  histo_GJets_40toInf_JESDown->Add(histo_GJets_400to600_JESDown);
  histo_GJets_40toInf_JESDown->Add(histo_GJets_600toInf_JESDown);
  histo_GJets_40toInf_PESUp->Add(histo_GJets_100to200_PESUp);
  histo_GJets_40toInf_PESUp->Add(histo_GJets_200to400_PESUp);
  histo_GJets_40toInf_PESUp->Add(histo_GJets_400to600_PESUp);
  histo_GJets_40toInf_PESUp->Add(histo_GJets_600toInf_PESUp);
  histo_GJets_40toInf_PESDown->Add(histo_GJets_100to200_PESDown);
  histo_GJets_40toInf_PESDown->Add(histo_GJets_200to400_PESDown);
  histo_GJets_40toInf_PESDown->Add(histo_GJets_400to600_PESDown);
  histo_GJets_40toInf_PESDown->Add(histo_GJets_600toInf_PESDown);
  // histo_GJets_40toInf->SetFillColor(kBlue-8);
  histo_GJets_40toInf->SetFillColor(kBlue-2);
  histo_vector.push_back(histo_GJets_40toInf);
  
  //DEBUG
  // cout<<"GJets got"<<endl;

  TFile *f_ZllG_130_under300 = new TFile("ZnnG_JESPES_ZLLGJets_130_under300.root");
  // TFile *f_ZllG_130_over300 = new TFile("ZnnG_JESPES_ZLLGJets_130_over300.root");
  TFile *f_ZllG_300 = new TFile("ZnnG_JESPES_ZLLGJets_300.root");
  TH1F* histo_ZllG_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname))->Clone("histo_ZllG_130_under300");
  TH1F* histo_ZllG_JESUp_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_JESUp))->Clone("histo_ZllG_JESUp_130_under300");
  TH1F* histo_ZllG_JESDown_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_JESDown))->Clone("histo_ZllG_JESDown_130_under300");
  TH1F* histo_ZllG_PESUp_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_PESUp))->Clone("histo_ZllG_PESUp_130_under300");
  TH1F* histo_ZllG_PESDown_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_PESDown))->Clone("histo_ZllG_PESDown_130_under300");
  TH1F* histo_ZllG_uncorrected_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_uncorrected))->Clone("histo_ZllG_uncorrected_130_under300");
  // TH1F* histo_ZllG_130_over300 = (TH1F*)((TH1F*)f_ZllG_130_over300->Get(histname))->Clone("histo_ZllG_130_over300");
  // TH1F* histo_ZllG_JESUp_130_over300 = (TH1F*)((TH1F*)f_ZllG_130_over300->Get(histname_JESUp))->Clone("histo_ZllG_JESUp_130_over300");
  // TH1F* histo_ZllG_JESDown_130_over300 = (TH1F*)((TH1F*)f_ZllG_130_over300->Get(histname_JESDown))->Clone("histo_ZllG_JESDown_130_over300");
  // TH1F* histo_ZllG_PESUp_130_over300 = (TH1F*)((TH1F*)f_ZllG_130_over300->Get(histname_PESUp))->Clone("histo_ZllG_PESUp_130_over300");
  // TH1F* histo_ZllG_PESDown_130_over300 = (TH1F*)((TH1F*)f_ZllG_130_over300->Get(histname_PESDown))->Clone("histo_ZllG_PESDown_130_over300");
  TH1F* histo_ZllG_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname))->Clone("histo_ZllG_300");
  TH1F* histo_ZllG_JESUp_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_JESUp))->Clone("histo_ZllG_JESUp_300");
  TH1F* histo_ZllG_JESDown_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_JESDown))->Clone("histo_ZllG_JESDown_300");
  TH1F* histo_ZllG_PESUp_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_PESUp))->Clone("histo_ZllG_PESUp_300");
  TH1F* histo_ZllG_PESDown_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_PESDown))->Clone("histo_ZllG_PESDown_300");
  TH1F* histo_ZllG_uncorrected_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_uncorrected))->Clone("histo_ZllG_uncorrected_300");
  histo_ZllG_130_under300->Scale(int_lumi*scale_factor*frac_below0p5*0.148/489423.0);
  histo_ZllG_JESUp_130_under300->Scale(int_lumi*scale_factor*frac_below0p5*0.148/489423.0);
  histo_ZllG_JESDown_130_under300->Scale(int_lumi*scale_factor*frac_below0p5*0.148/489423.0);
  histo_ZllG_PESUp_130_under300->Scale(int_lumi*scale_factor*frac_below0p5*0.148/489423.0);
  histo_ZllG_PESDown_130_under300->Scale(int_lumi*scale_factor*frac_below0p5*0.148/489423.0);
  histo_ZllG_uncorrected_130_under300->Scale(int_lumi*scale_factor*frac_below0p5*0.148/489423.0);
  histo_ZllG_300->Scale(int_lumi*scale_factor*frac_below0p5*0.0092/1707773.0);
  histo_ZllG_JESUp_300->Scale(int_lumi*scale_factor*frac_below0p5*0.0092/1707773.0);
  histo_ZllG_JESDown_300->Scale(int_lumi*scale_factor*frac_below0p5*0.0092/1707773.0);
  histo_ZllG_PESUp_300->Scale(int_lumi*scale_factor*frac_below0p5*0.0092/1707773.0);
  histo_ZllG_PESDown_300->Scale(int_lumi*scale_factor*frac_below0p5*0.0092/1707773.0);
  histo_ZllG_uncorrected_300->Scale(int_lumi*scale_factor*frac_below0p5*0.0092/1707773.0);
  TH1F* histo_ZllG_combined = (TH1F*)histo_ZllG_300->Clone("histo_ZllG_combined");
  TH1F* histo_ZllG_JESUp_combined = (TH1F*)histo_ZllG_JESUp_300->Clone("histo_ZllG_JESUp_combined");
  TH1F* histo_ZllG_JESDown_combined = (TH1F*)histo_ZllG_JESDown_300->Clone("histo_ZllG_JESDown_combined");
  TH1F* histo_ZllG_PESUp_combined = (TH1F*)histo_ZllG_PESUp_300->Clone("histo_ZllG_PESUp_combined");
  TH1F* histo_ZllG_PESDown_combined = (TH1F*)histo_ZllG_PESDown_300->Clone("histo_ZllG_PESDown_combined");
  TH1F* histo_ZllG_uncorrected_combined = (TH1F*)histo_ZllG_uncorrected_300->Clone("histo_ZllG_uncorrected_combined");
  histo_ZllG_combined->SetStats(0);
  histo_ZllG_combined->Add(histo_ZllG_130_under300);
  histo_ZllG_JESUp_combined->Add(histo_ZllG_JESUp_130_under300);
  histo_ZllG_JESDown_combined->Add(histo_ZllG_JESDown_130_under300);
  histo_ZllG_PESUp_combined->Add(histo_ZllG_PESUp_130_under300);
  histo_ZllG_PESDown_combined->Add(histo_ZllG_PESDown_130_under300);
  histo_ZllG_uncorrected_combined->Add(histo_ZllG_uncorrected_130_under300);
  Float_t int_ZllG_130_under300 = histo_ZllG_130_under300->Integral()+histo_ZllG_130_under300->GetBinContent(0)+histo_ZllG_130_under300->GetBinContent(nBins+1);
  Float_t int_ZllG_300 = histo_ZllG_300->Integral()+histo_ZllG_300->GetBinContent(0)+histo_ZllG_300->GetBinContent(nBins+1);
  Float_t int_ZllG = histo_ZllG_combined->Integral()+histo_ZllG_combined->GetBinContent(0)+histo_ZllG_combined->GetBinContent(nBins+1);
  Float_t int_ZllG_JESUp = histo_ZllG_JESUp_combined->Integral()+histo_ZllG_JESUp_combined->GetBinContent(0)+histo_ZllG_JESUp_combined->GetBinContent(nBins+1);
  Float_t int_ZllG_JESDown = histo_ZllG_JESDown_combined->Integral()+histo_ZllG_JESDown_combined->GetBinContent(0)+histo_ZllG_JESDown_combined->GetBinContent(nBins+1);
  Float_t int_ZllG_PESUp = histo_ZllG_PESUp_combined->Integral()+histo_ZllG_PESUp_combined->GetBinContent(0)+histo_ZllG_PESUp_combined->GetBinContent(nBins+1);
  Float_t int_ZllG_PESDown = histo_ZllG_PESDown_combined->Integral()+histo_ZllG_PESDown_combined->GetBinContent(0)+histo_ZllG_PESDown_combined->GetBinContent(nBins+1);
  Float_t int_ZllG_uncorrected = histo_ZllG_uncorrected_combined->Integral()+histo_ZllG_uncorrected_combined->GetBinContent(0)+histo_ZllG_uncorrected_combined->GetBinContent(nBins+1);
  double jeserr_ZllG = (fabs(int_ZllG_JESUp-int_ZllG)+fabs(int_ZllG_JESDown-int_ZllG))/2.0;
  double peserr_ZllG = (fabs(int_ZllG_PESUp-int_ZllG)+fabs(int_ZllG_PESDown-int_ZllG))/2.0;
  double xsecerr_ZllG = fabs(int_ZllG_uncorrected-int_ZllG);
  TH1F* histo_ZllG_RenUp_130_under300;
  TH1F* histo_ZllG_RenDown_130_under300;
  TH1F* histo_ZllG_FacUp_130_under300;
  TH1F* histo_ZllG_FacDown_130_under300;
  TH1F* histo_ZllG_PDFUp_130_under300;
  TH1F* histo_ZllG_PDFDown_130_under300;
  // TH1F* histo_ZllG_RenUp_130_over300;
  // TH1F* histo_ZllG_RenDown_130_over300;
  // TH1F* histo_ZllG_FacUp_130_over300;
  // TH1F* histo_ZllG_FacDown_130_over300;
  // TH1F* histo_ZllG_PDFUp_130_over300;
  // TH1F* histo_ZllG_PDFDown_130_over300;
  TH1F* histo_ZllG_RenUp_300;
  TH1F* histo_ZllG_RenDown_300;
  TH1F* histo_ZllG_FacUp_300;
  TH1F* histo_ZllG_FacDown_300;
  TH1F* histo_ZllG_PDFUp_300;
  TH1F* histo_ZllG_PDFDown_300;
  double renerr_ZllG = 0.0;
  double facerr_ZllG = 0.0;
  double pdferr_ZllG = 0.0;
  if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3"){
    // These should each already be scaled to the appropriate luminosity, scale factor, cross section, and total number of events
    string xaxis_variable = "XAXIS_VARIABLE_NOT_SET";
    if (histname == "Photon_Et_range_3") xaxis_variable = "Pt";
    else if (histname == "pfMET_3") xaxis_variable = "MET";
    else if (histname == "Mt_3") xaxis_variable = "Mt";
    TString filename = TString("histos_"+xaxis_variable+"_ZnnG_pdfscale_ZLLGJets_130_under300.root");
    TString filename_ext = TString("histos_"+xaxis_variable+"_ZnnG_pdfscale_ZLLGJets_300.root");
    TFile* f_pdfscale_ZllG_130_under300 = new TFile(filename);
    histo_ZllG_RenUp_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZnnG_pdfscale_ZLLGJets_130_under300_renUp");
    histo_ZllG_RenDown_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZnnG_pdfscale_ZLLGJets_130_under300_renDown");
    histo_ZllG_FacUp_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZnnG_pdfscale_ZLLGJets_130_under300_facUp");
    histo_ZllG_FacDown_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZnnG_pdfscale_ZLLGJets_130_under300_facDown");
    histo_ZllG_PDFUp_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZnnG_pdfscale_ZLLGJets_130_under300_pdfUp");
    histo_ZllG_PDFDown_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZnnG_pdfscale_ZLLGJets_130_under300_pdfDown");
    TFile* f_pdfscale_ZllG_300 = new TFile(filename_ext);
    histo_ZllG_RenUp_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZnnG_pdfscale_ZLLGJets_300_renUp");
    histo_ZllG_RenDown_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZnnG_pdfscale_ZLLGJets_300_renDown");
    histo_ZllG_FacUp_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZnnG_pdfscale_ZLLGJets_300_facUp");
    histo_ZllG_FacDown_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZnnG_pdfscale_ZLLGJets_300_facDown");
    histo_ZllG_PDFUp_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZnnG_pdfscale_ZLLGJets_300_pdfUp");
    histo_ZllG_PDFDown_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZnnG_pdfscale_ZLLGJets_300_pdfDown");
    // TFile* f_pdfscale_ZllG_130_over300 = new TFile("histos_ZnnG_pdfscale_ZLLGJets_130_over300.root");
    // histo_ZllG_RenUp_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZnnG_pdfscale_ZLLGJets_130_over300_renUp");
    // histo_ZllG_RenDown_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZnnG_pdfscale_ZLLGJets_130_over300_renDown");
    // histo_ZllG_FacUp_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZnnG_pdfscale_ZLLGJets_130_over300_facUp");
    // histo_ZllG_FacDown_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZnnG_pdfscale_ZLLGJets_130_over300_facDown");
    // histo_ZllG_PDFUp_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZnnG_pdfscale_ZLLGJets_130_over300_pdfUp");
    // histo_ZllG_PDFDown_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZnnG_pdfscale_ZLLGJets_130_over300_pdfDown");
    histo_ZllG_RenUp_300->Add(histo_ZllG_RenUp_130_under300);
    histo_ZllG_RenDown_300->Add(histo_ZllG_RenDown_130_under300);
    histo_ZllG_FacUp_300->Add(histo_ZllG_FacUp_130_under300);
    histo_ZllG_FacDown_300->Add(histo_ZllG_FacDown_130_under300);
    histo_ZllG_PDFUp_300->Add(histo_ZllG_PDFUp_130_under300);
    histo_ZllG_PDFDown_300->Add(histo_ZllG_PDFDown_130_under300);
    histo_ZllG_RenUp_300->Scale(frac_below0p5);
    histo_ZllG_RenDown_300->Scale(frac_below0p5);
    histo_ZllG_FacUp_300->Scale(frac_below0p5);
    histo_ZllG_FacDown_300->Scale(frac_below0p5);
    histo_ZllG_PDFUp_300->Scale(frac_below0p5);
    histo_ZllG_PDFDown_300->Scale(frac_below0p5);
    Float_t int_ZllG_RenUp = histo_ZllG_RenUp_300->Integral()+histo_ZllG_RenUp_300->GetBinContent(0)+histo_ZllG_RenUp_300->GetBinContent(nBins+1);
    Float_t int_ZllG_RenDown = histo_ZllG_RenDown_300->Integral()+histo_ZllG_RenDown_300->GetBinContent(0)+histo_ZllG_RenDown_300->GetBinContent(nBins+1);
    Float_t int_ZllG_FacUp = histo_ZllG_FacUp_300->Integral()+histo_ZllG_FacUp_300->GetBinContent(0)+histo_ZllG_FacUp_300->GetBinContent(nBins+1);
    Float_t int_ZllG_FacDown = histo_ZllG_FacDown_300->Integral()+histo_ZllG_FacDown_300->GetBinContent(0)+histo_ZllG_FacDown_300->GetBinContent(nBins+1);
    Float_t int_ZllG_PDFUp = histo_ZllG_PDFUp_300->Integral()+histo_ZllG_PDFUp_300->GetBinContent(0)+histo_ZllG_PDFUp_300->GetBinContent(nBins+1);
    Float_t int_ZllG_PDFDown = histo_ZllG_PDFDown_300->Integral()+histo_ZllG_PDFDown_300->GetBinContent(0)+histo_ZllG_PDFDown_300->GetBinContent(nBins+1);
    renerr_ZllG = (fabs(int_ZllG_RenUp-int_ZllG)+fabs(int_ZllG_RenDown-int_ZllG))/2.0;
    facerr_ZllG = (fabs(int_ZllG_FacUp-int_ZllG)+fabs(int_ZllG_FacDown-int_ZllG))/2.0;
    pdferr_ZllG = (fabs(int_ZllG_PDFUp-int_ZllG)+fabs(int_ZllG_PDFDown-int_ZllG))/2.0;
  }
  for(int i = 1; i <= nBins; i++){
    double int_bin_130_under300 = histo_ZllG_130_under300->GetBinContent(i);
    double int_bin_300 = histo_ZllG_300->GetBinContent(i);
    double int_bin = histo_ZllG_combined->GetBinContent(i);
    double jesup = histo_ZllG_JESUp_combined->GetBinContent(i);
    double jesdown = histo_ZllG_JESDown_combined->GetBinContent(i);
    double pesup = histo_ZllG_PESUp_combined->GetBinContent(i);
    double pesdown = histo_ZllG_PESDown_combined->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    double uncorrected = histo_ZllG_uncorrected_combined->GetBinContent(i);
    xsec_shift_ZLLG[i-1] = fabs(uncorrected-int_bin);
    if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3"){
      double renup = histo_ZllG_RenUp_300->GetBinContent(i);
      double rendown = histo_ZllG_RenDown_300->GetBinContent(i);
      double facup = histo_ZllG_FacUp_300->GetBinContent(i);
      double facdown = histo_ZllG_FacDown_300->GetBinContent(i);
      double pdfup = histo_ZllG_PDFUp_300->GetBinContent(i);
      double pdfdown = histo_ZllG_PDFDown_300->GetBinContent(i);
      renup_shift[i-1] += renup-int_bin;
      rendown_shift[i-1] += rendown-int_bin;
      facup_shift[i-1] += facup-int_bin;
      facdown_shift[i-1] += facdown-int_bin;
      pdfup_shift[i-1] += pdfup-int_bin;
      pdfdown_shift[i-1] += pdfdown-int_bin;
    }
    // cout<<"ZllG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = sqrt(int_bin_130_under300*(0.148*int_lumi*scale_factor*frac_below0p5-int_bin_130_under300)/489423.0 + int_bin_300*(0.0092*int_lumi*scale_factor*frac_below0p5-int_bin_300)/1707773.0);
    histo_ZllG_combined->SetBinError(i,err_bin);
  }
  Float_t err_ZllG = 0.0;
  if(int_ZllG > 0.0)
    err_ZllG = sqrt(int_ZllG_130_under300*int_ZllG_130_under300*((photon_scale_factor_unc*photon_scale_factor_unc)+(0.148*int_lumi*scale_factor*frac_below0p5-int_ZllG_130_under300)/(489423.0*int_ZllG_130_under300))+int_ZllG_300*int_ZllG_300*((photon_scale_factor_unc*photon_scale_factor_unc)+(0.0092*int_lumi*scale_factor*frac_below0p5-int_ZllG_300)/(1707773.0*int_ZllG_300))+(xsecerr_ZllG*xsecerr_ZllG)+(jeserr_ZllG*jeserr_ZllG)+(peserr_ZllG*peserr_ZllG)+(renerr_ZllG*renerr_ZllG)+(facerr_ZllG*facerr_ZllG)+(pdferr_ZllG*pdferr_ZllG));
  total_background += int_ZllG;
  background_unc_sumsquares += err_ZllG*err_ZllG;
  // histo_ZllG_130_over300->SetFillColor(kSpring-9);
  histo_ZllG_combined->SetFillColor(kRed-10);
  histo_ZllG_combined->SetLineColor(kRed-10);
  histo_vector.push_back(histo_ZllG_combined);
  
  //DEBUG
  // cout<<"ZllG got"<<endl;
    
  TFile *f_TTG = new TFile("ZnnG_JESPES_TTGJets.root");
  TH1F* histo_TTG = (TH1F*)((TH1F*)f_TTG->Get(histname))->Clone("histo_TTG");
  TH1F* histo_TTG_JESUp = (TH1F*)((TH1F*)f_TTG->Get(histname_JESUp))->Clone("histo_TTG_JESUp");
  TH1F* histo_TTG_JESDown = (TH1F*)((TH1F*)f_TTG->Get(histname_JESDown))->Clone("histo_TTG_JESDown");
  TH1F* histo_TTG_PESUp = (TH1F*)((TH1F*)f_TTG->Get(histname_PESUp))->Clone("histo_TTG_PESUp");
  TH1F* histo_TTG_PESDown = (TH1F*)((TH1F*)f_TTG->Get(histname_PESDown))->Clone("histo_TTG_PESDown");
  histo_TTG->SetStats(0);
  histo_TTG->Scale(int_lumi*scale_factor*frac_below0p5*3.697/3170400.0);
  histo_TTG_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*3.697/3170400.0);
  histo_TTG_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*3.697/3170400.0);
  histo_TTG_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*3.697/3170400.0);
  histo_TTG_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*3.697/3170400.0);
  Float_t int_TTG = histo_TTG->Integral()+histo_TTG->GetBinContent(0)+histo_TTG->GetBinContent(nBins+1);
  Float_t int_TTG_JESUp = histo_TTG_JESUp->Integral()+histo_TTG_JESUp->GetBinContent(0)+histo_TTG_JESUp->GetBinContent(nBins+1);
  Float_t int_TTG_JESDown = histo_TTG_JESDown->Integral()+histo_TTG_JESDown->GetBinContent(0)+histo_TTG_JESDown->GetBinContent(nBins+1);
  Float_t int_TTG_PESUp = histo_TTG_PESUp->Integral()+histo_TTG_PESUp->GetBinContent(0)+histo_TTG_PESUp->GetBinContent(nBins+1);
  Float_t int_TTG_PESDown = histo_TTG_PESDown->Integral()+histo_TTG_PESDown->GetBinContent(0)+histo_TTG_PESDown->GetBinContent(nBins+1);
  double jeserr_TTG = (fabs(int_TTG_JESUp-int_TTG)+fabs(int_TTG_JESDown-int_TTG))/2.0;
  double peserr_TTG = (fabs(int_TTG_PESUp-int_TTG)+fabs(int_TTG_PESDown-int_TTG))/2.0;
  Float_t err_TTG = 0.0;
  if(int_TTG > 0.0)
    err_TTG = sqrt(int_TTG*int_TTG*((photon_scale_factor_unc*photon_scale_factor_unc)+(3.697*int_lumi*scale_factor*frac_below0p5-int_TTG)/(3170400.0*int_TTG))+(jeserr_TTG*jeserr_TTG)+(peserr_TTG*peserr_TTG));
  total_background += int_TTG;
  background_unc_sumsquares += err_TTG*err_TTG;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_TTG->GetBinContent(i);
    double jesup = histo_TTG_JESUp->GetBinContent(i);
    double jesdown = histo_TTG_JESDown->GetBinContent(i);
    double pesup = histo_TTG_PESUp->GetBinContent(i);
    double pesdown = histo_TTG_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"TTG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((3.697*int_lumi*scale_factor*frac_below0p5-int_bin)/(3170400.0*int_bin));
    histo_TTG->SetBinError(i,err_bin);
  }
  // histo_TTG->SetFillColor(kOrange-5);
  histo_TTG->SetFillColor(kYellow+1);
  histo_vector.push_back(histo_TTG);
  
  //DEBUG
  // cout<<"TTG got"<<endl;
  
  TFile *f_TG = new TFile("ZnnG_JESPES_TGJets.root");
  TH1F* histo_TG = (TH1F*)((TH1F*)f_TG->Get(histname))->Clone("histo_TG");
  TH1F* histo_TG_JESUp = (TH1F*)((TH1F*)f_TG->Get(histname_JESUp))->Clone("histo_TG_JESUp");
  TH1F* histo_TG_JESDown = (TH1F*)((TH1F*)f_TG->Get(histname_JESDown))->Clone("histo_TG_JESDown");
  TH1F* histo_TG_PESUp = (TH1F*)((TH1F*)f_TG->Get(histname_PESUp))->Clone("histo_TG_PESUp");
  TH1F* histo_TG_PESDown = (TH1F*)((TH1F*)f_TG->Get(histname_PESDown))->Clone("histo_TG_PESDown");
  histo_TG->SetStats(0);
  histo_TG->Scale(int_lumi*scale_factor*frac_below0p5*2.967/310437.0);
  histo_TG_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*2.967/310437.0);
  histo_TG_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*2.967/310437.0);
  histo_TG_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*2.967/310437.0);
  histo_TG_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*2.967/310437.0);
  Float_t int_TG = histo_TG->Integral()+histo_TG->GetBinContent(0)+histo_TG->GetBinContent(nBins+1);
  Float_t int_TG_JESUp = histo_TG_JESUp->Integral()+histo_TG_JESUp->GetBinContent(0)+histo_TG_JESUp->GetBinContent(nBins+1);
  Float_t int_TG_JESDown = histo_TG_JESDown->Integral()+histo_TG_JESDown->GetBinContent(0)+histo_TG_JESDown->GetBinContent(nBins+1);
  Float_t int_TG_PESUp = histo_TG_PESUp->Integral()+histo_TG_PESUp->GetBinContent(0)+histo_TG_PESUp->GetBinContent(nBins+1);
  Float_t int_TG_PESDown = histo_TG_PESDown->Integral()+histo_TG_PESDown->GetBinContent(0)+histo_TG_PESDown->GetBinContent(nBins+1);
  double jeserr_TG = (fabs(int_TG_JESUp-int_TG)+fabs(int_TG_JESDown-int_TG))/2.0;
  double peserr_TG = (fabs(int_TG_PESUp-int_TG)+fabs(int_TG_PESDown-int_TG))/2.0;
  Float_t err_TG = 0.0;
  if(int_TG > 0.0)
    err_TG = sqrt(int_TG*int_TG*((photon_scale_factor_unc*photon_scale_factor_unc)+(2.967*int_lumi*scale_factor*frac_below0p5-int_TG)/(310437.0*int_TG))+(jeserr_TG*jeserr_TG)+(peserr_TG*peserr_TG));
  total_background += int_TG;
  background_unc_sumsquares += err_TG*err_TG;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_TG->GetBinContent(i);
    double jesup = histo_TG_JESUp->GetBinContent(i);
    double jesdown = histo_TG_JESDown->GetBinContent(i);
    double pesup = histo_TG_PESUp->GetBinContent(i);
    double pesdown = histo_TG_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"TG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((2.967*int_lumi*scale_factor*frac_below0p5-int_bin)/(310437.0*int_bin));
    histo_TG->SetBinError(i,err_bin);
  }
  // histo_TG->SetFillColor(kRed-7);
  histo_TG->SetFillColor(kRed-10);
  histo_TG->SetLineColor(kRed-10);
  histo_vector.push_back(histo_TG);
  
  //DEBUG
  // cout<<"TG got"<<endl;
  
  TFile* f_WWG = new TFile("ZnnG_JESPES_WWG.root");
  TH1F* histo_WWG = (TH1F*)((TH1F*)f_WWG->Get(histname))->Clone("histo_WWG");
  TH1F* histo_WWG_JESUp = (TH1F*)((TH1F*)f_WWG->Get(histname_JESUp))->Clone("histo_WWG_JESUp");
  TH1F* histo_WWG_JESDown = (TH1F*)((TH1F*)f_WWG->Get(histname_JESDown))->Clone("histo_WWG_JESDown");
  TH1F* histo_WWG_PESUp = (TH1F*)((TH1F*)f_WWG->Get(histname_PESUp))->Clone("histo_WWG_PESUp");
  TH1F* histo_WWG_PESDown = (TH1F*)((TH1F*)f_WWG->Get(histname_PESDown))->Clone("histo_WWG_PESDown");
  histo_WWG->SetStats(0);
  histo_WWG->Scale(int_lumi*scale_factor*frac_below0p5*0.2147/827604.0);
  histo_WWG_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.2147/827604.0);
  histo_WWG_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.2147/827604.0);
  histo_WWG_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.2147/827604.0);
  histo_WWG_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.2147/827604.0);
  Float_t int_WWG = histo_WWG->Integral()+histo_WWG->GetBinContent(0)+histo_WWG->GetBinContent(nBins+1);
  Float_t int_WWG_JESUp = histo_WWG_JESUp->Integral()+histo_WWG_JESUp->GetBinContent(0)+histo_WWG_JESUp->GetBinContent(nBins+1);
  Float_t int_WWG_JESDown = histo_WWG_JESDown->Integral()+histo_WWG_JESDown->GetBinContent(0)+histo_WWG_JESDown->GetBinContent(nBins+1);
  Float_t int_WWG_PESUp = histo_WWG_PESUp->Integral()+histo_WWG_PESUp->GetBinContent(0)+histo_WWG_PESUp->GetBinContent(nBins+1);
  Float_t int_WWG_PESDown = histo_WWG_PESDown->Integral()+histo_WWG_PESDown->GetBinContent(0)+histo_WWG_PESDown->GetBinContent(nBins+1);
  double jeserr_WWG = (fabs(int_WWG_JESUp-int_WWG)+fabs(int_WWG_JESDown-int_WWG))/2.0;
  double peserr_WWG = (fabs(int_WWG_PESUp-int_WWG)+fabs(int_WWG_PESDown-int_WWG))/2.0;
  Float_t err_WWG = 0.0;
  if(int_WWG > 0.0)
    err_WWG = sqrt(int_WWG*int_WWG*((photon_scale_factor_unc*photon_scale_factor_unc)+(0.2147*int_lumi*scale_factor*frac_below0p5-int_WWG)/(827604.0*int_WWG))+(jeserr_WWG*jeserr_WWG)+(peserr_WWG*peserr_WWG));
  total_background += int_WWG;
  background_unc_sumsquares += err_WWG*err_WWG;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_WWG->GetBinContent(i);
    double jesup = histo_WWG_JESUp->GetBinContent(i);
    double jesdown = histo_WWG_JESDown->GetBinContent(i);
    double pesup = histo_WWG_PESUp->GetBinContent(i);
    double pesdown = histo_WWG_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"WWG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((0.2147*int_lumi*scale_factor*frac_below0p5-int_bin)/(827604.0*int_bin));
    histo_WWG->SetBinError(i,err_bin);
  }
  // histo_WWG->SetFillColor(kTeal+3);
  histo_WWG->SetFillColor(kRed-10);
  histo_WWG->SetLineColor(kRed-10);
  histo_vector.push_back(histo_WWG);
  
  //DEBUG
  // cout<<"WWG got"<<endl;
  
  TFile* f_diphoton = new TFile("ZnnG_JESPES_Diphoton.root");
  TH1F* histo_diphoton = (TH1F*)((TH1F*)f_diphoton->Get(histname))->Clone("histo_diphoton");
  TH1F* histo_diphoton_JESUp = (TH1F*)((TH1F*)f_diphoton->Get(histname_JESUp))->Clone("histo_diphoton_JESUp");
  TH1F* histo_diphoton_JESDown = (TH1F*)((TH1F*)f_diphoton->Get(histname_JESDown))->Clone("histo_diphoton_JESDown");
  TH1F* histo_diphoton_PESUp = (TH1F*)((TH1F*)f_diphoton->Get(histname_PESUp))->Clone("histo_diphoton_PESUp");
  TH1F* histo_diphoton_PESDown = (TH1F*)((TH1F*)f_diphoton->Get(histname_PESDown))->Clone("histo_diphoton_PESDown");
  histo_diphoton->SetStats(0);
  histo_diphoton->Scale(int_lumi*scale_factor*frac_below0p5*135.1/18633300.0);
  histo_diphoton_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*135.1/18633300.0);
  histo_diphoton_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*135.1/18633300.0);
  histo_diphoton_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*135.1/18633300.0);
  histo_diphoton_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*135.1/18633300.0);
  Float_t int_diphoton = histo_diphoton->Integral()+histo_diphoton->GetBinContent(0)+histo_diphoton->GetBinContent(nBins+1);
  Float_t int_diphoton_JESUp = histo_diphoton_JESUp->Integral()+histo_diphoton_JESUp->GetBinContent(0)+histo_diphoton_JESUp->GetBinContent(nBins+1);
  Float_t int_diphoton_JESDown = histo_diphoton_JESDown->Integral()+histo_diphoton_JESDown->GetBinContent(0)+histo_diphoton_JESDown->GetBinContent(nBins+1);
  Float_t int_diphoton_PESUp = histo_diphoton_PESUp->Integral()+histo_diphoton_PESUp->GetBinContent(0)+histo_diphoton_PESUp->GetBinContent(nBins+1);
  Float_t int_diphoton_PESDown = histo_diphoton_PESDown->Integral()+histo_diphoton_PESDown->GetBinContent(0)+histo_diphoton_PESDown->GetBinContent(nBins+1);
  double jeserr_diphoton = (fabs(int_diphoton_JESUp-int_diphoton)+fabs(int_diphoton_JESDown-int_diphoton))/2.0;
  double peserr_diphoton = (fabs(int_diphoton_PESUp-int_diphoton)+fabs(int_diphoton_PESDown-int_diphoton))/2.0;
  Float_t err_diphoton = 0.0;
  if(int_diphoton > 0.0)
    err_diphoton = sqrt(int_diphoton*int_diphoton*((photon_scale_factor_unc*photon_scale_factor_unc)+(135.1*int_lumi*scale_factor*frac_below0p5-int_diphoton)/(18633300.0*int_diphoton))+(jeserr_diphoton*jeserr_diphoton)+(peserr_diphoton*peserr_diphoton));
  total_background += int_diphoton;
  background_unc_sumsquares += err_diphoton*err_diphoton;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_diphoton->GetBinContent(i);
    double jesup = histo_diphoton_JESUp->GetBinContent(i);
    double jesdown = histo_diphoton_JESDown->GetBinContent(i);
    double pesup = histo_diphoton_PESUp->GetBinContent(i);
    double pesdown = histo_diphoton_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"diphoton: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((135.1*int_lumi*scale_factor*frac_below0p5-int_bin)/(18633300.0*int_bin));
    histo_diphoton->SetBinError(i,err_bin);
  }
  // histo_diphoton->SetFillColor(kOrange+10);
  histo_diphoton->SetFillColor(28);
  histo_vector.push_back(histo_diphoton);
  
  //DEBUG
  // cout<<"diphoton got"<<endl;
  
  TFile *f_WZ = new TFile("ZnnG_JESPES_WZ.root");
  TH1F* histo_WZ = (TH1F*)((TH1F*)f_WZ->Get(histname))->Clone("histo_WZ");
  TH1F* histo_WZ_JESUp = (TH1F*)((TH1F*)f_WZ->Get(histname_JESUp))->Clone("histo_WZ_JESUp");
  TH1F* histo_WZ_JESDown = (TH1F*)((TH1F*)f_WZ->Get(histname_JESDown))->Clone("histo_WZ_JESDown");
  TH1F* histo_WZ_PESUp = (TH1F*)((TH1F*)f_WZ->Get(histname_PESUp))->Clone("histo_WZ_PESUp");
  TH1F* histo_WZ_PESDown = (TH1F*)((TH1F*)f_WZ->Get(histname_PESDown))->Clone("histo_WZ_PESDown");
  histo_WZ->SetStats(0);
  histo_WZ->Scale(int_lumi*scale_factor*frac_below0p5*66.1/2995783.0);
  histo_WZ_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*66.1/2995783.0);
  histo_WZ_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*66.1/2995783.0);
  histo_WZ_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*66.1/2995783.0);
  histo_WZ_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*66.1/2995783.0);
  Float_t int_WZ = histo_WZ->Integral()+histo_WZ->GetBinContent(0)+histo_WZ->GetBinContent(nBins+1);
  Float_t int_WZ_JESUp = histo_WZ_JESUp->Integral()+histo_WZ_JESUp->GetBinContent(0)+histo_WZ_JESUp->GetBinContent(nBins+1);
  Float_t int_WZ_JESDown = histo_WZ_JESDown->Integral()+histo_WZ_JESDown->GetBinContent(0)+histo_WZ_JESDown->GetBinContent(nBins+1);
  Float_t int_WZ_PESUp = histo_WZ_PESUp->Integral()+histo_WZ_PESUp->GetBinContent(0)+histo_WZ_PESUp->GetBinContent(nBins+1);
  Float_t int_WZ_PESDown = histo_WZ_PESDown->Integral()+histo_WZ_PESDown->GetBinContent(0)+histo_WZ_PESDown->GetBinContent(nBins+1);
  double jeserr_WZ = (fabs(int_WZ_JESUp-int_WZ)+fabs(int_WZ_JESDown-int_WZ))/2.0;
  double peserr_WZ = (fabs(int_WZ_PESUp-int_WZ)+fabs(int_WZ_PESDown-int_WZ))/2.0;
  Float_t err_WZ = 0.0;
  if(int_WZ > 0.0)
    err_WZ = sqrt(int_WZ*int_WZ*((photon_scale_factor_unc*photon_scale_factor_unc)+(66.1*int_lumi*scale_factor*frac_below0p5-int_WZ)/(2995783.0*int_WZ))+(jeserr_WZ*jeserr_WZ)+(peserr_WZ*peserr_WZ));
  total_background += int_WZ;
  background_unc_sumsquares += err_WZ*err_WZ;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_WZ->GetBinContent(i);
    double jesup = histo_WZ_JESUp->GetBinContent(i);
    double jesdown = histo_WZ_JESDown->GetBinContent(i);
    double pesup = histo_WZ_PESUp->GetBinContent(i);
    double pesdown = histo_WZ_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"WZ: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((66.1*int_lumi*scale_factor*frac_below0p5-int_bin)/(2995783.0*int_bin));
    histo_WZ->SetBinError(i,err_bin);
  }
  // histo_WZ->SetFillColor(kRed-5);
  histo_WZ->SetFillColor(kRed-10);
  histo_vector.push_back(histo_WZ);
  
  //DEBUG
  // cout<<"WZ got"<<endl;
  
  TFile *f_ZZ = new TFile("ZnnG_JESPES_ZZ.root");
  TH1F* histo_ZZ = (TH1F*)((TH1F*)f_ZZ->Get(histname))->Clone("histo_ZZ");
  TH1F* histo_ZZ_JESUp = (TH1F*)((TH1F*)f_ZZ->Get(histname_JESUp))->Clone("histo_ZZ_JESUp");
  TH1F* histo_ZZ_JESDown = (TH1F*)((TH1F*)f_ZZ->Get(histname_JESDown))->Clone("histo_ZZ_JESDown");
  TH1F* histo_ZZ_PESUp = (TH1F*)((TH1F*)f_ZZ->Get(histname_PESUp))->Clone("histo_ZZ_PESUp");
  TH1F* histo_ZZ_PESDown = (TH1F*)((TH1F*)f_ZZ->Get(histname_PESDown))->Clone("histo_ZZ_PESDown");
  histo_ZZ->SetStats(0);
  histo_ZZ->Scale(int_lumi*scale_factor*frac_below0p5*16.52/922975.0);
  histo_ZZ_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*16.52/922975.0);
  histo_ZZ_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*16.52/922975.0);
  histo_ZZ_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*16.52/922975.0);
  histo_ZZ_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*16.52/922975.0);
  Float_t int_ZZ = histo_ZZ->Integral()+histo_ZZ->GetBinContent(0)+histo_ZZ->GetBinContent(nBins+1);
  Float_t int_ZZ_JESUp = histo_ZZ_JESUp->Integral()+histo_ZZ_JESUp->GetBinContent(0)+histo_ZZ_JESUp->GetBinContent(nBins+1);
  Float_t int_ZZ_JESDown = histo_ZZ_JESDown->Integral()+histo_ZZ_JESDown->GetBinContent(0)+histo_ZZ_JESDown->GetBinContent(nBins+1);
  Float_t int_ZZ_PESUp = histo_ZZ_PESUp->Integral()+histo_ZZ_PESUp->GetBinContent(0)+histo_ZZ_PESUp->GetBinContent(nBins+1);
  Float_t int_ZZ_PESDown = histo_ZZ_PESDown->Integral()+histo_ZZ_PESDown->GetBinContent(0)+histo_ZZ_PESDown->GetBinContent(nBins+1);
  double jeserr_ZZ = (fabs(int_ZZ_JESUp-int_ZZ)+fabs(int_ZZ_JESDown-int_ZZ))/2.0;
  double peserr_ZZ = (fabs(int_ZZ_PESUp-int_ZZ)+fabs(int_ZZ_PESDown-int_ZZ))/2.0;
  Float_t err_ZZ = 0.0;
  if(int_ZZ > 0.0)
    err_ZZ = sqrt(int_ZZ*int_ZZ*((16.52*int_lumi*scale_factor*frac_below0p5-int_ZZ)/(922975.0*int_ZZ))+(jeserr_ZZ*jeserr_ZZ)+(peserr_ZZ*peserr_ZZ));
  total_background += int_ZZ;
  background_unc_sumsquares += err_ZZ*err_ZZ;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_ZZ->GetBinContent(i);
    double jesup = histo_ZZ_JESUp->GetBinContent(i);
    double jesdown = histo_ZZ_JESDown->GetBinContent(i);
    double pesup = histo_ZZ_PESUp->GetBinContent(i);
    double pesdown = histo_ZZ_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"ZZ: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((16.52*int_lumi*scale_factor*frac_below0p5-int_bin)/(922975.0*int_bin));
    histo_ZZ->SetBinError(i,err_bin);
  }
  // histo_ZZ->SetFillColor(kTeal-1);
  histo_ZZ->SetFillColor(kRed-10);
  histo_ZZ->SetLineColor(kRed-10);
  histo_vector.push_back(histo_ZZ);
  
  //DEBUG
  // cout<<"ZZ got"<<endl;
  
  TFile *f_WMuNu = new TFile("ZnnG_JESPES_WToMuNu.root");
  TH1F* histo_WMuNu = (TH1F*)((TH1F*)f_WMuNu->Get(histname))->Clone("histo_WMuNu");
  TH1F* histo_WMuNu_JESUp = (TH1F*)((TH1F*)f_WMuNu->Get(histname_JESUp))->Clone("histo_WMuNu_JESUp");
  TH1F* histo_WMuNu_JESDown = (TH1F*)((TH1F*)f_WMuNu->Get(histname_JESDown))->Clone("histo_WMuNu_JESDown");
  TH1F* histo_WMuNu_PESUp = (TH1F*)((TH1F*)f_WMuNu->Get(histname_PESUp))->Clone("histo_WMuNu_PESUp");
  TH1F* histo_WMuNu_PESDown = (TH1F*)((TH1F*)f_WMuNu->Get(histname_PESDown))->Clone("histo_WMuNu_PESDown");
  histo_WMuNu->SetStats(0);
  histo_WMuNu->Scale(int_lumi*scale_factor*frac_below0p5*174.0/1985332.0);
  histo_WMuNu_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*174.0/1985332.0);
  histo_WMuNu_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*174.0/1985332.0);
  histo_WMuNu_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*174.0/1985332.0);
  histo_WMuNu_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*174.0/1985332.0);
  Float_t int_WMuNu = histo_WMuNu->Integral()+histo_WMuNu->GetBinContent(0)+histo_WMuNu->GetBinContent(nBins+1);
  Float_t int_WMuNu_JESUp = histo_WMuNu_JESUp->Integral()+histo_WMuNu_JESUp->GetBinContent(0)+histo_WMuNu_JESUp->GetBinContent(nBins+1);
  Float_t int_WMuNu_JESDown = histo_WMuNu_JESDown->Integral()+histo_WMuNu_JESDown->GetBinContent(0)+histo_WMuNu_JESDown->GetBinContent(nBins+1);
  Float_t int_WMuNu_PESUp = histo_WMuNu_PESUp->Integral()+histo_WMuNu_PESUp->GetBinContent(0)+histo_WMuNu_PESUp->GetBinContent(nBins+1);
  Float_t int_WMuNu_PESDown = histo_WMuNu_PESDown->Integral()+histo_WMuNu_PESDown->GetBinContent(0)+histo_WMuNu_PESDown->GetBinContent(nBins+1);
  double jeserr_WMuNu = (fabs(int_WMuNu_JESUp-int_WMuNu)+fabs(int_WMuNu_JESDown-int_WMuNu))/2.0;
  double peserr_WMuNu = (fabs(int_WMuNu_PESUp-int_WMuNu)+fabs(int_WMuNu_PESDown-int_WMuNu))/2.0;
  Float_t err_WMuNu = 0.0;
  if(int_WMuNu > 0.0)
    err_WMuNu = sqrt(int_WMuNu*int_WMuNu*((174.0*int_lumi*scale_factor*frac_below0p5-int_WMuNu)/(1985332.0*int_WMuNu))+(jeserr_WMuNu*jeserr_WMuNu)+(peserr_WMuNu*peserr_WMuNu));
  total_background += int_WMuNu;
  background_unc_sumsquares += err_WMuNu*err_WMuNu;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_WMuNu->GetBinContent(i);
    double jesup = histo_WMuNu_JESUp->GetBinContent(i);
    double jesdown = histo_WMuNu_JESDown->GetBinContent(i);
    double pesup = histo_WMuNu_PESUp->GetBinContent(i);
    double pesdown = histo_WMuNu_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"WMuNu: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((174.0*int_lumi*scale_factor*frac_below0p5-int_bin)/(1985332.0*int_bin));
    histo_WMuNu->SetBinError(i,err_bin);
  }
  // histo_WMuNu->SetFillColor(12);
  histo_WMuNu->SetFillColor(TColor::GetColor("#FF6633"));
  histo_vector.push_back(histo_WMuNu);
  
  //DEBUG
  // cout<<"WMuNu got"<<endl;
  
  TFile *f_WTauNu = new TFile("ZnnG_JESPES_WToTauNu.root");
  TH1F* histo_WTauNu = (TH1F*)((TH1F*)f_WTauNu->Get(histname))->Clone("histo_WTauNu");
  TH1F* histo_WTauNu_JESUp = (TH1F*)((TH1F*)f_WTauNu->Get(histname_JESUp))->Clone("histo_WTauNu_JESUp");
  TH1F* histo_WTauNu_JESDown = (TH1F*)((TH1F*)f_WTauNu->Get(histname_JESDown))->Clone("histo_WTauNu_JESDown");
  TH1F* histo_WTauNu_PESUp = (TH1F*)((TH1F*)f_WTauNu->Get(histname_PESUp))->Clone("histo_WTauNu_PESUp");
  TH1F* histo_WTauNu_PESDown = (TH1F*)((TH1F*)f_WTauNu->Get(histname_PESDown))->Clone("histo_WTauNu_PESDown");
  histo_WTauNu->SetStats(0);
  histo_WTauNu->Scale(int_lumi*scale_factor*frac_below0p5*165.0/1371381.0);
  histo_WTauNu_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*165.0/1371381.0);
  histo_WTauNu_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*165.0/1371381.0);
  histo_WTauNu_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*165.0/1371381.0);
  histo_WTauNu_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*165.0/1371381.0);
  Float_t int_WTauNu = histo_WTauNu->Integral()+histo_WTauNu->GetBinContent(0)+histo_WTauNu->GetBinContent(nBins+1);
  Float_t int_WTauNu_JESUp = histo_WTauNu_JESUp->Integral()+histo_WTauNu_JESUp->GetBinContent(0)+histo_WTauNu_JESUp->GetBinContent(nBins+1);
  Float_t int_WTauNu_JESDown = histo_WTauNu_JESDown->Integral()+histo_WTauNu_JESDown->GetBinContent(0)+histo_WTauNu_JESDown->GetBinContent(nBins+1);
  Float_t int_WTauNu_PESUp = histo_WTauNu_PESUp->Integral()+histo_WTauNu_PESUp->GetBinContent(0)+histo_WTauNu_PESUp->GetBinContent(nBins+1);
  Float_t int_WTauNu_PESDown = histo_WTauNu_PESDown->Integral()+histo_WTauNu_PESDown->GetBinContent(0)+histo_WTauNu_PESDown->GetBinContent(nBins+1);
  double jeserr_WTauNu = (fabs(int_WTauNu_JESUp-int_WTauNu)+fabs(int_WTauNu_JESDown-int_WTauNu))/2.0;
  double peserr_WTauNu = (fabs(int_WTauNu_PESUp-int_WTauNu)+fabs(int_WTauNu_PESDown-int_WTauNu))/2.0;
  Float_t err_WTauNu = 0.0;
  if(int_WTauNu > 0.0)
    err_WTauNu = sqrt(int_WTauNu*int_WTauNu*((165.0*int_lumi*scale_factor*frac_below0p5-int_WTauNu)/(1371381.0*int_WTauNu))+(jeserr_WTauNu*jeserr_WTauNu)+(peserr_WTauNu*peserr_WTauNu));
  total_background += int_WTauNu;
  background_unc_sumsquares += err_WTauNu*err_WTauNu;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_WTauNu->GetBinContent(i);
    double jesup = histo_WTauNu_JESUp->GetBinContent(i);
    double jesdown = histo_WTauNu_JESDown->GetBinContent(i);
    double pesup = histo_WTauNu_PESUp->GetBinContent(i);
    double pesdown = histo_WTauNu_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"WTauNu: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((165.0*int_lumi*scale_factor*frac_below0p5-int_bin)/(1371381.0*int_bin));
    histo_WTauNu->SetBinError(i,err_bin);
  }
  // histo_WTauNu->SetFillColor(18);
  histo_WTauNu->SetFillColor(kMagenta-5);
  histo_vector.push_back(histo_WTauNu);
  
  TFile *f_WW = new TFile("ZnnG_JESPES_WW.root");
  TH1F* histo_WW = (TH1F*)((TH1F*)f_WW->Get(histname))->Clone("histo_WW");
  TH1F* histo_WW_JESUp = (TH1F*)((TH1F*)f_WW->Get(histname_JESUp))->Clone("histo_WW_JESUp");
  TH1F* histo_WW_JESDown = (TH1F*)((TH1F*)f_WW->Get(histname_JESDown))->Clone("histo_WW_JESDown");
  TH1F* histo_WW_PESUp = (TH1F*)((TH1F*)f_WW->Get(histname_PESUp))->Clone("histo_WW_PESUp");
  TH1F* histo_WW_PESDown = (TH1F*)((TH1F*)f_WW->Get(histname_PESDown))->Clone("histo_WW_PESDown");
  histo_WW->SetStats(0);
  histo_WW->Scale(int_lumi*scale_factor*frac_below0p5*64.21/6658858.0);
  histo_WW_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*64.21/6658858.0);
  histo_WW_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*64.21/6658858.0);
  histo_WW_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*64.21/6658858.0);
  histo_WW_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*64.21/6658858.0);
  Float_t int_WW = histo_WW->Integral()+histo_WW->GetBinContent(0)+histo_WW->GetBinContent(nBins+1);
  Float_t int_WW_JESUp = histo_WW_JESUp->Integral()+histo_WW_JESUp->GetBinContent(0)+histo_WW_JESUp->GetBinContent(nBins+1);
  Float_t int_WW_JESDown = histo_WW_JESDown->Integral()+histo_WW_JESDown->GetBinContent(0)+histo_WW_JESDown->GetBinContent(nBins+1);
  Float_t int_WW_PESUp = histo_WW_PESUp->Integral()+histo_WW_PESUp->GetBinContent(0)+histo_WW_PESUp->GetBinContent(nBins+1);
  Float_t int_WW_PESDown = histo_WW_PESDown->Integral()+histo_WW_PESDown->GetBinContent(0)+histo_WW_PESDown->GetBinContent(nBins+1);
  double jeserr_WW = (fabs(int_WW_JESUp-int_WW)+fabs(int_WW_JESDown-int_WW))/2.0;
  double peserr_WW = (fabs(int_WW_PESUp-int_WW)+fabs(int_WW_PESDown-int_WW))/2.0;
  Float_t err_WW = 0.0;
  if(int_WW > 0.0)
    err_WW = sqrt(int_WW*int_WW*((64.21*int_lumi*scale_factor*frac_below0p5-int_WW)/(6658858.0*int_WW))+(jeserr_WW*jeserr_WW)+(peserr_WW*peserr_WW));
  total_background += int_WW;
  background_unc_sumsquares += err_WW*err_WW;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_WW->GetBinContent(i);
    double jesup = histo_WW_JESUp->GetBinContent(i);
    double jesdown = histo_WW_JESDown->GetBinContent(i);
    double pesup = histo_WW_PESUp->GetBinContent(i);
    double pesdown = histo_WW_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"WW: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((64.21*int_lumi*scale_factor*frac_below0p5-int_bin)/(6658858.0*int_bin));
    histo_WW->SetBinError(i,err_bin);
  }
  // histo_WW->SetFillColor(18);
  histo_WW->SetFillColor(kRed-10);
  histo_vector.push_back(histo_WW);
  
  //DEBUG
  // cout<<"WW got"<<endl;

  TFile *f_WZG = new TFile("ZnnG_JESPES_WZG.root");
  TH1F* histo_WZG = (TH1F*)((TH1F*)f_WZG->Get(histname))->Clone("histo_WZG");
  TH1F* histo_WZG_JESUp = (TH1F*)((TH1F*)f_WZG->Get(histname_JESUp))->Clone("histo_WZG_JESUp");
  TH1F* histo_WZG_JESDown = (TH1F*)((TH1F*)f_WZG->Get(histname_JESDown))->Clone("histo_WZG_JESDown");
  TH1F* histo_WZG_PESUp = (TH1F*)((TH1F*)f_WZG->Get(histname_PESUp))->Clone("histo_WZG_PESUp");
  TH1F* histo_WZG_PESDown = (TH1F*)((TH1F*)f_WZG->Get(histname_PESDown))->Clone("histo_WZG_PESDown");
  histo_WZG->SetStats(0);
  histo_WZG->Scale(int_lumi*scale_factor*frac_below0p5*0.04123/844824.0);
  histo_WZG_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.04123/844824.0);
  histo_WZG_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.04123/844824.0);
  histo_WZG_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.04123/844824.0);
  histo_WZG_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.04123/844824.0);
  Float_t int_WZG = histo_WZG->Integral()+histo_WZG->GetBinContent(0)+histo_WZG->GetBinContent(nBins+1);
  Float_t int_WZG_JESUp = histo_WZG_JESUp->Integral()+histo_WZG_JESUp->GetBinContent(0)+histo_WZG_JESUp->GetBinContent(nBins+1);
  Float_t int_WZG_JESDown = histo_WZG_JESDown->Integral()+histo_WZG_JESDown->GetBinContent(0)+histo_WZG_JESDown->GetBinContent(nBins+1);
  Float_t int_WZG_PESUp = histo_WZG_PESUp->Integral()+histo_WZG_PESUp->GetBinContent(0)+histo_WZG_PESUp->GetBinContent(nBins+1);
  Float_t int_WZG_PESDown = histo_WZG_PESDown->Integral()+histo_WZG_PESDown->GetBinContent(0)+histo_WZG_PESDown->GetBinContent(nBins+1);
  double jeserr_WZG = (fabs(int_WZG_JESUp-int_WZG)+fabs(int_WZG_JESDown-int_WZG))/2.0;
  double peserr_WZG = (fabs(int_WZG_PESUp-int_WZG)+fabs(int_WZG_PESDown-int_WZG))/2.0;
  Float_t err_WZG = 0.0;
  if(int_WZG > 0.0)
    err_WZG = sqrt(int_WZG*int_WZG*((0.04123*int_lumi*scale_factor*frac_below0p5-int_WZG)/(844824.0*int_WZG))+(jeserr_WZG*jeserr_WZG)+(peserr_WZG*peserr_WZG));
  total_background += int_WZG;
  background_unc_sumsquares += err_WZG*err_WZG;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_WZG->GetBinContent(i);
    double jesup = histo_WZG_JESUp->GetBinContent(i);
    double jesdown = histo_WZG_JESDown->GetBinContent(i);
    double pesup = histo_WZG_PESUp->GetBinContent(i);
    double pesdown = histo_WZG_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"WZG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((0.04123*int_lumi*scale_factor*frac_below0p5-int_bin)/(844824.0*int_bin));
    histo_WZG->SetBinError(i,err_bin);
  }
  // histo_WZG->SetFillColor(18);
  histo_WZG->SetFillColor(kRed-10);
  histo_vector.push_back(histo_WZG);
  
  //DEBUG
  // cout<<"WZG got"<<endl;

  TFile *f_WGG = new TFile("ZnnG_JESPES_WGGJets.root");
  TH1F* histo_WGG = (TH1F*)((TH1F*)f_WGG->Get(histname))->Clone("histo_WGG");
  TH1F* histo_WGG_JESUp = (TH1F*)((TH1F*)f_WGG->Get(histname_JESUp))->Clone("histo_WGG_JESUp");
  TH1F* histo_WGG_JESDown = (TH1F*)((TH1F*)f_WGG->Get(histname_JESDown))->Clone("histo_WGG_JESDown");
  TH1F* histo_WGG_PESUp = (TH1F*)((TH1F*)f_WGG->Get(histname_PESUp))->Clone("histo_WGG_PESUp");
  TH1F* histo_WGG_PESDown = (TH1F*)((TH1F*)f_WGG->Get(histname_PESDown))->Clone("histo_WGG_PESDown");
  histo_WGG->SetStats(0);
  histo_WGG->Scale(int_lumi*scale_factor*frac_below0p5*1.711/428254.0);
  histo_WGG_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*1.711/428254.0);
  histo_WGG_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*1.711/428254.0);
  histo_WGG_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*1.711/428254.0);
  histo_WGG_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*1.711/428254.0);
  Float_t int_WGG = histo_WGG->Integral()+histo_WGG->GetBinContent(0)+histo_WGG->GetBinContent(nBins+1);
  Float_t int_WGG_JESUp = histo_WGG_JESUp->Integral()+histo_WGG_JESUp->GetBinContent(0)+histo_WGG_JESUp->GetBinContent(nBins+1);
  Float_t int_WGG_JESDown = histo_WGG_JESDown->Integral()+histo_WGG_JESDown->GetBinContent(0)+histo_WGG_JESDown->GetBinContent(nBins+1);
  Float_t int_WGG_PESUp = histo_WGG_PESUp->Integral()+histo_WGG_PESUp->GetBinContent(0)+histo_WGG_PESUp->GetBinContent(nBins+1);
  Float_t int_WGG_PESDown = histo_WGG_PESDown->Integral()+histo_WGG_PESDown->GetBinContent(0)+histo_WGG_PESDown->GetBinContent(nBins+1);
  double jeserr_WGG = (fabs(int_WGG_JESUp-int_WGG)+fabs(int_WGG_JESDown-int_WGG))/2.0;
  double peserr_WGG = (fabs(int_WGG_PESUp-int_WGG)+fabs(int_WGG_PESDown-int_WGG))/2.0;
  Float_t err_WGG = 0.0;
  if(int_WGG > 0.0)
    err_WGG = sqrt(int_WGG*int_WGG*((1.711*int_lumi*scale_factor*frac_below0p5-int_WGG)/(428254.0*int_WGG))+(jeserr_WGG*jeserr_WGG)+(peserr_WGG*peserr_WGG));
  total_background += int_WGG;
  background_unc_sumsquares += err_WGG*err_WGG;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_WGG->GetBinContent(i);
    double jesup = histo_WGG_JESUp->GetBinContent(i);
    double jesdown = histo_WGG_JESDown->GetBinContent(i);
    double pesup = histo_WGG_PESUp->GetBinContent(i);
    double pesdown = histo_WGG_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"WGG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((1.711*int_lumi*scale_factor*frac_below0p5-int_bin)/(428254.0*int_bin));
    histo_WGG->SetBinError(i,err_bin);
  }
  // histo_WGG->SetFillColor(18);
  histo_WGG->SetFillColor(kRed-10);
  histo_vector.push_back(histo_WGG);
  
  //DEBUG
  // cout<<"WGG got"<<endl;

  TFile *f_ZGGToNuNuGG = new TFile("ZnnG_JESPES_ZGGToNuNuGG.root");
  TH1F* histo_ZGGToNuNuGG = (TH1F*)((TH1F*)f_ZGGToNuNuGG->Get(histname))->Clone("histo_ZGGToNuNuGG");
  TH1F* histo_ZGGToNuNuGG_JESUp = (TH1F*)((TH1F*)f_ZGGToNuNuGG->Get(histname_JESUp))->Clone("histo_ZGGToNuNuGG_JESUp");
  TH1F* histo_ZGGToNuNuGG_JESDown = (TH1F*)((TH1F*)f_ZGGToNuNuGG->Get(histname_JESDown))->Clone("histo_ZGGToNuNuGG_JESDown");
  TH1F* histo_ZGGToNuNuGG_PESUp = (TH1F*)((TH1F*)f_ZGGToNuNuGG->Get(histname_PESUp))->Clone("histo_ZGGToNuNuGG_PESUp");
  TH1F* histo_ZGGToNuNuGG_PESDown = (TH1F*)((TH1F*)f_ZGGToNuNuGG->Get(histname_PESDown))->Clone("histo_ZGGToNuNuGG_PESDown");
  histo_ZGGToNuNuGG->SetStats(0);
  histo_ZGGToNuNuGG->Scale(int_lumi*scale_factor*frac_below0p5*0.07477/764438.0);
  histo_ZGGToNuNuGG_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.07477/764438.0);
  histo_ZGGToNuNuGG_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.07477/764438.0);
  histo_ZGGToNuNuGG_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.07477/764438.0);
  histo_ZGGToNuNuGG_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.07477/764438.0);
  Float_t int_ZGGToNuNuGG = histo_ZGGToNuNuGG->Integral()+histo_ZGGToNuNuGG->GetBinContent(0)+histo_ZGGToNuNuGG->GetBinContent(nBins+1);
  Float_t int_ZGGToNuNuGG_JESUp = histo_ZGGToNuNuGG_JESUp->Integral()+histo_ZGGToNuNuGG_JESUp->GetBinContent(0)+histo_ZGGToNuNuGG_JESUp->GetBinContent(nBins+1);
  Float_t int_ZGGToNuNuGG_JESDown = histo_ZGGToNuNuGG_JESDown->Integral()+histo_ZGGToNuNuGG_JESDown->GetBinContent(0)+histo_ZGGToNuNuGG_JESDown->GetBinContent(nBins+1);
  Float_t int_ZGGToNuNuGG_PESUp = histo_ZGGToNuNuGG_PESUp->Integral()+histo_ZGGToNuNuGG_PESUp->GetBinContent(0)+histo_ZGGToNuNuGG_PESUp->GetBinContent(nBins+1);
  Float_t int_ZGGToNuNuGG_PESDown = histo_ZGGToNuNuGG_PESDown->Integral()+histo_ZGGToNuNuGG_PESDown->GetBinContent(0)+histo_ZGGToNuNuGG_PESDown->GetBinContent(nBins+1);
  double jeserr_ZGGToNuNuGG = (fabs(int_ZGGToNuNuGG_JESUp-int_ZGGToNuNuGG)+fabs(int_ZGGToNuNuGG_JESDown-int_ZGGToNuNuGG))/2.0;
  double peserr_ZGGToNuNuGG = (fabs(int_ZGGToNuNuGG_PESUp-int_ZGGToNuNuGG)+fabs(int_ZGGToNuNuGG_PESDown-int_ZGGToNuNuGG))/2.0;
  Float_t err_ZGGToNuNuGG = 0.0;
  if(int_ZGGToNuNuGG > 0.0)
    err_ZGGToNuNuGG = sqrt(int_ZGGToNuNuGG*int_ZGGToNuNuGG*((0.07477*int_lumi*scale_factor*frac_below0p5-int_ZGGToNuNuGG)/(764438.0*int_ZGGToNuNuGG))+(jeserr_ZGGToNuNuGG*jeserr_ZGGToNuNuGG)+(peserr_ZGGToNuNuGG*peserr_ZGGToNuNuGG));
  total_background += int_ZGGToNuNuGG;
  background_unc_sumsquares += err_ZGGToNuNuGG*err_ZGGToNuNuGG;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_ZGGToNuNuGG->GetBinContent(i);
    double jesup = histo_ZGGToNuNuGG_JESUp->GetBinContent(i);
    double jesdown = histo_ZGGToNuNuGG_JESDown->GetBinContent(i);
    double pesup = histo_ZGGToNuNuGG_PESUp->GetBinContent(i);
    double pesdown = histo_ZGGToNuNuGG_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"ZGGToNuNuGG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((0.07477*int_lumi*scale_factor*frac_below0p5-int_bin)/(764438.0*int_bin));
    histo_ZGGToNuNuGG->SetBinError(i,err_bin);
  }
  // histo_ZGGToNuNuGG->SetFillColor(18);
  histo_ZGGToNuNuGG->SetFillColor(kRed-10);
  histo_vector.push_back(histo_ZGGToNuNuGG);
  
  //DEBUG
  // cout<<"ZGGToNuNuGG got"<<endl;

  TFile *f_DMVMx1Mv1000 = new TFile("ZnnG_JESPES_DM_V_Mx-1_Mv-1000.root");
  TH1F* histo_DMVMx1Mv1000 = (TH1F*)((TH1F*)f_DMVMx1Mv1000->Get(histname))->Clone("histo_DMVMx1Mv1000");
  TH1F* histo_DMVMx1Mv1000_JESUp = (TH1F*)((TH1F*)f_DMVMx1Mv1000->Get(histname_JESUp))->Clone("histo_DMVMx1Mv1000_JESUp");
  TH1F* histo_DMVMx1Mv1000_JESDown = (TH1F*)((TH1F*)f_DMVMx1Mv1000->Get(histname_JESDown))->Clone("histo_DMVMx1Mv1000_JESDown");
  TH1F* histo_DMVMx1Mv1000_PESUp = (TH1F*)((TH1F*)f_DMVMx1Mv1000->Get(histname_PESUp))->Clone("histo_DMVMx1Mv1000_PESUp");
  TH1F* histo_DMVMx1Mv1000_PESDown = (TH1F*)((TH1F*)f_DMVMx1Mv1000->Get(histname_PESDown))->Clone("histo_DMVMx1Mv1000_PESDown");
  histo_DMVMx1Mv1000->SetStats(0);
  histo_DMVMx1Mv1000->Scale(int_lumi*scale_factor*frac_below0p5*0.01446/1371381.0);
  histo_DMVMx1Mv1000_JESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.01446/1371381.0);
  histo_DMVMx1Mv1000_JESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.01446/1371381.0);
  histo_DMVMx1Mv1000_PESUp->Scale(int_lumi*scale_factor*frac_below0p5*0.01446/1371381.0);
  histo_DMVMx1Mv1000_PESDown->Scale(int_lumi*scale_factor*frac_below0p5*0.01446/1371381.0);
  Float_t int_DMVMx1Mv1000 = histo_DMVMx1Mv1000->Integral()+histo_DMVMx1Mv1000->GetBinContent(0)+histo_DMVMx1Mv1000->GetBinContent(nBins+1);
  Float_t int_DMVMx1Mv1000_JESUp = histo_DMVMx1Mv1000_JESUp->Integral()+histo_DMVMx1Mv1000_JESUp->GetBinContent(0)+histo_DMVMx1Mv1000_JESUp->GetBinContent(nBins+1);
  Float_t int_DMVMx1Mv1000_JESDown = histo_DMVMx1Mv1000_JESDown->Integral()+histo_DMVMx1Mv1000_JESDown->GetBinContent(0)+histo_DMVMx1Mv1000_JESDown->GetBinContent(nBins+1);
  Float_t int_DMVMx1Mv1000_PESUp = histo_DMVMx1Mv1000_PESUp->Integral()+histo_DMVMx1Mv1000_PESUp->GetBinContent(0)+histo_DMVMx1Mv1000_PESUp->GetBinContent(nBins+1);
  Float_t int_DMVMx1Mv1000_PESDown = histo_DMVMx1Mv1000_PESDown->Integral()+histo_DMVMx1Mv1000_PESDown->GetBinContent(0)+histo_DMVMx1Mv1000_PESDown->GetBinContent(nBins+1);
  double jeserr_DMVMx1Mv1000 = (fabs(int_DMVMx1Mv1000_JESUp-int_DMVMx1Mv1000)+fabs(int_DMVMx1Mv1000_JESDown-int_DMVMx1Mv1000))/2.0;
  double peserr_DMVMx1Mv1000 = (fabs(int_DMVMx1Mv1000_PESUp-int_DMVMx1Mv1000)+fabs(int_DMVMx1Mv1000_PESDown-int_DMVMx1Mv1000))/2.0;
  Float_t err_DMVMx1Mv1000 = 0.0;
  if(int_DMVMx1Mv1000 > 0.0)
    err_DMVMx1Mv1000 = sqrt(int_DMVMx1Mv1000*int_DMVMx1Mv1000*((0.01446*int_lumi*scale_factor*frac_below0p5-int_DMVMx1Mv1000)/(1371381.0*int_DMVMx1Mv1000))+(jeserr_DMVMx1Mv1000*jeserr_DMVMx1Mv1000)+(peserr_DMVMx1Mv1000*peserr_DMVMx1Mv1000));
  total_background += int_DMVMx1Mv1000;
  background_unc_sumsquares += err_DMVMx1Mv1000*err_DMVMx1Mv1000;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_DMVMx1Mv1000->GetBinContent(i);
    double jesup = histo_DMVMx1Mv1000_JESUp->GetBinContent(i);
    double jesdown = histo_DMVMx1Mv1000_JESDown->GetBinContent(i);
    double pesup = histo_DMVMx1Mv1000_PESUp->GetBinContent(i);
    double pesdown = histo_DMVMx1Mv1000_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"DMVMx1Mv1000: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((0.01446*int_lumi*scale_factor*frac_below0p5-int_bin)/(1371381.0*int_bin));
    histo_DMVMx1Mv1000->SetBinError(i,err_bin);
  }
  
  // TH1F* histo_ZNuNuG_phoSFUp = (TH1F*)histo_ZNuNuG->Clone("histo_ZNuNuG_phoSFUp");
  // TH1F* histo_ZNuNuG_phoSFDown = (TH1F*)histo_ZNuNuG->Clone("histo_ZNuNuG_phoSFDown");
  // TH1F* histo_WG_phoSFUp = (TH1F*)histo_WG->Clone("histo_WG_phoSFUp");
  // TH1F* histo_WG_phoSFDown = (TH1F*)histo_WG->Clone("histo_WG_phoSFDown");
  // TH1F* histo_GJets_40toInf_phoSFUp = (TH1F*)histo_GJets_40toInf->Clone("histo_GJets_40toInf_phoSFUp");
  // TH1F* histo_GJets_40toInf_phoSFDown = (TH1F*)histo_GJets_40toInf->Clone("histo_GJets_40toInf_phoSFDown");
  // TH1F* histo_ZllG_combined_phoSFUp = (TH1F*)histo_ZllG_combined->Clone("histo_ZllG_combined_phoSFUp");
  // TH1F* histo_ZllG_combined_phoSFDown = (TH1F*)histo_ZllG_combined->Clone("histo_ZllG_combined_phoSFDown");
  // TH1F* histo_TTG_phoSFUp = (TH1F*)histo_TTG->Clone("histo_TTG_phoSFUp");
  // TH1F* histo_TTG_phoSFDown = (TH1F*)histo_TTG->Clone("histo_TTG_phoSFDown");
  // TH1F* histo_TG_phoSFUp = (TH1F*)histo_TG->Clone("histo_TG_phoSFUp");
  // TH1F* histo_TG_phoSFDown = (TH1F*)histo_TG->Clone("histo_TG_phoSFDown");
  // TH1F* histo_WWG_phoSFUp = (TH1F*)histo_WWG->Clone("histo_WWG_phoSFUp");
  // TH1F* histo_WWG_phoSFDown = (TH1F*)histo_WWG->Clone("histo_WWG_phoSFDown");
  // TH1F* histo_diphoton_phoSFUp = (TH1F*)histo_diphoton->Clone("histo_diphoton_phoSFUp");
  // TH1F* histo_diphoton_phoSFDown = (TH1F*)histo_diphoton->Clone("histo_diphoton_phoSFDown");
  // TH1F* histo_WZ_phoSFUp = (TH1F*)histo_WZ->Clone("histo_WZ_phoSFUp");
  // TH1F* histo_WZ_phoSFDown = (TH1F*)histo_WZ->Clone("histo_WZ_phoSFDown");
  // TH1F* histo_ZZ_phoSFUp = (TH1F*)histo_ZZ->Clone("histo_ZZ_phoSFUp");
  // TH1F* histo_ZZ_phoSFDown = (TH1F*)histo_ZZ->Clone("histo_ZZ_phoSFDown");
  // TH1F* histo_WMuNu_phoSFUp = (TH1F*)histo_WMuNu->Clone("histo_WMuNu_phoSFUp");
  // TH1F* histo_WMuNu_phoSFDown = (TH1F*)histo_WMuNu->Clone("histo_WMuNu_phoSFDown");
  // TH1F* histo_WTauNu_phoSFUp = (TH1F*)histo_WTauNu->Clone("histo_WTauNu_phoSFUp");
  // TH1F* histo_WTauNu_phoSFDown = (TH1F*)histo_WTauNu->Clone("histo_WTauNu_phoSFDown");
  // TH1F* histo_WW_phoSFUp = (TH1F*)histo_WW->Clone("histo_WW_phoSFUp");
  // TH1F* histo_WW_phoSFDown = (TH1F*)histo_WW->Clone("histo_WW_phoSFDown");
  // TH1F* histo_WZG_phoSFUp = (TH1F*)histo_WZG->Clone("histo_WZG_phoSFUp");
  // TH1F* histo_WZG_phoSFDown = (TH1F*)histo_WZG->Clone("histo_WZG_phoSFDown");
  // TH1F* histo_WGG_phoSFUp = (TH1F*)histo_WGG->Clone("histo_WGG_phoSFUp");
  // TH1F* histo_WGG_phoSFDown = (TH1F*)histo_WGG->Clone("histo_WGG_phoSFDown");
  // TH1F* histo_ZGGToNuNuGG_phoSFUp = (TH1F*)histo_ZGGToNuNuGG->Clone("histo_ZGGToNuNuGG_phoSFUp");
  // TH1F* histo_ZGGToNuNuGG_phoSFDown = (TH1F*)histo_ZGGToNuNuGG->Clone("histo_ZGGToNuNuGG_phoSFDown");
  // for(int i = 1; i <= nBins; i++){
  //   histo_ZNuNuG->Write();
  //   histo_ZNuNuG_JESUp->Write();
  //   histo_ZNuNuG_JESDown->Write();
  //   histo_ZNuNuG_PESUp->Write();
  //   histo_ZNuNuG_PESDown->Write();
  //   histo_ZNuNuG_uncorrected->Write();
  //   histo_WG->Write();
  //   histo_WG_JESUp->Write();
  //   histo_WG_JESDown->Write();
  //   histo_WG_PESUp->Write();
  //   histo_WG_PESDown->Write();
  //   histo_WG_uncorrected->Write();
  //   histo_GJets_40toInf->Write();
  //   histo_GJets_40toInf_JESUp->Write();
  //   histo_GJets_40toInf_JESDown->Write();
  //   histo_GJets_40toInf_PESUp->Write();
  //   histo_GJets_40toInf_PESDown->Write();
  //   histo_ZllG_combined->Write();
  //   histo_ZllG_JESUp_combined->Write();
  //   histo_ZllG_JESDown_combined->Write();
  //   histo_ZllG_PESUp_combined->Write();
  //   histo_ZllG_PESDown_combined->Write();
  //   histo_ZllG_uncorrected_combined->Write();
  //   histo_TTG->Write();
  //   histo_TTG_JESUp->Write();
  //   histo_TTG_JESDown->Write();
  //   histo_TTG_PESUp->Write();
  //   histo_TTG_PESDown->Write();
  //   histo_TG->Write();
  //   histo_TG_JESUp->Write();
  //   histo_TG_JESDown->Write();
  //   histo_TG_PESUp->Write();
  //   histo_TG_PESDown->Write();
  //   histo_WWG->Write();
  //   histo_WWG_JESUp->Write();
  //   histo_WWG_JESDown->Write();
  //   histo_WWG_PESUp->Write();
  //   histo_WWG_PESDown->Write();
  //   histo_diphoton->Write();
  //   histo_diphoton_JESUp->Write();
  //   histo_diphoton_JESDown->Write();
  //   histo_diphoton_PESUp->Write();
  //   histo_diphoton_PESDown->Write();
  //   histo_WZ->Write();
  //   histo_WZ_JESUp->Write();
  //   histo_WZ_JESDown->Write();
  //   histo_WZ_PESUp->Write();
  //   histo_WZ_PESDown->Write();
  //   histo_ZZ->Write();
  //   histo_ZZ_JESUp->Write();
  //   histo_ZZ_JESDown->Write();
  //   histo_ZZ_PESUp->Write();
  //   histo_ZZ_PESDown->Write();
  //   histo_WMuNu->Write();
  //   histo_WMuNu_JESUp->Write();
  //   histo_WMuNu_JESDown->Write();
  //   histo_WMuNu_PESUp->Write();
  //   histo_WMuNu_PESDown->Write();
  //   histo_WTauNu->Write();
  //   histo_WTauNu_JESUp->Write();
  //   histo_WTauNu_JESDown->Write();
  //   histo_WTauNu_PESUp->Write();
  //   histo_WTauNu_PESDown->Write();
  //   histo_WW->Write();
  //   histo_WW_JESUp->Write();
  //   histo_WW_JESDown->Write();
  //   histo_WW_PESUp->Write();
  //   histo_WW_PESDown->Write();
  //   histo_WZG->Write();
  //   histo_WZG_JESUp->Write();
  //   histo_WZG_JESDown->Write();
  //   histo_WZG_PESUp->Write();
  //   histo_WZG_PESDown->Write();
  //   histo_WGG->Write();
  //   histo_WGG_JESUp->Write();
  //   histo_WGG_JESDown->Write();
  //   histo_WGG_PESUp->Write();
  //   histo_WGG_PESDown->Write();
  //   histo_ZGGToNuNuGG->Write();
  //   histo_ZGGToNuNuGG_JESUp->Write();
  //   histo_ZGGToNuNuGG_JESDown->Write();
  //   histo_ZGGToNuNuGG_PESUp->Write();
  //   histo_ZGGToNuNuGG_PESDown->Write();
  // }
  
  if (histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3")
  {
    TFile* f_ZnnG_histos;
    if(histname == "Photon_Et_range_3")
      f_ZnnG_histos = new TFile("ZnnG_histos_below0p5_Pt.root","RECREATE");
    else if(histname == "pfMET_3")
      f_ZnnG_histos = new TFile("ZnnG_histos_below0p5_MET.root","RECREATE");
    else if(histname == "Mt_3")
      f_ZnnG_histos = new TFile("ZnnG_histos_below0p5_Mt.root","RECREATE");
    f_ZnnG_histos->cd();
    histo_data->Write();
    histo_jetfake->Write();
    histo_jetfake_errUp->Write();
    histo_jetfake_errDown->Write();
    histo_spikes->Write();
    histo_elefake->Write();
    histo_elefake_errUp->Write();
    histo_elefake_errDown->Write();
    histo_bhalo->Write();
    histo_ZNuNuG->Write();
    histo_ZNuNuG_JESUp->Write();
    histo_ZNuNuG_JESDown->Write();
    histo_ZNuNuG_PESUp->Write();
    histo_ZNuNuG_PESDown->Write();
    histo_ZNuNuG_uncorrected->Write();
    histo_WG->Write();
    histo_WG_JESUp->Write();
    histo_WG_JESDown->Write();
    histo_WG_PESUp->Write();
    histo_WG_PESDown->Write();
    histo_WG_uncorrected->Write();
    histo_GJets_40toInf->Write();
    histo_GJets_40toInf_JESUp->Write();
    histo_GJets_40toInf_JESDown->Write();
    histo_GJets_40toInf_PESUp->Write();
    histo_GJets_40toInf_PESDown->Write();
    histo_ZllG_combined->Write();
    histo_ZllG_JESUp_combined->Write();
    histo_ZllG_JESDown_combined->Write();
    histo_ZllG_PESUp_combined->Write();
    histo_ZllG_PESDown_combined->Write();
    histo_ZllG_uncorrected_combined->Write();
    histo_TTG->Write();
    histo_TTG_JESUp->Write();
    histo_TTG_JESDown->Write();
    histo_TTG_PESUp->Write();
    histo_TTG_PESDown->Write();
    histo_TG->Write();
    histo_TG_JESUp->Write();
    histo_TG_JESDown->Write();
    histo_TG_PESUp->Write();
    histo_TG_PESDown->Write();
    histo_WWG->Write();
    histo_WWG_JESUp->Write();
    histo_WWG_JESDown->Write();
    histo_WWG_PESUp->Write();
    histo_WWG_PESDown->Write();
    histo_diphoton->Write();
    histo_diphoton_JESUp->Write();
    histo_diphoton_JESDown->Write();
    histo_diphoton_PESUp->Write();
    histo_diphoton_PESDown->Write();
    histo_WZ->Write();
    histo_WZ_JESUp->Write();
    histo_WZ_JESDown->Write();
    histo_WZ_PESUp->Write();
    histo_WZ_PESDown->Write();
    histo_ZZ->Write();
    histo_ZZ_JESUp->Write();
    histo_ZZ_JESDown->Write();
    histo_ZZ_PESUp->Write();
    histo_ZZ_PESDown->Write();
    histo_WMuNu->Write();
    histo_WMuNu_JESUp->Write();
    histo_WMuNu_JESDown->Write();
    histo_WMuNu_PESUp->Write();
    histo_WMuNu_PESDown->Write();
    histo_WTauNu->Write();
    histo_WTauNu_JESUp->Write();
    histo_WTauNu_JESDown->Write();
    histo_WTauNu_PESUp->Write();
    histo_WTauNu_PESDown->Write();
    histo_WW->Write();
    histo_WW_JESUp->Write();
    histo_WW_JESDown->Write();
    histo_WW_PESUp->Write();
    histo_WW_PESDown->Write();
    histo_WZG->Write();
    histo_WZG_JESUp->Write();
    histo_WZG_JESDown->Write();
    histo_WZG_PESUp->Write();
    histo_WZG_PESDown->Write();
    histo_WGG->Write();
    histo_WGG_JESUp->Write();
    histo_WGG_JESDown->Write();
    histo_WGG_PESUp->Write();
    histo_WGG_PESDown->Write();
    histo_ZGGToNuNuGG->Write();
    histo_ZGGToNuNuGG_JESUp->Write();
    histo_ZGGToNuNuGG_JESDown->Write();
    histo_ZGGToNuNuGG_PESUp->Write();
    histo_ZGGToNuNuGG_PESDown->Write();
    histo_DMVMx1Mv1000->Write();
    histo_DMVMx1Mv1000_JESUp->Write();
    histo_DMVMx1Mv1000_JESDown->Write();
    histo_DMVMx1Mv1000_PESUp->Write();
    histo_DMVMx1Mv1000_PESDown->Write();
    f_ZnnG_histos->Close();
  }
    
  if (histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "h_min_dphijetmet_2"){
    for(int i = 1; i <= nBins; i++){
      double binWidth = histo_data->GetBinWidth(i);
      if (histname == "h_min_dphijetmet_2")
        binWidth = binWidth*(1.0/0.25);
      histo_data->SetBinContent(i,histo_data->GetBinContent(i)/binWidth);
      histo_data->SetBinError(i,histo_data->GetBinError(i)/binWidth);
      histo_jetfake->SetBinContent(i,histo_jetfake->GetBinContent(i)/binWidth);
      histo_jetfake->SetBinError(i,histo_jetfake->GetBinError(i)/binWidth);
      histo_spikes->SetBinContent(i,histo_spikes->GetBinContent(i)/binWidth);
      histo_spikes->SetBinError(i,histo_spikes->GetBinError(i)/binWidth);
      histo_elefake->SetBinContent(i,histo_elefake->GetBinContent(i)/binWidth);
      histo_elefake->SetBinError(i,histo_elefake->GetBinError(i)/binWidth);
      histo_bhalo->SetBinContent(i,histo_bhalo->GetBinContent(i)/binWidth);
      histo_bhalo->SetBinError(i,histo_bhalo->GetBinError(i)/binWidth);
      histo_ZNuNuG->SetBinContent(i,histo_ZNuNuG->GetBinContent(i)/binWidth);
      histo_ZNuNuG->SetBinError(i,histo_ZNuNuG->GetBinError(i)/binWidth);
      histo_WG->SetBinContent(i,histo_WG->GetBinContent(i)/binWidth);
      histo_WG->SetBinError(i,histo_WG->GetBinError(i)/binWidth);
      histo_GJets_40toInf->SetBinContent(i,histo_GJets_40toInf->GetBinContent(i)/binWidth);
      histo_GJets_40toInf->SetBinError(i,histo_GJets_40toInf->GetBinError(i)/binWidth);
      histo_ZllG_combined->SetBinContent(i,histo_ZllG_combined->GetBinContent(i)/binWidth);
      histo_ZllG_combined->SetBinError(i,histo_ZllG_combined->GetBinError(i)/binWidth);
      histo_TTG->SetBinContent(i,histo_TTG->GetBinContent(i)/binWidth);
      histo_TTG->SetBinError(i,histo_TTG->GetBinError(i)/binWidth);
      histo_TG->SetBinContent(i,histo_TG->GetBinContent(i)/binWidth);
      histo_TG->SetBinError(i,histo_TG->GetBinError(i)/binWidth);
      histo_WWG->SetBinContent(i,histo_WWG->GetBinContent(i)/binWidth);
      histo_WWG->SetBinError(i,histo_WWG->GetBinError(i)/binWidth);
      histo_diphoton->SetBinContent(i,histo_diphoton->GetBinContent(i)/binWidth);
      histo_diphoton->SetBinError(i,histo_diphoton->GetBinError(i)/binWidth);
      histo_WZ->SetBinContent(i,histo_WZ->GetBinContent(i)/binWidth);
      histo_WZ->SetBinError(i,histo_WZ->GetBinError(i)/binWidth);
      histo_ZZ->SetBinContent(i,histo_ZZ->GetBinContent(i)/binWidth);
      histo_ZZ->SetBinError(i,histo_ZZ->GetBinError(i)/binWidth);
      histo_WMuNu->SetBinContent(i,histo_WMuNu->GetBinContent(i)/binWidth);
      histo_WMuNu->SetBinError(i,histo_WMuNu->GetBinError(i)/binWidth);
      histo_WTauNu->SetBinContent(i,histo_WTauNu->GetBinContent(i)/binWidth);
      histo_WTauNu->SetBinError(i,histo_WTauNu->GetBinError(i)/binWidth);
      histo_WW->SetBinContent(i,histo_WW->GetBinContent(i)/binWidth);
      histo_WW->SetBinError(i,histo_WW->GetBinError(i)/binWidth);
      histo_WZG->SetBinContent(i,histo_WZG->GetBinContent(i)/binWidth);
      histo_WZG->SetBinError(i,histo_WZG->GetBinError(i)/binWidth);
      histo_WGG->SetBinContent(i,histo_WGG->GetBinContent(i)/binWidth);
      histo_WGG->SetBinError(i,histo_WGG->GetBinError(i)/binWidth);
      histo_ZGGToNuNuGG->SetBinContent(i,histo_ZGGToNuNuGG->GetBinContent(i)/binWidth);
      histo_ZGGToNuNuGG->SetBinError(i,histo_ZGGToNuNuGG->GetBinError(i)/binWidth);
    }
  }
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.26,0.99,0.99);
  pad1->Draw(); pad1->cd();
  pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
  pad1->SetBottomMargin(0.);
  
  TH1F *histo_allbackgrounds = (TH1F*)histo_jetfake->Clone("histo_allbackgrounds");
  if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3")
    histo_allbackgrounds->Add(histo_spikes);
  histo_allbackgrounds->Add(histo_elefake);
  histo_allbackgrounds->Add(histo_bhalo);
  histo_allbackgrounds->Add(histo_ZNuNuG);
  histo_allbackgrounds->Add(histo_WG);
  histo_allbackgrounds->Add(histo_GJets_40toInf);
  histo_allbackgrounds->Add(histo_ZllG_combined);
  histo_allbackgrounds->Add(histo_TTG);
  histo_allbackgrounds->Add(histo_TG);
  histo_allbackgrounds->Add(histo_WWG);
  histo_allbackgrounds->Add(histo_diphoton);
  histo_allbackgrounds->Add(histo_WZ);
  histo_allbackgrounds->Add(histo_ZZ);
  histo_allbackgrounds->Add(histo_WMuNu);
  histo_allbackgrounds->Add(histo_WTauNu);
  histo_allbackgrounds->Add(histo_WW);
  histo_allbackgrounds->Add(histo_WZG);
  histo_allbackgrounds->Add(histo_WGG);
  histo_allbackgrounds->Add(histo_ZGGToNuNuGG);
  for(int i = 1; i <= nBins; i++){
    double background = histo_allbackgrounds->GetBinContent(i);
    // Add statistical errors
    double sum_binerrors_squared = 0.0;
    sum_binerrors_squared += pow(histo_jetfake->GetBinError(i),2);
    if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3")
      sum_binerrors_squared += pow(histo_spikes->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_elefake->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_bhalo->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_ZNuNuG->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WG->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_GJets_40toInf->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_ZllG_combined->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_TTG->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_TG->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WWG->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_diphoton->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WZ->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_ZZ->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WMuNu->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WTauNu->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WW->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WZG->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WGG->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_ZGGToNuNuGG->GetBinError(i),2);
    double binerror = sqrt(sum_binerrors_squared); // Include just the statistical error
    double jeserr = (fabs(jesup_shift[i-1])+fabs(jesdown_shift[i-1]))/2.0;
    double peserr = (fabs(pesup_shift[i-1])+fabs(pesdown_shift[i-1]))/2.0;
    double renerr = (fabs(renup_shift[i-1])+fabs(rendown_shift[i-1]))/2.0;
    double facerr = (fabs(facup_shift[i-1])+fabs(facdown_shift[i-1]))/2.0;
    double pdferr = (fabs(pdfup_shift[i-1])+fabs(pdfdown_shift[i-1]))/2.0;
    double xsecerr_ZNuNuG = fabs(xsec_shift_ZNuNuG[i-1]);
    double xsecerr_WG = fabs(xsec_shift_WG[i-1]);
    double xsecerr_ZLLG = fabs(xsec_shift_ZLLG[i-1]);
    double jetfakeerr = (fabs(syst_shiftUp_jetfake[i-1])+fabs(syst_shiftDown_jetfake[i-1]))/2.0;
    double elefakeerr = (fabs(syst_shiftUp_elefake[i-1])+fabs(syst_shiftDown_elefake[i-1]))/2.0;
    if (histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "h_min_dphijetmet_2"){
      double binWidth = histo_data->GetBinWidth(i);
      if (histname == "h_min_dphijetmet_2")
        binWidth = binWidth*(1.0/0.25);
      jeserr /= binWidth;
      peserr /= binWidth;
      renerr /= binWidth;
      facerr /= binWidth;
      pdferr /= binWidth;
      xsecerr_ZNuNuG /= binWidth;
      xsecerr_WG /= binWidth;
      xsecerr_ZLLG /= binWidth;
      jetfakeerr /= binWidth;
      elefakeerr /= binWidth;
    }
    binerror = sqrt(sum_binerrors_squared+pow(background*photon_scale_factor_unc,2)+pow(jeserr,2)+pow(peserr,2)+pow(renerr,2)+pow(facerr,2)+pow(pdferr,2)+pow(xsecerr_ZNuNuG,2)+pow(xsecerr_WG,2)+pow(xsecerr_ZLLG,2)+pow(jetfakeerr,2)+pow(elefakeerr,2));
    histo_allbackgrounds->SetBinError(i,binerror);
  }
  histo_allbackgrounds->SetFillColorAlpha(kGray+1,0.6);
  histo_vector.push_back(histo_allbackgrounds);
  
  TH1F *histo_allbackgrounds_outline = (TH1F*)histo_allbackgrounds->Clone("histo_allbackgrounds_outline");
  histo_allbackgrounds_outline->SetFillColorAlpha(kWhite,0.0);
  histo_allbackgrounds_outline->SetLineWidth(1);
  histo_vector.push_back(histo_allbackgrounds_outline);
  
  // Only necessary until unblinding
  if (histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3")
  {
    TH1F *histo_expectedData = (TH1F*)histo_allbackgrounds->Clone("data_obs");
    if (histname == "Photon_Et_range_3" || histname == "pfMET_3"){
      for(int i = 1; i <= nBins; i++){
        double binWidth = histo_expectedData->GetBinWidth(i);
        Float_t adjusted_bin_content = histo_expectedData->GetBinContent(i)*binWidth;
        histo_expectedData->SetBinContent(i,adjusted_bin_content);
        histo_expectedData->SetBinError(i,sqrt(adjusted_bin_content));
      }
    }
    else{
      for(int i = 1; i <= nBins; i++){
        histo_expectedData->SetBinError(i,sqrt(histo_expectedData->GetBinContent(i)));
      }
    }
    TFile* f_ZnnG_histos;
    if(histname == "Photon_Et_range_3")
      f_ZnnG_histos = new TFile("ZnnG_histos_below0p5_Pt.root","UPDATE");
    else if(histname == "pfMET_3")
      f_ZnnG_histos = new TFile("ZnnG_histos_below0p5_MET.root","UPDATE");
    else if(histname == "Mt_3")
      f_ZnnG_histos = new TFile("ZnnG_histos_below0p5_Mt.root","UPDATE");
    f_ZnnG_histos->cd();
    histo_expectedData->Write();
    f_ZnnG_histos->Close();
  }
  
  histo_WZ->Add(histo_ZllG_combined);
  histo_WZ->Add(histo_TG);
  histo_WZ->Add(histo_WWG);
  histo_WZ->Add(histo_ZZ);
  histo_WZ->Add(histo_WW);
  histo_WZ->Add(histo_WZG);
  histo_WZ->Add(histo_WGG);
  histo_WZ->Add(histo_ZGGToNuNuGG);
    
  THStack *stackHisto = new THStack("stackHisto","Title");
  // stackHisto->Add(histo_ZZ);  // Already in WZ
  // stackHisto->Add(histo_WWG);  // Already in WZ
  // stackHisto->Add(histo_TG);  // Already in WZ
  // stackHisto->Add(histo_ZllG_combined);  // Already in WZ
  stackHisto->Add(histo_WZ);
  stackHisto->Add(histo_bhalo);
  stackHisto->Add(histo_TTG);
  stackHisto->Add(histo_WTauNu);
  if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3")
    stackHisto->Add(histo_spikes);
  stackHisto->Add(histo_diphoton);
  stackHisto->Add(histo_WMuNu);
  stackHisto->Add(histo_jetfake);
  stackHisto->Add(histo_GJets_40toInf);
  stackHisto->Add(histo_elefake);
  stackHisto->Add(histo_WG);
  stackHisto->Add(histo_ZNuNuG);
  stackHisto->SetTitle("");
  
  for(int i = 0; i < int(histo_vector.size()); i++){
    histo_vector[i]->SetStats(0);
    histo_vector[i]->SetTitle("");
    histo_vector[i]->SetLineColor(kBlack);
    histo_vector[i]->GetXaxis()->SetTitle(xaxis_title);
    histo_vector[i]->GetXaxis()->SetLabelFont(42);
    histo_vector[i]->GetXaxis()->SetLabelSize(0.06);
    histo_vector[i]->GetXaxis()->SetTitleFont(42);
    histo_vector[i]->GetXaxis()->SetTitleSize(0.06);
    histo_vector[i]->GetYaxis()->SetTitle("Events / bin");
    if (histname == "Photon_Et_range_3" || histname == "pfMET_3")
      histo_vector[i]->GetYaxis()->SetTitle("Events / GeV");
    else if (histname == "h_min_dphijetmet_2")
      histo_vector[i]->GetYaxis()->SetTitle("Events / 0.25 radians");
    histo_vector[i]->GetYaxis()->SetLabelFont(42);
    histo_vector[i]->GetYaxis()->SetLabelSize(0.06);
    histo_vector[i]->GetYaxis()->SetTitleFont(42);
    histo_vector[i]->GetYaxis()->SetTitleSize(0.06);
    histo_vector[i]->GetYaxis()->SetTitleOffset(0.9);
  }
  
  //Accommodate both the data and background plots
  double ymax_data = 0.0;
  double ymax_background = 0.0;
  for(int i = 1; i <= nBins; i++){
    double y_data = histo_data->GetBinContent(i);
    double y_error_data = histo_data->GetBinError(i);
    double y_high_data = y_data+y_error_data;
    if(y_high_data > ymax_data)
      ymax_data = y_high_data;
    double y_background = histo_allbackgrounds->GetBinContent(i);
    double y_error_background = histo_allbackgrounds->GetBinError(i);
    double y_high_background = y_background+y_error_background;
    if(y_high_background > ymax_background)
      ymax_background = y_high_background;
  }
  
  double ymin = 0.0003;
  double ymax = 1.7*ymax_data;
  if(ymax_background > ymax_data)
    ymax = 1.7*ymax_background;
  if (histname == "Photon_Et_range_3" || histname == "pfMET_3"){
    pad1->SetLogy();
    ymax *= 5;
  }
  histo_data->GetYaxis()->SetRangeUser(ymin,ymax);
  if (histname == "nJet_3")
    histo_data->GetXaxis()->SetRangeUser(0,10);
  histo_data->Draw();
  stackHisto->Draw("HIST SAME");
  histo_allbackgrounds->Draw("E2 SAME");
  histo_allbackgrounds_outline->Draw("HIST SAME");
  histo_data->SetLineColor(kBlack);
  histo_data->SetMarkerColor(kBlack);
  histo_data->Draw("E0 P0 SAME");
  gPad->RedrawAxis();
  
  //Central location of leg defined to be location of leg in phoPt plot
  TLegend* leg = new TLegend(0.40+leg_xoffset,0.46075+leg_yoffset,0.855387+leg_xoffset,0.862969+leg_yoffset,"");
  leg->AddEntry(histo_data,"Data");
  leg->AddEntry(histo_ZNuNuG,"Z(#nu#nu)#gamma","F");
  leg->AddEntry(histo_WG,"W#gamma#rightarrowl#nu#gamma","F");
  leg->AddEntry(histo_elefake,"e#rightarrow#gamma MisID","F");
  leg->AddEntry(histo_GJets_40toInf,"#gamma+jet","F");
  leg->AddEntry(histo_jetfake,"jet#rightarrow#gamma MisID","F");
  leg->AddEntry(histo_WMuNu,"W(#mu#nu)","F");
  leg->AddEntry(histo_diphoton,"#gamma#gamma","F");
  if(histname == "Photon_Et_range_3" || histname == "pfMET_3" || histname == "Mt_3")
    leg->AddEntry(histo_spikes,"Spikes","F");
  leg->AddEntry(histo_WTauNu,"W(#tau#nu)","F");
  leg->AddEntry(histo_TTG,"tt#gamma","F");
  leg->AddEntry(histo_bhalo,"Beam halo","F");
  leg->AddEntry(histo_WZ,"Diboson,Triboson,t#gamma","F");
  // leg->AddEntry(histo_WZ,"WZ","F");
  // leg->AddEntry(histo_ZllG_combined,"Z(ll)#gamma","F");
  // leg->AddEntry(histo_TG,"t#gamma","F");
  // leg->AddEntry(histo_WWG,"WW#gamma","F");
  // leg->AddEntry(histo_ZZ,"ZZ","F");
  leg->SetNColumns(2);
  leg->SetFillColor(kWhite);
  leg->SetShadowColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.040);
  leg->Draw();
  
  float lumiTextSize = 0.6;
  float lumiTextOffset = 0.2;
  float cmsTextSize = 0.75;
  TLatex *texS = new TLatex(0.60023,0.917173,"36.8 fb^{-1} (13 TeV)");
  texS->SetNDC();
  texS->SetTextFont(42);
  texS->SetTextSize(lumiTextSize*t_m);
  texS->Draw();
  TLatex *texS1 = new TLatex(0.13592,0.817173,"#bf{CMS} #it{Preliminary}");
  texS1->SetNDC();
  texS1->SetTextFont(42);
  texS1->SetTextSize(cmsTextSize*t_m);
  texS1->Draw();
  
  c->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.26);
  pad2->Draw(); pad2->cd();
  pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.35);
  
  double max_ratio = 2.5;
  
  TH1F* Ratio = (TH1F*)histo_data->Clone("Ratio");
  TH1F* Ratio_background = (TH1F*)histo_allbackgrounds->Clone("Ratio_background");
  for(int i = 1; i <= nBins; i++){
    double y_data = histo_data->GetBinContent(i);
    double y_error_data = histo_data->GetBinError(i);
    double y_background = histo_allbackgrounds->GetBinContent(i);
    double y_error_background = histo_allbackgrounds->GetBinError(i);
    double Ratiocontent = 0.0;
    double Ratioerror = max_ratio;
    double Ratioerror_background = max_ratio;
    if(y_background > 0.){
      Ratiocontent = y_data/y_background;
      Ratioerror_background = y_error_background/y_background;
      if(y_error_data > 0.)
        Ratioerror = y_error_data/y_background;
    }
    else if(y_data > 0.){
      Ratiocontent = 3.*max_ratio;
    }
    Ratio->SetBinContent(i,Ratiocontent);
    Ratio->SetBinError(i,Ratioerror);
    Ratio_background->SetBinContent(i,1);
    Ratio_background->SetBinError(i,Ratioerror_background);
  }
  
  Ratio_background->GetYaxis()->SetRangeUser(0.0,max_ratio-0.01);
  Ratio_background->GetYaxis()->SetTitle("Data/SM");
  Ratio_background->GetYaxis()->CenterTitle();
  Ratio_background->GetYaxis()->SetLabelSize(0.14);
  Ratio_background->GetYaxis()->SetTitleSize(0.15);
  Ratio_background->GetYaxis()->SetLabelFont(42);
  Ratio_background->GetYaxis()->SetTitleFont(42);
  Ratio_background->GetYaxis()->SetTitleOffset(0.30);
  Ratio_background->GetYaxis()->SetNdivisions(305);
  Ratio_background->GetXaxis()->SetTitle(xaxis_title);
  Ratio_background->GetXaxis()->SetLabelSize(0.16);
  Ratio_background->GetXaxis()->SetTitleSize(0.18);
  Ratio_background->GetXaxis()->SetLabelFont(42);
  Ratio_background->GetXaxis()->SetTitleFont(42);
  Ratio_background->GetXaxis()->SetTitleOffset(0.9);
  Ratio_background->GetXaxis()->SetTickLength(0.05);
  Ratio_background->SetStats(0);
  Ratio->SetMarkerStyle(0);
  double xmin = histo_data->GetXaxis()->GetBinLowEdge(1);
  double xmax = histo_data->GetXaxis()->GetBinUpEdge(nBins);
  TLine* line = new TLine(xmin,1.,xmax,1.);
  line->SetLineStyle(2);
  line->SetLineColor(kBlack);
  gStyle->SetLineStyleString(11,"3 12");
  TLine* line1 = new TLine(xmin,1.5,xmax,1.5);
  line1->SetLineStyle(11);
  line1->SetLineColor(kBlack);
  TLine* line2 = new TLine(xmin,2.0,xmax,2.0);
  line2->SetLineStyle(11);
  line2->SetLineColor(kBlack);
  TLine* line3 = new TLine(xmin,0.5,xmax,0.5);
  line3->SetLineStyle(11);
  line3->SetLineColor(kBlack);
  TLine* line4 = new TLine(xmin,2.5,xmax,2.5);
  line4->SetLineStyle(11);
  line4->SetLineColor(kBlack);
  TLine* line5 = new TLine(xmin,3.0,xmax,3.0);
  line5->SetLineStyle(11);
  line5->SetLineColor(kBlack);
  Ratio_background->Draw("E2");
  line->Draw("SAME");
  line1->Draw("SAME");
  line2->Draw("SAME");
  line3->Draw("SAME");
  line4->Draw("SAME");
  line5->Draw("SAME");
  Ratio->Draw("E0 P0 SAME");

  double background_unc = sqrt(background_unc_sumsquares);
  if(histname == "Photon_Et_range_3"){
    cout<<"ZnnG region (above0p5)"<<endl;
    cout<<"------------------------------------"<<endl;
    cout<<"Jet faking photon: "<<int_jetfake<<" +- "<<err_jetfake<<endl;
    cout<<"Spikes: "<<int_spikes<<" +- "<<err_spikes<<endl;
    cout<<"Electron faking photon: "<<int_elefake<<" +- "<<err_elefake<<endl;
    cout<<"Beam halo: "<<int_bhalo<<" +- "<<err_bhalo<<endl;
    cout<<"ZNuNu+gamma: "<<int_ZNuNuG<<" +- "<<err_ZNuNuG<<endl;
    cout<<"W+gamma: "<<int_WG<<" +- "<<err_WG<<endl;
    cout<<"GJets: "<<int_sum_GJets<<" +- "<<err_sum_GJets<<endl;
    cout<<"Z(ll)+Gamma: "<<int_ZllG<<" +- "<<err_ZllG<<endl;
    cout<<"tt+Gamma: "<<int_TTG<<" +- "<<err_TTG<<endl;
    cout<<"t+Gamma: "<<int_TG<<" +- "<<err_TG<<endl;
    cout<<"WWG: "<<int_WWG<<" +- "<<err_WWG<<endl;
    cout<<"Diphoton: "<<int_diphoton<<" +- "<<err_diphoton<<endl;
    cout<<"WZ: "<<int_WZ<<" +- "<<err_WZ<<endl;
    cout<<"ZZ: "<<int_ZZ<<" +- "<<err_ZZ<<endl;
    cout<<"WMuNu: "<<int_WMuNu<<" +- "<<err_WMuNu<<endl;
    cout<<"WTauNu: "<<int_WTauNu<<" +- "<<err_WTauNu<<endl;
    cout<<"WW: "<<int_WW<<" +- "<<err_WW<<endl;
    cout<<"WZG: "<<int_WZG<<" +- "<<err_WZG<<endl;
    cout<<"WGG: "<<int_WGG<<" +- "<<err_WGG<<endl;
    cout<<"ZGGToNuNuGG: "<<int_ZGGToNuNuGG<<" +- "<<err_ZGGToNuNuGG<<endl;
    cout<<"Total background: "<<total_background<<" +- "<<background_unc<<endl;
    cout<<"------------------------------------"<<endl;
    cout<<"Data: "<<int_data<<endl;
    cout<<"------------------------------------"<<endl;
  }

  c->SaveAs(TString("znng_below0p5_"+plotname+".png"));
  c->SaveAs(TString("znng_below0p5_"+plotname+".pdf"));
  delete(c);
}

void znng_below0p5_plotter()
{
  std::vector<string> histnames;
  histnames.clear();
  std::vector<Double_t> leg_xoffsets;
  leg_xoffsets.clear();
  std::vector<Double_t> leg_yoffsets;
  leg_yoffsets.clear();
  std::vector<TString> xaxis_titles;
  xaxis_titles.clear();
  std::vector<TString> plotnames;
  plotnames.clear();

//  histnames.push_back(TString("_4"));
//  leg_xoffsets.push_back(0.);
//  leg_yoffsets.push_back(0.);
//  xaxis_titles.push_back(TString(""));
//  plotnames.push_back(TString(""));

  histnames.push_back("Photon_Et_range");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("phoPt"));

  histnames.push_back("Photon_SCeta");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #eta"));
  plotnames.push_back(TString("phoEta"));

  histnames.push_back("Photon_SCphi");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #phi"));
  plotnames.push_back(TString("phoPhi"));

  histnames.push_back("pfMET");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("pfMET [GeV]"));
  plotnames.push_back(TString("pfMET"));

  histnames.push_back("nJet");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Number of Jets"));
  plotnames.push_back(TString("nJet"));

  histnames.push_back("h_dPhi");
  leg_xoffsets.push_back(-0.2);
  leg_yoffsets.push_back(-0.1);
  xaxis_titles.push_back(TString("#Delta#phi(Photon,MET)"));
  plotnames.push_back(TString("dPhiPhoMET"));

  histnames.push_back("PTMET");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} / MET"));
  plotnames.push_back(TString("phoPtOverMET"));
  
  histnames.push_back("METoverSqrtSumEt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("MET / #sqrt{#SigmaE_{T}}"));
  plotnames.push_back(TString("METoverSumET"));

  histnames.push_back("Mt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon-MET M_{T} [GeV]"));
  plotnames.push_back(TString("phoMETmT"));

  histnames.push_back("h_min_dphijetmet");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Min. #Delta#phi(jets,MET)"));
  plotnames.push_back(TString("dPhiJetsMET"));

  for(int i = 0; i < histnames.size(); i++){
    plot(histnames[i],leg_xoffsets[i],leg_yoffsets[i],xaxis_titles[i],plotnames[i]);
  }
  
  TFile* f_DM_NLO_histos_Pt = new TFile("DM_NLO_histos_below0p5_Pt.root", "RECREATE");
  f_DM_NLO_histos_Pt->Close();
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-1000_Mv-1000", 3.954e-06, 20933, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-10_Mv-10000", 3.312e-07, 20411, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-10_Mv-100", 4.708e+00, 26007, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-10_Mv-50", 1.174e+01, 23291, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-140_Mv-300", 1.095e-01, 20166, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-150_Mv-10000", 2.777e-07, 19489, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-150_Mv-1000", 4.016e-02, 20755, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-150_Mv-500", 2.122e-01, 22092, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-1_Mv-10000", 3.297e-07, 20859, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-1_Mv-1000", 4.362e-02, 20910, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-1_Mv-100", 4.845e+00, 22224, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-1_Mv-2000", 2.900e-03, 19762, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-1_Mv-200", 1.957e+00, 23430, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-1_Mv-300", 9.612e-01, 22381, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-1_Mv-500", 3.112e-01, 22192, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-1_Mv-50", 1.200e+01, 25786, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-500_Mv-10000", 1.194e-07, 19460, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-500_Mv-2000", 2.082e-03, 20529, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-500_Mv-500", 1.248e-04, 21015, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-50_Mv-10000", 3.180e-07, 18047, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-50_Mv-10", 5.034e-02, 15429, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-50_Mv-200", 1.555e+00, 23682, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-50_Mv-300", 8.902e-01, 22688, true);
  plot_DM("Photon_Et_range", "DM_NLO_Axial_Mx-990_Mv-2000", 5.658e-05, 19406, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-0_Mv-20", 4.822e-01, 24769, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1000_Mv-1000", 1.716e-05, 21196, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1000_Mv-10", 1.142e-05, 22629, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-10_Mv-10000", 3.269e-07, 20208, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-10_Mv-100", 4.017e+00, 23231, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-10_Mv-10", 3.615e-01, 23455, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-10_Mv-50", 5.298e+00, 23448, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-140_Mv-300", 6.212e-01, 23056, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-150_Mv-10000", 3.082e-07, 21289, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-150_Mv-1000", 4.282e-02, 22100, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-150_Mv-10", 1.098e-02, 19876, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-150_Mv-200", 1.904e-02, 20851, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-150_Mv-500", 2.870e-01, 24272, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1_Mv-10000", 3.197e-07, 21054, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1_Mv-1000", 4.321e-02, 20151, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1_Mv-100", 4.021e+00, 21924, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1_Mv-2000", 2.852e-03, 19199, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1_Mv-200", 1.839e+00, 20045, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1_Mv-20", 4.861e-01, 22678, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1_Mv-300", 9.479e-01, 22845, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1_Mv-500", 2.999e-01, 22704, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-1_Mv-50", 5.294e+00, 24229, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-40_Mv-100", 3.537e+00, 23175, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-490_Mv-1000", 1.723e-02, 21617, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-500_Mv-10000", 1.892e-07, 20064, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-500_Mv-10", 2.885e-04, 20860, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-500_Mv-2000", 2.727e-03, 21807, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-500_Mv-500", 4.098e-04, 25445, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-50_Mv-10000", 3.151e-07, 21446, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-50_Mv-10", 8.874e-02, 23590, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-50_Mv-200", 1.830e+00, 22112, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-50_Mv-300", 9.407e-01, 22053, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-50_Mv-50", 1.074e-01, 23703, true);
  plot_DM("Photon_Et_range", "DM_NLO_Vector_Mx-990_Mv-2000", 8.433e-04, 19573, true);
  
  TFile* f_DM_NLO_histos_MET = new TFile("DM_NLO_histos_below0p5_MET.root", "RECREATE");
  f_DM_NLO_histos_MET->Close();
  plot_DM("pfMET", "DM_NLO_Axial_Mx-1000_Mv-1000", 3.954e-06, 20933, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-10_Mv-10000", 3.312e-07, 20411, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-10_Mv-100", 4.708e+00, 26007, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-10_Mv-50", 1.174e+01, 23291, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-140_Mv-300", 1.095e-01, 20166, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-150_Mv-10000", 2.777e-07, 19489, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-150_Mv-1000", 4.016e-02, 20755, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-150_Mv-500", 2.122e-01, 22092, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-1_Mv-10000", 3.297e-07, 20859, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-1_Mv-1000", 4.362e-02, 20910, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-1_Mv-100", 4.845e+00, 22224, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-1_Mv-2000", 2.900e-03, 19762, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-1_Mv-200", 1.957e+00, 23430, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-1_Mv-300", 9.612e-01, 22381, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-1_Mv-500", 3.112e-01, 22192, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-1_Mv-50", 1.200e+01, 25786, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-500_Mv-10000", 1.194e-07, 19460, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-500_Mv-2000", 2.082e-03, 20529, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-500_Mv-500", 1.248e-04, 21015, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-50_Mv-10000", 3.180e-07, 18047, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-50_Mv-10", 5.034e-02, 15429, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-50_Mv-200", 1.555e+00, 23682, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-50_Mv-300", 8.902e-01, 22688, true);
  plot_DM("pfMET", "DM_NLO_Axial_Mx-990_Mv-2000", 5.658e-05, 19406, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-0_Mv-20", 4.822e-01, 24769, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1000_Mv-1000", 1.716e-05, 21196, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1000_Mv-10", 1.142e-05, 22629, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-10_Mv-10000", 3.269e-07, 20208, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-10_Mv-100", 4.017e+00, 23231, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-10_Mv-10", 3.615e-01, 23455, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-10_Mv-50", 5.298e+00, 23448, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-140_Mv-300", 6.212e-01, 23056, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-150_Mv-10000", 3.082e-07, 21289, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-150_Mv-1000", 4.282e-02, 22100, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-150_Mv-10", 1.098e-02, 19876, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-150_Mv-200", 1.904e-02, 20851, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-150_Mv-500", 2.870e-01, 24272, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1_Mv-10000", 3.197e-07, 21054, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1_Mv-1000", 4.321e-02, 20151, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1_Mv-100", 4.021e+00, 21924, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1_Mv-2000", 2.852e-03, 19199, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1_Mv-200", 1.839e+00, 20045, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1_Mv-20", 4.861e-01, 22678, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1_Mv-300", 9.479e-01, 22845, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1_Mv-500", 2.999e-01, 22704, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-1_Mv-50", 5.294e+00, 24229, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-40_Mv-100", 3.537e+00, 23175, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-490_Mv-1000", 1.723e-02, 21617, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-500_Mv-10000", 1.892e-07, 20064, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-500_Mv-10", 2.885e-04, 20860, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-500_Mv-2000", 2.727e-03, 21807, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-500_Mv-500", 4.098e-04, 25445, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-50_Mv-10000", 3.151e-07, 21446, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-50_Mv-10", 8.874e-02, 23590, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-50_Mv-200", 1.830e+00, 22112, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-50_Mv-300", 9.407e-01, 22053, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-50_Mv-50", 1.074e-01, 23703, true);
  plot_DM("pfMET", "DM_NLO_Vector_Mx-990_Mv-2000", 8.433e-04, 19573, true);

  TFile* f_DM_NLO_histos_Mt = new TFile("DM_NLO_histos_below0p5_Mt.root", "RECREATE");
  f_DM_NLO_histos_Mt->Close();
  plot_DM("Mt", "DM_NLO_Axial_Mx-1000_Mv-1000", 3.954e-06, 20933, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-10_Mv-10000", 3.312e-07, 20411, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-10_Mv-100", 4.708e+00, 26007, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-10_Mv-50", 1.174e+01, 23291, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-140_Mv-300", 1.095e-01, 20166, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-150_Mv-10000", 2.777e-07, 19489, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-150_Mv-1000", 4.016e-02, 20755, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-150_Mv-500", 2.122e-01, 22092, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-1_Mv-10000", 3.297e-07, 20859, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-1_Mv-1000", 4.362e-02, 20910, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-1_Mv-100", 4.845e+00, 22224, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-1_Mv-2000", 2.900e-03, 19762, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-1_Mv-200", 1.957e+00, 23430, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-1_Mv-300", 9.612e-01, 22381, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-1_Mv-500", 3.112e-01, 22192, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-1_Mv-50", 1.200e+01, 25786, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-500_Mv-10000", 1.194e-07, 19460, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-500_Mv-2000", 2.082e-03, 20529, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-500_Mv-500", 1.248e-04, 21015, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-50_Mv-10000", 3.180e-07, 18047, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-50_Mv-10", 5.034e-02, 15429, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-50_Mv-200", 1.555e+00, 23682, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-50_Mv-300", 8.902e-01, 22688, true);
  plot_DM("Mt", "DM_NLO_Axial_Mx-990_Mv-2000", 5.658e-05, 19406, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-0_Mv-20", 4.822e-01, 24769, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1000_Mv-1000", 1.716e-05, 21196, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1000_Mv-10", 1.142e-05, 22629, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-10_Mv-10000", 3.269e-07, 20208, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-10_Mv-100", 4.017e+00, 23231, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-10_Mv-10", 3.615e-01, 23455, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-10_Mv-50", 5.298e+00, 23448, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-140_Mv-300", 6.212e-01, 23056, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-150_Mv-10000", 3.082e-07, 21289, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-150_Mv-1000", 4.282e-02, 22100, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-150_Mv-10", 1.098e-02, 19876, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-150_Mv-200", 1.904e-02, 20851, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-150_Mv-500", 2.870e-01, 24272, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1_Mv-10000", 3.197e-07, 21054, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1_Mv-1000", 4.321e-02, 20151, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1_Mv-100", 4.021e+00, 21924, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1_Mv-2000", 2.852e-03, 19199, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1_Mv-200", 1.839e+00, 20045, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1_Mv-20", 4.861e-01, 22678, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1_Mv-300", 9.479e-01, 22845, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1_Mv-500", 2.999e-01, 22704, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-1_Mv-50", 5.294e+00, 24229, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-40_Mv-100", 3.537e+00, 23175, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-490_Mv-1000", 1.723e-02, 21617, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-500_Mv-10000", 1.892e-07, 20064, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-500_Mv-10", 2.885e-04, 20860, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-500_Mv-2000", 2.727e-03, 21807, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-500_Mv-500", 4.098e-04, 25445, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-50_Mv-10000", 3.151e-07, 21446, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-50_Mv-10", 8.874e-02, 23590, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-50_Mv-200", 1.830e+00, 22112, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-50_Mv-300", 9.407e-01, 22053, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-50_Mv-50", 1.074e-01, 23703, true);
  plot_DM("Mt", "DM_NLO_Vector_Mx-990_Mv-2000", 8.433e-04, 19573, true);
  
  TFile* f_DM_LO_histos_Pt = new TFile("DM_LO_histos_below0p5_Pt.root", "RECREATE");
  f_DM_LO_histos_Pt->Close();
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1000_Mv-10000", 1.572e-08, 49997, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1000_Mv-1000", 1.830e-06, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1000_Mv-10", 1.344e-06, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1000_Mv-1995", 1.649e-05, 14198, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-10_Mv-10000", 1.217e-07, 49998, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-10_Mv-100", 3.776e-01, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-10_Mv-10", 2.898e-02, 800, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-10_Mv-15", 3.174e-02, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-10_Mv-50", 3.952e-01, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-150_Mv-10000", 1.076e-07, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-150_Mv-1000", 1.349e-02, 26800, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-150_Mv-10", 1.284e-03, 49998, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-150_Mv-200", 1.800e-03, 49200, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-150_Mv-295", 4.456e-03, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-150_Mv-500", 5.312e-02, 49997, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1_Mv-10000", 1.218e-07, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1_Mv-1000", 1.465e-02, 49997, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1_Mv-100", 3.876e-01, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1_Mv-10", 5.202e-01, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1_Mv-2000", 1.213e-03, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1_Mv-200", 2.681e-01, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1_Mv-20", 4.870e-01, 49198, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1_Mv-300", 1.783e-01, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-1_Mv-50", 4.468e-01, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-500_Mv-10000", 5.267e-08, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-500_Mv-10", 3.674e-05, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-500_Mv-2000", 8.975e-04, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-500_Mv-500", 4.725e-05, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-500_Mv-995", 2.678e-04, 48799, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-50_Mv-10000", 1.198e-07, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-50_Mv-10", 7.734e-03, 2400, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-50_Mv-200", 2.144e-01, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-50_Mv-300", 1.631e-01, 48798, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-50_Mv-50", 8.656e-03, 49998, false);
  plot_DM("Photon_Et_range", "DM_LO_AV_Mx-50_Mv-95", 1.482e-02, 49198, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1000_Mv-10000", 3.499e-08, 48399, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1000_Mv-1000", 7.679e-06, 49196, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1000_Mv-10", 5.262e-06, 4000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1000_Mv-1995", 2.533e-04, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-10_Mv-10000", 1.216e-07, 49998, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-10_Mv-100", 3.885e-01, 49998, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-10_Mv-10", 3.781e-02, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-10_Mv-15", 4.527e-02, 49199, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-10_Mv-50", 4.422e-01, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-150_Mv-10000", 1.195e-07, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-150_Mv-1000", 1.437e-02, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-150_Mv-10", 2.668e-03, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-150_Mv-200", 4.398e-03, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-150_Mv-295", 2.901e-02, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-150_Mv-500", 7.160e-02, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-10000", 1.222e-07, 49399, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-1000", 1.446e-02, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-100", 3.877e-01, 41000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-10", 5.114e-01, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-2000", 1.207e-03, 47599, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-200", 2.683e-01, 49997, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-20", 4.752e-01, 49998, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-300", 1.785e-01, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-500", 7.464e-02, 49998, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-1_Mv-50", 4.451e-01, 49400, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-500_Mv-10000", 8.307e-08, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-500_Mv-10", 1.087e-04, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-500_Mv-2000", 1.152e-03, 50000, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-500_Mv-500", 1.510e-04, 48400, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-500_Mv-995", 3.136e-03, 47600, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-50_Mv-10000", 1.223e-07, 49998, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-50_Mv-200", 2.646e-01, 49999, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-50_Mv-300", 1.776e-01, 49997, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-50_Mv-50", 1.464e-02, 49998, false);
  plot_DM("Photon_Et_range", "DM_LO_V_Mx-50_Mv-95", 4.502e-02, 50000, false);

  TFile* f_DM_LO_histos_MET = new TFile("DM_LO_histos_below0p5_MET.root", "RECREATE");
  f_DM_LO_histos_MET->Close();
  plot_DM("pfMET", "DM_LO_AV_Mx-1000_Mv-10000", 1.572e-08, 49997, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1000_Mv-1000", 1.830e-06, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1000_Mv-10", 1.344e-06, 49999, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1000_Mv-1995", 1.649e-05, 14198, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-10_Mv-10000", 1.217e-07, 49998, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-10_Mv-100", 3.776e-01, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-10_Mv-10", 2.898e-02, 800, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-10_Mv-15", 3.174e-02, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-10_Mv-50", 3.952e-01, 49999, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-150_Mv-10000", 1.076e-07, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-150_Mv-1000", 1.349e-02, 26800, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-150_Mv-10", 1.284e-03, 49998, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-150_Mv-200", 1.800e-03, 49200, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-150_Mv-295", 4.456e-03, 49999, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-150_Mv-500", 5.312e-02, 49997, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1_Mv-10000", 1.218e-07, 49999, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1_Mv-1000", 1.465e-02, 49997, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1_Mv-100", 3.876e-01, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1_Mv-10", 5.202e-01, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1_Mv-2000", 1.213e-03, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1_Mv-200", 2.681e-01, 49999, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1_Mv-20", 4.870e-01, 49198, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1_Mv-300", 1.783e-01, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-1_Mv-50", 4.468e-01, 49999, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-500_Mv-10000", 5.267e-08, 49999, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-500_Mv-10", 3.674e-05, 49999, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-500_Mv-2000", 8.975e-04, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-500_Mv-500", 4.725e-05, 49999, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-500_Mv-995", 2.678e-04, 48799, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-50_Mv-10000", 1.198e-07, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-50_Mv-10", 7.734e-03, 2400, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-50_Mv-200", 2.144e-01, 50000, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-50_Mv-300", 1.631e-01, 48798, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-50_Mv-50", 8.656e-03, 49998, false);
  plot_DM("pfMET", "DM_LO_AV_Mx-50_Mv-95", 1.482e-02, 49198, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1000_Mv-10000", 3.499e-08, 48399, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1000_Mv-1000", 7.679e-06, 49196, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1000_Mv-10", 5.262e-06, 4000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1000_Mv-1995", 2.533e-04, 50000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-10_Mv-10000", 1.216e-07, 49998, false);
  plot_DM("pfMET", "DM_LO_V_Mx-10_Mv-100", 3.885e-01, 49998, false);
  plot_DM("pfMET", "DM_LO_V_Mx-10_Mv-10", 3.781e-02, 50000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-10_Mv-15", 4.527e-02, 49199, false);
  plot_DM("pfMET", "DM_LO_V_Mx-10_Mv-50", 4.422e-01, 49999, false);
  plot_DM("pfMET", "DM_LO_V_Mx-150_Mv-10000", 1.195e-07, 49999, false);
  plot_DM("pfMET", "DM_LO_V_Mx-150_Mv-1000", 1.437e-02, 49999, false);
  plot_DM("pfMET", "DM_LO_V_Mx-150_Mv-10", 2.668e-03, 50000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-150_Mv-200", 4.398e-03, 49999, false);
  plot_DM("pfMET", "DM_LO_V_Mx-150_Mv-295", 2.901e-02, 49999, false);
  plot_DM("pfMET", "DM_LO_V_Mx-150_Mv-500", 7.160e-02, 49999, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-10000", 1.222e-07, 49399, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-1000", 1.446e-02, 49999, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-100", 3.877e-01, 41000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-10", 5.114e-01, 50000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-2000", 1.207e-03, 47599, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-200", 2.683e-01, 49997, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-20", 4.752e-01, 49998, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-300", 1.785e-01, 50000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-500", 7.464e-02, 49998, false);
  plot_DM("pfMET", "DM_LO_V_Mx-1_Mv-50", 4.451e-01, 49400, false);
  plot_DM("pfMET", "DM_LO_V_Mx-500_Mv-10000", 8.307e-08, 50000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-500_Mv-10", 1.087e-04, 50000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-500_Mv-2000", 1.152e-03, 50000, false);
  plot_DM("pfMET", "DM_LO_V_Mx-500_Mv-500", 1.510e-04, 48400, false);
  plot_DM("pfMET", "DM_LO_V_Mx-500_Mv-995", 3.136e-03, 47600, false);
  plot_DM("pfMET", "DM_LO_V_Mx-50_Mv-10000", 1.223e-07, 49998, false);
  plot_DM("pfMET", "DM_LO_V_Mx-50_Mv-200", 2.646e-01, 49999, false);
  plot_DM("pfMET", "DM_LO_V_Mx-50_Mv-300", 1.776e-01, 49997, false);
  plot_DM("pfMET", "DM_LO_V_Mx-50_Mv-50", 1.464e-02, 49998, false);
  plot_DM("pfMET", "DM_LO_V_Mx-50_Mv-95", 4.502e-02, 50000, false);

  TFile* f_DM_LO_histos_Mt = new TFile("DM_LO_histos_below0p5_Mt.root", "RECREATE");
  f_DM_LO_histos_Mt->Close();
  plot_DM("Mt", "DM_LO_AV_Mx-1000_Mv-10000", 1.572e-08, 49997, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1000_Mv-1000", 1.830e-06, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1000_Mv-10", 1.344e-06, 49999, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1000_Mv-1995", 1.649e-05, 14198, false);
  plot_DM("Mt", "DM_LO_AV_Mx-10_Mv-10000", 1.217e-07, 49998, false);
  plot_DM("Mt", "DM_LO_AV_Mx-10_Mv-100", 3.776e-01, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-10_Mv-10", 2.898e-02, 800, false);
  plot_DM("Mt", "DM_LO_AV_Mx-10_Mv-15", 3.174e-02, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-10_Mv-50", 3.952e-01, 49999, false);
  plot_DM("Mt", "DM_LO_AV_Mx-150_Mv-10000", 1.076e-07, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-150_Mv-1000", 1.349e-02, 26800, false);
  plot_DM("Mt", "DM_LO_AV_Mx-150_Mv-10", 1.284e-03, 49998, false);
  plot_DM("Mt", "DM_LO_AV_Mx-150_Mv-200", 1.800e-03, 49200, false);
  plot_DM("Mt", "DM_LO_AV_Mx-150_Mv-295", 4.456e-03, 49999, false);
  plot_DM("Mt", "DM_LO_AV_Mx-150_Mv-500", 5.312e-02, 49997, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1_Mv-10000", 1.218e-07, 49999, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1_Mv-1000", 1.465e-02, 49997, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1_Mv-100", 3.876e-01, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1_Mv-10", 5.202e-01, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1_Mv-2000", 1.213e-03, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1_Mv-200", 2.681e-01, 49999, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1_Mv-20", 4.870e-01, 49198, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1_Mv-300", 1.783e-01, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-1_Mv-50", 4.468e-01, 49999, false);
  plot_DM("Mt", "DM_LO_AV_Mx-500_Mv-10000", 5.267e-08, 49999, false);
  plot_DM("Mt", "DM_LO_AV_Mx-500_Mv-10", 3.674e-05, 49999, false);
  plot_DM("Mt", "DM_LO_AV_Mx-500_Mv-2000", 8.975e-04, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-500_Mv-500", 4.725e-05, 49999, false);
  plot_DM("Mt", "DM_LO_AV_Mx-500_Mv-995", 2.678e-04, 48799, false);
  plot_DM("Mt", "DM_LO_AV_Mx-50_Mv-10000", 1.198e-07, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-50_Mv-10", 7.734e-03, 2400, false);
  plot_DM("Mt", "DM_LO_AV_Mx-50_Mv-200", 2.144e-01, 50000, false);
  plot_DM("Mt", "DM_LO_AV_Mx-50_Mv-300", 1.631e-01, 48798, false);
  plot_DM("Mt", "DM_LO_AV_Mx-50_Mv-50", 8.656e-03, 49998, false);
  plot_DM("Mt", "DM_LO_AV_Mx-50_Mv-95", 1.482e-02, 49198, false);
  plot_DM("Mt", "DM_LO_V_Mx-1000_Mv-10000", 3.499e-08, 48399, false);
  plot_DM("Mt", "DM_LO_V_Mx-1000_Mv-1000", 7.679e-06, 49196, false);
  plot_DM("Mt", "DM_LO_V_Mx-1000_Mv-10", 5.262e-06, 4000, false);
  plot_DM("Mt", "DM_LO_V_Mx-1000_Mv-1995", 2.533e-04, 50000, false);
  plot_DM("Mt", "DM_LO_V_Mx-10_Mv-10000", 1.216e-07, 49998, false);
  plot_DM("Mt", "DM_LO_V_Mx-10_Mv-100", 3.885e-01, 49998, false);
  plot_DM("Mt", "DM_LO_V_Mx-10_Mv-10", 3.781e-02, 50000, false);
  plot_DM("Mt", "DM_LO_V_Mx-10_Mv-15", 4.527e-02, 49199, false);
  plot_DM("Mt", "DM_LO_V_Mx-10_Mv-50", 4.422e-01, 49999, false);
  plot_DM("Mt", "DM_LO_V_Mx-150_Mv-10000", 1.195e-07, 49999, false);
  plot_DM("Mt", "DM_LO_V_Mx-150_Mv-1000", 1.437e-02, 49999, false);
  plot_DM("Mt", "DM_LO_V_Mx-150_Mv-10", 2.668e-03, 50000, false);
  plot_DM("Mt", "DM_LO_V_Mx-150_Mv-200", 4.398e-03, 49999, false);
  plot_DM("Mt", "DM_LO_V_Mx-150_Mv-295", 2.901e-02, 49999, false);
  plot_DM("Mt", "DM_LO_V_Mx-150_Mv-500", 7.160e-02, 49999, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-10000", 1.222e-07, 49399, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-1000", 1.446e-02, 49999, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-100", 3.877e-01, 41000, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-10", 5.114e-01, 50000, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-2000", 1.207e-03, 47599, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-200", 2.683e-01, 49997, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-20", 4.752e-01, 49998, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-300", 1.785e-01, 50000, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-500", 7.464e-02, 49998, false);
  plot_DM("Mt", "DM_LO_V_Mx-1_Mv-50", 4.451e-01, 49400, false);
  plot_DM("Mt", "DM_LO_V_Mx-500_Mv-10000", 8.307e-08, 50000, false);
  plot_DM("Mt", "DM_LO_V_Mx-500_Mv-10", 1.087e-04, 50000, false);
  plot_DM("Mt", "DM_LO_V_Mx-500_Mv-2000", 1.152e-03, 50000, false);
  plot_DM("Mt", "DM_LO_V_Mx-500_Mv-500", 1.510e-04, 48400, false);
  plot_DM("Mt", "DM_LO_V_Mx-500_Mv-995", 3.136e-03, 47600, false);
  plot_DM("Mt", "DM_LO_V_Mx-50_Mv-10000", 1.223e-07, 49998, false);
  plot_DM("Mt", "DM_LO_V_Mx-50_Mv-200", 2.646e-01, 49999, false);
  plot_DM("Mt", "DM_LO_V_Mx-50_Mv-300", 1.776e-01, 49997, false);
  plot_DM("Mt", "DM_LO_V_Mx-50_Mv-50", 1.464e-02, 49998, false);
  plot_DM("Mt", "DM_LO_V_Mx-50_Mv-95", 4.502e-02, 50000, false);
}
