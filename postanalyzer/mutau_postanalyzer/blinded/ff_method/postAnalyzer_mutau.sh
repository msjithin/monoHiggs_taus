
./Make.sh postAnalyzer_mutau.C



./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/blinded/fakerate/Data_2017.root files_fakebkg/Data.root data_obs data_obs 0 

#################### MC ###############################3
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DYJetsToLL_LO.root files_fakebkg/DYJetsToLL_LO.root DY_LO DY_LO 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DYJetsToLL_LOext.root files_fakebkg/DYJetsToLL_LOext.root DY_LO DY_LO 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/TTTo2L2Nu.root files_fakebkg/TTTo2L2Nu.root TTTo2L2Nu TTTo2L2Nu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/TTToHadronic.root files_fakebkg/TTToHadronic.root TTToHadronic TTToHadronic 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/TTToSemiLeptonic.root files_fakebkg/TTToSemiLeptonic.root TTToSemiLeptonic TTToSemiLeptonic 0

#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DYJetsToLL_LO.root files_fakebkg/DYJetsToLL_LO_inc.root ZTTinc ZTTinc 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DYJetsToLL_LOext.root files_fakebkg/DYJetsToLL_LOext_inc.root ZTTinc ZTTinc 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DY1JetsToLL.root files_fakebkg/DY1JetsToLL.root ZTT1jet ZTT1jet 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DY2JetsToLL.root files_fakebkg/DY2JetsToLL.root ZTT2jet ZTT2jet 0 
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DY3JetsToLL.root files_fakebkg/DY3JetsToLL.root ZTT3jet ZTT3jet 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DY4JetsToLL.root files_fakebkg/DY4JetsToLL.root ZTT4jet ZTT4jet 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DYJetsToLL_0.root files_fakebkg/DYJetsToLL_0.root ZTTinc ZTTinc 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DYJetsToLL_1.root files_fakebkg/DYJetsToLL_1.root ZTTinc ZTTinc 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/DYJetsToLL_2.root files_fakebkg/DYJetsToLL_2.root ZTTinc ZTTinc 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/EWKWMinus2Jets_WToLNu.root files_fakebkg/EWKWMinus2Jets_WToLNu.root EWKWMinus2Jets EWKWMinus2Jets 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/EWKWPlus2Jets_WToLNu.root files_fakebkg/EWKWPlus2Jets_WToLNu.root EWKWPlus2Jets EWKWPlus2Jets 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/EWKZ2Jets_ZToLL.root files_fakebkg/EWKZ2Jets_ZToLL.root EWKZ2Jets_ZToLL EWKZ2Jets_ZToLL 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/EWKZ2Jets_ZToNuNu.root files_fakebkg/EWKZ2Jets_ZToNuNu.root EWKZ2Jets_ZToNuNu EWKZ2Jets_ZToNuNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/GluGluHToTauTau.root files_fakebkg/GluGluHToTauTau.root GluGluHToTauTau GluGluHToTauTau 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ST_t-channel_antitop.root files_fakebkg/ST_t-channel_antitop.root ST_t-channel_antitop ST_t-channel_antitop 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ST_t-channel_top.root files_fakebkg/ST_t-channel_top.root ST_t-channel_top ST_t-channel_top 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ST_tW_antitop.root files_fakebkg/ST_tW_antitop.root ST_tW_antitop ST_tW_antitop 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ST_tW_top.root files_fakebkg/ST_tW_top.root ST_tW_top ST_tW_top 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/TTJets_0.root files_fakebkg/TTJets_0.root TTJets TTJets 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/TTJets_1.root files_fakebkg/TTJets_1.root TTJets TTJets 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/TTJets_2.root files_fakebkg/TTJets_2.root TTJets TTJets 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/VBFHToTauTau.root files_fakebkg/VBFHToTauTau.root VBFHToTauTau VBFHToTauTau 0

#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WJetsToLNu.root files_fakebkg/WJetsToLNu.root WJetsToLNu WJetsToLNu 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W_cr/Wjet_shape/WJetsToLNu_ext.root files_fakebkg/WJetsToLNu_ext.root WJetsToLNu WJetsToLNu 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W_cr/Wjet_shape/WJetsToLNu_HT-100To200.root files_fakebkg/WJetsToLNu_HT100To200.root WJetsToLNu_HT100To200 WJetsToLNu_HT100To200 0 
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W_cr/Wjet_shape/WJetsToLNu_HT-1200To2500.root files_fakebkg/WJetsToLNu_HT1200To2500.root WJetsToLNu_HT1200To2500 WJetsToLNu_HT1200To2500 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W_cr/Wjet_shape/WJetsToLNu_HT-200To400.root files_fakebkg/WJetsToLNu_HT200To400.root WJetsToLNu_HT200To400 WJetsToLNu_HT200To400 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W_cr/Wjet_shape/WJetsToLNu_HT-2500ToInf.root files_fakebkg/WJetsToLNu_HT2500ToInf.root WJetsToLNu_HT2500ToInf WJetsToLNu_HT2500ToInf 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W_cr/Wjet_shape/WJetsToLNu_HT-400To600.root files_fakebkg/WJetsToLNu_HT400To600.root WJetsToLNu_HT400To600 WJetsToLNu_HT400To600 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W_cr/Wjet_shape/WJetsToLNu_HT-600To800.root files_fakebkg/WJetsToLNu_HT600To800.root WJetsToLNu_HT600To800 WJetsToLNu_HT600To800 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W_cr/Wjet_shape/WJetsToLNu_HT-800To1200.root files_fakebkg/WJetsToLNu_HT800To1200.root WJetsToLNu_HT800To1200 WJetsToLNu_HT800To1200 0

#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WJetsToLNu.root files_fakebkg/WJetsToLNu_inc.root WJetsToLNu_inc WJetsToLNu_inc 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WJetsToLNu_ext.root files_fakebkg/WJetsToLNu_ext_inc.root WJetsToLNu_inc WJetsToLNu_inc 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/W1JetsToLNu.root files_fakebkg/W1JetsToLNu.root W1JetsToLNu W1JetsToLNu 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/W2JetsToLNu.root files_fakebkg/W2JetsToLNu.root W2JetsToLNu W2JetsToLNu 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/W3JetsToLNu.root files_fakebkg/W3JetsToLNu.root W3JetsToLNu W3JetsToLNu 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/W4JetsToLNu.root files_fakebkg/W4JetsToLNu.root W4JetsToLNu W4JetsToLNu 0

./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WWTo2L2Nu.root files_fakebkg/WWTo2L2Nu.root WWTo2L2Nu VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WWToLNuQQ.root files_fakebkg/WWToLNuQQ.root WWToLNuQQ VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WWW.root files_fakebkg/WWW.root WWW VVV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WWZ.root files_fakebkg/WWZ.root WWZ VVV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WZTo1L1Nu2Q.root files_fakebkg/WZTo1L1Nu2Q.root WZTo1L1Nu2Q VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WZTo2L2Q.root files_fakebkg/WZTo2L2Q.root WZTo2L2Q VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WZTo3LNu.root files_fakebkg/WZTo3LNu.root WZTo3LNu VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WZZ.root files_fakebkg/WZZ.root WZZ VVV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WminusHToTauTau.root files_fakebkg/WminusHToTauTau.root WMinusHToTauTau WMinusHToTauTau 0
#./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WpWpJJ_EWK-QCD.root files_fakebkg/WpWpJJ_EWK-QCD.root WpWpJJ_EWK_QCD WpWpJJ_EWK_QCD 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WplusHToTauTau.root files_fakebkg/WplusHToTauTau.root WPlusHToTauTau WPlusHToTauTau 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZHToTauTau.root files_fakebkg/ZHToTauTau.root ZHToTauTau ZHToTauTau 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZJetsToNuNu_HT-100To200.root files_fakebkg/ZJetsToNuNu_HT-100To200.root ZJetsToNuNu_HT100To200 ZJetsToNuNu 0 
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZJetsToNuNu_HT-1200To2500.root files_fakebkg/ZJetsToNuNu_HT-1200To2500.root ZJetsToNuNu_HT1200To2500 ZJetsToNuNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZJetsToNuNu_HT-200To400.root files_fakebkg/ZJetsToNuNu_HT-200To400.root ZJetsToNuNu_HT200To400 ZJetsToNuNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZJetsToNuNu_HT-2500ToInf.root files_fakebkg/ZJetsToNuNu_HT-2500ToInf.root ZJetsToNuNu_HT2500ToInf ZJetsToNuNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZJetsToNuNu_HT-400To600.root files_fakebkg/ZJetsToNuNu_HT-400To600.root ZJetsToNuNu_HT400To600 ZJetsToNuNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZJetsToNuNu_HT-600To800.root files_fakebkg/ZJetsToNuNu_HT-600To800.root ZJetsToNuNu_HT600To800 ZJetsToNuNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZJetsToNuNu_HT-800To1200.root files_fakebkg/ZJetsToNuNu_HT-800To1200.root ZJetsToNuNu_HT800To1200 ZJetsToNuNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZZTo2L2Q.root files_fakebkg/ZZTo2L2Q.root ZZTo2L2Q VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZZTo4L_0.root files_fakebkg/ZZTo4L_0.root ZZTo4L VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZZTo4L_1.root files_fakebkg/ZZTo4L_1.root ZZTo4L VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZZZ.root files_fakebkg/ZZZ.root ZZZ VVV 0


./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/VBFHToWWTo2L2Nu.root files_fakebkg/VBFHToWWTo2L2Nu.root VBFHToWWTo2L2Nu VBFHToWWTo2L2Nu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WWTo2L2Nu_DoubleScattering.root files_fakebkg/WWTo2L2Nu_DoubleScattering.root WWTo2L2Nu_DoubleScattering VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WZTo1L3Nu.root files_fakebkg/WZTo1L3Nu.root WZTo1L3Nu VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WpWpJJ_EWK.root files_fakebkg/WpWpJJ_EWK.root WpWpJJ_EWK WpWpJJ_EWK_QCD 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/WpWpJJ_QCD.root files_fakebkg/WpWpJJ_QCD.root WpWpJJ_QCD WpWpJJ_EWK_QCD 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZZTo2Q2Nu_0.root files_fakebkg/ZZTo2Q2Nu_0.root ZZTo2Q2Nu VV 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff_1/muTau/fakerate/MC_2017_1/ZZTo2Q2Nu_1.root files_fakebkg/ZZTo2Q2Nu_1.root ZZTo2Q2Nu VV 0

hadd -f f_mutau_fakebkg.root files_fakebkg/*.root

echo "*************** root file made ***************"

#sh do_plots_mutau.sh

echo "*************** plots made ***************"
