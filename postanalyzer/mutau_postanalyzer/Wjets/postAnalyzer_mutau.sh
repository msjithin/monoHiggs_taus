
./Make.sh postAnalyzer_mutau.C




./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W1JetsToLNu.root files_initial/W1JetsToLNu.root W1JetsToLNu W1JetsToLNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W2JetsToLNu.root files_initial/W2JetsToLNu.root W2JetsToLNu W2JetsToLNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W3JetsToLNu.root files_initial/W3JetsToLNu.root W3JetsToLNu W3JetsToLNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/W4JetsToLNu.root files_initial/W4JetsToLNu.root W4JetsToLNu W4JetsToLNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu.root files_initial/WJetsToLNu.root WJetsToLNu WJetsToLNu 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_HT-100To200.root  files_initial/WJetsToLNu_HT100To200.root WJetsToLNu_HT100To200 WJetsToLNu_HT100To200 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_HT-1200To2500.root files_initial/WJetsToLNu_HT1200To2500.root WJetsToLNu_HT1200To2500 WJetsToLNu_HT1200To2500 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_HT-200To400.root files_initial/WJetsToLNu_HT200To400.root WJetsToLNu_HT200To400 WJetsToLNu_HT200To400 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_HT-2500ToInf.root files_initial/WJetsToLNu_HT2500ToInf.root WJetsToLNu_HT2500ToInf WJetsToLNu_HT2500ToInf 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_HT-400To600.root files_initial/WJetsToLNu_HT400To600.root WJetsToLNu_HT400To600 WJetsToLNu_HT400To600 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_HT-600To800.root files_initial/WJetsToLNu_HT600To800.root WJetsToLNu_HT600To800 WJetsToLNu_HT600To800 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_HT-800To1200.root files_initial/WJetsToLNu_HT800To1200.root WJetsToLNu_HT800To1200 WJetsToLNu_HT800To1200 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_ext.root files_initial/WJetsToLNu_ext.root WJetsToLNu_inc WJetsToLNu_inc 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_inc.root files_initial/WJetsToLNu_inc.root WJetsToLNu_inc WJetsToLNu_inc 0
./postAnalyzer_mutau.exe /nfs_scratch/jmadhusu/CMSSW_9_4_9_cand2/src/monoHiggs_w_ff/muTau/WJetsToLNu_inc_ext.root files_initial/WJetsToLNu_inc_ext.root WJetsToLNu_inc WJetsToLNu_inc 0


hadd -f f_mutau_initial.root files_initial/*.root
