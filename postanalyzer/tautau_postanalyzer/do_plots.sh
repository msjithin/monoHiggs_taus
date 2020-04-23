#!/bin/bash
#rm plots/png/estimate.csv

histNumber="0"
#for ((i=1;i<=6;i+=1)
for i in {2..3}

do
    echo $i
    histNumber=$i
#    python makeplot_tautauh.py Tau_Pt_ $histNumber pT tau 
    python makeplot_tautauh.py Tau_Pt_ $histNumber 'pT tau2' 'GeV' 
    python makeplot_tautauh.py Muon_Pt_ $histNumber 'pT tau1' 'GeV' 
    python makeplot_tautauh.py leadingJetPt_ $histNumber 'pT jet' 'GeV' 
    python makeplot_tautauh.py h_dPhi_ $histNumber '#Delta' 'phi'  
    python makeplot_tautauh.py VisibleMass_ $histNumber 'm_{vis}' '_' 
    python makeplot_tautauh.py pfMET_ $histNumber MET '_' 
    python makeplot_tautauh.py HiggsPt_ $histNumber pT 'vector sum' 
#python makeplot_tautauh.py pTsum pT 'scalar sum' 
#python makeplot_tautauh.py mttot_1 total 'transverse mass' 
    python makeplot_tautauh.py Mt_ $histNumber 'ele-met' ' transverse mass GeV'  
#python makeplot_tautauh.py deltaR '#Delta' R 
#python makeplot_tautauh.py mjj mjj 
#python makeplot_tautauh.py tauIso tau iso 
#python makeplot_tautauh.py taudecay tau decay_mode 
#python makeplot_tautauh.py iso_1 ele iso 
    python makeplot_tautauh.py Tau_eta_ $histNumber tau2 '#eta' 
    python makeplot_tautauh.py Muon_eta_ $histNumber tau1 '#eta' 
    python makeplot_tautauh.py Tau_phi_ $histNumber tau2 '#phi' 
    python makeplot_tautauh.py Muon_phi_ $histNumber tau1 '#phi' 
#python makeplot_tautauh.py dphi_eMet '#Delta phi' 'ele-Met' 
#python makeplot_tautauh.py dphi_tauMet '#Delta phi' 'tau-Met' 
# python makeplot_unroll.py visMmet visible mass 
    python makeplot_tautauh.py nVtx_ $histNumber npv " " 
#python makeplot_tautauh.py nJet_ njets " " 


done

python makeplot_tautauh.py Events_level_ _ events eachlevel -logYaxis logy
python makeplot_tautauh.py Cutflow _ events eachlevel -logYaxis logy
#python makeplot_tautauh.py MET _0  MET '_' 
#python makeplot_tautauh.py MET _1  MET '_' 
#python makeplot_tautauh.py MET _2  MET '_' 
#python makeplot_tautauh.py MET _3  MET '_' 
#python makeplot_tautauh.py MET _4  MET '_' 


