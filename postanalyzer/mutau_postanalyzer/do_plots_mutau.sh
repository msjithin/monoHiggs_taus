#!/bin/bash
#rm plots/png/estimate.csv

histNumber="0"
#for ((i=1;i<=6;i+=1)
for i in {3..3}

do
    echo $i
    histNumber=$i
#    python makeplot_mutauh.py Tau_Pt_ $histNumber pT tau  
    python makeplot_mutauh.py Tau_Pt_ $histNumber pT tau  
    python makeplot_mutauh.py Muon_Pt_ $histNumber pT mu  
    python makeplot_mutauh.py leadingJetPt_ $histNumber pT jet  
    python makeplot_mutauh.py h_dPhi_ $histNumber '#Delta' 'phi'   
    python makeplot_mutauh.py VisibleMass_ $histNumber 'm_{vis}' '_'  
    python makeplot_mutauh.py pfMET_ $histNumber MET '_'  
    python makeplot_mutauh.py HiggsPt_ $histNumber pT 'vector sum'  
#python makeplot_mutauh.py pTsum pT 'scalar sum'  
#python makeplot_mutauh.py mttot_1 total 'transverse mass'  
    python makeplot_mutauh.py Mt_ $histNumber 'mu-met' ' transverse mass'   
#python makeplot_mutauh.py deltaR '#Delta' R  
#python makeplot_mutauh.py mjj mjj  
    python makeplot_mutauh.py Tau_iso_ $histNumber Tau isolation  
    python makeplot_mutauh.py Tau_Decay_Mode_ $histNumber Tau decay_mode
    python makeplot_mutauh.py Tau_mass_ $histNumber tau mass
    python makeplot_mutauh.py Tau_eta_ $histNumber tau '#eta'  
    python makeplot_mutauh.py Muon_eta_ $histNumber mu '#eta'  
    python makeplot_mutauh.py Tau_phi_ $histNumber tau '#phi'  
    python makeplot_mutauh.py Muon_phi_ $histNumber mu '#phi'  
#python makeplot_mutauh.py dphi_eMet '#Delta phi' 'ele-Met'  
#python makeplot_mutauh.py dphi_tauMet '#Delta phi' 'tau-Met'  
# python makeplot_unroll.py visMmet visible mass  
    python makeplot_mutauh.py nVtx_ $histNumber npv " "  
#python makeplot_mutauh.py nJet_ njets " " 


done

python makeplot_mutauh.py Events_level_ _ events eachlevel -logYaxis logy
python makeplot_mutauh.py Cutflow _ events eachlevel -logYaxis logy
python makeplot_mutauh.py MET _0  MET '_'  
python makeplot_mutauh.py MET _1  MET '_'  
python makeplot_mutauh.py MET _2  MET '_'  
python makeplot_mutauh.py MET _3  MET '_'  
python makeplot_mutauh.py MET _4  MET '_'  


