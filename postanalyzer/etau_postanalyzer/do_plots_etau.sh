#!/bin/bash
#rm plots/png/estimate.csv

histNumber="0"
#for ((i=1;i<=6;i+=1)
for i in {3..3}

do
    echo $i
    histNumber=$i
#    python makeplot_etauh.py Tau_Pt_ $histNumber pT tau  
    python makeplot_etauh.py Tau_Pt_ $histNumber pT tau  
    python makeplot_etauh.py Electron_Pt_ $histNumber pT ele  
    python makeplot_etauh.py leadingJetPt_ $histNumber pT jet  
    python makeplot_etauh.py h_dPhi_ $histNumber '#Delta' 'phi'   
    python makeplot_etauh.py VisibleMass_ $histNumber 'm_{vis}' '_'  
    python makeplot_etauh.py pfMET_ $histNumber MET '_'  
    python makeplot_etauh.py HiggsPt_ $histNumber pT 'vector sum'  
#python makeplot_etauh.py pTsum pT 'scalar sum'  
#python makeplot_etauh.py mttot_1 total 'transverse mass'  
    python makeplot_etauh.py Mt_ $histNumber 'mu-met' ' transverse mass'   
#python makeplot_etauh.py deltaR '#Delta' R  
#python makeplot_etauh.py mjj mjj  
    python makeplot_etauh.py Tau_iso_ $histNumber Tau isolation  
    python makeplot_etauh.py Tau_Decay_Mode_ $histNumber Tau decay_mode
    python makeplot_etauh.py Tau_mass_ $histNumber tau mass
    python makeplot_etauh.py Tau_eta_ $histNumber tau '#eta'  
    python makeplot_etauh.py Electron_eta_ $histNumber ele '#eta'  
    python makeplot_etauh.py Tau_phi_ $histNumber tau '#phi'  
    python makeplot_etauh.py Electron_phi_ $histNumber ele '#phi'  
#python makeplot_etauh.py dphi_eMet '#Delta phi' 'ele-Met'  
#python makeplot_etauh.py dphi_tauMet '#Delta phi' 'tau-Met'  
# python makeplot_unroll.py visMmet visible mass  
    python makeplot_etauh.py nVtx_ $histNumber npv " "  
#python makeplot_etauh.py nJet_ njets " " 


done

python makeplot_etauh.py Events_level_ _ events eachlevel -logYaxis logy
python makeplot_etauh.py Cutflow _ events eachlevel -logYaxis logy
python makeplot_etauh.py MET _0  MET '_'  
python makeplot_etauh.py MET _1  MET '_'  
python makeplot_etauh.py MET _2  MET '_'  
python makeplot_etauh.py MET _3  MET '_'  
python makeplot_etauh.py MET _4  MET '_'  


