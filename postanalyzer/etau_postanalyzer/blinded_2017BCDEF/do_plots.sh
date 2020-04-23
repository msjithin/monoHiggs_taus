#!/bin/bash
#rm plots/png/estimate.csv

histNumber="0"
#for ((i=1;i<=6;i+=1)
for i in {3..5}

do
    echo $i
    histNumber=$i
#    python makeplot_etauh.py Tau_Pt_ $histNumber pT tau  
    python makeplot_etauh.py Tau_Pt_ $histNumber 'pT tau' 'GeV'  
    python makeplot_etauh.py Electron_Pt_ $histNumber 'pT ele' 'GeV'  
    python makeplot_etauh.py leadingJetPt_ $histNumber 'pT jet' 'GeV'  
    python makeplot_etauh.py h_dPhi_ $histNumber '#Delta' 'phi'   
    python makeplot_etauh.py VisibleMass_ $histNumber 'm_{vis}' '_'  
    python makeplot_etauh.py pfMET_ $histNumber MET '_'  
    python makeplot_etauh.py HiggsPt_ $histNumber pT 'vector sum'  
#python makeplot_etauh.py pTsum pT 'scalar sum'  
    python makeplot_etauh.py TotalTrMass_ $histNumber 'total' ' transverse mass GeV'  
    python makeplot_etauh.py Mt_ $histNumber 'ele-met' ' transverse mass GeV'   
#python makeplot_etauh.py deltaR '#Delta' R  
#python makeplot_etauh.py mjj mjj  
#python makeplot_etauh.py tauIso tau iso  
#python makeplot_etauh.py taudecay tau decay_mode  
#python makeplot_etauh.py iso_1 ele iso  
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

python makeplot_etauh.py Tau_Pt_ 6 'pT tau' 'GeV' 
python makeplot_etauh.py Electron_Pt_ 6 'pT ele' 'GeV' 
python makeplot_etauh.py leadingJetPt_ 6 'pT jet' 'GeV' 
python makeplot_etauh.py h_dPhi_ 6 '#Delta' 'phi' 
python makeplot_etauh.py VisibleMass_ 6 'm_{vis}' '_' 
python makeplot_etauh.py pfMET_ 6 MET '_' 
python makeplot_etauh.py HiggsPt_ 6 pT 'vector sum' 
python makeplot_etauh.py Mt_ 6 'ele-met' ' transverse mass GeV' 
python makeplot_etauh.py TotalTrMass_ 6 'total' ' transverse mass GeV'

python makeplot_etauh.py Tau_eta_ 6 tau '#eta' 
python makeplot_etauh.py Electron_eta_ 6 ele '#eta' 
python makeplot_etauh.py Tau_phi_ 6 tau '#phi' 
python makeplot_etauh.py Electron_phi_ 6 ele '#phi' 
python makeplot_etauh.py nVtx_ 6 npv " " 

python makeplot_etauh.py Events_level_ _ events eachlevel -logYaxis logy 
python makeplot_etauh.py Cutflow _ events eachlevel -logYaxis logy   
#python makeplot_etauh.py MET _0  MET '_' -logYaxis logy 
#python makeplot_etauh.py MET _1  MET '_' -logYaxis logy
#python makeplot_etauh.py MET _2  MET '_' -logYaxis logy 
#python makeplot_etauh.py MET _3  MET '_' -logYaxis logy  
#python makeplot_etauh.py MET _4  MET '_' -logYaxis logy


