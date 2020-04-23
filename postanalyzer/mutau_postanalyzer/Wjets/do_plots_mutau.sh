#!/bin/bash
#rm plots/png/estimate.csv

histNumber="0"
#for ((i=1;i<=6;i+=1)
for i in {0..3}

do
    echo $i
    histNumber=$i
#    python makeplot_mutauh.py Tau_Pt_ $histNumber pT tau  
    python makeplot_mutauh.py Tau_Pt_ $histNumber pT_tau GeV

done

