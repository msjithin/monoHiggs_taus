
echo "Copy plot_etau.C? Enter 1 for yes, 0 for no. "
read isCopyETau

if [[ $isCopyETau = "1" ]]; then
    echo "Copying postAnalyzer_etau.C ..... "
    scp jmadhusu@lxplus.cern.ch:/afs/cern.ch/user/j/jmadhusu/CMSSW_9_4_4/src/monoHiggs/analyzer/postAnalyzer_etau.C .
else
    echo "Files copied from   ~/CMSSW_9_4_4/src/monoHiggs/analyzer/ "

    echo "Type the name of files seprated with commas and no spaces :"
    read files

    if [[ $files = *","* ]]; then
	echo "More than one files copying"
	scp jmadhusu@lxplus.cern.ch:/afs/cern.ch/user/j/jmadhusu/CMSSW_9_4_4/src/monoHiggs/analyzer/\{$files\} .
    else
	echo "One file copying"
	scp jmadhusu@lxplus.cern.ch:/afs/cern.ch/user/j/jmadhusu/CMSSW_9_4_4/src/monoHiggs/analyzer/$files .
    fi
fi