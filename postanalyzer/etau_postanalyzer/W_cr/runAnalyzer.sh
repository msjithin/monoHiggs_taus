
set -e

./rootcom Analyzer_etau analyze
if [ $? -eq 0 ]; then
  echo Compiling successful
else
  echo Compiling failed;
  exit 1
fi

#To run, assuming this is compiled to an executable named 'analyze':
./analyze /hdfs/store/user/jmadhusu/MonoHiggs_MC2017/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/crab_ZZZ/180603_185329/0000/ /afs/hep.wisc.edu/user/ms/CMSSW_9_4_4/src/2017_analysis/analyzer/output.root -1 10000