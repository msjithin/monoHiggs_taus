#!/bin/sh


export CMSSW_RELEASE_BASE=/cvmfs/cms.cern.ch/slc7_amd64_gcc630/cms/cmssw/CMSSW_9_4_1/

cat>Job_${6}.sh<<EOF
#!/bin/sh
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd $CMSSW_RELEASE_BASE
eval `scramv1 runtime -sh`
cd \${_CONDOR_SCRATCH_DIR}
./${1} ${2} ${3} ${4} ${5} ${6}
EOF

chmod 775 Job_${6}.sh

cat>condor_${6}<<EOF
x509userproxy = /tmp/x509up_u10034
universe = vanilla
Executable = Job_${6}.sh
Notification         = never
WhenToTransferOutput = On_Exit
ShouldTransferFiles  = yes
Requirements = (TARGET.UidDomain == "hep.wisc.edu" && TARGET.HAS_CMS_HDFS && OpSysAndVer == "CENTOS7" && TARGET.Arch == "X86_64" && (MY.RequiresSharedFS=!=true || TARGET.HasAFS_OSG) && (TARGET.OSG_major =!= undefined || TARGET.IS_GLIDEIN=?=true) && IsSlowSlot=!=true)
on_exit_remove       = (ExitBySignal == FALSE && (ExitCode == 0 || ExitCode == 42 || NumJobStarts>3))
+IsFastQueueJob      = True
getenv               = True
request_memory       = 1992
request_disk         = 2048000
Transfer_Input_Files = ${1}, /nfs_scratch/tost/CMSSW_9_4_9/src/tauTau/tauTriggerEfficiencies2017.root , /nfs_scratch/tost/CMSSW_9_4_9/src/tauTau/tauTriggerEfficiencies2017.root , /nfs_scratch/tost/CMSSW_9_4_9/src/tauTau/L1PrefiringMaps_new.root
output               = \$(Cluster)_\$(Process)_${6}.out
error                = \$(Cluster)_\$(Process)_${6}.err
Log                  = \$(Cluster)_\$(Process)_${6}.log
Queue
EOF

condor_submit condor_${6}
