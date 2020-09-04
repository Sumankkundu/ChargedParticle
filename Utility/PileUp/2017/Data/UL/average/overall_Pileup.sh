pileupCalc.py -i Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 PileUp_17UL_ture_overall.root
pileupCalc.py -i Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode observed --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 PileUp_17UL_Observer_overall.root

