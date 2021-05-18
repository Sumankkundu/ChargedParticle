brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_DiPFJetAve40_v*" -o output_40.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_DiPFJetAve60_v*" -o output_60.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_DiPFJetAve80_v*" -o output_80.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_DiPFJetAve140_v*" -o output_140.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_DiPFJetAve200_v*" -o output_200.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_DiPFJetAve260_v*" -o output_260.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_DiPFJetAve320_v*" -o output_320.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_DiPFJetAve400_v*" -o output_400.csv
brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --hltpath "HLT_DiPFJetAve500_v*" -o output_500.csv

pileupReCalc_HLTpaths.py -i output_40.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_40.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_60.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_60.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_80.csv --inputLumiJSON  /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_80.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_140.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_140.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_200.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_200.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_260.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_260.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_320.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_320.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_400.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_400.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_500.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_500.txt --runperiod Run2

pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_40.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_40.root
pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_60.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_60.root
pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_80.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_80.root
pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_140.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_140.root
pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_200.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_200.root
pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_260.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_260.root
pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON HLT_PileupJSON_320.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_320.root
pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt  --inputLumiJSON HLT_PileupJSON_400.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_400.root
pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt  --inputLumiJSON HLT_PileupJSON_500.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_500.root

