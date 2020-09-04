#brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i HLT_diPFJetAve_2017_ReReco_v1.txt --hltpath HLT_DiPFJetAve60_v12 -o output_60.csv
#brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i HLT_diPFJetAve_2017_ReReco_v1.txt --hltpath HLT_DiPFJetAve80_v11 -o output_80.csv
#brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i HLT_diPFJetAve_2017_ReReco_v1.txt --hltpath HLT_DiPFJetAve140_v11 -o output_140.csv
#brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i HLT_diPFJetAve_2017_ReReco_v1.txt --hltpath HLT_DiPFJetAve200_v11 -o output_200.csv
#brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i HLT_diPFJetAve_2017_ReReco_v1.txt --hltpath HLT_DiPFJetAve260_v12 -o output_260.csv
#brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i HLT_diPFJetAve_2017_ReReco_v1.txt --hltpath HLT_DiPFJetAve320_v12 -o output_320.csv
#brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i HLT_diPFJetAve_2017_ReReco_v1.txt --hltpath HLT_DiPFJetAve400_v12 -o output_400.csv
#brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i HLT_diPFJetAve_2017_ReReco_v1.txt --hltpath HLT_DiPFJetAve500_v12 -o output_500.csv

pileupReCalc_HLTpaths.py -i output_60.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_60.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_80.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_80.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_140.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_140.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_200.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_200.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_260.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_260.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_320.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_320.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_400.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_400.txt --runperiod Run2
pileupReCalc_HLTpaths.py -i output_500.csv --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt -o HLT_PileupJSON_500.txt --runperiod Run2

pileupCalc.py -i HLT_PileupJSON_60.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_60.root
pileupCalc.py -i HLT_PileupJSON_80.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_80.root
pileupCalc.py -i HLT_PileupJSON_140.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_140.root
pileupCalc.py -i HLT_PileupJSON_200.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_200.root
pileupCalc.py -i HLT_PileupJSON_260.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_260.root
pileupCalc.py -i HLT_PileupJSON_320.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_320.root
pileupCalc.py -i HLT_PileupJSON_400.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_400.root
pileupCalc.py -i HLT_PileupJSON_500.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 HLT_Histogram_500.root

