make bin
./testunfoldmain2c
root -l -b -q coorel_Lcurve.C
root -l -b -q DataMC_Hist_Jet.C
root -l -b -q DataMC_Hist_Charged.C
root -l -b -q DataMC_Hist_Charged_NoReg.C
root -l -b -q DataMC_Hist_Jet_NoReg.C

pdfjam Esv_dataMC_Jets.pdf Esv_dataMC_charged.pdf --nup 8x5 --landscape --outfile ESV_DataMC.pdf
pdfjam Esv_dataMC_Jets_NoReg.pdf Esv_dataMC_charged_NoReg.pdf --nup 8x5 --landscape --outfile ESV_DataMC_NoReg.pdf
pdfjam Lcurve_corr_Jet.pdf Lcurve_corr_Charge.pdf --nup 8x5 --landscape --outfile Correl_Lcurve.pdf
pdfjam Lcurve_Jet.pdf Lcurve_Charged.pdf --nup 8x5 --landscape --outfile Lcurve.pdf
pdfjam NoReg_corr_Jet.pdf NoReg_corr_Charge.pdf --nup 8x5 --landscape --outfile NoReg_corr.pdf

pdftk ESV_DataMC.pdf Correl_Lcurve.pdf Lcurve.pdf ESV_DataMC_NoReg.pdf NoReg_corr.pdf cat output TUnfolding_all.pdf
