make bin
./testunfoldmain2c
root -l -b -q coorel_Lcurve.C
root -l -b -q DataMC_Hist_Charged_binned.C
root -l -b -q DataMC_Hist_Jet_binned.C
root -l -b -q DataMC_Jet_binned_NoReg.C
root -l -b -q DataMC_Charged_binned_NoReg.C

pdfjam Esv_dataMC_Jets.pdf Esv_dataMC_charged.pdf --nup 8x5 --landscape --outfile ESV_DataMC.pdf
pdfjam Lcurve_corr_Jet.pdf Lcurve_corr_Charge.pdf --nup 8x5 --landscape --outfile Correl_Lcurve.pdf
pdfjam NoReg_corr_Jet.pdf NoReg_corr_Charge.pdf  --nup 8x5 --landscape --outfile Correl_NoReg.pdf
pdfjam Lcurve_Jet.pdf Lcurve_Charged.pdf --nup 8x5 --landscape --outfile Lcurve.pdf
pdfjam Esv_Jets_Binned_NoReg.pdf Esv_charged_Binned_NoReg.pdf --nup 8x5 --landscape --outfile ESV_Binned_NoReg.pdf
pdfjam Response_Jet_binned.pdf Response_char_binned.pdf --nup 8x5 --landscape --outfile Response_binned.pdf

pdftk ESV_DataMC.pdf Correl_Lcurve.pdf Lcurve.pdf ESV_Binned_NoReg.pdf Correl_NoReg.pdf Response_binned.pdf cat output Binned_Lcurve_ALL.pdf
