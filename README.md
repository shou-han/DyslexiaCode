# DyslexiaCode
The codes in this folder is to collect and anaylse the data for the paper Stefanac, N., Zhou, SH., O'Connell, R., Spencer-Smith, M., Castles, A., Bellgrove, M. "Insufficient Evidence Accumulation in Dyslexia." 

The data is collected using a random dots paradigm built using psychtoolbox http://psychtoolbox.org/. The code for the random dots paradigm is given in RandomContinuousDotsParadigm.m. 

The analysis files are split into three folders: PreprocessHAPPE, analysisDyslexia and R_analysis.
The PreprocessHAPPE uses the Harvard Automated Processing Pipeline for EEG (HAPPE) to preprocess the EEG data, as described in the paper. The file to use for this is runafew_HAPPE.m. After the initial preprocessing, another process is made to extract the relevant ERPs and the SSVEPs using the files runafew_after.m and runafew_after_SSVEP.m respectively. 

After the preprocess, the data is analysed in analysisDyslexia. The main analysis is performed in erp_analysis_fn.m to extract the CPP waveform and its features (e.g. slope, amplitude ). Similarly, the SSVEP and the fast/slow bin analysis are performed in erp_analysis_fn_SSVEP.m and erp_analysis_fn_bins.m respectively. 

Finally, after the analysed data plots can be visualised in DylexiaGraphs.m, DyslexiaGraphs_SSVEP.m and DyslexiaGraphs_bins.m. 

Statistics of the extracted features of the CPP/SSVEP waveforms within hand, hemifield and bins are determiend in R_analysis. The main statistical analysis is performed in R_analysis_Paper.Rmd. The statistical analysis for SSVEP and Bins are performed in R_analysis_Paper_Bins.Rmd and R_analysis_Paper_SSVEP.Rmd respectively. 

Further details of the preprocessing, analysis, the visualised plots and the statistics are found in the paper. 

The HAPPE protocol can be found in 
Gabard-Durnam, L. J., Mendez Leal, A. S., Wilkinson, C. L., & Levin, A. R. (2018). The Harvard Automated Processing Pipeline for Electroencephalography (HAPPE): standardized processing software for developmental and high-artifact data. Frontiers in neuroscience, 12, 97.
