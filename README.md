This github repository contains the phyloseq object "CellarBats_ps_min250reads.rds" saved as RDS file and the skipt "CellarBats_Analysis.R" containing all the critical analyses steps performed in R for the manuscript Schmid et al. (unpublished). 

The work describes an experiment performed to a) assess whether diurnal microbial rhythms are detectable in captive bats and b) to see whether shifting feeding times alters/shifts microbial responses/oscillations

the phyloseq object's meta data contains 12 columns: 
sample ID = unique identifier
replicate = whether the same was replicated 
sample_plaze = which is the plate the sample was on
day = the day the sample was collected
treatment = whether the sample was from the control (c) or delayed feeding (t) group
time = day time in 24 hours
time_point = day time in 24 hours but as factor
weight_mg = fecal sample weight in mg
volume_ml = water added in ml for measuring pH
pH_value = measured pH
consistency = samples physical consistency (s = solid; ss = semi-solid, d = dried liquid)
sample_or_control = column used during filtering to subset extraction blanks 

Detailed information to experimental design and analysis can be found in the main text of the publication:
