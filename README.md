# WhiteMatterDevelopment
DTI developmental trajectories in infants

1. WhiteMatterDevelopment_DataPrepHarmonization.html: applies longitudinal combat to the DTI regions of interest.

2. WhiteMatterDevelopment_Analysis_WholeBrain_v3.Rmd: fits GAMMs with the REML based smoothing.
 a.  DTI metrics ~ s(corrected age)
 b. DTI metrics ~ s(gestational age, chronological age)
 c. DTI metrics ~ s(corrected age, by=sex)

3. WhiteMatterDevelopment_Analysis_Mahalanobis.Rmd: calculate distance between infant DTI metrics and adult DTI, adult DTI calculated from 40 HCP participants  (compare distribution between infant and adult DTI)
   
4. WhiteMatterDevelopment_Analysis_TractSpecific.Rmd : 11 tracts analyses

5. WhiteMatterDevelopment_Analysis_TractSpecificGestationalAge.Rmd tensor splines fit to DTI ~ s(gestational age, chronological age)
