# Tumor-mutational-patterns-in-CMMRD

The code in this repository was used to perform analyses and create figures for the manuscript "Unraveling mutagenic processes influencing the tumor mutational patterns of individuals with Constitutional Mismatch Repair Deficiency" by Weijers et al.

The order in and purpose for which the different scripts were used is described below. A supplementary file containing the mutational matrices used in the analyses is provided, as well as a file containing the conversion between the manuscript IDs and the IDs in the datasets uploaded to EGA.
  
<br />  
**SOMATIC CALLING**
  
For samples sequenced with Twist v1, run:  
generateMutect2CallingJobs.WGS-WXS.IntersectTwist.sh  

For all other samples, run:  
generateMutect2CallingJobs.WGS-WXS.Intersect.sh  
  
Both are followed by:  
generateMutect2FilterJobs.WXS.noPASS.Intersect.noMNP.sh  
  
<br />  
**VCF FILTERING**

For all samples:  
FilterVCFs.R  
  
<br />  
**MUTATIONAL SIGNATURES**

First, create the mutational matrices:  
Create_MutationalMatrices.R  
  
Extract Mutational Signatures using a combination of:  
generateDeNovoExtractionJob.sh and denovoextraction.R (SNV)  
generateDeNovoExtractionJobIndel.sh and denovoextraction_indel.R (indel)  
  
Perform further downstream signature analyses using:  
SignatureAnalysesRefits.R  
  
<br />  
**PLOT MANUSCRIPT FIGURES**

PlottingCMMRDFigures.R  
  
<br />  
**PERFORM STATISTICS**

StatisticsCMMRD.R  
  
<br />  
**SUPPLEMENTARY FILES**

Mutational Matrices: NCOMMS-23-52197_MutationalMatrices.xlsx  
Conversion file manuscript ID to EGA ID: NCOMMS-23-52197_Conversion_ManuscriptID_EGA-ID.xlsx  

