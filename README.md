# Tumor-mutational-patterns-in-CMMRD

The code in this repository was used to perform analyses and create figures for the manuscript "**Unraveling mutagenic processes influencing the tumor mutational patterns of individuals with Constitutional Mismatch Repair Deficiency**" by *Weijers et al*.

The order in and purpose for which the different scripts were used is described below. Further instructions for running the code are included in the code as comments. All unfiltered VCF files containing the somatic mutations of all individuals and a supplementary file containing the mutational matrices used in the analyses is provided, as well as a file containing the conversion between the manuscript IDs and the IDs in the datasets uploaded to EGA.

No separate demo data was provided, because (part of) the original data can be used to demo the code.

The following software and tools are used in the code provided:
```Burrows-Wheeler Aligner v0.7.13
Genome Analysis Toolkit v4.2.0.0 
Java/1.8.0_60
Ensembl VEP v92 https://github.com/Ensembl/ensembl-vep/tree/release/92
samtools/1.9
R v4.2.1
R v4.1.2 (for the scripts under "Extract Mutational Signatures" specifically)
MutationalPatterns v3.10.0  https://bioconductor.org/packages/release/bioc/html/MutationalPatterns.html
BSgenome.Hsapiens.UCSC.hg38 v1.4.5
python3
SigProfilerMatrixGenerator v1.1.28 https://github.com/AlexandrovLab/SigProfilerMatrixGenerator
SigProfilerExtractor v1.1.0 https://github.com/AlexandrovLab/SigProfilerExtractor
ggplot2 v3.4.2 https://ggplot2.tidyverse.org/
stats v4.2.1 https://rdrr.io/r/stats/stats-package.html
FSA v0.9.5 https://fishr-core-team.github.io/FSA/
contingencytables v3.0.1 https://github.com/ocbe-uio/contingencytables
VariantAnnotation v1.50.0
GenomeInfoDb 1.40.1 ```


The code has not been tested on other versions of these software and tools.
Only previously developed software and tools were used.
Instructions for installing these software and tools and, where relevant, the typical install time can be consulted on the developers webpage or repository.

Part of the code has been run on a high-performance cluster, these sections are indicated below.
  
<br />  

**SOMATIC CALLING**
  
For samples sequenced with Twist v1, run:  
```generateMutect2CallingJobs.WGS-WXS.IntersectTwist.sh```  

For all other samples, run:  
```generateMutect2CallingJobs.WGS-WXS.Intersect.sh```  
  
Both are followed by:  
```generateMutect2FilterJobs.WXS.noPASS.Intersect.noMNP.sh```  

The code above has been run on a high-performance cluster.
The expected output is a VCF file containing all somatic mutations with Ensembl VEP annotations.
  
<br />  

**VCF FILTERING**

For all samples:  
```FilterVCFs.R```  

The expected output is a filtered VCF file.
  
<br />  

**MUTATIONAL SIGNATURES**

First, create the mutational matrices:  
```Create_MutationalMatrices.R```  
The expected output is a mutational matrix containing the numbers of mutations per tumor per mutation type.
  
Extract Mutational Signatures using a combination of:  
```generateDeNovoExtractionJob.sh and denovoextraction.R (SNV)```  
```generateDeNovoExtractionJobIndel.sh and denovoextraction_indel.R (indel)```  
The code for "Extract Mutational Signatures" has been run on a high-performance cluster.
The expected output is an Rdata file containing the extracted signatures.
  
Perform further downstream signature analyses using:  
```SignatureAnalysesRefits.R```  
Expected output is a variety of R objects containing specific analyses.
  
<br />  

**SigProfiler**

Create the mutational matrix using SigProfiler:
```Part_I_create_mutational_matrices.R```
```Part_II_Mutational_Matrices.txt```
The expected output is a mutational matrix containing the numbers of mutations per tumor per mutation type.

Signature extraction using SigProfiler:
```Mutational_Signatures.txt```
The expected output is an object containing the extracted signatures.

<br />

**PLOT MANUSCRIPT FIGURES**

```PlottingCMMRDFigures.R```
```MS_contributions.R```  
```Supplemental_Fig5C.R```
The expected output is a variety of figures and graphs.
  
<br />  

**PERFORM STATISTICS**

```StatisticsCMMRD.R```
The expected output are values associated with these statistical tests, among which are p-values.
  
<br />  

**SUPPLEMENTARY FILES**

Mutational Matrices: ```NCOMMS-23-52197_MutationalMatrices.xlsx```  
VCF files: 
Conversion file manuscript ID to EGA ID: ```NCOMMS-23-52197_Conversion_ManuscriptID_EGA-ID.xlsx```  

