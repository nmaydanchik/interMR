# Preliminary data analysis using InterMR to identify risk factors with sex-specific effects on ADHD 
This repository contains coding samples from work I did for a biostatistics laboratory at UChicago as a Temporary Research Coordinator where I did statistical analysis using our proposed new Mendelian randomization method, InterMR, to **find exposures with sex-specific effects on ADHD**.

Project Dates: Jul-Aug 2024 (summer before Senior year of High School)

This project highlights my ability to:

* Review existing literature and collect and process data
* Work with datasets with millions of values in R
* Conduct Mendelian randomization analysis
* Create presentations

I presented my findings 8/16/24 to a lab group with a professor and Ph.D. students. The presentation is available as AugPresentation.pdf. This was part of a greater research effort by a Ph.D. student that resulted in a paper published in PLOS Genetics, for which I am a co-author. This paper can be accessed here: https://doi.org/10.1371/journal.pgen.1011819.

## Project and File Description

The objective of this project was to apply our proposed method Mendelian randomization data to real data and showcase its ability to detect group-specific effects. Additionally, we wished to showcase the increased power of the method to detect group-specific effects when multiple outcome GWAS datasets are integrated. One of the two ways we applied was in ADHD data I found that had information of the proportion of males and females in the study, which was my task.

We considered 54 possible exposures, which I analyzed in ADHD_analysis.R. Additionally, I further investigated exposures such as Birth weight, Autism spectrum disorder, and Cannabis use disorder by using a variety of datasets to see if results would be consistent when using GWAS data from different sources. I processed this data in process_data.R and did exposure-specific analyses in exposure_analysis.R.

## Bonus File

In July 2025, when we published the paper along with a package that implements the proposed method (now renamed int2MR), I was asked to help rewrite the package documentation and a vignette that users could refer to as a tutorial to make the package more accessible and increase the impact of the research. The vignette is available as int2MR.pdf.
