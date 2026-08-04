# NMA-TVHR 
Pairwise (PWMA) and network meta-analysis (NMA) of time to event outcomes assuming non-proportional hazards using a time varying hazard ratio approach.

# Citation
Please cite: Owen, R.K., Abrams, K.R. Meta-analysis and network meta-analysis of time-to-event outcomes with non-proportional hazards: a Bayesian time-varying hazard ratio approach. BMC Med Res Methodol (2026). https://doi.org/10.1186/s12874-026-02963-5

Bibtex citation:
@article{Owen2026,
  author  = {Owen, Rhiannon K. and Abrams, Keith R.},
  title   = {Meta-analysis and network meta-analysis of time-to-event outcomes with non-proportional hazards: a {B}ayesian time-varying hazard ratio approach},
  journal = {BMC Medical Research Methodology},
  year    = {2026},
  volume  = {26},
  publisher = {BioMed Central},
  doi     = {10.1186/s12874-026-02963-5},
  url     = {https://doi.org/10.1186/s12874-026-02963-5},
  issn    = {1471-2288}
}

# Example Data 
Pairwise meta-analysis example data: https://search.r-project.org/CRAN/refmans/surrosurv/html/gastadv.html 

Network meta-analysis example data: [melanoma_ipd.csv](melanoma_ipd.csv)

# Pairwise meta-analysis
## Descriptive analyses
Pairwise meta-analysis descriptive analyses: [descriptive_pw_RScript.R](descriptive_pw_RScript.R)

## Constant and TVHR
Pairwise meta-analysis constant and TVHR: [tvhr_pw_RE_RScript.R](tvhr_pw_RE_RScript.R)

# Network meta-analysis 
## Descriptive analyses
Network meta-analysis descriptive analyses: [descriptive_nma_RScript.R](descriptive_nma_RScript.R)

## Constant HR  
Network meta-analysis constant HR - Fixed effect model: [ph_nma_FE_RScript.R](ph_nma_FE_RScript.R) 

Network meta-analysis constant HR - Random effects model: [ph_nma_RE_RScript.R](ph_nma_RE_RScript.R)

## TVHR
Network meta-analysis TVHR - Fixed effect model: [tvhr_nma_FE_RScript.R](tvhr_nma_FE_RScript.R)

Network meta-analysis TVHR - Random effects model: [tvhr_nma_RE_RScript.R](tvhr_nma_RE_RScript.R)

