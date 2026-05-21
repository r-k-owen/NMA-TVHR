# NMA-TVHR 
Pairwise (PWMA) and network meta-analysis (NMA) of time to event outcomes assuming non-proportional hazards using a time varying hazard ratio approach.

# Citation
Please cite: https://doi.org/10.48550/arXiv.2605.20567.

@article{owen2026tvhr,
      title={Meta-analysis and network meta-analysis of time-to-event outcomes with non-proportional hazards: a Bayesian time-varying hazard ratio approach}, 
      author={Rhiannon K Owen and Keith R Abrams},
      year={2026},
      eprint={2605.20567},
      archivePrefix={arXiv},
      primaryClass={stat.ME},
      url={https://arxiv.org/abs/2605.20567}, 
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

