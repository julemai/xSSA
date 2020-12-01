# Simultaneously Determining Global Sensitivities of Model Parameters and Model Structure
*by Juliane Mai, James R Craig, and Bryan A Tolson (University of Waterloo, Canada)*

## Abstract
Model structure uncertainty is known to be one of the three main sources of hydrologic model uncertainty along with input and parameter uncertainty. Some recent hydrological modeling frameworks address model structure uncertainty by supporting multiple options for representing hydrological processes. It is, however, still unclear how best to analyze structural sensitivity using these frameworks. In this work, we apply the Extended Sobol’ Sensitivity Analysis (xSSA) method that operates on grouped parameters rather than individual parameters. The method can estimate not only traditional model parameter sensitivities but is also able to provide measures of the sensitivities of process options (e.g., linear vs. non-linear storage) and sensitivities of model processes (e.g., infiltration vs. baseflow) with respect to a model output. Key to the xSSA method’s applicability to process option and process sensitivity is the novel introduction of process option weights in the Raven hydrological modeling framework. The method is applied to both artificial benchmark models and a watershed model built with the Raven framework. The results show that: 
1. The xSSA method provides sensitivity estimates consistent with those derived analytically for individual as well as grouped parameters linked to model structure.
2. The xSSA method with process weighting is computationally less expensive than the alternative aggregate sensitivity analysis approach performed for the exhaustive set of structural model configurations, with savings of 81.9% for the benchmark model and 98.6% for the watershed case study.
3. The xSSA method applied to the hydrologic case study analyzing simulated streamflow showed that model parameters adjusting forcing functions were responsible for 42.1% of the overall model variability while surface processes cause 38.5% of the overall model variability in a mountainous catchment; such information may readily inform model calibration and uncertainty analysis.
4. The analysis of time dependent process sensitivities regarding simulated streamflow is a helpful tool to understand model internal dynamics over the course of the year.

## Examples
We provide the examples that were used in our publication:
- Benchmark model (see [here](https://github.com/julemai/xSSA/wiki/Examples#benchmark-model))
- Salmon River catchment using RAVEN modeling framework (see [here](https://github.com/julemai/xSSA/wiki/Examples#salmon-river-catchment-using-raven-hydrologic-modeling-framework))

## Citation

### Journal Publication
Mai, J., Craig, J. R., and Tolson, B. A. (2020).<br>
Simultaneously Determining Global Sensitivities of Model Parameters and Model Structure. <br>
*Hydrol. Earth Syst. Sci.*, accepted. <br>
https://doi.org/10.5194/hess-2020-215

### Code Publication
Mai, J., Craig, J. R., and Tolson, B. A. (2020).<br>
Simultaneously Determining Global Sensitivities of Model Parameters and Model Structure. <br>
*Zenodo*<br>
https://doi.org/10.5281/zenodo.???

