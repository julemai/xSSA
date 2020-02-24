# Extended Sobol' Sensitivity Analysis xSSA: Simultaneous Global Sensitivity Analysis of Model Parameters and Model Structure
*by Juliane Mai, James R Craig, and Bryan A Tolson (University of Waterloo, Canada)*

## Abstract
Model structure uncertainty is known to be one of the three main sources of hydrologic model uncertainty along with input and parameter uncertainty. Some recent hydrological modeling frameworks address model structure uncertainty by supporting multiple options for representing hydrological processes. It is, however, still unclear how best to analyze structural sensitivity using these frameworks. In this work, we propose an *Extended Sobol’ Sensitivity Analysis (xSSA)* method that operates on grouped parameters rather than indi- vidual parameters. The method can estimate not only traditional model parameter sensitivities but is also able to provide measures of the sensitivities of process options (e.g., linear vs. non-linear storage) and sensitivities of model processes (e.g., in ltration vs. base ow) with respect to a model output. The method is applied to both am arti cial benchmark model and a watershed model built with the Raven framework. The results show that: 
1. The xSSA method provides sensitivity estimates consistent with those derived analytically. 
2. The xSSA method with process weighting is computationally less expensive than a sensitivity analysis performed for individual models, with savings of 82% for the benchmark model and 98% for the watershed case study. 
3. The xSSA method applied to the case study showed that surface processes were responsible for 80.3% of the overall model variability in a mountainous catchment; such information may readily inform model calibration. 
4. The analysis of time dependent process sensitivities is a helpful tool to understand model internal dynamics over the course of the year.

## Step-by-step tutorial
The step-by-step tutorial describes all the steps to estimate sensitivity indexes for a model that allows for multiple process options such as the hydrologic modeling framework [RAVEN](http://raven.uwaterloo.ca). It explains step-by-step on how to setup and adjust the provided codes. Details can be found [here](https://github.com/julemai/xSSA/wiki/Step-by-step-tutorial).

## Examples
We provide the examples that were used in our publication:
- Benchmark model (see [here](https://github.com/julemai/xSSA/wiki/Examples#benchmark-model))
- Salmon River catchment using RAVEN modeling framework (see [here](https://github.com/julemai/xSSA/wiki/Examples#raven-hydrologic-modeling-framework-Salmon-River))
- MOPEX 12 using RAVEN modeling framework (see [here](https://github.com/julemai/xSSA/wiki/Examples#raven-hydrologic-modeling-framework-MOPEX-12))

## Setup your own model
This section provides you with a guideline on how to setup your own model and run a sensitivity analysis for your model. Details can be found [here](https://github.com/julemai/xSSA/wiki/Setup-your-own-model).

## Citation
J Mai, JR Craig, and BA Tolson (2019). <br>
Simultaneous Global Sensitivity Analysis of Model Variables and Model Structures. <br>
*Water Resources Research*, ??, ??.<br>
https://doi.org/??/??
