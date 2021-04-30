<p align="center">
  <a href="https://github.com/niupuhua1234/TRIMER"> <img width="2500px" src="https://github.com/niupuhua1234/TRIMER/blob/main/logo.png"></a> 
  <br />
  <br />
  <a href="https://www.mathworks.com/"><img alt="Python Version" src="https://img.shields.io/badge/MATLAB-%3E2016-brightgreen" /></a>
  <a href="https://www.mathworks.com/"><img alt="Python Version" src="https://img.shields.io/badge/R%20version-%3E3.0-orange" /></a>
  <a href="https://github.com//niupuhua1234/TRIMER/blob/main/LICENSE"><img alt="MIT License" src="https://img.shields.io/badge/license-MIT-blue.svg" /></a>
</p>
<!-- <div align="center">
     <a href="https://github.com/niupuhua1234/TRIMER"> <img width="2500px" src="https://github.com/niupuhua1234/TRIMER/blob/main/logo.png"></a> 
</div> -->

--------------------------------------------------------------------------------

# Welcome to TRIMER Library

**TRIMER** is a package for building genome-scale integrated metabolic–regulatory models based on Bayesian network. The integrated model can be used for tasks such as knockout phenotypes and  knockout flux prediction.


## Table of Contents

- [Introduction](#Introduction)
- [Dependency](#Dependency)
- [Usage](#usage)
- [Citation](#citation)

## Introduction
Transcriptional regulation plays a key role in controlling metabolism and a forefront challenge in modeling organisms today is to build integrated models of regulation and metabolism. Predicting the effect of transcriptional regulations on the metabolic network can lead to accurate predictions on how genetic mutations and perturbations are translated into flux responses at the metabolic level. TRIMER enables the quantitative integration of regulatory and metabolic networks to build genome-scale integrated metabolic–regulatory models. TRIMER is a Bayesian extension of PROM (Probabilistic Regulation of Metabolism).In TRIMER, transcriptional reguation on the metabolic network is represented by BN(Bayesian Network).

The construction of an integrated metabolic-regulatory network using TRIMER requires the following: 1) the reconstructed genome scale metabolic network 2) regulatory network structure, consisting of transcription factors (TF) and their targets 3) gene expression data.  We used TRIMER to build genome-scale models for various model organisms and showed that TRIMER can detect drug targets, identify gene knockout phenotypes with accuracies as high as 95% and predict microbial growth-rates of transcription factor knockout strains quantitatively with correlation of 0.96.

## Structure of the code
<p align="center">
  <a href="https://github.com/niupuhua1234/TRIMER"> <img width="500px" src="https://github.com/niupuhua1234/TRIMER/blob/main/code_structure.png"></a> 
</p>


## Prerequisites

1. [__CPLEX__](https://www.ibm.com/analytics/cplex-optimizer)
2. [__GLPK__](https://www.ibm.com/analytics/cplex-optimizer)
```

## Usage 
1. The raw gene expression data and binarized gene expression data  can be found in the google drive:
 https://drive.google.com/file/d/197DwrvBz8IMmwi3nTO64TV02_23gjIaT/view?usp=sharing, https://drive.google.com/file/d/1n0MDIhO17n7_Jy158euCU10kuySncbsp/view?usp=sharing .The metabolic model ***iAF260*** and interaction list are already saved  found in folder source_data for convenience.
 
3. ***(Under Any R environment)***  Run R script shown below for Bayesian network learning before estimating the regulatory bound for flux prediction. The BN learning is seperated from other part of the package as the process is time-consuming. The input are binarized gene expression data and interaction list which serve as prior knowledge for structure leaning. The learned BN is saved in **.bif** format which can be read by function read.bif in **bnlearn** package.
    ```text
    bnlearn.R
    ```
2. ***(Under Matlab Environment)*** The two matlab script shown below are demo codes of  knock-out flux predicton for indole and biomass.

    ```text
    flux_indole.m
    flix_biomass.m
    ```
    The following matlab script is  a demo code for phenotypes prediction.Phenotypes label, boolean network(***imc1010***) and tiger-trimer model  found in used in the demo code are saved in source_data folder.
    ```text
    prediction_phenotype.m
     ```
## Citation

If you find our library useful, please considering citing our paper in your publications. We provide a BibTeX entry below.

```bibtex
@article{niu2021trimer,
  title={TRIMER: Transcription Regulation Integrated with MEtabolic Regulation},
  author={Niu, Puhua and Soto, Maria J and Yoon, Byung-Jun and Dougherty, Edward R and Alexander, Francis J and Blaby, Ian and Qian, Xiaoning},
  journal={bioRxiv},
  year={2021},
  publisher={Cold Spring Harbor Laboratory}
}
