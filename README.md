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

- [Security](#security)
- [Introduction](#Introduction)
- [Install](#install)
- [Usage](#usage)
- [API](#api)
- [Contributing](#contributing)

## Introduction
Transcriptional regulation plays a key role in controlling metabolism and a forefront challenge in modeling organisms today is to build integrated models of regulation and metabolism. Predicting the effect of transcriptional regulations on the metabolic network can lead to accurate predictions on how genetic mutations and perturbations are translated into flux responses at the metabolic level. TRIMER enables the quantitative integration of regulatory and metabolic networks to build genome-scale integrated metabolic–regulatory models. TRIMER is a Bayesian extension of PROM (Probabilistic Regulation of Metabolism).In TRIMER, transcriptional reguation on the metabolic network is represented by BN(Bayesian Network).

The construction of an integrated metabolic-regulatory network using TRIMER requires the following: 1) the reconstructed genome scale metabolic network 2) regulatory network structure, consisting of transcription factors (TF) and their targets 3) gene expression data.  We used TRIMER to build genome-scale models for various model organisms and showed that TRIMER can detect drug targets, identify gene knockout phenotypes with accuracies as high as 95% and predict microbial growth-rates of transcription factor knockout strains quantitatively with correlation of 0.96.

## Structure of the code
<p align="center">
  <a href="https://github.com/niupuhua1234/TRIMER"> <img width="500px" src="https://github.com/niupuhua1234/TRIMER/blob/main/code_structure.png"></a> 
</p>


## Install


```
```


## Usage 
1. The raw gene expression data and binarized gene expression data  can be found in the google drive:
 https://drive.google.com/file/d/197DwrvBz8IMmwi3nTO64TV02_23gjIaT/view?usp=sharing, https://drive.google.com/file/d/1n0MDIhO17n7_Jy158euCU10kuySncbsp/view?usp=sharing
3. ***(Under Any R environment)***  Run R script shown below for Bayesian network learning before estimating the regulatory bound for flux prediction. The BN learning is seperated from other part of the package as the process is time-consuming. The Input is binarized gene expression data and interaction list as prior knowledge for structure leaning. The learned BN is saved in **.bif** format which can be read by function read.bif in **bnlearn** package.
    ```text
    bnlearn.R
    ```
2. ***(Under Matlab Environment)*** The two maltba script shown below are demo codes of  Knock-out flux predicton for indole and biomass.

    ```text
    $ python Extract-Raw-Data-Into-Matlab-Files.py
    ```

3. Preprocessed the Dataset via the Matlab and save the data into the Excel files (training_set, training_label, test_set, and test_label) via [these scripts](https://github.com/SuperBruceJia/EEG-DL/tree/master/Preprocess_EEG_Data) with regards to different models. FYI, every lines of the Excel file is a sample, and the columns can be regarded as features, e.g., 4096 columns mean 64 channels X 64 time points. Later, the models will reshape 4096 columns into a Matrix with the shape 64 channels X 64 time points. You should can change the number of columns to fit your own needs, e.g., the real dimension of your own Dataset.

4. ***(Prerequsites)*** Train and test deep learning models **under the Python 3.6 Environment (Highly Recommended)** for EEG signals / tasks classification via [the EEG-DL library](https://github.com/SuperBruceJia/EEG-DL/tree/master/Models), which provides multiple SOTA DL models.

    ```text
    Python Version: Python 3.6 (Recommended)
    TensorFlow Version: TensorFlow 1.13.1
    ```

    Use the below command to install TensorFlow GPU Version 1.13.1:

    ```python
    $ pip install --upgrade --force-reinstall tensorflow-gpu==1.13.1 --user
    ```

5. Read evaluation criterias (through iterations) via the [Tensorboard](https://www.tensorflow.org/tensorboard). You can follow [this tutorial](https://www.guru99.com/tensorboard-tutorial.html). When you finished training the model, you will find the "events.out.tfevents.***" in the folder, e.g., "/Users/shuyuej/Desktop/trained_model/". You can use the following command in your terminal:

    ```python
    $ tensorboard --logdir="/Users/shuyuej/Desktop/trained_model/" --host=127.0.0.1
    ```

    You can open the website in the [Google Chrome](https://www.google.com/chrome/) (Highly Recommended). 
    
    ```html
    http://127.0.0.1:6006/
    ```

    Then you can read and save the criterias into Excel .csv files.

6. Finally, draw beautiful paper photograph using Matlab or Python. Please follow [these scripts](https://github.com/SuperBruceJia/EEG-DL/tree/master/Draw_Photos).


Note: The `license` badge image link at the top of this file should be updated with the correct `:user` and `:repo`.

### Any optional sections

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
