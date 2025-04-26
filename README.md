# Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information*

This repository contains all the codes needed to replicate the results included in the paper.

**Authors:** Juan David Herreño and Carlos Rondón-Moreno.

Comments and suggestions are welcome at: crondon[at]bcentral[dot]cl.

## Table of Contents

- [Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information\*](#overborrowing-and-systemic-externalities-in-the-business-cycle-under-imperfect-information)
  - [Table of Contents](#table-of-contents)
  - [0. Folder Structure](#0-folder-structure)
  - [1. System Requirements](#1-system-requirements)
  - [2. Instructions for Replication (Main Results)](#2-instructions-for-replication-main-results)
    - [2.1 Step0\_Parameters.m](#21-step0_parametersm)
    - [2.2 Step1\_transition\_matrix.m](#22-step1_transition_matrixm)
    - [2.3 Step2\_dctrlzd\_eqm\_imp\_info.m](#23-step2_dctrlzd_eqm_imp_infom)
    - [2.4 Step3\_planner\_imp\_info.m](#24-step3_planner_imp_infom)
    - [2.5 Step4\_dctrlzd\_eqm\_pfct\_info.m](#25-step4_dctrlzd_eqm_pfct_infom)
    - [2.6 Step5\_planner\_pfct\_info.m](#26-step5_planner_pfct_infom)
    - [2.7 Step6\_Welfare\_Costs.m](#27-step6_welfare_costsm)
  - [3. Figures](#3-figures)
  - [4. Tables](#4-tables)
  - [5. Calibration](#5-calibration)

## 0. Folder Structure

The folder structure is as follows:

```
├── README.md
├── Codes: 
├──├──Figures
├──├──├──Appendix
├──├──Tables
├──├──Replication
├──├──├──Calibration
├──├──├──Data
├──├──├──Step0_Parameters.m
├──├──├──Step1_transition_matrix.m
├──├──├──Step2_dctrlzd_eqm_imp_info.m
├──├──├──Step3_planner_imp_info.m
├──├──├──Step4_dctrlzd_eqm_pfct_info.m
├──├──├──Step5_planner_pfct_info.m
├──├──├──Step6_Welfare_Costs.m
```

## 1. System Requirements

The model results were estimated on the University of Notre Dame's computer cluster. Executing the codes with the baseline calibration's grid sizes requires approximately 1 TB of RAM.

**IMPORTANT**: If you wish to bypass model solving, please refer to points 4 and 5 below for instructions on downloading our simulated data.

The codes are for MATLAB.

## 2. Instructions for Replication (Main Results)

To replicate the full set of results, run `main_file.m`, which will create the necessary parameter and transition structures for solving the model and replicating our baseline findings. 

Modifying parameters in `Step0_Parameters.m` is required for the sensitivity analysis.

**IMPORTANT:** Refer to points 3 and 4 below for instructions on replicating figures and tables. 

All results from the subsequent steps should be saved in:

```
./Replication/Data
```

### 2.1 Step0_Parameters.m

This file contains all parameters for the model solution, including calibration parameters and constants used in simulations. Note that reproducing the sensitivity analysis requires modifying the parameters within this file.

### 2.2 Step1_transition_matrix.m

This file discretizes the VAR process used for the calibration. The result is two *.mat* files. The first one is the matrix used for the full information model. The second one is used for the model with imperfect information.

For this file, you need two ancillary files:

```
- tpm.m: Stephanie Smith-Grohé and Martín Uribe's routine to discretize VAR(1) processes.
- KF_transition_matrix.m: Adjusts the imperfect information transition matrix to reflect that the agent only observes aggregate income and not its components.
```

### 2.3 Step2_dctrlzd_eqm_imp_info.m

In this code, we use Policy Function iteration to solve the decentralized equilibrium model with imperfect information. The program also produces the simulations of the model.

### 2.4 Step3_planner_imp_info.m

In this code, we use Value Function iteration to solve the decentralized equilibrium model with imperfect information. The program also produces the simulations of the model.

### 2.5 Step4_dctrlzd_eqm_pfct_info.m

In this code, we use Policy Function iteration to solve the decentralized equilibrium model with perfect information. The program also produces the simulations of the model.

### 2.6 Step5_planner_pfct_info.m

In this code, we use Value Function iteration to solve the decentralized equilibrium model with perfect information. The program also produces the simulations of the model.

### 2.7 Step6_Welfare_Costs.m

With the .mat files containing the solutions from the previously solved models, this code simulates the value functions and the welfare costs for the perfect and imperfect information models.

## 3. Figures

To replicate the figures, you have two options. The first one is to run the file `main_file.m`, which will generate all the `.mat` files needed. The second option is to download our data from:

[http://tiny.cc/HR2025rep](http://tiny.cc/HR2025rep)

Download the data and place it in the folder:

```
./Replication/Data
```

Then, you can run run the individual files located in the Figures folder. Subfolder Appendix contains the codes to replicate the figures in the appendix.

## 4. Tables

The Tables folder contains the files that generate each table in the paper. Similar to the figures, you can either run `main_file.m` to generate the `.mat` files or download our data from the [link above](http://tiny.cc/HR2025rep) and place it in the ``Data`` folder.

## 5. Calibration

The files needed to replicate the calibration are located in:

```
./Replication/Calibration
```

See the README file in that folder for detailed instructions.
