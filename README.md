# GEMcompare

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Content](#content)
- [Contact](#contact)
- [License](#license)

---

## Overview

**GEMcompare** is a two-stage computational pipeline designed for the preparation and comparison of **genome-scale metabolic models (GEMs)**, with a focus on *Bacillus subtilis*.  
This repository also functions as a centralized collection of curated GEMs for this organism.

This repository serves three main purposes:

1. To centralize most of the available B. subtilis genome-scale metabolic models (GEMs). 2. A pipeline "1_Preparation" that can be used to prepare GEMs by standardizing so they can be compared. 3. A pipeline "2_Comparison" to compare GEMs.
1. **Model repository** – It compiles most publicly available *B. subtilis* GEMs in one place.  
2. **Preparation pipeline (`1_Preparation`)** – Standardizes GEMs to ensure they are comparable. 
3. **Comparison pipeline (`2_Comparison`)** – Performs different comparative analysis on GEMs.

The pipeline has been tested with the following *B. subtilis* GEMs:

- **iYO844**  
- **iBB1018**  
- **iBsu1147R**

To modify the parameters in the pipeline, it is recommended to use the `config.yaml` file

---

## Installation

Clone the repository:

```bash
git clone https://github.com/<felipeg210303>/GEMcompare.git
cd GEMcompare

Install dependencies:
pip install -r requirements.txt
```

---

## Content

### **InitialModels/**
Contains the raw GEMs exactly as obtained from their original sources (various formats, naming conventions, and annotation styles). These files are kept unmodified.

### **Models/**
Includes the curated models converted to **SBML** and standardized to the **BiGG ID** naming system.  
Note: Only *iYO844*, *iBB1018*, and *iBsu1147R* have undergone extensive manual curation; other models may still contain inconsistent identifiers.

### **Pipeline/**
The computational workflow of GEMcompare, divided into two stages:

### 1. Preparation
The preparation stage ensures that all models share a common structure, naming system, and simulation environment.

#### **1.1 Translation**
Applies unified BiGG-style identifiers for reactions, metabolites, and genes using `mergem`. Note: Should be skipped if the model is already with BiGG IDs, otherwise it will break them. 

#### **1.2 Validation**
Evaluates model integrity by checking internal consistency, SBML compatibility and syntaxis.

#### **1.3 Quality Control**
Runs `memote`, a metabolic model testing package which returns a quality report with a quality score. This quality reports are in the foler QualityReports.

#### **1.4 Set Media**
Applies a standardized growth medium used for all models, ensuring that comparisons are performed under identical environmental conditions.

### 2. Comparison
The comparison stage performs different analyses across the curated GEMs. Results are shown in graphs.

#### **2.1 Structural Comparison**
Gives genes, metabolites and reactions for all models.

#### **2.2 Functional**
Runs FBA simulations to evaluate flux distributions and pathway activity in central carbon metabolism across models.

#### **2.3 Different carbon sources**
Performs FBA with biomass as objective function and gives growth rates under different carbon sources.

#### **2.4 Gene essetiality analysis**
Performs single gene knockouts to find which genes are essential for growth. Also compares with a dataset of known essential genes to determine true positives (TN), false positives (FP), true negatives (TP) and false negatives (FN). Using those data metrics such as precision, recall, Matthews Correlation Coefficient (MCC), etc. are calculated.

#### **2.4 Carbon utilization analysis**
Alternates the carbon source and performs FBA to find which carbon sources can support growth. Also compares with a dataset of known carbon source utilization and determines true positives (TN), false positives (FP), true negatives (TP) and false negatives (FN).


### **Quality Reports/**
This folder contains the quality reports generated with `memote`.

### **Outputs/**
Contains the outputs from each of the analysis ran.

---

## Contact
The associated undergraduate thesis project can be found at:
<***>

If you need help or have any questions, feel free to reach out:
<gr.gabriel@javeriana.edu.co>

---

## License
The software is licensed under the MIT license.