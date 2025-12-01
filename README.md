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

1. **Model repository** – It compiles most publicly available *B. subtilis* GEMs in one place.  
2. **Preparation pipeline (`1_Preparation`)** – Standardizes, validates, and formats GEMs to ensure they are comparable across tools and formats.  
3. **Comparison pipeline (`2_Comparison`)** – Performs structural and functional comparisons, flux balance analysis, carbon-source simulations, and gene essentiality predictions.

The pipeline has been tested with the following *B. subtilis* GEMs:

- **iYO844**  
- **iBB1018**  
- **iBsu1147R**

---

## Installation

Clone the repository:

```bash
git clone https://github.com/<your-username>/GEMcompare.git
cd GEMcompare

Install dependencies:
pip install -r requirements.txt
```
## Content

## Contact
The associated undergraduate thesis project can be found at:
<insert-link-here>

If you need help or have any questions about the pipeline, feel free to reach out:
📧 gr.gabriel@javeriana.edu.co

## License
The software licensed under the MIT license.