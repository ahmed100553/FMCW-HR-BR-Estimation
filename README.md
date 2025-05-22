# FCMW-HR-BR-Estimation
I present an end-to-end framework for non-contact monitoring of heart rate (HR) and breath- ing rate (BR) in children using a compact 60 GHz FMCW radar. Two parallel signal-processing branches are investigated

## Overview
This repository contains all code, data stubs, and documentation for the project  
**“Contact-less Heart- and Breath-Rate Estimation in Children with a 60 GHz FMCW Radar.”**  
The project delivers an **end-to-end MATLAB framework** that extracts heart-rate (HR) and breathing-rate (BR) traces from raw radar IQ data, compares two alternative back-ends (fixed Butterworth band-passes vs. adaptive MODWT separation), and reproduces the figures and metrics reported in the accompanying technical report.

Key contributions    
* A child-friendly data-acquisition protocol and an open 250-minute dataset (50 children, 7–13 yr).  
* A modular eight-stage MATLAB pipeline with shared radar front-end and two interchangeable physiological back-ends  
* Quantitative evaluation on **25 000** radar frames with instantaneous error traces and rank-ordered error curves
## 📦 Accessing the Child Vital-Sign Dataset

The radar IQ files and synchronized reference traces used in this project are **publicly available on Figshare** as part of the dataset that accompanies the open-access paper by Yoo *et al.* (2021):contentReference[oaicite:0]{index=0}.

| Resource | Link |
|----------|------|
| Paper | <https://pmc.ncbi.nlm.nih.gov/articles/PMC8036835/> |
| Dataset landing page | <https://figshare.com/articles/dataset/13515977> |

> **Licence** The dataset is distributed under the Creative Commons CC-BY 4.0 licence. Please cite the original work when you use these data.

---

### Step-by-step download

1. Open the dataset page on Figshare.  
2. Click **“Download all”** (this fetches a single ZIP archive containing three folders: `Rawdata`, `Vital sign`, `Participants`).  
3. Unzip the archive into `data/raw/` inside this repository. Your tree should now look like

   ```text
   data/
     └── raw/
         ├── Rawdata/
         ├── Vital sign/
         └── Participants/
