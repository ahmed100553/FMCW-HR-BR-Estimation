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
