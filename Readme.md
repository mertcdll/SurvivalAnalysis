# Survival Analysis of Seed Germination Experiments

This repository contains R scripts used for the statistical analysis of seed germination experiments conducted in two experimental phases. The analyses focus on germination dynamics under different moisture treatments and provenances, using classical survival analysis and accelerated failure time (AFT) models.

The workflow integrates descriptive summaries, exploratory visualization, correlation analysis, Kaplan–Meier estimators, and parametric survival models.

---

## Overview of the Study Design

* **Experiment phases**

  * **Phase 1**: Germination monitored for up to 71 days
  * **Phase 2**: Germination monitored for up to 103 days (with a short 71-day subset for comparison)

* **Species**

  * *Abies alba*
  * *Abies nordmanniana*
  * *Fagus sylvatica*
  * *Fagus orientalis*

* **Treatments**

  * Moist vs. dry moisture treatments
  * Multiple provenances per species

* **Experimental unit**

  * Individual seeds (100 seeds per tray)
  * Germination time recorded as time-to-event data
  * Non-germinated seeds treated as right-censored observations


---

## Required R Packages

The analysis relies on the following R packages:

```r
library(survival)
library(ggplot2)
library(flexsurv)
library(dplyr)
library(tibble)
library(survminer)
library(PerformanceAnalytics)
```

All packages are available from CRAN.

---

## Analysis Workflow

### 1. Data Import and Preprocessing

* Phase 1 and Phase 2 datasets are imported separately.
* Tray identifiers are made unique across phases.
* Metadata and cumulative germination counts are separated from individual seed observations.

### 2. Germination Overview

* Final germination counts and percentages are extracted per provenance.
* Moist and dry treatments are combined where appropriate.
* Seed dry weight and post-stratification moisture content are added to the summary tables.

### 3. Exploratory Visualization

* Time-series plots of cumulative germination counts by species, provenance, and treatment.
* Visual inspection of germination dynamics across experimental phases.

### 4. Correlation Analysis

* Pearson correlation analysis between:

  * Final germination counts
  * Seed dry weight
  * Moisture content
* Results indicate no strong linear association between final germination success and seed weight or moisture content.

### 5. Survival Data Construction

* Individual seed-level datasets are constructed.
* Germination time is treated as the event time.
* Non-germinated seeds are right-censored at:

  * Day 72 (Phase 1 / short Phase 2)
  * Day 104 (long Phase 2)

### 6. Kaplan–Meier Survival Analysis

* Survival curves (probability of not germinating) are estimated:

  * Separately for *Abies* and *Fagus*
  * Stratified by provenance
* Results are visualized using Kaplan–Meier plots.

### 7. Accelerated Failure Time (AFT) Models

* Parametric AFT models are fitted using:

  * Exponential
  * Weibull
  * Log-normal
  * Log-logistic
  * Gaussian distributions
* Model selection is based on AIC.
* Log-normal models provide the best fit for both *Abies* and *Fagus*.
* Final models include provenance and moisture treatment as covariates.

### 8. Phase 2 Comparison

* Short (71-day) and long (103-day) Phase 2 observations are compared.
* AFT models are fitted separately to assess the impact of observation window length on parameter estimates.

---

## Outputs

* Cleaned and combined data tables for downstream analysis
* Publication-ready figures for germination dynamics and survival curves
* Parametric survival model summaries with effect estimates and confidence intervals

---

## Notes

* Provenances with no or extremely low germination are excluded from parametric modeling to ensure model stability.
* The scripts are written as a complete, linear analysis pipeline and are intended to be run top-to-bottom.

---

## Author

Mert Çelik
PhD-level analysis script prepared for academic research in plant ecology and seed biology.

---
