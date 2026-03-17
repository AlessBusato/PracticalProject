# PracticalProject

This repository implements a **univariate twin modeling analysis (ACE model)** on **Nested Spectral Partitioning (NSP)** measures of **Functional Connectivity** using data from the **Human Connectome Project (HCP)**.

## Overview
The goal of this project is to investigate the **genetic and environmental contributions** to functional brain connectivity. Using **Nested Spectral Partitioning (NSP)** metrics derived from functional connectivity networks, we apply a **classical twin design (ACE model)** to decompose variance into:

- **A (Additive Genetic Effects)** – heritable influences shared more strongly by monozygotic twins.
- **C (Common/Shared Environment)** – environmental factors shared between twins.
- **E (Unique Environment)** – individual-specific environmental influences and measurement error.

## Data
The analysis uses neuroimaging data from the **Human Connectome Project (HCP)**, which includes monozygotic (MZ) and dizygotic (DZ) twin pairs. Integration and Segregation measures are computed and summarized using **Nested Spectral Partitioning (NSP)**. 

## Methods
1. Extract **NSP-derived metrics** from FC matrices.
2. Prepare excel file with NSP measures for **twin pairs (MZ and DZ)**.
3. Apply a **univariate ACE twin model** to estimate genetic and environmental variance components.
