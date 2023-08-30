# The overlooked effect of amplitude on within-speaker vowel variation
**Code and analysis repository**  
**NZILBB Vowel Clustering Project**

## Overview

This repository contains analysis and code for the paper "The overlooked effect
of amplitude on within-speaker vowel variation". The bulk of this repository
is made up of four R markdown files corresponding to each stage of the analysis
in the paper and found in the `supplementary_material` folder along with knitted
html versions which can be read without rerunning the code yourself.

## Directory structure

- `supplementary material`: contains the primary R markdown files running through the four stages of our analysis:
  1. preprocessing (hosted [here](https://nzilbb.github.io/amp_f1_public/supplementary_material/SM1_preprocessing.html))
  2. interval representation (hosted [here](https://nzilbb.github.io/amp_f1_public/supplementary_material/SM2_interval_representation.html)),
  3. corpus level PCA analysis (hosted [here](https://nzilbb.github.io/amp_f1_public/supplementary_material/SM3_corpus_pca.html)), and
  4. GAMM models (hosted [here](https://nzilbb.github.io/amp_f1_public/supplementary_material/SM4_models.html)).
- `labbcat_data`: Folder for anonymised raw data from LaBB-CAT. **Empty directory. To get data go to <https://osf.io/m8nkh/> and download the file `merged_data.rds` and place it in this directory.**
- `scripts`: Contains R scripts for running permutation tests referred to in the supplementary materials.
- `processed_data`: Data output at various stages of the analysis by the scripts of supplementary materials.
- `plots`: Plots output by scripts or supplementaries.
- `models`: Folder for models reported in paper. **Empty directory. To get data go to <https://osf.io/m8nkh/> and download the files in the `models` directory.**

## How to use this repository

The above files are sufficient to repeat the analysis we have carried out. Some
models take a long time to fit. You can download our models
[here](https://www.dropbox.com/sh/vb167ng4v7vkfiv/AACbyEr-KeJGaZTCe_3G2BOfa?dl=0).
Place these in the project directory in a directory named `models`.

It is best to open this repository as an Rproject. One way to do this in RStudio
is to go to 'File -> New project', select 'Existing directory', and navigate to
the downloaded repository.
wherever
