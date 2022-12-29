# HGCGANLDA
Predicting lncRNA-disease associations based on heterogeneous graph convolutional generative adversarial network.

Author: Zhonghao Lu.


## Overview
This repository contains codes necessary to run the HGCGANLDA algorithm. 

## Running Environment
* Windows environment, Python 3.7
* PyTorch >= 1.10.1

## Datasets
All datasets are available at Data.

The details of the data are shown in the table.

|              | lncRNA | Disease | miRNA | LDA  | MDA   | LMA   |
|--------------|--------|---------|-------|------|-------|-------|
| **Dataset1** | 861    | 432     | 673   | 4518 | 4189  | 2105  |
| **Dataset2** | 1363   | 501     | 1190  | 5338 | 6763  | 2291  |
| **Dataset3** | 240    | 412     | 495   | 2697 | 13562 | 1002  |
| **Dataset4** | 1723   | 236     | 675   | 1151 | 4634  | 10102 |

# Model
* model.py: the core model proposed in the paper.
* main.py: the main program in the project. Run the entire project by running main.py.
* get_adj.py: run the program to obtain the correlation matrix.
* get_k_mer_feat.py: run the program to obtain the sequence features of lncRNA.
* hetero_graph.py: construct and save heterogeneous graphs.
* train.py: training Model.

## Question
+ If you have any problems or find mistakes in this code, please contact with me: Zhonghao Lu: lzhao6688@163.com

