# eNet: an algorithm to build enhancer networks based on scATAC-seq and scRNA-seq data
## Introduction
eNet is an algorithm designed to integrate single-cell chromatin accessibility and gene expression profiles and build enhancer networks, delineating how multiple enhancers interact with each other in gene regulation. 

## Workflow
### Step1. Preparing input matrix (Input)
Two matrices are needed for the input of eNet. 1) scATAC-seq matrix (peak-cell); 2) scRNA-seq matrix (gene-cell).
### Step2. Identifying the putative enhancer cluster (Node)
Node <- FindNode(dis=200000)
### Step3. Identifying the predicted enhancer interactions (Edge)
Edge <- FindEdge()
### Step4. Building enhancer networks (Network)
NetworkInfo <- BuildNetwork(cutoff=0.1)
### Step5. Calculating network complexity (Network complexity)
NetworkComplexity <- NetworkComlexity()
### Step6. Classification of enhancer networks (Mode)
NetworkMode <- NetworkMode(Dcutoff=1, Ncutoff=5)

## How to cite eNet
1. Danni Hong#, Hongli Lin#, Lifang Liu, Muya Shu, Jianwu Dai, Falong Lu, Jialiang Huang*. Complexity of enhancer networks predicts cell identity and disease genes. (Submitted)
2. Muya Shu#, Danni Hong#, Hongli Lin, Jixiang Zhang, Zhengnan Luo, Yi Du, Zheng Sun, Man Yin, Yanyun Yin, Shilai Bao, Zhiyong Liu, Falong Lu*, Jialiang Huang*, Jianwu Dai*. Enhancer networks driving mouse spinal cord development revealed by single-cell multi-omics analysis. (Submitted)
