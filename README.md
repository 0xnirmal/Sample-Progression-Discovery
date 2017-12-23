
# Advanced Topics in Genomic Data Analysis- Mini Project 2 #

# Part II: Reproducibility #

## Running Instructions ##
```
git clone https://github.com/nkrishn9/Sample-Progression-Discovery.git
cd Sample-Progression-Discovery/
chmod -x run.sh
source run.sh
```
This code is intended to run on the JHU CS ugrad cluster. Figures 1 and 2 were generated from the Jupyter notebook provided; however, the raw data used to generate figure 2 is printed when the py file is run. 

If you would like to see my implementation, I recommend taking a peek at the spd.ipynb (using the jupyter notebook command), not the spd.py because the notebook contains significant markdown comments. 

## Paper Selection ##
I chose to replicate sample progression discovery (SPD) method from ["Discovering Biological Progression Underlying Microarray Samples"](https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/spd.PDF) paper by Qiu et al. Given that I will soon need to run this process using their Matlab software package provided through their [website](http://pengqiu.gatech.edu/software/SPD/index.html) on the iPSC data, I thought this would be a useful paper to select in order to broaden my understanding of how the process works. 

The goal of SPD is "to discover biological progression from gene expression microarray data." It accomplishes this by four steps: 
1) cluster genes into modules of co-expressed genes
2) construct minimum spanning tree (MST) for each module
3) select modules that support common MSTs
4) reconstruct an overall MST based on all the genes of all the selected modules

### Figure 0: SPD Methods ###
![figure0]

Since this is a fairly long, multifaceted process and we were only given four hours- I chose to replicate steps 4 and 1 and compare with the results from the provided software package. The software package includes a simulated data file on cell progression, which I used in all of my methods. 

## Replicating Step 4- Reconstructing an overall MST based on all the genes of all the selected modules ##

I began with step 4 because it is fairly independent of all the other steps. Once the modules have been selected, according to the paper, you simply translate the reduced data matrix [feature selected genes x samples] to a graph representation with nodes as samples and edge weights as Euclidean distance between samples which you use to produce a MST. 

To create the modules, I used the software package on the simulated data. I selected 6 modules (using the process described in the paper- visual inspection), which provided me with a list of 630 genes to use. I outputted this list of genes from the Matlab software package to my Jupyter notebook and created the reduced matrix accordingly. Using the reduced matrix, I translated this to the graph representation fairly easily. 

In order to create the overall MST using the graph representation, the paper uses Boruvka’s algorithm, a standard greedy algorithm. Using the software package with the selected modules, the resulting progression graph can be seen below.

### Figure 1: Overall MST using Matlab Software Package ###
![figure1]

As you can see, the software package is fairly successful at recovering the true progression of the cells. It only breaks continuity where there are edges from 1 to 2 and from 1 to 3 instead of from 1 to 2 and from 2 to 3 as well as edges from 14 to 16 and from 16 to 15 instead of eges from 14 to 15 and from 15 to 16. Other than these instances, the resulting graph completely captures the effect of time on progression. 

### Figure 2: Overall MST using my code ####
![figure2]

My results replicate the celluar progression well--there is full continuity from time step 0 to 14 with a jump to 16 after 14 and then a reversion to 15. This is nearly identical to the results in the provided code.  I spent a significant amount of time trying to understand why the MSTs generated were different, but couldn't seem to identify the difference between the two methods--both my code and the provided software package use Boruvska's algorithm with Euclidean distance as the edge weight to generate the MST. I decided to move on and try to replicate step 1 as well, since I did not want to use all four hours on debugging this.

## Replicating Step 1- Cluster genes into modules of co-expressed genes ##
SPD uses a augmentation of k-means as a clustering method for co-expressed described directly in the paper: 
> Given an N by M gene expression data matrix, we perform the k-means algorithm L times, with random initialization, to cluster the N genes into k = 2 clusters. Clustering results are arranged into an N by L matrix, where the (i,j) element is the cluster assignment of gene i in the j’th run of k-means. In order to draw the consensus of the L runs of k-means, we apply k-means again based on the N by L matrix, the collection of clustering results of the L runs, to divide genes into two clusters. For each of the two clusters, the coherence is computed as the average Pearson correlation between each gene in the cluster and the cluster mean. If the coherence of a cluster is higher than a pre-specified threshold c1, this cluster is considered to be a coherent gene module. Otherwise, this cluster is further partitioned by iterating the algorithm. After the iterative process ends, we examine the resulting coherent modules pairwisely. If the Pearson correlation of two modules’ centers is higher than a pre-specified threshold c2, these two modules are merged. This step iterates until no module-pair shares correlation higher than c2. The stopping criterion of cluster coherence guarantees that all resulting modules satisfy the pre-specified coherence threshold c1. Modules that share correlation higher than c2 are merged, so that the resulting gene modules are not highly correlated with each other. We typically set the algorithm parameters to the following values: L~200, c1 ~0:7, c2 ~0:9:

Using the software package, there were 54 resulting co-expressed gene modules generated. I attempted to implement the code in the exact way described above; however, my results yielded over 600. I think this may be due to the ambiguity with how the paper worded the recombination process of modules: "If the Pearson correlation of two modules’ centers is higher than a pre-specified threshold c2, these two modules are merged." I assumed this meant to take the means of modules at each time step and then find the correlation among the means of the two modules, but this is slightly ambiguous. It is also possible I have a bug in my code in this step; however, my four hours were up after spending 20 minutes attempting to debug this.


## Reflection on Replication ##
It was challenging to reproduce the SPD results on the simulated data even given the source code used by the researchers. Given only four hours, it was difficult to understand why my results and the provided results did not match, primarily because of the disorganization of the provided code base and lacking documentation. However, I think given more time it is reasonable that I could have better understood the provided source code and found the disparities in my own code accordingly. 

I think this exercise has taught me a lot about code structuring and documentation--if these are lacking, it makes reproducing your study very difficult with limited time constraints. 

# Part III: Final Project Proposal #
For my final project, I would like to evaluate the predictive accuracy of Gaussian Processes (GP) on gene expression data for the iPSC dataset. GPs have quickly gained acclaim in the scientific community for their flexibility and predictive accuracy; however, I am curious as to how they compare to traditional models in this domain. Specifically, the iPSC dataset contains daily (t=0 to t=15) gene expression data for iPSCs undergoing differentiation into cardiomycte cells. In order to determine how well this model can predict expression, for each gene and cell line, I will hold a single time step out and train a model using all N-1 data points between time (the x) and expression (the y). Then, using this model, I will predict expression of this gene at that time step and repeat for all time steps. I will do this for all genes x cell lines and then compute error statistics accordingly. This same process will be repeated with traditional models like linear regression and the results will be compared.

I plan on varying my kernel and hyperparameters for the GPs and testing across a number of configurations. However, the real goal of this project is to determine whether the flexibility of the GP enables great predictive accuracy, or whether it dramatically encourages overfitting. 

In figure 3, we can see the flexibility the GP has to fit the data. Figure 3 contains a fully fitted model for both a GP (RBF kernel) and a basic linear regression. As we can see, the GP perfectly fits all the points while the linear regression in this situation struggles to reconsole the non-linear looking progression. 

### Figure 3 ###
![figure3]

However, in figure 4, we now hold a sample out and train on all other samples and then predict on the held out sample. In this situation, it appears the performance between the GP and linear model is significantly closer, indicating that GP may be overfitting and a linear model may perform better in this situation.

### Figure 4 ###
![figure4]

This study aims to take a deeper look at this question in this domain. 

# Part I: Rewriting Methods #
Original excerpt from ["RSEM: accurate transcript quantification from
RNA-Seq data with or without a reference
genome"](https://d1b10bmlvqabco.cloudfront.net/attach/j6zot2yz1ti44r/j783o2sgiadfq/j9e774zkefp6/RSEM_Transcript_Quantification.pdf) by Li et al:
> Gibbs sampling: 

> In addition to computing ML estimates, RSEM uses a Bayesian version of its model to compute PME and 95% CIs of abundances. In the Bayesian model, the θ parameters are treated as latent random variables with a Dirichlet prior distribution. The parameters of the Dirichlet distribution (a) are set to one, which makes the prior equivalent to a uniform distribution and the maximum a posteriori estimates of θ equal to the ML estimates. 

> RSEM computes PMEs and 95% CIs with a two-stage sampling process. First, a standard application of the collapsed Gibbs sampling algorithm is used to obtain a sampled set of count vectors, where each vector represents the number of fragments that are mapped to each transcript. During each round of the Gibbs sampling algorithm, the true mapping of each fragment is resampled given the current mappings of all other fragments. The initial mapping of each fragment is sampled according to the ML parameters computed by the EM algorithm. The algorithm is run to sample 1000 count vectors. 

Rewritten excerpt:
> Gibbs sampling:

> In the standard maximum likelihood (ML) approach for parameter estimation, typically only point estimates of the θ parameters are available. Therefore, RSEM also leverages a Bayesian approach in order to compute 95% credible intervals (CI) around the posterior mean estimate (PME). The prior distribution for the the θ parameters in the Bayesian approach are all Dirichlet prior distributions set to one. A Dirichlet distribution set to one is equivalent to a uniform distribution, meaning our prior is non-informative. As a consequence of this non-informative prior, the maximum a posteriori estimates of θ are equal to the ML estimates, with the added benefit of 95% CIs around the posterior mean. 

> RSEM computes the PMEs and 95% CIs for the Bayesian approach with a two-stage sampling process. First, a standard application of the collapsed Gibbs sampling algorithm is used to sample from the posterior distribution of count vectors. A sampled count vector represents the number of fragments that are mapped to each transcript. This allows us to easily construct CIs and a PME because we can simply compute quantiles of the sampled set of count vectors. For each iteration in the Gibbs sampling algorithm, the true mapping of each fragment is resampled given the current mappings of all other fragments. The initial mapping of each fragment is sampled according to the ML parameters computed by the EM algorithm. The algorithm is run with 1000 iterations. 

[figure0]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/figure_0.png
[figure1]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/figure_1.png
[figure2]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/figure_2.png
[figure3]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/gp_1.png
[figure4]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/gp_2.png
