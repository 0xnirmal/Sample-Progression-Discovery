
# Advanced Topics in Genomic Data Analysis- Mini Project 2 #

# Running Instructions #
```
git clone https://github.com/nkrishn9/Sample-Progression-Discovery.git
chmod -x run.sh
source run.sh
```
This code is intended to run on the JHU CS ugrad cluster. Figures 1 and 2 were generated from the Jupyter notebook provided; however, the raw data used to generate figure 2 is printed when the py file is run. 

# Part II: Reproducibility #

I recommend taking a peek at the spd.ipynb (using the jupyter notebook command) to gain an understanding of the process, not the spd.py because the notebook contains significant markdown comments. 

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



## Discussion ##
Based on the results in figures 1A and 2, the ridge regression model using thyroid cis-eQTL data is the best at predicting the age of the patient. In figure 1A, we see that across training-increment configurations, thyroid is the lowest line in the graph, indicating the lowest error on prediction in the testing set. As expected, across tissues the RMSE is minimal when the number of examples used in training is maximial, indicating that after 200 examples are seen, a good approximation of the true distribution is reached. In figure 2, we can see that in both configurations, using all samples and using the smallest common number of training examples across tissues, thyroid performs optimally. 

In terms of intepreting the results, we cannot take the RMSE to indicate the actual average error of prediction in terms of real-valued age (i.e just because our RMSE for thyroid is ~7, we cannot say that on average our prediction is off by 7 years). This is because of how the labels were encoded. Consider a case where a patient is aged 69. We place them in the 60 bucket, and then upon prediction, we predict their age using thyroid data to be 53. The error recorded is only 7, however, in actual biological space, the error is 16. Therefore, RMSE can be only used to compare the predictive accuracies of tissues, not to evaluate thyroid as a predictor of age itself. We could obviously improve our analysis here if the actual ages were provided to us upon sequencing, as there would no longer be noise associated with our labels. 

Given a longer timeframe, there are significant improvements that could be made. The GTEx consortium reports covariates for each of the patients, so we could "regress" these out of our gene expression matrices by obtaining the residuals of a linear model trained using the covariates on each gene. This could help remove batch effects and other covariates that typically overpower the signal in genomic data sets. Additionally, instead of using the top 1000 genes, we could instead cross-validate the number of genes we should use, which would likely improve our model in testing. Another idea I had was to do a polynomial feature mapping:


This feature mapping would allow you to capture squared effects of features, while also allowing you to capture feature interaction. It is possible that two genes have little signal by themselves, but multiplied, produce a better predictor of age. 

In terms of the biological interpretation of thyroid being the best predictor, what this indicates to us is that the genes expressed in the thyroid vary the most (linearly) across age. That is, the genes expressed in the thyroid either change in response to or drive the aging process. In the other tissues, our model does a poorer job of picking up this age-varying process of gene expression. This could be a result of thyroid in fact being the best predictor, or that our model is not sufficient to pick up the signal in the other tissues (i.e. non-linear signal). Consequently, a natural extension of this project is to use other non-linear models and see how this affects the RMSE across tissues. 

[figure0]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/figure_0.png
[figure1]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/figure_1.png
[figure2]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/figure_2.png
