
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
I chose to replicate sample progression discovery (SPD) method from "Discovering Biological Progression Underlying Microarray Samples" paper by Qiu et al. Given that I will soon need to run this process using the code provided from their [website](http://pengqiu.gatech.edu/software/SPD/index.html) on the iPSC data, I thought this would be a useful paper to reproduce to broaden my understanding of the process works. 

The goal of SPD is "to discover biological progression from gene expression microarray data." It accomplishes this by four steps: 
1) cluster genes into modules of co-expressed genes
2) construct minimum spanning tree (MST) for each module
3) select modules that supports common MSTs
4) reconstruct an overall MST based on all the genes of all the selected modules



### Cross-Validation Hyperparameter Tuning ###
The training data for each tissue was split using leave-one-out cross-validation in order to evaluate regularization parameters (alphas). We evaluated alphas across different scales (0.1, 1.0, 10.0, 100), and used the alpha for each tissue that resulted in the highest root mean-squared error. 

### Training/ Prediction ###
 
![figure1]


![figure2]


## Discussion ##
Based on the results in figures 1A and 2, the ridge regression model using thyroid cis-eQTL data is the best at predicting the age of the patient. In figure 1A, we see that across training-increment configurations, thyroid is the lowest line in the graph, indicating the lowest error on prediction in the testing set. As expected, across tissues the RMSE is minimal when the number of examples used in training is maximial, indicating that after 200 examples are seen, a good approximation of the true distribution is reached. In figure 2, we can see that in both configurations, using all samples and using the smallest common number of training examples across tissues, thyroid performs optimally. 

In terms of intepreting the results, we cannot take the RMSE to indicate the actual average error of prediction in terms of real-valued age (i.e just because our RMSE for thyroid is ~7, we cannot say that on average our prediction is off by 7 years). This is because of how the labels were encoded. Consider a case where a patient is aged 69. We place them in the 60 bucket, and then upon prediction, we predict their age using thyroid data to be 53. The error recorded is only 7, however, in actual biological space, the error is 16. Therefore, RMSE can be only used to compare the predictive accuracies of tissues, not to evaluate thyroid as a predictor of age itself. We could obviously improve our analysis here if the actual ages were provided to us upon sequencing, as there would no longer be noise associated with our labels. 

Given a longer timeframe, there are significant improvements that could be made. The GTEx consortium reports covariates for each of the patients, so we could "regress" these out of our gene expression matrices by obtaining the residuals of a linear model trained using the covariates on each gene. This could help remove batch effects and other covariates that typically overpower the signal in genomic data sets. Additionally, instead of using the top 1000 genes, we could instead cross-validate the number of genes we should use, which would likely improve our model in testing. Another idea I had was to do a polynomial feature mapping:


This feature mapping would allow you to capture squared effects of features, while also allowing you to capture feature interaction. It is possible that two genes have little signal by themselves, but multiplied, produce a better predictor of age. 

In terms of the biological interpretation of thyroid being the best predictor, what this indicates to us is that the genes expressed in the thyroid vary the most (linearly) across age. That is, the genes expressed in the thyroid either change in response to or drive the aging process. In the other tissues, our model does a poorer job of picking up this age-varying process of gene expression. This could be a result of thyroid in fact being the best predictor, or that our model is not sufficient to pick up the signal in the other tissues (i.e. non-linear signal). Consequently, a natural extension of this project is to use other non-linear models and see how this affects the RMSE across tissues. 

[figure1]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/figure_1.png
[figure2]: https://github.com/nkrishn9/Sample-Progression-Discovery/blob/master/figures/figure_2.png
