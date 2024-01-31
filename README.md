# Learning Theory from First Principles
Python code for the figures of the book "Learning Theory from First Principles" by Francis Bach


See available draft of the book [here](https://www.di.ens.fr/%7Efbach/ltfp_book.pdf).

**Contributors** : Francis Bach, Maria Bastiani, Gabriel Fiastre, Shane Hoeberichts, Camille Leempoels, Berné Nortier

## Table of Contents
[To Do](#to-do) | [Contribution guidelines](#guidelines) | [Python code](#python-code) | [Matlab code](#matlab-code) 
<br>

## To Do 

**Figures still to be done in python** :

- [ ] Figure 5.3
- [ ] Figure 7.3
- [ ] Figure 8.2
- [ ] Figures 9.1, 9.2

<br>

## Guidelines
<br>

## Python Code
The Python code is organized into individual notebooks for each chapter :
- Chapter 1: [Mathematical preliminaries](python/1_mathematical_preliminaries.ipynb)
- Chapter 2: [Introduction to supervised learning](python/2_introduction_supervised_learning.ipynb)
- Chapter 3: [Linear least-squares regres](python/3_least_squares.ipynb)
- Chapter 4: [Empirical risk minimization](python/4_empirical_risk_minimization.ipynb)
- Chapter 5: [Optimization](python/5_optimization.ipynb)
- Chapter 6: [Local averaging methods](python/6_local_averaging.ipynb)
- Chapter 7: [Kernel methods](python/7_kernels.ipynb)
- Chapter 8: [Sparse methods](python/8_model_selection.ipynb)
- Chapter 9: [Neural networks](python/9_neural_networks.ipynb)

<br>

## Matlab code
### Generic helper functions
- [affine_fit.m](/matlab/affine_fit.m)
- [sq_dist.m](/matlab/sq_dist.m)


### Chapter 1: Mathematical preliminaries

- [Figure 1.1](matlab/expectation_of_max.m) (expectation of maximum of Gaussian random variables)


### Chapter 2: Introduction to supervised learning

- [Figure 2.1](matlab/intro_supervised_learning/polynomial_regression.m) (polynomial regression with increasing orders - predictions)
- [Figure 2.2](matlab/polynomial_regression_with_replications.m) (polynomial regression with increasing orders - errors)


### Chapter 3: Linear least-squares regression

- [Figure 3.1](matlab/least_squares/OLS_polynomals_plots.m) (polynomial regression with varying number of observations)
- [Figure 3.2](matlab/least_squares/OLS_polynomals_rates.m) (convergence rate for polynomial regression)
- [Figure 3.3](matlab/least_squares/ridge_regression.m) (polynomial ridge regression)


### Chapter 4: Empirical risk minimization

- [Figure 4.1](matlab/empirical_risk_minimization/plot_losses_theory_class.m) (convex surrogates)
- [Figure 4.2](matlab/empirical_risk_minimization/plot_binary_classification_scores.m) (optimal score functions for Gaussian class-conditional densities)


### Chapter 5: Optimization

- [Figure 5.1](matlab/optimization/grad_descent_comparison.m) (gradient descent on two least-squares problems)
- [Figure 5.2](matlab/optimization/hinge_sgd.m) (comparison of step-sizes for SGD for the support vector machine)
- [Figure 5.3](matlab/optimization/logistic_sgd_saga.m) (comparison of step-sizes for SGD for logistic regression)


### Chapter 6: Local averaging methods

- [Figure 6.2](matlab/local_averaging/regressogram.m) (regressogram in one dimension)
- [Figure 6.3](matlab/local_averaging/knn.m) (k-nearest neighbor in one dimension)
- [Figure 6.4](matlab/local_averaging/nadaraya.m) (Nadaraya-Watson in one dimension)
- [Figure 6.5](matlab/local_averaging/all_learning_curves.m) (learning curves for local averaging)
- [Figure 6.6](matlab/local_averaging/regressogram_poly.m) (locally linear partitioning estimate)


### Chapter 7: Kernel methods

- [Figure 7.2](matlab/kernels/interpolate_kernel.m) (minimum norm interpolator)
- [Figure 7.3](matlab/kernels/kernel_simulations_1d_single_plot.m) (comparison of kernels)


### Chapter 8: Sparse methods

- [Figure 8.1](matlab/model_selection/path_lasso.m) (regularization path)
- [Figure 8.2](matlab/model_selection/model_selection.m) (comparison of estimators) + [script_model_selection.m](matlab/model_selection/script_model_selection.m) +  [script_model_selectionROT.m](matlab/model_selection/script_model_selectionROT.m)


### Chapter 9: Neural networks

- [Figure 9.1](matlab/neural_networks/neural_networks_1d_testerrors.m) (global convergence for different numbers of neurons) + [launch_training_relu_nn.m](matlab/neural_networks/launch_training_relu_nn.m)
- [Figure 9.2](matlab/neural_networks/random_features_interpolation.m) (random features - kernels)
- [Figure 9.3](matlab/neural_networks/neural_networks_1d.m) (neural networks fitting)
