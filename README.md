# Sample Splitting in Bayesian Approaches to Estimate Conditional Average Treatment Effects
## Spring 2022, STA640 Final Project

Authors: Xige Huang, Youran Wu

### Abstract

To obtain valid and accurate estimation in treatment effects, two strategies of reducing model sensitivity are often considered, namely, balancing the covariates in the design stage and specifying a flexible outcome model. A Double-Robust (DR) approach in the frequentist context combines such two strategies and thus usually involves specification of both propensity score model and outcome model; recent literature has provided numerous flexible machine learning models as the outcome model. Additionally, sampling splitting is a key procedure in obtaining good estimates of the causal effects when using machine learning methods. Kennedy (2020) has proven that sample-splitting improves accuracy in estimating the Conditional Average Treatment Effect (CATE) using a DR-estimator. 

The question of interest is whether sample splitting and cross fitting helps in the Bayesian paradigm. To answer this question, we design a simulation study with different data generating processes and adopt three methods in the Bayesian paradigm, among which two involve two-step estimation are can be seen as analogue to a DR-approach. Specifically, we will consider Bayesian Additive Regression Trees (BART), BART with propensity score as a covariates (BART-ps), and Bayesian Causal Forest (BCF); the latter two include estimating the propensity score and building the potential outcome model. We consider several different but related data generating processes: homogeneous and heterogeneous treatment effects, low and high dimensional covariates, sparse and non-sparse covariates matrix, randomized and non-randomized treatment assignment.

### Ackownledgement

  Edward H Kennedy. (2019). Optimal doubly robust estimation of heterogeneous causal effects. ArXiv preprint. https://arxiv.org/abs/2004.14497

Our data generation process is inspired by and modified from codes from https://github.com/QuantLet/DataGenerationForCausalInference
