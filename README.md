# Scaling-of-whole-brain-resting-state-dynamics

-PRG_function:

This function implements the Phenomenological Renormalization Group
(PRG) method, introduced by Meshulam et al. (2018, 2019). It also can
implement the extended version of PRG based on connectivity.

Within this method, the collective activity is iteratively coarse-grained 
by grouping maximally correlated variables (or maximally coupled variables 
in the connectivity-based case). At each coarse-graining step 
k=0,1,…,kmax, clusters of size K=2^k are built, resulting in a system of 
N/K coarse-grained variables and successively ignoring degrees of
freedom. This code computes several observables of the coarse-grained 
variables as a function of K.   
